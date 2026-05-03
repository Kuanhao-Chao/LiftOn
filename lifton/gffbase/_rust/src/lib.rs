// ---------------------------------------------------------------------------
// Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
// ---------------------------------------------------------------------------
//! gffbase native core — PyO3 entry point.
//!
//! Exposes a single `parse_file(path, force_dialect_check=False, checklines=10)`
//! callable plus a `parse_bytes(data, ...)` callable. Both yield Python tuples
//! with the canonical 11-tuple shape that `gffbase.parser` consumes:
//!
//!     (seqid, source, featuretype, start, end, score, strand, frame,
//!      attributes_blob, attributes_pairs, extra)
//!
//! `attributes_pairs` is a list[(key, value, idx)]; `idx` preserves multi-value
//! ordering. `start`/`end` are int or None (`.` becomes None). `attributes_blob`
//! is the raw col-9 bytes for byte-faithful round-trip.

use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList, PyTuple};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::create_exception;

mod dialect;
mod attributes;
mod escape;
mod parser;
mod validate;

use parser::{ParseOptions, RecordIter, FileSource};
use validate::GffError;

// Phase 16: descriptive parser errors. Subclassing `PyValueError` keeps
// legacy `pytest.raises(ValueError)` calls working while letting users
// catch the more specific `GFFFormatError` for rich line-numbered context.
create_exception!(_native, GFFFormatError, PyValueError);

/// Convert a structured `GffError` into a Python `GFFFormatError` with
/// `.line_no`, `.kind`, and `.message` attributes attached on the
/// exception instance.
fn gff_error_to_py(py: Python<'_>, e: GffError) -> PyErr {
    let msg = format!("line {}: {}", e.line_no, e.message);
    let err = GFFFormatError::new_err(msg);
    // Attach structured fields onto the exception's `value` object.
    if let Ok(value_obj) = err.value_bound(py).clone().into_any().getattr("__class__") {
        // mypy / runtime path — set on the bound value, not the type.
        let _ = value_obj;
    }
    if let Some(inst) = err.value_bound(py).extract::<Bound<'_, PyAny>>().ok() {
        let _ = inst.setattr("line_no", e.line_no);
        let _ = inst.setattr("kind", e.kind.as_str());
        let _ = inst.setattr("message", e.message.clone());
    }
    err
}

/// Parse a path (plain text or .gz). Yields one tuple per feature.
#[pyfunction]
#[pyo3(signature = (path, checklines=10, force_dialect_check=false, force_gff=false, strict=true))]
fn parse_file(
    py: Python<'_>,
    path: &str,
    checklines: usize,
    force_dialect_check: bool,
    force_gff: bool,
    strict: bool,
) -> PyResult<PyObject> {
    let opts = ParseOptions {
        checklines,
        force_dialect_check,
        force_gff,
        strict,
    };
    let source = FileSource::open(path)
        .map_err(|e| PyIOError::new_err(format!("could not open {}: {}", path, e)))?;
    let iter = RecordIter::new(source, opts)
        .map_err(|e| PyValueError::new_err(format!("parser error: {}", e)))?;
    let py_iter = PyRecordIterator { inner: Some(iter) };
    Ok(Py::new(py, py_iter)?.into_any())
}

/// Parse an in-memory byte buffer. Yields one tuple per feature.
#[pyfunction]
#[pyo3(signature = (data, checklines=10, force_dialect_check=false, force_gff=false, strict=true))]
fn parse_bytes(
    py: Python<'_>,
    data: &[u8],
    checklines: usize,
    force_dialect_check: bool,
    force_gff: bool,
    strict: bool,
) -> PyResult<PyObject> {
    let opts = ParseOptions {
        checklines,
        force_dialect_check,
        force_gff,
        strict,
    };
    let source = FileSource::from_bytes(data.to_vec());
    let iter = RecordIter::new(source, opts)
        .map_err(|e| PyValueError::new_err(format!("parser error: {}", e)))?;
    let py_iter = PyRecordIterator { inner: Some(iter) };
    Ok(Py::new(py, py_iter)?.into_any())
}

/// Detect the dialect of an input by reading the first `checklines` features
/// and return it as a Python dict. Useful for tests and for the API layer.
#[pyfunction]
#[pyo3(signature = (path, checklines=10))]
fn detect_dialect(py: Python<'_>, path: &str, checklines: usize) -> PyResult<PyObject> {
    let source = FileSource::open(path)
        .map_err(|e| PyIOError::new_err(format!("could not open {}: {}", path, e)))?;
    let opts = ParseOptions {
        checklines,
        force_dialect_check: false,
        force_gff: false,
        // Dialect detection is non-strict by design: malformed lines in
        // the first `checklines` get skipped without poisoning detection.
        strict: false,
    };
    let iter = RecordIter::new(source, opts)
        .map_err(|e| PyValueError::new_err(format!("parser error: {}", e)))?;
    Ok(dialect_to_pydict(py, iter.dialect())?)
}

#[pyclass]
struct PyRecordIterator {
    inner: Option<RecordIter>,
}

#[pymethods]
impl PyRecordIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>, py: Python<'_>) -> PyResult<Option<PyObject>> {
        let iter = match slf.inner.as_mut() {
            Some(i) => i,
            None => return Ok(None),
        };
        match iter.next() {
            Some(Ok(rec)) => Ok(Some(record_to_pytuple(py, &rec)?)),
            // Structured GffError → Python GFFFormatError with .line_no /
            // .kind / .message attributes.
            Some(Err(e)) => Err(gff_error_to_py(py, e)),
            // Do NOT null out `inner` here: callers expect to be able to read
            // `.dialect()` and `.directives()` after the iterator is exhausted.
            None => Ok(None),
        }
    }

    /// List of `GffError`s collected during a non-strict run. Each item
    /// is a dict with keys `line_no`, `kind`, and `message`. Empty when
    /// the iterator was created with `strict=True` (the default), since
    /// errors propagate via `__next__` instead in that mode.
    fn warnings(&self, py: Python<'_>) -> PyResult<PyObject> {
        let list = PyList::empty_bound(py);
        if let Some(it) = &self.inner {
            for w in it.warnings() {
                let d = PyDict::new_bound(py);
                d.set_item("line_no", w.line_no)?;
                d.set_item("kind", w.kind.as_str())?;
                d.set_item("message", w.message.as_str())?;
                list.append(d)?;
            }
        }
        Ok(list.into_any().unbind())
    }

    /// Return the inferred dialect dict. Available after iteration begins or
    /// after the peek phase completes.
    fn dialect(&self, py: Python<'_>) -> PyResult<PyObject> {
        match &self.inner {
            Some(it) => dialect_to_pydict(py, it.dialect()),
            None => Ok(py.None()),
        }
    }

    /// Return the list of `##` directives encountered so far.
    fn directives(&self, py: Python<'_>) -> PyResult<PyObject> {
        let list = PyList::empty_bound(py);
        if let Some(it) = &self.inner {
            for d in it.directives() {
                list.append(d.as_str())?;
            }
        }
        Ok(list.into_any().unbind())
    }
}

fn record_to_pytuple(py: Python<'_>, rec: &parser::Record) -> PyResult<PyObject> {
    let start_obj: PyObject = match rec.start {
        Some(v) => v.into_py(py),
        None => py.None(),
    };
    let end_obj: PyObject = match rec.end {
        Some(v) => v.into_py(py),
        None => py.None(),
    };
    let attrs_blob = PyBytes::new_bound(py, &rec.attributes_blob);
    let pairs = PyList::empty_bound(py);
    for (k, v, idx) in &rec.attributes_pairs {
        let t = PyTuple::new_bound(py, &[k.into_py(py), v.into_py(py), (*idx as i64).into_py(py)]);
        pairs.append(t)?;
    }
    let extra_list = PyList::empty_bound(py);
    for e in &rec.extra {
        extra_list.append(e.as_str())?;
    }
    let tup = PyTuple::new_bound(
        py,
        &[
            rec.seqid.clone().into_py(py),
            rec.source.clone().into_py(py),
            rec.featuretype.clone().into_py(py),
            start_obj,
            end_obj,
            rec.score.clone().into_py(py),
            rec.strand.clone().into_py(py),
            rec.frame.clone().into_py(py),
            attrs_blob.into_any().unbind(),
            pairs.into_any().unbind(),
            extra_list.into_any().unbind(),
        ],
    );
    Ok(tup.into_any().unbind())
}

fn dialect_to_pydict(py: Python<'_>, d: &dialect::Dialect) -> PyResult<PyObject> {
    let dict = pyo3::types::PyDict::new_bound(py);
    dict.set_item("fmt", d.fmt_str())?;
    dict.set_item("field separator", d.field_separator.as_str())?;
    dict.set_item("keyval separator", d.keyval_separator.to_string())?;
    dict.set_item("multival separator", d.multival_separator.to_string())?;
    dict.set_item("leading semicolon", d.leading_semicolon)?;
    dict.set_item("trailing semicolon", d.trailing_semicolon)?;
    dict.set_item("quoted GFF2 values", d.quoted_gff2_values)?;
    dict.set_item("repeated keys", d.repeated_keys)?;
    dict.set_item("semicolon in quotes", d.semicolon_in_quotes)?;
    let order = PyList::empty_bound(py);
    for k in &d.order {
        order.append(k.as_str())?;
    }
    dict.set_item("order", order)?;
    Ok(dict.into_any().unbind())
}

#[pymodule]
fn _native(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_file, m)?)?;
    m.add_function(wrap_pyfunction!(parse_bytes, m)?)?;
    m.add_function(wrap_pyfunction!(detect_dialect, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    // Phase 16 — expose the descriptive Python exception type.
    m.add("GFFFormatError", py.get_type_bound::<GFFFormatError>())?;
    Ok(())
}
