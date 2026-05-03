// ---------------------------------------------------------------------------
// Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
// ---------------------------------------------------------------------------
//! Streaming GFF3/GTF record parser.
//!
//! Reads from a file (mmap'd plain text or gzip'd) or from an in-memory byte
//! slice and yields `Record`s lazily. Splits lines and tab fields with `memchr`
//! for SIMD-accelerated scanning. Stops feature emission at the GFF3 `##FASTA`
//! sentinel.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use memchr::memchr;

use crate::attributes::parse_attributes;
use crate::dialect::{self, Dialect};
use crate::validate::{validate_attributes_pairs, validate_fields, ErrorKind, GffError};

#[derive(Debug, Clone)]
pub struct Record {
    pub seqid: String,
    pub source: String,
    pub featuretype: String,
    pub start: Option<i64>,
    pub end: Option<i64>,
    pub score: String,
    pub strand: String,
    pub frame: String,
    pub attributes_blob: Vec<u8>,
    pub attributes_pairs: Vec<(String, String, u16)>,
    pub extra: Vec<String>,
}

pub struct ParseOptions {
    pub checklines: usize,
    pub force_dialect_check: bool,
    pub force_gff: bool,
    /// When true (default), the iterator yields `Err(GffError)` on the
    /// first malformed line and the caller is expected to surface it as
    /// a Python `GFFFormatError`. When false, the iterator silently
    /// drops malformed lines and pushes a `GffError` onto its
    /// `warnings` vector — the caller can inspect them via
    /// `RecordIter::warnings()`.
    pub strict: bool,
}

/// A simple text source: either a Vec<u8> we own or a Read trait object.
pub enum FileSource {
    Bytes(Vec<u8>),
    Stream(Box<dyn Read + Send>),
}

impl FileSource {
    pub fn open<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let path = path.as_ref();
        let f = File::open(path)?;
        let is_gz = path.extension().map(|e| e == "gz").unwrap_or(false);
        if is_gz {
            Ok(FileSource::Stream(Box::new(BufReader::new(MultiGzDecoder::new(f)))))
        } else {
            // For plain text we read fully into memory. mmap could be a future
            // optimization but adds platform complexity; the parser is already
            // fast enough that the bottleneck moves elsewhere.
            let mut buf = Vec::new();
            BufReader::new(f).read_to_end(&mut buf)?;
            Ok(FileSource::Bytes(buf))
        }
    }

    pub fn from_bytes(b: Vec<u8>) -> Self {
        FileSource::Bytes(b)
    }
}

pub struct RecordIter {
    buf: Vec<u8>,            // entire input materialized
    pos: usize,
    dialect: Dialect,
    directives: Vec<String>,
    fasta_reached: bool,
    line_no: usize,
    strict: bool,
    warnings: Vec<GffError>,
}

impl RecordIter {
    pub fn new(source: FileSource, opts: ParseOptions) -> Result<Self, String> {
        let buf = match source {
            FileSource::Bytes(b) => b,
            FileSource::Stream(mut r) => {
                let mut v = Vec::new();
                r.read_to_end(&mut v).map_err(|e| e.to_string())?;
                v
            }
        };
        let mut iter = RecordIter {
            buf,
            pos: 0,
            dialect: Dialect::default(),
            directives: Vec::new(),
            fasta_reached: false,
            line_no: 0,
            strict: opts.strict,
            warnings: Vec::new(),
        };
        iter.peek_dialect(&opts);
        Ok(iter)
    }

    /// Errors collected when running with `strict=false`. Empty in
    /// strict mode (errors propagate via `Iterator::next` instead).
    pub fn warnings(&self) -> &[GffError] {
        &self.warnings
    }

    pub fn dialect(&self) -> &Dialect {
        &self.dialect
    }

    pub fn directives(&self) -> &[String] {
        &self.directives
    }

    /// Peek up to `checklines` features (without consuming them) to compute
    /// the dialect. We snapshot `pos`, walk forward, then reset.
    fn peek_dialect(&mut self, opts: &ParseOptions) {
        let saved_pos = self.pos;
        let saved_line = self.line_no;
        let saved_directives = self.directives.clone();
        let saved_fasta = self.fasta_reached;

        let mut samples: Vec<Dialect> = Vec::new();
        let limit = if opts.force_dialect_check { usize::MAX } else { opts.checklines };
        while samples.len() < limit {
            match self.next_raw_record() {
                Some(Ok((_, _, _, _, _, _, _, _, blob, _))) => {
                    let blob_str = std::str::from_utf8(&blob).unwrap_or("");
                    let (_pairs, obs) = parse_attributes(blob_str);
                    samples.push(obs);
                }
                Some(Err(_)) => break,
                None => break,
            }
        }
        self.dialect = dialect::choose(&samples);

        // Reset.
        self.pos = saved_pos;
        self.line_no = saved_line;
        self.directives = saved_directives;
        self.fasta_reached = saved_fasta;

        // force_gff overrides format detection.
        if opts.force_gff {
            self.dialect.fmt = crate::dialect::Format::Gff3;
            self.dialect.keyval_separator = '=';
        }
    }

    /// Read the next non-comment, non-blank line as 9-or-more tab fields.
    /// Returns the raw fields plus the attributes blob (col 9 raw bytes).
    /// Errors are structured `GffError` values carrying the offending
    /// line number; `lib.rs` converts them to Python `GFFFormatError`s.
    #[allow(clippy::type_complexity)]
    fn next_raw_record(
        &mut self,
    ) -> Option<Result<
        (
            String, // seqid
            String, // source
            String, // featuretype
            Option<i64>, // start
            Option<i64>, // end
            String, // score
            String, // strand
            String, // frame
            Vec<u8>, // blob
            Vec<String>, // extra
        ),
        GffError,
    >> {
        loop {
            if self.fasta_reached || self.pos >= self.buf.len() {
                return None;
            }
            // We materialize the line into an owned Vec to drop the borrow on
            // self.buf before we touch other &self fields below.
            let line_owned: Vec<u8> = match self.read_line() {
                Some(slice) => slice.to_vec(),
                None => return None,
            };
            let cur_line_no = self.line_no;
            if line_owned.is_empty() {
                continue;
            }
            let line = line_owned.as_slice();
            // Directive / comment handling.
            if line.starts_with(b"##") {
                let s = std::str::from_utf8(line).unwrap_or("").to_string();
                if s.starts_with("##FASTA") {
                    self.fasta_reached = true;
                    return None;
                }
                self.directives.push(s);
                continue;
            }
            if line.starts_with(b"#") {
                continue;
            }
            // Trim trailing \r (Windows files).
            let line = trim_cr(line);
            // Tab split.
            let fields = split_tabs(line);
            if fields.len() < 9 {
                return Some(Err(GffError::new(
                    cur_line_no,
                    ErrorKind::TooFewFields,
                    format!(
                        "expected at least 9 tab-separated fields, found {}",
                        fields.len()
                    ),
                )));
            }
            let seqid = bytes_to_string(fields[0]);
            let source = bytes_to_string(fields[1]);
            let featuretype = bytes_to_string(fields[2]);
            // Coordinate parsing distinguishes "valid `.` / empty"
            // (yielding `None`) from "non-numeric trash" (a hard error).
            let start = match parse_coord_strict(fields[3]) {
                Ok(v) => v,
                Err(_) => {
                    return Some(Err(GffError::new(
                        cur_line_no,
                        ErrorKind::InvalidCoordinate,
                        format!(
                            "start coordinate is not an integer: {:?}",
                            std::str::from_utf8(fields[3]).unwrap_or("<non-utf8>")
                        ),
                    )));
                }
            };
            let end = match parse_coord_strict(fields[4]) {
                Ok(v) => v,
                Err(_) => {
                    return Some(Err(GffError::new(
                        cur_line_no,
                        ErrorKind::InvalidCoordinate,
                        format!(
                            "end coordinate is not an integer: {:?}",
                            std::str::from_utf8(fields[4]).unwrap_or("<non-utf8>")
                        ),
                    )));
                }
            };
            let score = bytes_to_string(fields[5]);
            let strand = bytes_to_string(fields[6]);
            let frame = bytes_to_string(fields[7]);
            let blob = fields[8].to_vec();
            let mut extra: Vec<String> = Vec::new();
            for ex in fields.iter().skip(9) {
                extra.push(bytes_to_string(ex));
            }
            return Some(Ok((seqid, source, featuretype, start, end, score, strand, frame, blob, extra)));
        }
    }

    /// Return a slice covering the next line (without the terminating LF).
    /// Increments `line_no` internally so the caller can keep borrowing the
    /// returned slice without conflicting with `&mut self`.
    fn read_line(&mut self) -> Option<&[u8]> {
        if self.pos >= self.buf.len() {
            return None;
        }
        let start = self.pos;
        let nl = memchr(b'\n', &self.buf[start..]);
        let (end, advance) = match nl {
            Some(i) => (start + i, i + 1),
            None => (self.buf.len(), self.buf.len() - start),
        };
        self.pos += advance;
        self.line_no += 1;
        Some(&self.buf[start..end])
    }
}

impl Iterator for RecordIter {
    type Item = Result<Record, GffError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line_no_for_err = self.line_no.saturating_add(1);
            let raw = match self.next_raw_record() {
                Some(Ok(r)) => r,
                Some(Err(e)) => {
                    if self.strict {
                        return Some(Err(e));
                    }
                    self.warnings.push(e);
                    // Try the next line.
                    continue;
                }
                None => return None,
            };
            let (seqid, source, featuretype, start, end, score, strand, frame, blob, extra) = raw;

            let is_gtf = matches!(self.dialect.fmt, crate::dialect::Format::Gtf);
            if let Err(e) = validate_fields(
                self.line_no,           // line we just consumed
                &seqid,
                &featuretype,
                start,
                end,
                &score,
                &strand,
                &frame,
                &blob,
                is_gtf,
            ) {
                if self.strict {
                    return Some(Err(e));
                }
                self.warnings.push(e);
                continue;
            }

            let blob_str = std::str::from_utf8(&blob).unwrap_or("");
            let (pairs, _obs) = parse_attributes(blob_str);
            // Post-parse attribute structure check (handles both GFF3 and
            // GTF correctly because it inspects what the parser produced).
            if let Err(e) = validate_attributes_pairs(
                self.line_no, pairs.len(), &blob, is_gtf,
            ) {
                if self.strict {
                    return Some(Err(e));
                }
                self.warnings.push(e);
                continue;
            }
            // Suppress unused-warning under release builds (line_no_for_err is a
            // defensive snapshot — not used on the happy path).
            let _ = line_no_for_err;
            return Some(Ok(Record {
                seqid,
                source,
                featuretype,
                start,
                end,
                score,
                strand,
                frame,
                attributes_blob: blob,
                attributes_pairs: pairs,
                extra,
            }));
        }
    }
}

fn split_tabs(line: &[u8]) -> Vec<&[u8]> {
    let mut out: Vec<&[u8]> = Vec::with_capacity(9);
    let mut start = 0usize;
    let mut i = 0usize;
    while i < line.len() {
        if line[i] == b'\t' {
            out.push(&line[start..i]);
            start = i + 1;
        }
        i += 1;
    }
    out.push(&line[start..]);
    out
}

fn trim_cr(line: &[u8]) -> &[u8] {
    if let Some((&b'\r', rest)) = line.split_last().map(|(b, r)| (b, r)) {
        rest
    } else {
        line
    }
}

fn bytes_to_string(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

fn parse_coord(b: &[u8]) -> Option<i64> {
    if b == b"." || b.is_empty() {
        return None;
    }
    let s = std::str::from_utf8(b).ok()?;
    s.parse::<i64>().ok()
}

/// Like `parse_coord` but distinguishes "valid `.` / empty" from
/// "non-numeric trash". Returns `Ok(None)` for `.` / empty, `Ok(Some(n))`
/// for a real integer, and `Err(())` for anything else. The caller turns
/// the `Err` into a structured `GffError`.
fn parse_coord_strict(b: &[u8]) -> Result<Option<i64>, ()> {
    if b == b"." || b.is_empty() {
        return Ok(None);
    }
    let s = std::str::from_utf8(b).map_err(|_| ())?;
    s.parse::<i64>().map(Some).map_err(|_| ())
}
