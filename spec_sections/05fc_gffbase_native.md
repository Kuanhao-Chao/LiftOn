### 5.6.C gffbase native Rust extension & pure-Python fallback

This subsection specifies the **column-9 / record-line parsing layer** of gffbase: the compiled Rust extension `gffbase._native` (sources under `lifton/gffbase/_rust/src/`), the pure-Python oracle/fallback `gffbase._pyfallback` (`parser.py`, `attributes.py`), and the selection shim in `gffbase/parser.py`. The two implementations are required by the test suite to be **bit-for-bit equivalent** in everything they yield (record tuples, dialect dicts, directive lists, warning dicts, raised errors). A reimplementer must be able to build **either** path from this text alone and have them agree on every input.

**Build status in this checkout:** No compiled `_native*.so` exists in `lifton/gffbase/`. The `from . import _native` import in `gffbase/parser.py:17` and `gffbase/__init__.py:24` therefore raises `ImportError`, `_NATIVE == False`, and **the pure-Python fallback is the live path**. The Rust sources are nonetheless the *normative reference* for behaviour; the fallback is documented as "mirror of `<rust file>`" throughout. The Cargo package is `gffbase-core` version `0.1.0` (`_rust/Cargo.toml`); `_native.__version__` is `env!("CARGO_PKG_VERSION")` = `"0.1.0"` (`lib.rs:252`).

#### 5.6.C.1 Engine selection: native vs fallback

Three call sites decide which engine runs.

**(a) Module-load probe** — `gffbase/parser.py:16-22`:

```python
try:
    from . import _native as _rust   # the compiled .so
    _NATIVE = True
except ImportError:
    _rust = None
    _NATIVE = False
```

`native_available()` (`parser.py:25-27`) returns `_NATIVE`. The probe runs once at import; there is no re-probe.

**(b) Per-call engine resolution** — `_resolve_engine(engine)` (`parser.py:74-83`):

| `engine` argument | Result | Condition |
|---|---|---|
| `None` or `"auto"` (default) | `"rust"` if `_NATIVE` else `"python"` | always |
| `"rust"` | `"rust"` | raises `RuntimeError("Rust extension not built. Run `maturin develop --release` or pass engine='python'.")` if `_NATIVE` is False |
| `"python"` | `"python"` | always |
| anything else | — | raises `ValueError(f"unknown engine '{engine}'")` |

The three public functions `parse_gff` (`parser.py:86`), `parse_bytes` (`parser.py:127`), `detect_dialect` (`parser.py:156`) each call `_resolve_engine` then branch: on `"rust"` they call `_rust.parse_file` / `_rust.parse_bytes` / `_rust.detect_dialect` and wrap the result in `_Iterator(it, native=True)`; on `"python"` they call `_pyparser.parse_file` / `parse_bytes` / `detect_dialect` and wrap in `_Iterator(it, native=False)`. `detect_dialect` returns a plain dict from either engine, not an iterator.

**(c) Exception-class rebinding** — `gffbase/__init__.py:23-26`:

```python
try:
    from ._native import GFFFormatError       # Rust PyO3 exception type
except ImportError:
    GFFFormatError = _PyGFFFormatError         # pure-Python class
```

So `gffbase.GFFFormatError` is *whichever class the live engine will actually raise*. Both inherit from `ValueError` (Rust via `create_exception!(_native, GFFFormatError, PyValueError)` at `lib.rs:34`; Python via `class GFFFormatError(ValueError)` at `exceptions.py:13`), so `except ValueError` catches both. **Gotcha:** `_pyfallback/parser.py` must *not* import `GFFFormatError` at module top, because at the time that submodule loads the `__init__` rebinding has not happened. It resolves the class lazily through `_gff_format_error_class()` (`_pyfallback/parser.py:23-38`), which tries `from .._native import GFFFormatError` first and falls back to `from ..exceptions import GFFFormatError`. This guarantees a caller's `except` clause catches the *same* class regardless of which path raised.

#### 5.6.C.2 The `_Iterator` adapter

`gffbase/parser.py:30-71` wraps either underlying iterator and normalises the surface. `__slots__ = ("_inner", "_native")`.

| Method / property | Behaviour |
|---|---|
| `__iter__()` | returns `self` |
| `__next__()` | `item = next(self._inner)`. If `self._native`: return `ParsedFeature.from_tuple(item)` (Rust yields 11-tuples). Else: return `item` directly (fallback already yields `ParsedFeature`). |
| `dialect()` | delegates to `self._inner.dialect()` → dict |
| `directives()` | `list(self._inner.directives())` |
| `warnings` (property) | if `self._inner` has a `warnings` attribute, return `list(w() if callable(w) else w)` — handles the Rust method form and the Python property form uniformly; else `[]` |

So the *only* externally visible difference between engines is the tuple-vs-object yield, which `from_tuple` erases. Downstream code (`ingest.from_file`, `Annotation`) sees identical `ParsedFeature` streams.

#### 5.6.C.3 The native extension's responsibilities and Python-facing contract

`gffbase._native` (PyO3 module, `lib.rs:247-256`) exposes exactly three free functions, one iterator class, and one exception type.

**Module-level exports** (`lib.rs:248-255`):

| Name | Kind | Signature (Python view) |
|---|---|---|
| `parse_file` | function | `parse_file(path, checklines=10, force_dialect_check=False, force_gff=False, strict=True) -> PyRecordIterator` |
| `parse_bytes` | function | `parse_bytes(data, checklines=10, force_dialect_check=False, force_gff=False, strict=True) -> PyRecordIterator` |
| `detect_dialect` | function | `detect_dialect(path, checklines=10) -> dict` |
| `__version__` | str | `"0.1.0"` |
| `GFFFormatError` | type | subclass of `ValueError`, with `.line_no`, `.kind`, `.message` |

**`PyRecordIterator`** (`lib.rs:124-188`) is the object returned by `parse_file` / `parse_bytes`. Its Python-facing methods:

| Method | Returns | Notes |
|---|---|---|
| `__iter__()` | self | |
| `__next__()` | the **11-tuple** (or raises `GFFFormatError`, or `StopIteration` at exhaustion) | `lib.rs:135-149` |
| `warnings()` | `list[dict]` each `{line_no, kind, message}` | non-empty only when `strict=False`; `lib.rs:155-167` |
| `dialect()` | dict | the inferred dialect; `lib.rs:171-176` |
| `directives()` | `list[str]` | the `##…` lines seen so far; `lib.rs:179-187` |

**Gotcha (`lib.rs:146-147`):** `__next__` does **not** null out the internal iterator on `None`. Callers are expected to read `.dialect()` and `.directives()` *after* exhaustion. The fallback replicates this (its generator's `directives` list is mutated in place; see §5.6.C.6).

**The canonical 11-tuple** yielded per feature (`lib.rs:190-226`, built by `record_to_pytuple`), in order:

| # | Field | Python type | Meaning |
|---|---|---|---|
| 0 | `seqid` | str | column 1 |
| 1 | `source` | str | column 2 |
| 2 | `featuretype` | str | column 3 |
| 3 | `start` | int or `None` | column 4; `.`/empty → `None` |
| 4 | `end` | int or `None` | column 5; `.`/empty → `None` |
| 5 | `score` | str | column 6 (kept as string, incl. `"."`) |
| 6 | `strand` | str | column 7 |
| 7 | `frame` | str | column 8 (a.k.a. phase) |
| 8 | `attributes_blob` | bytes | raw column-9 bytes, UTF-8, **verbatim** for byte-faithful round-trip |
| 9 | `attributes_pairs` | list[tuple(str,str,int)] | `(key, value, multivalue_index)` triples |
| 10 | `extra` | list[str] | tab-separated columns past column 9 (rare) |

`ParsedFeature.from_tuple` (`feature.py:58-86`) unpacks this exact order, coercing `attributes_blob` to `bytes`, `attributes_pairs` to `[(k, v, int(i))]`, and `extra` to `list`. The dataclass field order (`feature.py:23-39`) **is** the tuple order. A reimplementer's native module MUST emit the tuple in precisely this order.

The Rust `Record` struct (`parser.rs:22-35`) mirrors these fields; note `attributes_pairs: Vec<(String, String, u16)>` — the multivalue index is a `u16` (capped at 65535 multi-values per key); the Python side casts via `(*idx as i64)` at `lib.rs:202`.

**Error conversion** — `gff_error_to_py` (`lib.rs:39-53`): a structured Rust `GffError {line_no, kind, message}` becomes a `GFFFormatError` whose message is `format!("line {}: {}", e.line_no, e.message)` and which has `.line_no` (int), `.kind` (str, from `ErrorKind::as_str()`), `.message` (the bare message without the "line N:" prefix) attached as instance attributes. **Gotcha:** the displayed message string includes the `line N:` prefix, but the `.message` attribute does **not** — it holds the bare message only.

The `detect_dialect` function (`lib.rs:106-122`) builds a `RecordIter` with `force_dialect_check=False`, `force_gff=False`, **`strict=False`** (malformed lines in the peek window are skipped, never poison detection) and returns `dialect_to_pydict(iter.dialect())`. It reads only the peek window; it does not iterate records.

`dialect_to_pydict` (`lib.rs:228-245`) emits a dict with keys in this exact spelling (note the spaces): `"fmt"`, `"field separator"`, `"keyval separator"`, `"multival separator"`, `"leading semicolon"`, `"trailing semicolon"`, `"quoted GFF2 values"`, `"repeated keys"`, `"semicolon in quotes"`, `"order"`. `keyval_separator` and `multival_separator` are single chars stringified.

#### 5.6.C.4 Native parsing algorithm (`parser.rs`)

The Rust parser is the normative reference. The fallback must reproduce it line-for-line.

**Input sourcing** — `FileSource` (`parser.rs:51-76`):
- `FileSource::open(path)`: if the path extension is `gz`, wrap in `MultiGzDecoder` (multi-member gzip) over a `BufReader`; otherwise read the **entire file into memory** as a `Vec<u8>` (no mmap; comment at `parser.rs:65-66` notes this is deliberate).
- `FileSource::from_bytes(b)`: owns the byte vector directly.
- `RecordIter::new` (`parser.rs:90-111`) fully materialises any stream into `buf: Vec<u8>` so the whole input is in memory before iteration. It then immediately calls `peek_dialect(&opts)`.

**`RecordIter` state** (`parser.rs:78-87`): `buf` (whole input), `pos` (cursor), `dialect`, `directives: Vec<String>`, `fasta_reached: bool`, `line_no`, `strict`, `warnings: Vec<GffError>`.

**Line reading** — `read_line` (`parser.rs:273-286`): from `pos`, find the next `\n` via `memchr`; the line slice excludes the LF; advance `pos` past the LF (or to EOF if none); increment `line_no`. Returns `None` at EOF.

**`next_raw_record`** (`parser.rs:168-268`) — the per-line tokenizer, looping until a feature line or EOF:
1. If `fasta_reached` or `pos >= buf.len()` → return `None`.
2. Read one line (owned copy via `.to_vec()` to drop the `buf` borrow). Snapshot `cur_line_no = self.line_no`.
3. If the line is empty → `continue` (skip blank lines).
4. If line starts with `b"##"`: decode to UTF-8 string. If it starts with `"##FASTA"` → set `fasta_reached = true`, return `None` (stops all further feature emission). Otherwise push the directive string onto `directives` and `continue`.
5. If line starts with a single `b"#"` (a plain comment, not `##`) → `continue` (not recorded).
6. Trim a single trailing `\r` (`trim_cr`, `parser.rs:377-383`) for Windows CRLF files.
7. Split on **tab only** (`split_tabs`, `parser.rs:362-375`) — empty fields are preserved; the final field captures the remainder.
8. If `fields.len() < 9` → return `Err(GffError{kind: TooFewFields, msg: "expected at least 9 tab-separated fields, found {n}"})` (line `cur_line_no`).
9. Columns 0,1,2 (`seqid`, `source`, `featuretype`) → `String::from_utf8_lossy` (`bytes_to_string`, `parser.rs:385-387`).
10. Columns 3,4 (`start`, `end`) → `parse_coord_strict` (`parser.rs:401-407`): `.` or empty → `Ok(None)`; a valid i64 → `Ok(Some(n))`; anything else → `Err(())` which the caller turns into `GffError{kind: InvalidCoordinate, msg: "start/end coordinate is not an integer: {…}"}`.
11. Columns 5,6,7 (`score`, `strand`, `frame`) → `bytes_to_string`.
12. Column 8 → raw bytes `Vec<u8>` (the `blob`).
13. Columns ≥9 → `extra: Vec<String>` (each lossy-decoded).
14. Return the 10-element raw tuple `(seqid, source, featuretype, start, end, score, strand, frame, blob, extra)`.

**`Iterator::next`** (`parser.rs:289-360`) — full record assembly + validation, looping (so non-strict skips advance to the next line):
1. Call `next_raw_record`. On `Some(Err(e))`: if `strict` → return `Some(Err(e))`; else push `e` onto `warnings` and `continue`. On `None` → return `None`.
2. `is_gtf = matches!(dialect.fmt, Format::Gtf)`.
3. `validate_fields(line_no, …, is_gtf)` (see §5.6.C.7). On error: strict → return; else push warning + `continue`.
4. `parse_attributes(blob_str)` where `blob_str = from_utf8(&blob).unwrap_or("")` — this is the *second* call to `parse_attributes` on this line (the peek phase already called it once; this is the authoritative parse).
5. `validate_attributes_pairs(line_no, pairs.len(), &blob, is_gtf)`. On error: strict → return; else warning + `continue`.
6. Build and return `Record` with `attributes_blob = blob` (the raw bytes, **not** re-serialised) and `attributes_pairs = pairs`.

**`peek_dialect`** (`parser.rs:129-161`) — the two-pass dialect inference:
1. Snapshot `pos`, `line_no`, `directives` (cloned), `fasta_reached`.
2. `limit = usize::MAX` if `force_dialect_check` else `checklines` (default 10).
3. While `samples.len() < limit`: call `next_raw_record`; on `Some(Ok(raw))` run `parse_attributes(blob_str)` and push the per-line `obs` dialect onto `samples`; on `Some(Err(_))` → **break**; on `None` → break.
4. `dialect = dialect::choose(&samples)`.
5. **Restore** the snapshot of `pos`/`line_no`/`directives`/`fasta_reached` — the peek does not consume records; iteration starts over from the top.
6. If `force_gff`: override `dialect.fmt = Format::Gff3` and `dialect.keyval_separator = '='`.

**Gotcha (peek/iterate double-parse):** every feature line in the first `checklines` is parsed *twice* — once during `peek_dialect` (to sample dialect) and once during `next` (the authoritative parse). The two parses must be deterministic (they are — `parse_attributes` is pure). Validation runs **only** in `next`, never in `peek`, so a malformed line within the peek window does not raise during peek (it just breaks the sampling loop); it raises (or warns) when iteration reaches it.

**Gotcha (`##FASTA` terminates everything):** once `fasta_reached` is set, `next_raw_record` returns `None` forever. Any feature rows physically after a `##FASTA` directive are unreachable. The fallback replicates this exactly.

#### 5.6.C.5 Native attribute parser (`attributes.rs`)

`parse_attributes(blob: &str) -> (Vec<(String,String,u16)>, Dialect)` (`attributes.rs:18-107`). The fallback's `parse_attributes` (`_pyfallback/attributes.py:17-73`) must match byte-for-byte. Algorithm:

1. If `blob` is empty → return `([], Dialect::default())`.
2. **Leading/trailing semicolon flags:** `obs.leading_semicolon = blob.trim_start().starts_with(';')`; `obs.trailing_semicolon = blob.trim_end().ends_with(';')`.
3. **Field-separator inference** (first match wins): if `blob.contains("; ")` → `"; "`; else if `blob.contains(" ; ")` → `" ; "`; else → `";"`. **Gotcha:** the Python fallback checks `"; "` first too (`attributes.py:24-29`); the order matters because `"; "` is a substring test that would also be true when `" ; "` is present, so `"; "` always wins when both could match. Keep this exact order.
4. **Top-level semicolon split** — `split_top_level_semicolons` (`attributes.rs:110-129`, fallback `attributes.py:76-90`): walk chars; a `"` toggles `in_quotes`; a `;` at top level (`!in_quotes`) ends a segment; a `;` while `in_quotes` sets `obs.semicolon_in_quotes = true` (it is *not* a split point). After the loop, push the trailing segment `blob[start..]`. **Gotcha:** the Rust guard `if start <= blob.len()` (`attributes.rs:125`) is always true, so the trailing segment is always pushed — matching the fallback's unconditional `out.append(blob[start:])` (`attributes.py:89`).
5. For each segment (`for seg in segments`):
   a. `seg = seg.trim()`; if empty → `continue`.
   b. `split_keyval(seg)` (`attributes.rs:133-144`, fallback `attributes.py:93-101`): if `=` is present, split at the **first** `=`; key = `seg[:eq].trim()`, value = `seg[eq+1:].trim()`, separator `'='`. Else if a whitespace char is present (`is_whitespace()` in Rust; `' '` or `'\t'` in the fallback at `attributes.py:99`), split at the first whitespace; key/value trimmed; separator `' '`. Else key = `seg.trim()`, value = `""`, separator `'='`.
   c. If `key` is empty → `continue`.
   d. `local_fmt = Gff3` if separator is `'='` else `Gtf`. The **first** non-empty segment's `local_fmt` sets `detected_fmt` (subsequent differing segments do not change it — `attributes.rs:61-67`; fallback only sets it once at `attributes.py:46-47`).
   e. `strip_quotes(raw_val)` (`attributes.rs:146-152`, fallback `attributes.py:104-108`): trim, then if `len ≥ 2 and starts_with('"') and ends_with('"')` → strip both quotes, `was_quoted = true`. If `was_quoted` → `obs.quoted_gff2_values = true`.
   f. **Multi-value split:** if `local_fmt == Gff3` → `split_unquoted_commas(clean_val)` (split on `,` not inside quotes; `attributes.rs:155-172`, fallback `attributes.py:111-122`); else (GTF) → a single-element list `[clean_val]`.
   g. If `key` not already in `order` → append it.
   h. `counter = keys_seen.get(key, 0)`; if `counter > 0` → `obs.repeated_keys = true`.
   i. For each value `v` in the multi-value list: `decoded = unescape(v)` if GFF3 else `v` verbatim (GTF values are **not** percent-decoded); push `(key, decoded, counter)`; `counter += 1`. After the loop, `keys_seen[key] = counter`.
6. `obs.fmt = detected_fmt.unwrap_or(Gff3)`. `obs.keyval_separator = ' '` if fmt is GTF else `'='`. `obs.order = order`.
7. Return `(out, obs)`.

**Percent-decoding** — Rust `unescape` (`escape.rs:10-32`); the **fallback uses `urllib.parse.unquote`** (`attributes.py:65`). These must agree on real GFF3 inputs. The Rust rules (`escape.rs:17-31`):
- Fast path: if no `%` byte, return the input borrowed unchanged.
- Otherwise scan: at `%` where `i+2 < len` and both following bytes are hex (`hex_val`, `escape.rs:35-42`: `0-9`,`a-f`,`A-F`), emit `(h<<4)|l` and skip 3 bytes; else emit the literal byte and advance 1.
- After building the byte vector, attempt `String::from_utf8`; **on UTF-8 failure, return the *original* input unchanged** (`escape.rs:28-31`).

**Gotcha (`%` decoding edge):** the Rust condition is `i + 2 < bytes.len()` (strict less-than), meaning a `%` in the **last two byte positions** is left literal (e.g. `"100%"` → `"100%"`, `"%ZZ"` → `"%ZZ"`, both per the `escape.rs` tests). `urllib.parse.unquote` has the same "pass through malformed escapes" behaviour, which is why the fallback can use it. A reimplementer choosing a different decoder MUST verify: malformed escapes pass through literally, and a `%` needing 2 trailing hex chars at the very end of the string is left as a literal `%`.

**Gotcha (GTF values not decoded):** percent-unescaping is applied **only** when `local_fmt == Gff3`. GTF values are inserted verbatim. Both implementations agree (`attributes.rs:93-97`, `attributes.py:64-66`).

#### 5.6.C.6 Pure-Python fallback parser (`_pyfallback/parser.py`)

This module is the byte-compatible oracle and the live path when the `.so` is absent. Public surface re-exported via `_pyfallback/__init__.py:8`: `parse_file`, `parse_bytes`, `detect_dialect`.

**Line normalisation** (`_iter_lines`, `parser.py:76-82`): strips a single trailing `\n` then a single trailing `\r` per line. **Gotcha:** files are opened with `newline=""` (`_open`, `parser.py:70-73`) so Python does not do universal-newline translation; the strip logic handles CRLF explicitly to match the Rust `trim_cr`. `.gz` paths open via `gzip.open(..., "rt", encoding="utf-8", newline="")`.

**`_parse_line_into_feature(line, line_no)`** (`parser.py:184-219`) mirrors `next_raw_record` + the validation in `next`:
1. `fields = line.split("\t")`; if `< 9` → raise `TooFewFields` error.
2. Unpack `fields[:8]` and `blob = fields[8]`, `extra = fields[9:]`.
3. `pairs, _obs = parse_attributes(blob)`.
4. `start = _coord_or_error(start_s, …, "start")`, same for `end` (`parser.py:171-181`): `.`/empty → `None`; non-int → raise `InvalidCoordinate`.
5. `_validate(...)` (`parser.py:91-168`) — the unified field+attribute check (mirror of `validate_fields` + `validate_attributes_pairs`, see §5.6.C.7). If it returns a non-None error → raise it.
6. Return a `ParsedFeature` with `attributes_blob = blob.encode("utf-8")`, `attributes_pairs = pairs`, `extra = extra`.

**Streaming generator `_stream_features`** (`parser.py:222-304`) — replicates the peek/iterate split:
1. **Phase A (peek + buffer):** iterate lines counting `line_no` from 1. Skip blanks. On `##`: if `##FASTA` set `fasta_reached=True` and **break**; else append to `directives` and continue. On single `#`: continue. Otherwise parse the line into a feature via `_parse_line_into_feature`, catching the format-error class; on error call `_maybe_handle` (strict → re-raise; non-strict → append `{line_no,kind,message}` to `warnings` and `continue`). On success: re-call `parse_attributes` on the decoded blob to get `obs`, append `obs` to `samples` and the feature to `buffered`. If `not force_dialect_check and len(buffered) >= checklines` → **break** (stop peeking).
2. `dialect = merge_dialects(samples)` (the dialect-choose). If `force_gff`: set `dialect["fmt"]="gff3"` and `dialect["keyval separator"]="="`.
3. **Yield the buffered features** (each as `(feat, directives, dialect)`).
4. If `fasta_reached` → return (no more features).
5. **Phase B (stream the rest):** continue reading lines from where Phase A left off (the underlying stream is a single iterator, so position is preserved), applying the same directive/comment/FASTA/parse/validate rules, yielding each feature.

**Gotcha (peek double-parse, fallback edition):** in Phase A the fallback parses each buffered line via `_parse_line_into_feature` (which calls `parse_attributes` once and validates) **and then calls `parse_attributes` again** on the decoded blob to collect `obs` (`parser.py:270`). This is the same double-work as Rust but in a different order — the fallback validates *during* the first parse, whereas Rust validates only in `next`. Net effect is identical for valid inputs. For malformed inputs in the peek window: the fallback's first parse (inside `_parse_line_into_feature`) will raise/warn (because it validates), whereas Rust's peek just `break`s and re-encounters the line in `next`. Both ultimately raise/warn on the same line — only the dialect sample set may differ by the offending line in non-strict mode. **A reimplementer should preserve "validation happens once per feature, errors surface on the offending line, dialect samples exclude lines that error" as the contract, not the exact internal ordering.**

**`_FallbackIterator`** (`parser.py:307-370`) wraps the generator with the Rust-iterator surface:
- `__next__` pulls `(feat, directives, dialect)` from the generator, caches `self._dialect`, returns `feat`.
- `dialect()` — if `_dialect is None`, call `_drain_for_metadata()` first; return `self._dialect or {}`.
- `directives()` — if `_directives` empty, `_drain_for_metadata()`; return a *copy* of `_directives`.
- `warnings` property — copy of `_warnings`.
- `_drain_for_metadata()` (`parser.py:335-349`) — forces the generator to its first yield (running the whole peek phase, populating `directives`+`dialect`), then **stashes** that first record by re-chaining `[first]` ahead of the remaining generator (`_chain`, `parser.py:351-356`) so `__next__` still sees it. Idempotent via `_exhausted`. This is what makes `dialect()`/`directives()` work on directives-only or empty inputs *without* consuming the first feature. **Gotcha:** `directives` is passed into the generator and mutated *in place* (`parser.py:233-234`), so callers can read directives even when the file has zero feature rows.

**`detect_dialect(path, checklines=10)`** (`parser.py:395-404`): non-strict by design. Builds the iterator with `strict=False`, drains up to `checklines` records (swallowing `StopIteration`), then returns `it.dialect()`. Mirrors the Rust `detect_dialect`'s non-strict peek.

Default parameter values are identical across both engines and both wrappers: `checklines=10`, `force_dialect_check=False`, `force_gff=False`, `strict=True`.

#### 5.6.C.7 Validation rules (must match across engines)

The Rust validator lives in `validate.rs`; the fallback packs the same rules into `_validate` (`_pyfallback/parser.py:91-168`). Both produce a `GFFFormatError`-family error carrying `line_no`, `kind` (one of the `ErrorKind` strings), and `message`. **All checks fire in this exact order** (first failure wins); the message strings must match.

`ErrorKind` values (`validate.rs:19-46`): `TooFewFields`, `EmptySeqid`, `EmptyFeaturetype`, `InvalidFeaturetype`, `InvalidCoordinate`, `InvalidStrand`, `InvalidPhase`, `InvalidScore`, `InvalidAttribute`.

Field validation order (`validate_fields`, `validate.rs:82-177`; fallback `_validate`, `parser.py:91-168`):

| Order | Condition (error if true) | kind | message |
|---|---|---|---|
| 1 | `seqid` is empty | `EmptySeqid` | `seqid (col 1) is empty` |
| 2 | `featuretype` is empty | `EmptyFeaturetype` | `featuretype (col 3) is empty` |
| 3 | `featuretype` contains any whitespace char | `InvalidFeaturetype` | `featuretype contains whitespace: <repr>` |
| 4 | `start is not None and start < 1` | `InvalidCoordinate` | `start coordinate must be >= 1 (got N)` |
| 5 | `start` and `end` both present and `end < start` | `InvalidCoordinate` | `end < start (E < S)` |
| 6 | `strand not in {"+","-","?","."}` | `InvalidStrand` | `strand must be one of '+', '-', '?', '.'; got <repr>` |
| 7 | `frame not in {".","0","1","2"}` | `InvalidPhase` | `phase must be 0, 1, 2, or '.'; got <repr>` |
| 8 | `featuretype == "CDS" and frame == "."` | `InvalidPhase` | `CDS row missing required phase (must be 0, 1, or 2)` |
| 9 | `score not in {"","."}` and `float(score)` fails | `InvalidScore` | `score must be a float or '.'; got <repr>` |

**Gotcha (TooFewFields and coordinate-not-integer fire earlier):** the `<9 fields` check and the non-integer coordinate checks happen in `next_raw_record` / `_parse_line_into_feature` *before* `validate_fields`/`_validate` runs, so their `line_no` is the line being read. In the fallback, `_coord_or_error` raises `InvalidCoordinate` with message `start/end coordinate is not an integer: <repr>` before `_validate` is reached.

**Attribute-structure validation** (`validate_attributes_pairs`, `validate.rs:188-220`; folded into `_validate` at `parser.py:152-167`):
1. `trimmed = blob.strip()`. If `trimmed` is empty or `== "."` → OK (no error). Real-world files emit bare `.` column-9.
2. Else compute `has_eq = "=" in trimmed`, `has_quote = '"' in trimmed`, `has_pair = n_pairs > 0`.
3. `malformed = (not has_pair) or (not has_eq and not has_quote)`. If `malformed` → `InvalidAttribute` with message `attribute string did not parse into any key=value pair: <trimmed[:60] repr>`.

**Gotcha (the `is_gtf` flag is informational only):** `validate_attributes_pairs` ignores `is_gtf` (`validate.rs:207` `let _ = is_gtf;`). A `force_gff=True` override that mislabels GTF-attribute data as GFF3 must **not** trigger a spurious failure — accepting EITHER a `=` (GFF3) OR a `"` (GTF) covers both. The fallback condition at `parser.py:162` is `(not has_eq and not has_quote) or n_pairs == 0`, logically identical to the Rust `!has_pair || (!has_eq && !has_quote)`. **Gotcha:** raw garbage like `"justsomegarbage"` parses into a stub `(token, "")` pair (so `n_pairs==1`) but contains neither `=` nor `"` — hence it is correctly rejected as `InvalidAttribute`. This is the precise reason both `has_pair` *and* the `=`/`"` structural check are required.

**Strict vs non-strict semantics** (identical both engines):
- `strict=True` (default): the iterator yields/raises the first `GFFFormatError` and stops at that line; `warnings` stays empty.
- `strict=False`: malformed lines are *silently dropped*, each recorded as a dict `{line_no, kind, message}` appended to the iterator's `warnings`. Iteration continues to the next line.

#### 5.6.C.8 Dialect inference & the dialect dict

Per-line observations are reconciled by `dialect::choose` (Rust, `dialect.rs:55-97`) / `merge_dialects` (Python, `dialect.py:28-62`). Rules (must match):
1. If `samples` empty → `default_dialect()` (Rust `Dialect::default()`, `dialect.rs:27-42`).
2. **Format vote:** `n_gtf = count(fmt==gtf)`; `fmt = "gtf"` iff `n_gtf > n_gff3` (i.e. strict majority GTF; **ties go to GFF3**). `keyval separator = " "` for GTF else `"="`.
3. **Boolean flags** (`leading semicolon`, `trailing semicolon`, `quoted GFF2 values`, `repeated keys`, `semicolon in quotes`): logical **OR** across all samples (any-true → true).
4. **Field separator:** plurality vote (most common value among samples; Rust uses `max_by_key` on counts, Python uses `max(set(...), key=list.count)`). **Gotcha:** tie-break behaviour for the plurality vote is unspecified/hash-order-dependent in Rust's `HashMap` and `set`-based in Python; for byte-identity this is only safe when one separator strictly dominates the peek window, which holds for all real homogeneous annotation files. A reimplementer must not rely on a deterministic tie-break here.
5. **Order:** union of per-sample `order` lists by **first appearance** across samples (a `seen` set guards duplicates).

The `default_dialect()` template (Python `dialect.py:13-25`, Rust `Dialect::default()`):

| Key | Default |
|---|---|
| `fmt` | `"gff3"` |
| `field separator` | `";"` |
| `keyval separator` | `"="` |
| `multival separator` | `","` |
| `leading semicolon` | `False` |
| `trailing semicolon` | `False` |
| `quoted GFF2 values` | `False` |
| `repeated keys` | `False` |
| `semicolon in quotes` | `False` |
| `order` | `[]` |

This shape **deliberately mirrors `gffutils.constants.dialect`** (comments at `dialect.rs:5`, `dialect.py:4`) so downstream legacy consumers and the `Feature._format_attributes` round-trip (`feature.py:357-395`) can consume it untranslated.

#### 5.6.C.9 Equivalence checklist for a reimplementer

To keep the two paths byte-identical (the property the test suite pins), a reimplementer of either engine must guarantee:

1. **Identical 11-tuple shape and order** (§5.6.C.3) — `start`/`end` are `int | None` with `.`/empty → `None`; `attributes_blob` is the *raw verbatim* column-9 bytes (never re-serialised); `attributes_pairs` are `(key, value, idx)` with `idx` being the 0-based multi-value index per key.
2. **Same attribute split + decode rules** (§5.6.C.5) — first-`=` keyval split, first-whitespace GTF split, quote stripping only when both ends are `"`, comma split only outside quotes and only for GFF3, percent-decode only for GFF3 values, malformed-escape pass-through, GTF values verbatim.
3. **Same validation order, kinds, and message strings** (§5.6.C.7), and the same strict/non-strict demotion to `warnings`.
4. **Same dialect inference** (§5.6.C.8) — majority-GTF-else-GFF3 format vote, OR of boolean flags, plurality field-separator, first-appearance order, with the `force_gff` override forcing `fmt="gff3"` + `keyval="="`.
5. **`##FASTA` stops all feature emission; single-`#` comments are dropped and not recorded; `##` directives are collected; blank lines skipped.**
6. **Peek window of `checklines` (default 10) lines; `force_dialect_check` peeks the whole file; the peek does not consume records (records are re-read in iteration).**
7. **`directives()`/`dialect()` remain readable after exhaustion and on directives-only/empty inputs.**
8. **The raised exception class is `ValueError`-derived with `.line_no` / `.kind` / `.message` attributes, and the live class is whichever engine is loaded (rebound in `gffbase/__init__.py`).**
