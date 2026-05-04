// ---------------------------------------------------------------------------
// Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
// ---------------------------------------------------------------------------
//! NCBI GFF3 compliance validation.
//!
//! Enforces the rules at
//! <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/>:
//! 9 tab-separated columns, 1-based integer coordinates with `start <= end`,
//! strand ∈ {`+`,`-`,`?`,`.`}, phase ∈ {`0`,`1`,`2`,`.`} (mandatory `0/1/2`
//! for `CDS` rows), score ∈ {float, `.`}, non-empty `seqid` and
//! whitespace-free `featuretype`, attributes parseable as `key=value` pairs.
//!
//! Errors are STRUCTURED — every `GffError` carries the offending line
//! number, an `ErrorKind` enum, and a human-readable message ready to be
//! shown to a user grepping a 3 GB annotation file.

use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorKind {
    TooFewFields,
    EmptySeqid,
    EmptyFeaturetype,
    InvalidFeaturetype,
    InvalidCoordinate,
    InvalidStrand,
    InvalidPhase,
    InvalidScore,
    InvalidAttribute,
}

impl ErrorKind {
    pub fn as_str(&self) -> &'static str {
        match self {
            ErrorKind::TooFewFields       => "TooFewFields",
            ErrorKind::EmptySeqid         => "EmptySeqid",
            ErrorKind::EmptyFeaturetype   => "EmptyFeaturetype",
            ErrorKind::InvalidFeaturetype => "InvalidFeaturetype",
            ErrorKind::InvalidCoordinate  => "InvalidCoordinate",
            ErrorKind::InvalidStrand      => "InvalidStrand",
            ErrorKind::InvalidPhase       => "InvalidPhase",
            ErrorKind::InvalidScore       => "InvalidScore",
            ErrorKind::InvalidAttribute   => "InvalidAttribute",
        }
    }
}

#[derive(Debug, Clone)]
pub struct GffError {
    pub line_no: usize,
    pub kind: ErrorKind,
    pub message: String,
}

impl GffError {
    pub fn new(line_no: usize, kind: ErrorKind, message: impl Into<String>) -> Self {
        Self {
            line_no,
            kind,
            message: message.into(),
        }
    }
}

impl fmt::Display for GffError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "GFF3 format error on line {} [{}]: {}",
            self.line_no,
            self.kind.as_str(),
            self.message,
        )
    }
}

/// Validate the parsed 9-column GFF3 record. `start` and `end` are
/// `None` when the source file used `.` (some real-world files do this
/// for chromosome-level rows; we accept it but check ordering when both
/// are present).
#[allow(clippy::too_many_arguments)]
pub fn validate_fields(
    line_no: usize,
    seqid: &str,
    featuretype: &str,
    start: Option<i64>,
    end: Option<i64>,
    score: &str,
    strand: &str,
    frame: &str,
    attrs_blob: &[u8],
    is_gtf: bool,
) -> Result<(), GffError> {
    if seqid.is_empty() {
        return Err(GffError::new(
            line_no,
            ErrorKind::EmptySeqid,
            "seqid (col 1) is empty",
        ));
    }

    if featuretype.is_empty() {
        return Err(GffError::new(
            line_no,
            ErrorKind::EmptyFeaturetype,
            "featuretype (col 3) is empty",
        ));
    }
    if featuretype.chars().any(|c| c.is_whitespace()) {
        return Err(GffError::new(
            line_no,
            ErrorKind::InvalidFeaturetype,
            format!("featuretype contains whitespace: {:?}", featuretype),
        ));
    }

    if let Some(s) = start {
        if s < 1 {
            return Err(GffError::new(
                line_no,
                ErrorKind::InvalidCoordinate,
                format!("start coordinate must be >= 1 (got {})", s),
            ));
        }
    }
    if let (Some(s), Some(e)) = (start, end) {
        if e < s {
            return Err(GffError::new(
                line_no,
                ErrorKind::InvalidCoordinate,
                format!("end < start ({} < {})", e, s),
            ));
        }
    }

    if !matches!(strand, "+" | "-" | "?" | ".") {
        return Err(GffError::new(
            line_no,
            ErrorKind::InvalidStrand,
            format!(
                "strand must be one of '+', '-', '?', '.'; got {:?}",
                strand
            ),
        ));
    }

    if !matches!(frame, "." | "0" | "1" | "2") {
        return Err(GffError::new(
            line_no,
            ErrorKind::InvalidPhase,
            format!("phase must be 0, 1, 2, or '.'; got {:?}", frame),
        ));
    }
    if featuretype == "CDS" && frame == "." {
        return Err(GffError::new(
            line_no,
            ErrorKind::InvalidPhase,
            "CDS row missing required phase (must be 0, 1, or 2)",
        ));
    }

    if score != "." && !score.is_empty() && score.parse::<f64>().is_err() {
        return Err(GffError::new(
            line_no,
            ErrorKind::InvalidScore,
            format!("score must be a float or '.'; got {:?}", score),
        ));
    }

    // Attribute-string structure is validated AFTER parsing in
    // `parser.rs::Iterator::next` — see `validate_attributes_pairs`
    // below. Mixing the check in here would force this function to
    // duplicate the parser's GTF/GFF3 dispatch.
    let _ = (attrs_blob, is_gtf);

    Ok(())
}

/// Post-parse attribute check. If the parser yielded zero `(key, value)`
/// pairs from a non-empty col-9 string, the blob is structurally
/// malformed regardless of dialect. Returns `Ok(())` for the empty / `.`
/// blob (some real-world files emit these).
///
/// For GFF3 inputs we additionally require the raw blob to contain at
/// least one `=` somewhere — the attribute parser is permissive enough
/// to yield a `(key, "")` pair from raw garbage like `"justsomegarbage"`,
/// which the GFF3 spec forbids.
pub fn validate_attributes_pairs(
    line_no: usize,
    n_pairs: usize,
    attrs_blob: &[u8],
    is_gtf: bool,
) -> Result<(), GffError> {
    let trimmed = std::str::from_utf8(attrs_blob).unwrap_or("").trim();
    if trimmed.is_empty() || trimmed == "." {
        return Ok(());
    }
    let has_eq    = trimmed.contains('=');
    let has_quote = trimmed.contains('"');
    let has_pair  = n_pairs > 0;

    // Accept the blob if EITHER `=` (GFF3) OR `"` (GTF) appears AND the
    // parser produced at least one pair. The dialect flag is informational
    // — `force_gff=True` may falsely label a GTF-attribute line as GFF3,
    // and we don't want that user override to trigger a spurious
    // validation failure.
    let _ = is_gtf;
    let malformed = !has_pair || (!has_eq && !has_quote);
    if malformed {
        return Err(GffError::new(
            line_no,
            ErrorKind::InvalidAttribute,
            format!(
                "attribute string did not parse into any key=value pair: {:?}",
                if trimmed.len() > 60 { &trimmed[..60] } else { trimmed }
            ),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ok(seqid: &str, ft: &str, s: Option<i64>, e: Option<i64>,
          score: &str, strand: &str, frame: &str, attrs: &[u8]) {
        validate_fields(1, seqid, ft, s, e, score, strand, frame, attrs, false).unwrap();
    }
    fn err_kind(seqid: &str, ft: &str, s: Option<i64>, e: Option<i64>,
                score: &str, strand: &str, frame: &str, attrs: &[u8]) -> ErrorKind {
        validate_fields(1, seqid, ft, s, e, score, strand, frame, attrs, false).unwrap_err().kind
    }

    #[test]
    fn happy_path() {
        ok("chr1", "exon", Some(1), Some(100), ".", "+", ".", b"ID=x");
    }

    #[test]
    fn empty_seqid() {
        assert_eq!(
            err_kind("", "exon", Some(1), Some(10), ".", "+", ".", b"ID=x"),
            ErrorKind::EmptySeqid,
        );
    }

    #[test]
    fn whitespace_in_featuretype() {
        assert_eq!(
            err_kind("chr1", "exon foo", Some(1), Some(10), ".", "+", ".", b"ID=x"),
            ErrorKind::InvalidFeaturetype,
        );
    }

    #[test]
    fn coord_zero_or_negative() {
        assert_eq!(
            err_kind("chr1", "exon", Some(0), Some(10), ".", "+", ".", b"ID=x"),
            ErrorKind::InvalidCoordinate,
        );
        assert_eq!(
            err_kind("chr1", "exon", Some(-5), Some(10), ".", "+", ".", b"ID=x"),
            ErrorKind::InvalidCoordinate,
        );
    }

    #[test]
    fn end_less_than_start() {
        assert_eq!(
            err_kind("chr1", "exon", Some(100), Some(50), ".", "+", ".", b"ID=x"),
            ErrorKind::InvalidCoordinate,
        );
    }

    #[test]
    fn invalid_strand() {
        assert_eq!(
            err_kind("chr1", "exon", Some(1), Some(10), ".", "@", ".", b"ID=x"),
            ErrorKind::InvalidStrand,
        );
        assert_eq!(
            err_kind("chr1", "exon", Some(1), Some(10), ".", "+-", ".", b"ID=x"),
            ErrorKind::InvalidStrand,
        );
    }

    #[test]
    fn strand_question_and_dot_ok() {
        ok("chr1", "exon", Some(1), Some(10), ".", "?", ".", b"ID=x");
        ok("chr1", "exon", Some(1), Some(10), ".", ".", ".", b"ID=x");
    }

    #[test]
    fn invalid_phase() {
        assert_eq!(
            err_kind("chr1", "exon", Some(1), Some(10), ".", "+", "5", b"ID=x"),
            ErrorKind::InvalidPhase,
        );
    }

    #[test]
    fn cds_requires_concrete_phase() {
        assert_eq!(
            err_kind("chr1", "CDS", Some(1), Some(10), ".", "+", ".", b"ID=x"),
            ErrorKind::InvalidPhase,
        );
        ok("chr1", "CDS", Some(1), Some(10), ".", "+", "0", b"ID=x");
    }

    #[test]
    fn invalid_score() {
        assert_eq!(
            err_kind("chr1", "exon", Some(1), Some(10), "abc", "+", ".", b"ID=x"),
            ErrorKind::InvalidScore,
        );
        ok("chr1", "exon", Some(1), Some(10), "0.95", "+", ".", b"ID=x");
        ok("chr1", "exon", Some(1), Some(10), ".", "+", ".", b"ID=x");
    }

    #[test]
    fn validate_attributes_pairs_rules() {
        // Empty/dot accepted.
        validate_attributes_pairs(1, 0, b"", false).unwrap();
        validate_attributes_pairs(1, 0, b".", false).unwrap();
        validate_attributes_pairs(1, 0, b"   ", false).unwrap();
        // GFF3 garbage → error even if the parser produced one stub pair.
        let e = validate_attributes_pairs(7, 1, b"justsomegarbage", false).unwrap_err();
        assert_eq!(e.kind, ErrorKind::InvalidAttribute);
        assert_eq!(e.line_no, 7);
        // GFF3 well-formed → OK.
        validate_attributes_pairs(1, 1, b"ID=x", false).unwrap();
        // GTF well-formed → OK.
        validate_attributes_pairs(1, 1, br#"gene_id "ENSG";"#, true).unwrap();
        // GTF unquoted garbage (no `=`, no `"`) → error.
        let e = validate_attributes_pairs(1, 1, b"unquoted", true).unwrap_err();
        assert_eq!(e.kind, ErrorKind::InvalidAttribute);
        // GFF3 user forced over GTF data: blob has quotes, no `=`. Should
        // still be accepted because the structure is recognizable.
        validate_attributes_pairs(1, 1, br#"gene_id "ENSG";"#, false).unwrap();
    }
}
