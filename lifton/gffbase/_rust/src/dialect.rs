// ---------------------------------------------------------------------------
// Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
// ---------------------------------------------------------------------------
//! Dialect representation. Mirrors the `gffutils.constants.dialect` shape so
//! Python code can consume it without translation.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Gff3,
    Gtf,
}

#[derive(Debug, Clone)]
pub struct Dialect {
    pub fmt: Format,
    pub field_separator: String,        // ";", "; ", " ; "
    pub keyval_separator: char,         // '=' for GFF3, ' ' for GTF
    pub multival_separator: char,       // typically ','
    pub leading_semicolon: bool,
    pub trailing_semicolon: bool,
    pub quoted_gff2_values: bool,
    pub repeated_keys: bool,
    pub semicolon_in_quotes: bool,
    pub order: Vec<String>,
}

impl Default for Dialect {
    fn default() -> Self {
        Dialect {
            fmt: Format::Gff3,
            field_separator: ";".into(),
            keyval_separator: '=',
            multival_separator: ',',
            leading_semicolon: false,
            trailing_semicolon: false,
            quoted_gff2_values: false,
            repeated_keys: false,
            semicolon_in_quotes: false,
            order: Vec::new(),
        }
    }
}

impl Dialect {
    pub fn fmt_str(&self) -> &'static str {
        match self.fmt {
            Format::Gff3 => "gff3",
            Format::Gtf => "gtf",
        }
    }
}

/// Reconcile per-line observations into a single dialect. `samples` is one
/// Dialect per peeked line; we pick the dominant format and combine flags.
pub fn choose(samples: &[Dialect]) -> Dialect {
    if samples.is_empty() {
        return Dialect::default();
    }
    let n_gtf = samples.iter().filter(|d| d.fmt == Format::Gtf).count();
    let n_gff3 = samples.len() - n_gtf;
    let fmt = if n_gtf > n_gff3 { Format::Gtf } else { Format::Gff3 };

    // For each flag, take the OR — if *any* line had a trailing semicolon,
    // the file has trailing semicolons. Mirrors gffutils' weighted-vote outcome
    // closely enough for round-trip purposes.
    let mut chosen = Dialect::default();
    chosen.fmt = fmt;
    chosen.keyval_separator = if fmt == Format::Gtf { ' ' } else { '=' };
    for s in samples {
        chosen.leading_semicolon |= s.leading_semicolon;
        chosen.trailing_semicolon |= s.trailing_semicolon;
        chosen.quoted_gff2_values |= s.quoted_gff2_values;
        chosen.repeated_keys |= s.repeated_keys;
        chosen.semicolon_in_quotes |= s.semicolon_in_quotes;
    }
    // Field separator: pick the most common.
    let mut counts: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for s in samples {
        *counts.entry(s.field_separator.as_str()).or_insert(0) += 1;
    }
    chosen.field_separator = counts
        .into_iter()
        .max_by_key(|(_, c)| *c)
        .map(|(k, _)| k.to_string())
        .unwrap_or_else(|| ";".into());

    // Build attribute key order from first appearance across samples.
    let mut seen = std::collections::HashSet::new();
    for s in samples {
        for k in &s.order {
            if seen.insert(k.clone()) {
                chosen.order.push(k.clone());
            }
        }
    }
    chosen
}
