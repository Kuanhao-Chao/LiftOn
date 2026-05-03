// ---------------------------------------------------------------------------
// Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
// ---------------------------------------------------------------------------
//! Column-9 attribute parser. Hand-written state machine that handles:
//!  - GFF3:   `key=val;key=val,val2`
//!  - GTF:    `key "val"; key "val";`
//!  - Mixed:  quoted values containing `;` or `,`
//!  - Spaces around `;` (some vendors emit `; ` or ` ; `)
//!  - Trailing semicolons, leading semicolons.
//!  - Repeated keys (we add another (key, value, idx) row).

use crate::dialect::{Dialect, Format};
use crate::escape::unescape;

/// Parse a single feature's column-9 string and return:
///   - the list of (key, value, multivalue_index) triples
///   - per-line dialect observations to feed `dialect::choose`
pub fn parse_attributes(blob: &str) -> (Vec<(String, String, u16)>, Dialect) {
    let mut out: Vec<(String, String, u16)> = Vec::new();
    let mut obs = Dialect::default();

    if blob.is_empty() {
        return (out, obs);
    }

    // Detect leading / trailing semicolons.
    obs.leading_semicolon = blob.trim_start().starts_with(';');
    obs.trailing_semicolon = blob.trim_end().ends_with(';');
    if blob.contains("; ") {
        obs.field_separator = "; ".into();
    } else if blob.contains(" ; ") {
        obs.field_separator = " ; ".into();
    } else {
        obs.field_separator = ";".into();
    }

    // Decide GFF3 vs GTF: GFF3 has `key=value`; GTF has `key "value"` or
    // `key value` (space separator). The first non-empty record decides per-line.
    let segments = split_top_level_semicolons(blob, &mut obs);

    let mut keys_seen: std::collections::HashMap<String, u16> = std::collections::HashMap::new();
    let mut order: Vec<String> = Vec::new();
    let mut detected_fmt: Option<Format> = None;

    for seg in segments {
        let seg = seg.trim();
        if seg.is_empty() {
            continue;
        }

        // Determine local key/val split.
        let (key, raw_val, kv_sep) = split_keyval(seg);
        if key.is_empty() {
            continue;
        }
        let local_fmt = match kv_sep {
            '=' => Format::Gff3,
            ' ' => Format::Gtf,
            _ => Format::Gff3,
        };
        match detected_fmt {
            None => detected_fmt = Some(local_fmt),
            Some(prev) if prev != local_fmt => {
                // Mixed-format line — keep first decision but flag.
            }
            _ => {}
        }

        // Quoted GTF value handling.
        let (clean_val, was_quoted) = strip_quotes(raw_val);
        if was_quoted {
            obs.quoted_gff2_values = true;
        }

        // Multi-value handling.
        // GFF3: split on commas (canonical).
        // GTF:  rarely multi-value; we still split on comma for compatibility.
        let multi_values: Vec<&str> = if local_fmt == Format::Gff3 {
            split_unquoted_commas(clean_val)
        } else {
            vec![clean_val]
        };

        if !order.contains(&key.to_string()) {
            order.push(key.to_string());
        }
        let counter = keys_seen.entry(key.to_string()).or_insert(0);
        if *counter > 0 {
            obs.repeated_keys = true;
        }

        for v in multi_values {
            let decoded = if local_fmt == Format::Gff3 {
                unescape(v).into_owned()
            } else {
                v.to_string()
            };
            out.push((key.to_string(), decoded, *counter));
            *counter += 1;
        }
    }

    obs.fmt = detected_fmt.unwrap_or(Format::Gff3);
    obs.keyval_separator = if obs.fmt == Format::Gtf { ' ' } else { '=' };
    obs.order = order;
    (out, obs)
}

/// Split a column-9 string on top-level semicolons, honoring quoted ranges.
fn split_top_level_semicolons<'a>(blob: &'a str, obs: &mut Dialect) -> Vec<&'a str> {
    let bytes = blob.as_bytes();
    let mut out: Vec<&str> = Vec::new();
    let mut start = 0;
    let mut in_quotes = false;
    for (i, &b) in bytes.iter().enumerate() {
        if b == b'"' {
            in_quotes = !in_quotes;
        } else if b == b';' && !in_quotes {
            out.push(&blob[start..i]);
            start = i + 1;
        } else if b == b';' && in_quotes {
            obs.semicolon_in_quotes = true;
        }
    }
    if start <= blob.len() {
        out.push(&blob[start..]);
    }
    out
}

/// Split `key<sep>value` where sep is `=` (GFF3) or whitespace (GTF).
/// Returns (key, value, separator_char).
fn split_keyval(seg: &str) -> (&str, &str, char) {
    if let Some(eq) = seg.find('=') {
        let (k, v) = seg.split_at(eq);
        return (k.trim(), v[1..].trim(), '=');
    }
    // GTF-style: split on first whitespace.
    if let Some(ws) = seg.find(|c: char| c.is_whitespace()) {
        let (k, v) = seg.split_at(ws);
        return (k.trim(), v.trim_start().trim(), ' ');
    }
    (seg.trim(), "", '=')
}

fn strip_quotes(s: &str) -> (&str, bool) {
    let s = s.trim();
    if s.len() >= 2 && s.starts_with('"') && s.ends_with('"') {
        return (&s[1..s.len() - 1], true);
    }
    (s, false)
}

/// Split on commas that are not inside double quotes.
fn split_unquoted_commas(s: &str) -> Vec<&str> {
    let bytes = s.as_bytes();
    let mut out: Vec<&str> = Vec::new();
    let mut start = 0;
    let mut in_quotes = false;
    for (i, &b) in bytes.iter().enumerate() {
        match b {
            b'"' => in_quotes = !in_quotes,
            b',' if !in_quotes => {
                out.push(&s[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    out.push(&s[start..]);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gff3_simple() {
        let (pairs, dialect) = parse_attributes("ID=g1;Name=foo;Parent=p1");
        assert_eq!(pairs.len(), 3);
        assert_eq!(pairs[0], ("ID".into(), "g1".into(), 0));
        assert_eq!(dialect.fmt, Format::Gff3);
    }

    #[test]
    fn gff3_multivalue() {
        let (pairs, _) = parse_attributes("Parent=a,b,c");
        assert_eq!(pairs.len(), 3);
        assert_eq!(pairs[0].2, 0);
        assert_eq!(pairs[1].2, 1);
        assert_eq!(pairs[2].2, 2);
    }

    #[test]
    fn gtf_quoted() {
        let (pairs, dialect) = parse_attributes(r#"gene_id "ENSG"; transcript_id "ENST";"#);
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].0, "gene_id");
        assert_eq!(pairs[0].1, "ENSG");
        assert_eq!(dialect.fmt, Format::Gtf);
        assert!(dialect.quoted_gff2_values);
        assert!(dialect.trailing_semicolon);
    }

    #[test]
    fn semicolon_in_quotes() {
        let (pairs, dialect) = parse_attributes(r#"note "a;b";ID=g1"#);
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].1, "a;b");
        assert!(dialect.semicolon_in_quotes);
    }

    #[test]
    fn percent_escapes() {
        let (pairs, _) = parse_attributes("Note=hello%20world%2C%20you");
        assert_eq!(pairs[0].1, "hello world, you");
    }
}
