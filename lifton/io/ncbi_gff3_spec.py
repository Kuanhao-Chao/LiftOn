"""NCBI GFF3 spec invariants used at runtime by the strict validator.

Constants only — no logic. Source:
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/
"""

from __future__ import annotations

# Reserved characters that MUST be percent-encoded inside an attribute
# value (per NCBI § Attribute Specifications, which inherits the
# Sequence Ontology GFF3 specification).
RESERVED_CHARS = frozenset({";", "=", "&", ",", "\t", "\n", "\r"})

# Official GFF3 attribute names — capitalised initial letter (NCBI §
# Official GFF3 attributes).
OFFICIAL_ATTRS = frozenset({
    "ID", "Parent", "Name", "Alias", "Target", "Gap",
    "Derives_from", "Note", "Dbxref", "Ontology_term", "Is_circular",
})

# Attributes whose value is a comma-separated list (NCBI multi-value §).
MULTI_VALUE_ATTRS = frozenset({
    "Parent", "Alias", "Note", "Dbxref", "Ontology_term",
})

# Allowed values for col 7 (strand) per the GFF3 specification.
VALID_STRANDS = frozenset({"+", "-", ".", "?"})

# Allowed values for col 8 (phase). NCBI permits "." for non-CDS rows;
# CDS rows should carry 0/1/2 but the spec explicitly tolerates "." for
# pseudogenes and other features with internal frameshifts.
VALID_PHASES = frozenset({"0", "1", "2", "."})

# Directive prefixes recognised by NCBI.
DIRECTIVE_PREFIX = "##"
NCBI_DIRECTIVE_PREFIX = "#!"

# Mandatory first-line directive.
GFF_VERSION_DIRECTIVE = "##gff-version 3"
