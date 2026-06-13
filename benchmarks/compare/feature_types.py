"""Generate the Liftoff/LiftOn ``-f`` feature-types file for the "allfeat" mode.

The "genes" benchmark mode runs Liftoff/LiftOn with the default feature set
(``gene`` only). The "allfeat" mode lifts EVERY top-level annotation type present
in the reference, so pseudogenes, mobile_genetic_elements, ncRNA genes, regulatory
features, etc. are lifted too.

Top-level types are derived by scanning the subset reference GFF for records that
have NO ``Parent=`` attribute (i.e. roots) — NOT from the manifest's flat
``feature_counts`` census, which also counts child types (mRNA/miRNA/tRNA) that
must not be passed to Liftoff's ``-f`` (they would be lifted detached from their
genes, breaking the hierarchy).

Scope ("all annotations, no alignment evidence"): exclude the chromosome
``region`` directive and alignment-evidence pseudo-features (``match`` /
``cDNA_match`` / ``match_part``); keep regulatory features (enhancer, silencer,
promoter, …) and every biological annotation type. ``gene`` is always included —
LiftOn's ``get_parent_features_to_lift`` (lifton_utils.py:230-246) RESETS to the
file's contents when ``-f`` is given, so an omitted ``gene`` would drop genes.
"""
from __future__ import annotations

from pathlib import Path

from .evaluator import _SUBFEATURE_TYPES

# alignment-evidence / sequence-region pseudo-features — not annotations to lift
_NON_ANNOTATION = frozenset({"region", "match", "cDNA_match", "match_part"})


def collect_toplevel_types(ref_gff: str) -> list[str]:
    """Return the sorted top-level annotation feature types in ``ref_gff``.

    Top-level := at least one record of that type has no ``Parent=`` attribute.
    Excludes sub-features, the chromosome ``region``, and alignment evidence.
    Always includes ``gene``.
    """
    toplevel: set[str] = set()
    with open(ref_gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if "Parent=" in parts[8]:          # has a parent -> not top-level
                continue
            toplevel.add(parts[2])
    keep = {t for t in toplevel
            if t not in _NON_ANNOTATION and t not in _SUBFEATURE_TYPES}
    keep.add("gene")
    return sorted(keep)


def write_types_file(manifest: dict, work_dir: Path, log=print) -> Path:
    """Write ``work/<id>/subset/all_feature_types.txt`` (one type per line) and
    return its path. Idempotent."""
    out = Path(work_dir) / "subset" / "all_feature_types.txt"
    types = collect_toplevel_types(manifest["paths"]["ref_gff"])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(types) + "\n")
    log(f"  [allfeat] {len(types)} top-level feature type(s): {types}")
    return out
