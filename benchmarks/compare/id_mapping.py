"""Per-tool resolution of a lifted mRNA feature back to its reference mRNA id.

Liftoff / LiftOn preserve the reference transcript id (possibly with a ``_N``
copy suffix). miniprot emits its own ids (``MP...``) whose ``Target`` attribute
carries the reference protein accession (NP_/XP_) — or, for the MANE benchmark
where we headered reference proteins with ``rna-NM_`` ids, the reference mRNA id
directly.
"""
from __future__ import annotations

import re

_COPY_SUFFIX = re.compile(r"_(\d+)$")


def strip_copy_suffix(feature_id: str, ref_ids: set[str]) -> tuple[str, int]:
    """Return (base_ref_id_or_self, copy_index).

    A trailing ``_<int>`` is treated as a LiftOn/Liftoff copy suffix ONLY if the
    stripped base exists in ``ref_ids`` (mirrors lifton_utils.get_ID_base's
    safety check). copy_index is 0 for the primary, >=1 for extra copies.
    """
    if feature_id in ref_ids:
        return feature_id, 0
    m = _COPY_SUFFIX.search(feature_id)
    if m:
        base = feature_id[: m.start()]
        if base in ref_ids:
            return base, int(m.group(1))
    return feature_id, 0


def resolve_liftoff_lifton(mrna_feature, ref_ids: set[str]) -> tuple[str | None, int]:
    """Resolve a Liftoff/LiftOn mRNA feature to a reference mRNA id.

    Tries the feature id, then ``transcript_id`` / ``Name`` attributes, stripping
    copy suffixes against ``ref_ids``.
    """
    candidates = [mrna_feature.id]
    for key in ("transcript_id", "Name", "ID"):
        try:
            v = mrna_feature.attributes.get(key)
        except AttributeError:
            v = None
        if v:
            candidates.append(v[0] if isinstance(v, (list, tuple)) else v)
    for cand in candidates:
        if not cand:
            continue
        base, copy_idx = strip_copy_suffix(str(cand), ref_ids)
        if base in ref_ids:
            return base, copy_idx
    return None, 0


def build_miniprot_target_map(m_db) -> dict[str, str]:
    """miniprot mRNA id -> Target first token (protein acc or transcript id).

    Reuses the same parse as lifton_utils.miniprot_id_mapping.
    """
    out: dict[str, str] = {}
    for feature in m_db.features_of_type("mRNA"):
        try:
            tgt = str(feature.attributes["Target"][0]).split(" ")[0]
        except (KeyError, IndexError):
            continue
        out[feature.id] = tgt
    return out


def resolve_miniprot(mrna_feature, target_map: dict[str, str],
                     space: str, acc_to_mrna: dict[str, str],
                     ref_ids: set[str]) -> str | None:
    """Resolve a miniprot mRNA to a reference mRNA id.

    space == 'transcript': Target already equals the reference mRNA id.
    space == 'protein':    Target is a protein accession -> map via acc_to_mrna.
    """
    tgt = target_map.get(mrna_feature.id)
    if tgt is None:
        return None
    if space == "transcript":
        return tgt if tgt in ref_ids else None
    mrna = acc_to_mrna.get(tgt)
    if mrna and mrna in ref_ids:
        return mrna
    # Fallback: some refs key proteins to the mRNA directly.
    return tgt if tgt in ref_ids else None
