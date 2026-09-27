"""Repair duplicate-ID Liftoff parent links without changing its source file.

Liftoff can emit separate fragments of a trans-spliced gene with the same GFF3
ID. gffutils' ``create_unique`` keeps every root but attaches a child naming
that ID to the FIRST root, wherever the child lies -- on another fragment's
sequence (rice mitochondrial ``nad5``), or on the same sequence inside another
fragment (drosophila ``mod(mdg4)``: five fragments, all 31 transcripts bound to
the fragment containing none of them). The first root's span is then
normalised to its children, so that fragment is written at another's
coordinates and its own locus has no gene row.

A child its bound root does not contain is reassigned to the unique family
member on its sequence, of a compatible strand, that does. All other
relationships retain their original database behavior. When more than one
member qualifies the branch is left as the database bound it and counted: the
input cannot say which fragment it belongs to, and aborting the run would lose
every other gene over one that earlier releases lifted.
"""

from collections import defaultdict

from gffutils.exceptions import FeatureNotFoundError

from lifton import logger


def _attribute_id(feature):
    values = getattr(feature, "attributes", {}).get("ID", [])
    if isinstance(values, str):
        return values
    return str(values[0]) if values else str(feature.id)


def _same_strand(left, right):
    return (left == right or left in (".", "?") or right in (".", "?"))


def _contains(root, child):
    return (root.seqid == child.seqid
            and root.start <= child.start <= child.end <= root.end)


class SameSeqidParentOverlay:
    """Read-only ``children`` view correcting only proven family misbindings."""

    def __init__(self, database, root_types):
        self.database = database
        self.root_types = tuple(root_types)
        self._removed = defaultdict(set)
        self._added = defaultdict(list)
        self.ambiguous = []
        if "gene" not in root_types:
            return

        # ``create_unique`` changes the database key of the second root but
        # leaves its declared GFF3 ID in attributes. Keep only such collisions
        # in memory: retaining every root would scale with genome size.
        families = defaultdict(dict)
        renamed_roots = [root for root in database.features_of_type("gene")
                         if _attribute_id(root) != root.id]
        for root in renamed_roots:
            declared_id = _attribute_id(root)
            try:
                canonical = database[declared_id]
            except (KeyError, ValueError, FeatureNotFoundError):
                continue
            if (canonical.featuretype != "gene"
                    or _attribute_id(canonical) != declared_id):
                continue
            families[declared_id][canonical.id] = canonical
            families[declared_id][root.id] = root

        for members in families.values():
            family = list(members.values())
            for source in family:
                for child in database.children(source, level=1):
                    if _contains(source, child):
                        continue
                    candidates = [
                        root for root in family
                        if (_contains(root, child)
                            and _same_strand(root.strand, child.strand))
                    ]
                    if not candidates:
                        # A genuine trans-spliced branch can have no containing
                        # gene fragment; leave its relationship untouched.
                        continue
                    if len(candidates) != 1:
                        self.ambiguous.append(child.id)
                        continue
                    target = candidates[0]
                    branch = [child, *database.children(child)]
                    self._removed[source.id].update(feature.id for feature in branch)
                    self._added[target.id].append((child, branch))

        if self.ambiguous:
            logger.log_warning(
                f"{len(self.ambiguous)} transcript(s) lie outside the gene "
                "they are bound to, and more than one fragment of that "
                "duplicate-ID gene could hold them; left as the input binds "
                f"them (e.g. {self.ambiguous[0]!r}).")

    @property
    def repaired_branches(self):
        return sum(len(branches) for branches in self._added.values())

    def reopen_for_thread(self, dbfn):
        """Give a worker its own DB connection and the same corrected hierarchy."""
        return type(self)(type(self.database)(dbfn), self.root_types)

    def children(self, feature, *args, **kwargs):
        root_id = feature.id if hasattr(feature, "id") else str(feature)
        if root_id not in self._removed and root_id not in self._added:
            return self.database.children(feature, *args, **kwargs)

        children = [child for child in self.database.children(feature, *args, **kwargs)
                    if child.id not in self._removed[root_id]]
        level = kwargs.get("level")
        if level is None and args:
            level = args[0]
        for direct, branch in self._added[root_id]:
            children.extend([direct] if level == 1 else branch)
        featuretype = kwargs.get("featuretype")
        if featuretype is not None:
            allowed = {featuretype} if isinstance(featuretype, str) else set(featuretype)
            children = [child for child in children if child.featuretype in allowed]
        if kwargs.get("order_by") == "start":
            children.sort(key=lambda child: (child.start, child.end, child.id))
        return iter(children)

    def __getattr__(self, name):
        if name == "children_batched_features":
            batched = getattr(self.database, name)

            def corrected(anchors, **kwargs):
                result = batched(anchors, **kwargs)
                for anchor in anchors:
                    root_id = anchor.id if hasattr(anchor, "id") else str(anchor)
                    if root_id in self._removed or root_id in self._added:
                        result[root_id] = tuple(self.children(anchor, **kwargs))
                return result

            return corrected
        return getattr(self.database, name)

    def __getitem__(self, key):
        return self.database[key]


def bind_same_seqid_parents(database, root_types, counts=None):
    """Return the original DB when no correction is needed (the common case).

    ``counts``, when given, receives ``repaired`` and ``ambiguous`` so the
    caller can record both even when the original DB is returned.
    """
    overlay = SameSeqidParentOverlay(database, root_types)
    if counts is not None:
        counts["repaired"] = overlay.repaired_branches
        counts["ambiguous"] = len(overlay.ambiguous)
    return overlay if overlay.repaired_branches else database
