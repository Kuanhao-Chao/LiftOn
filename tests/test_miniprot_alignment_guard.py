"""Tier-1 audit fix — an unusable miniprot candidate must not drop a whole locus.

`align.lifton_parasail_align` returns ``None`` when it cannot build a protein
alignment (its ``protein_seq is None`` early return) — e.g. a miniprot mRNA with no
``CDS``/``stop_codon`` children, so the coding sequence is empty. That happens with a
truncated `miniprot.gff3` or any user-supplied ``-M``.

`LiftOn_miniprot_alignment` used to assign that ``None`` and immediately read
``m_lifton_aln.identity`` → ``AttributeError``. `process_locus` catches per-locus
exceptions, so ONE unusable miniprot candidate silently removed the ENTIRE Liftoff
gene from the output even though its DNA lift was fine.

It also has to stop advertising ``has_valid_miniprot`` when no alignment survived:
callers gate ``protein_maximization.chaining_algorithm(...)`` on that flag and would
dereference the ``None``.

The companion test pins the *ordering contract* of the whole-side fallbacks — see
`test_chaining_fallback_order.py`.
"""
from types import SimpleNamespace

from intervaltree import IntervalTree

from lifton import lifton_utils


class _FakeEntry:
    def __init__(self, fid, start, end, seqid="chr1", strand="+"):
        self.id = fid
        self.start = start
        self.end = end
        self.seqid = seqid
        self.strand = strand
        self.featuretype = "mRNA"
        self.attributes = {"ID": [fid]}


class _FakeMiniprotDb:
    """Minimal m_feature_db: __getitem__ + children()."""
    def __init__(self, entries):
        self._entries = entries

    def __getitem__(self, key):
        return self._entries[key]

    def children(self, *_a, **_k):
        return iter(())          # NO CDS children -> empty coding seq -> None aln


def _drive(monkeypatch, align_returns):
    """Run LiftOn_miniprot_alignment with one miniprot candidate whose parasail
    alignment yields `align_returns`."""
    m_entry = _FakeEntry("mp1", 100, 400)
    m_db = _FakeMiniprotDb({"mp1": m_entry})
    transcript = SimpleNamespace(id="tx1", start=100, end=400, seqid="chr1")

    tree = IntervalTree()
    tree.addi(100, 400, "gene1")
    tree_dict = {"chr1": tree}

    monkeypatch.setattr(lifton_utils.align, "lifton_parasail_align",
                        lambda *a, **k: align_returns)

    status = SimpleNamespace(miniprot=0.0)
    return lifton_utils.LiftOn_miniprot_alignment(
        "chr1", transcript, {"tx1": ["mp1"]}, m_db, tree_dict,
        None, {"tx1": "PROT"}, "tx1", status)


def test_unusable_miniprot_candidate_is_skipped_not_fatal(monkeypatch):
    """A None parasail alignment must be skipped, not dereferenced."""
    aln, has_valid = _drive(monkeypatch, None)
    assert aln is None
    # Must NOT advertise a valid miniprot alignment — the caller would pass this
    # None into chaining_algorithm and crash on `.cds_children`.
    assert has_valid is False


def test_usable_miniprot_candidate_still_selected(monkeypatch):
    """Control: a real alignment is still picked up and scored."""
    good = SimpleNamespace(identity=0.75, cds_children=[object()])
    aln, has_valid = _drive(monkeypatch, good)
    assert aln is good
    assert has_valid is True
