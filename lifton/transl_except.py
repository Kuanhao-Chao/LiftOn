"""Codons the reference annotation declares to translate differently.

RefSeq writes ``transl_except=(pos:complement(30445548..30445550),aa:Sec)`` on
the CDS rows of every selenoprotein, stop-readthrough isoform, non-AUG start
and mitochondrial stop completed by polyadenylation. LiftOn ignored it: the
reference and the lifted protein both carry ``*`` at a selenocysteine, every
identity cut-off read it as a premature stop, and on a same-species
GRCh38 -> CHM13 lift the 53 selenoprotein transcripts scored a mean 0.666 --
SEPHS2 lost its first 118 residues to a miniprot fragment. The carried
attribute also kept reference coordinates.

This module is the registry of those declarations, keyed by reference
transcript id and filled once, after Step 3 and before any thread or fork:

* :func:`readthrough` gives the scorers the reference residues a model may
  read through (``lifton_parasail_align``, the ORF search, chaining, stop
  completion). A transcript with no declaration gets None and every caller
  takes its original path.
* :func:`cds_values` gives the serializer what a transcript's CDS rows should
  carry: the declarations rewritten into the model's own (target) coordinates,
  or nothing -- never the reference's coordinates.

It imports only ``coding`` and ``logger`` at module level so the serializer can
import it; alignment is imported where it is used.
"""
import threading

from lifton import coding, logger

_REGISTRY = {}
_READTHROUGH = {}
_HANDLES = {}
#: What each output transcript's rows carry: ``{model id: (written, dropped)}``.
#: Keyed by model so a model rendered twice (staged, then written) or again
#: after its CDS moved is counted once, as it was last rendered.
_RENDERED = {}
_LOCK = threading.Lock()

_CODING_CHILDREN = ("start_codon", "CDS", "stop_codon")


def clear():
    """Forget every declaration and handle (start of each run)."""
    _REGISTRY.clear()
    _READTHROUGH.clear()
    _HANDLES.clear()
    _RENDERED.clear()


def install(mapping, ref_proteins=None, target_fasta=None):
    """Register ``{reference transcript id: (coding.TranslExcept, ...)}``.

    A read-through residue is kept only if the protein actually loaded agrees
    it is a recoded codon: ``*`` from genome translation, or the ``U``/``O``/``X``
    an NCBI protein FASTA (``-P``) writes there. Returns the number of
    transcripts registered.
    """
    clear()
    for ref_id, entries in (mapping or {}).items():
        if not entries:
            continue
        protein = None
        if ref_proteins is not None and ref_id in ref_proteins.keys():
            protein = str(ref_proteins[ref_id]).upper()
        checked = []
        for entry in entries:
            if entry.readthrough and protein is not None and not (
                    entry.residue < len(protein) and protein[entry.residue] in "*UOX"):
                logger.log_warning(
                    f"{ref_id}: {coding.TRANSL_EXCEPT_ATTRIBUTE} residue "
                    f"{entry.residue + 1} ({entry.aa}) is not a recoded codon in "
                    f"the reference protein; it is not read through.")
                entry = entry._replace(readthrough=False)
            checked.append(entry)
        _REGISTRY[ref_id] = tuple(checked)
        residues = frozenset(e.residue for e in checked if e.readthrough)
        if residues:
            _READTHROUGH[ref_id] = residues
    _HANDLES.update(ref_proteins=ref_proteins, target_fasta=target_fasta)
    return len(_REGISTRY)


def get(ref_id):
    """The declarations for this reference transcript, or None."""
    return _REGISTRY.get(ref_id) if ref_id is not None else None


def readthrough(ref_id):
    """Reference residues this transcript reads through, or None."""
    return _READTHROUGH.get(ref_id) if ref_id is not None else None


def counts():
    """``{"emitted": n, "dropped": n}`` for the values written so far."""
    with _LOCK:
        return {"emitted": sum(written for written, _ in _RENDERED.values()),
                "dropped": sum(dropped for _, dropped in _RENDERED.values())}


def scan_reference(ref_db, ref_fai):
    """Declarations for a run whose Step 3 was skipped (``-P`` and ``-T``).

    One materialised pass finds the CDS rows that carry the attribute; only
    then is each of their transcripts read, so nothing queries the handle while
    its scan is open.
    """
    from lifton import extract_sequence
    db = ref_db.db_connection
    carriers = [cds for cds in db.features_of_type("CDS")
                if coding.TRANSL_EXCEPT_ATTRIBUTE in (cds.attributes or {})]
    parents = []
    for cds in carriers:
        for parent in cds.attributes.get("Parent", []):
            if parent not in parents:
                parents.append(parent)
    mapping = {}
    for parent_id in parents:
        try:
            feature = db[parent_id]
            children = [child for child in db.children(feature, level=1, order_by="start")
                        if child.featuretype in _CODING_CHILDREN]
            entries = extract_sequence.transl_excepts_of(feature, ref_fai, children)
        except Exception as error:
            logger.log_warning(f"{parent_id}: {coding.TRANSL_EXCEPT_ATTRIBUTE} "
                               f"could not be placed ({error}).")
            continue
        if entries:
            mapping[parent_id] = entries
    return mapping


def _render_location(positions, strand):
    runs = []
    for position in sorted(positions):
        if runs and position == runs[-1][1] + 1:
            runs[-1][1] = position
        else:
            runs.append([position, position])
    parts = [f"{a}..{b}" if a != b else f"{a}" for a, b in runs]
    location = parts[0] if len(parts) == 1 else f"join({','.join(parts)})"
    return f"complement({location})" if strand == "-" else location


def _coding_positions(cds_entries, strand, fasta):
    """``(coding sequence, genomic position of each coding base)`` in
    transcript orientation, after the initial phase."""
    from Bio.Seq import Seq
    segments = sorted((entry.start, entry.end) for entry in cds_entries)
    if strand == "-":
        segments.reverse()
    seqid = cds_entries[0].seqid
    sequence, positions = [], []
    for start, end in segments:
        piece = str(fasta[seqid][start - 1:end]).upper()
        if strand == "-":
            sequence.append(str(Seq(piece).reverse_complement()))
            positions.extend(range(end, start - 1, -1))
        else:
            sequence.append(piece)
            positions.extend(range(start, end + 1))
    phase = coding.initial_phase(cds_entries, strand)
    return "".join(sequence)[phase:], positions[phase:]


def _remap(trans, cds_entries, entries):
    """The declarations rewritten onto this model, and how many were dropped."""
    fasta = _HANDLES.get("target_fasta")
    proteins = _HANDLES.get("ref_proteins")
    ref_id = trans.ref_tran_id
    strand = getattr(trans.entry, "strand", None) or cds_entries[0].strand
    if (fasta is None or proteins is None or ref_id not in proteins.keys()
            or len({entry.seqid for entry in cds_entries}) != 1
            or strand not in ("+", "-")):
        return (), len(entries)
    coding_seq, positions = _coding_positions(cds_entries, strand, fasta)
    table = trans.transl_table() if hasattr(trans, "transl_table") else coding.DEFAULT_TRANSL_TABLE
    protein = coding.translate(coding_seq, table)
    reference = str(proteins[ref_id]).upper()
    if not protein or not reference:
        return (), len(entries)
    from lifton import align
    traceback = align.parasail_align_protein_base(protein, reference).traceback
    target_of = {}
    query_index = ref_index = -1
    for query_char, ref_char in zip(traceback.query, traceback.ref):
        if query_char != "-":
            query_index += 1
        if ref_char != "-":
            ref_index += 1
        if query_char != "-" and ref_char != "-":
            target_of[ref_index] = query_index
    values, dropped = [], 0
    last = len(protein) - 1
    for entry in entries:
        codon = None
        if entry.residue is not None:
            if entry.aa == coding.TRANSL_EXCEPT_TERM and entry.length < 3:
                # A stop completed by polyadenylation: the model must end in
                # that partial codon, right after the residue before it.
                partial = len(coding_seq) - 3 * len(protein)
                if partial == entry.length and target_of.get(entry.residue - 1) == last:
                    codon = positions[3 * len(protein):3 * len(protein) + entry.length]
            else:
                target = target_of.get(entry.residue)
                if target is not None:
                    if entry.readthrough:
                        # Only where the model has the recoded stop inside it:
                        # Sec -> Cys needs no exception, and a model ending
                        # there stops there.
                        keep = protein[target] == "*" and target < last
                    elif entry.aa == "Met" and entry.residue == 0:
                        keep = target == 0
                    else:
                        # Any other declaration (an amino acid over a sense
                        # codon, a full-codon TERM) belongs to the reference
                        # codon: it holds on the target only where the codon is
                        # the same. It used to be written wherever it aligned
                        # -- aa:TERM onto a sense codon on dog -> cat.
                        keep = bool(entry.codon) and (
                            coding_seq[3 * target:3 * target + 3] == entry.codon)
                    if keep:
                        codon = positions[3 * target:3 * target + 3]
                        if len(codon) != 3:
                            codon = None
        if codon:
            values.append(f"(pos:{_render_location(codon, strand)},aa:{entry.aa})")
        else:
            dropped += 1
    return tuple(values), dropped


def cds_values(trans):
    """What the CDS rows of ``trans`` should carry as ``transl_except``.

    None: leave the rows exactly as they are (no declaration, nothing carried).
    ``()``: remove the attribute -- a carried value holds reference coordinates.
    A tuple of values: write those, in target coordinates.
    Computed once per model geometry and cached on the transcript, so the
    stage-then-write paths, which render twice, count once.
    """
    cds_entries = [exon.cds.entry for exon in getattr(trans, "exons", ())
                   if getattr(exon, "cds", None) is not None]
    if not cds_entries:
        return None
    entries = get(getattr(trans, "ref_tran_id", None))
    carried = any(coding.TRANSL_EXCEPT_ATTRIBUTE in (getattr(entry, "attributes", None) or {})
                  for entry in cds_entries)
    if not entries:
        return () if carried else None
    key = tuple((entry.seqid, entry.start, entry.end, entry.strand, entry.frame)
                for entry in cds_entries)
    cached = getattr(trans, "_transl_except_rendered", None)
    if cached is not None and cached[0] == key:
        return cached[1]
    try:
        values, dropped = _remap(trans, cds_entries, entries)
    except Exception as error:
        logger.log_warning(f"{getattr(trans.entry, 'id', '?')}: "
                           f"{coding.TRANSL_EXCEPT_ATTRIBUTE} not rewritten ({error}).")
        values, dropped = (), len(entries)
    with _LOCK:
        _RENDERED[getattr(trans.entry, "id", None) or id(trans)] = (len(values), dropped)
    trans._transl_except_rendered = (key, values)
    return values
