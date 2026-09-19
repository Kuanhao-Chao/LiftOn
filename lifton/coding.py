"""Coding semantics shared by extraction, scoring and model completion."""
import os

from lifton.exceptions import LiftOnInputError


def initial_phase(features, strand):
    """Initial incomplete-codon bases, counted once in transcript orientation.

    Later CDS phases describe codons spanning splice junctions; trimming every
    segment would delete coding bases. Older interval-only callers have no phase.
    """
    cds = [f for f in features if getattr(f, 'featuretype', 'CDS') == 'CDS']
    if not cds:
        return 0
    first = max(cds, key=lambda f: f.end) if strand == '-' else min(cds, key=lambda f: f.start)
    phase = str(getattr(first, 'frame', '0'))
    if phase == '.':
        return 0
    if phase not in ('0', '1', '2'):
        raise LiftOnInputError(f'Invalid initial CDS phase {phase!r} for {getattr(first, "id", "CDS")}')
    return int(phase)


def phase_adjusted_lengths(lengths, phase):
    """Lengths contributing to translated sequence, preserving segment count."""
    adjusted = []
    for length in lengths:
        skipped = min(length, phase)
        adjusted.append(length - skipped)
        phase -= skipped
    return adjusted


#: NCBI genetic code used when an annotation declares none. GFF3 has no default
#: of its own; NCBI's convention is that an absent ``transl_table`` means 1.
DEFAULT_TRANSL_TABLE = 1

#: The GFF3 attribute NCBI writes the genetic code into, on the CDS rows.
TRANSL_TABLE_ATTRIBUTE = 'transl_table'

#: Codon the ORF search treats as a start. Several tables admit more (table 2
#: adds ATT/ATC/ATA/GTG, table 11 adds TTG/CTG/...), but widening the scan would
#: change the table-1 result too, which is a separate decision with its own
#: evidence. Every table LiftOn can be given accepts ATG, so the search finds
#: the same ORFs it finds today and stops where the declared table says to stop.
ORF_START_CODON = 'ATG'

_TABLE_CACHE = {}


def _codon_table(table_id):
    """The Biopython table for this NCBI id, or LiftOnInputError."""
    cached = _TABLE_CACHE.get(table_id)
    if cached is not None:
        return cached
    from Bio.Data import CodonTable
    try:
        table = CodonTable.unambiguous_dna_by_id[table_id]
    except KeyError:
        known = ', '.join(str(i) for i in sorted(
            CodonTable.unambiguous_dna_by_id))
        raise LiftOnInputError(
            f'Unknown NCBI translation table {table_id!r}. '
            f'Known tables: {known}.') from None
    _TABLE_CACHE[table_id] = table
    return table


def parse_transl_table(value, context=None):
    """The integer NCBI table id in ``value``, validated against Biopython."""
    try:
        table_id = int(str(value).strip())
    except (TypeError, ValueError):
        where = f' on {context}' if context else ''
        raise LiftOnInputError(
            f'Invalid {TRANSL_TABLE_ATTRIBUTE}={value!r}{where}; '
            f'it must be an NCBI translation table number.') from None
    try:
        _codon_table(table_id)
    except LiftOnInputError as error:
        where = f' on {context}' if context else ''
        raise LiftOnInputError(f'{error}{where}') from None
    return table_id


def table_from_attributes(attributes, context=None):
    """The declared table in one feature's attributes, or None."""
    if not attributes:
        return None
    value = attributes.get(TRANSL_TABLE_ATTRIBUTE)
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        value = value[0]
    return parse_transl_table(value, context)


def resolve_transl_table(features, context=None, default=DEFAULT_TRANSL_TABLE):
    """The genetic code this transcript's CDS rows declare.

    The CDS rows of one transcript are segments of a single discontinuous
    feature and so declare one code; the first in the given order wins, and a
    disagreement is reported once rather than being resolved silently. Anything
    that declares nothing -- an annotation without the attribute, a
    miniprot-derived model -- gets ``default``, which reproduces every result
    LiftOn produced before the attribute was read.
    """
    resolved = None
    for feature in features or ():
        if getattr(feature, 'featuretype', 'CDS') != 'CDS':
            continue
        table_id = table_from_attributes(
            getattr(feature, 'attributes', None),
            context or getattr(feature, 'id', None))
        if table_id is None:
            continue
        if resolved is None:
            resolved = table_id
        elif table_id != resolved:
            from lifton import logger
            logger.log_warning(
                f'{context or "transcript"}: CDS rows declare more than one '
                f'{TRANSL_TABLE_ATTRIBUTE} ({resolved} and {table_id}); '
                f'using {resolved}, the first in transcript order.')
            break
    return default if resolved is None else resolved


def stop_codons(table_id=DEFAULT_TRANSL_TABLE):
    """The stop codons of this table. Table 2 reads TGA as tryptophan and AGA
    and AGG as stops, so scanning with table 1's set both invents premature
    stops and runs through real ones."""
    return frozenset(_codon_table(table_id).stop_codons)


_FAST_TABLES = {}


def _codon_map(table_id):
    """``{codon: amino acid}`` over all 64 unambiguous codons, stops as ``*``.

    Biopython reaches its table through ``CodonTable.__getitem__``, a
    Python-level call made once per codon. Profiling Step 7 on whole genomes
    (notes/step7_profile_2026-09.md) counted **62 million** of them on rice and
    **257 million** on dog to cat -- 15 s and 60 s of pure interpreter
    overhead, inside a translation block that is 8-9 % of the phase. A plain
    dict answers the same question at C speed.
    """
    cached = _FAST_TABLES.get(table_id)
    if cached is None:
        table = _codon_table(table_id)
        cached = dict(table.forward_table)
        for stop in table.stop_codons:
            # setdefault, NOT assignment. Tables 27, 28 and 31 (Karyorelict,
            # Condylostoma and Blastocrithidia nuclear) list codons that are
            # BOTH a stop and an amino acid, and Biopython translates those as
            # the amino acid. Overwriting with '*' disagreed with it on TGA
            # under table 27 and TAA under 28 and 31 -- caught by translating
            # all 64 codons under all 25 tables, which is the only reason this
            # is a fixed bug rather than a shipped one.
            cached.setdefault(stop, '*')
        _FAST_TABLES[table_id] = cached
    return cached


def _legacy_translate():
    return os.environ.get("LIFTON_LEGACY_TRANSLATE", "").strip().lower() \
        not in ("", "0", "false", "no", "off")


def translate(dna, table_id=DEFAULT_TRANSL_TABLE):
    """Translate ``dna`` with this table, dropping a trailing partial codon.

    Biopython has historically truncated an incomplete terminal codon while
    warning that the behaviour may change. Every LiftOn translation depends on
    that result, so it is made explicit -- and warning-free -- in one place.

    The codon-map path answers only what a 64-entry dict can answer exactly.
    Anything else -- an ``N``, an ambiguity code, a gap -- raises ``KeyError``
    and the whole sequence goes to Biopython untouched, so the two agree by
    construction rather than by argument. ``LIFTON_LEGACY_TRANSLATE=1`` forces
    the Biopython path.
    """
    dna = str(dna)
    complete_length = len(dna) - (len(dna) % 3)
    if not complete_length:
        return ''
    dna = dna[:complete_length]
    if not _legacy_translate():
        codons = _codon_map(table_id)
        try:
            return "".join([codons[dna[i:i + 3]]
                            for i in range(0, complete_length, 3)])
        except KeyError:
            pass
        try:
            upper = dna.upper()
            return "".join([codons[upper[i:i + 3]]
                            for i in range(0, complete_length, 3)])
        except KeyError:
            pass
    from Bio.Seq import Seq
    return str(Seq(dna).translate(table=table_id))
