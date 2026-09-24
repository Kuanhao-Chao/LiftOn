"""Coding semantics shared by extraction, scoring and model completion."""
import os
import re
from typing import NamedTuple, Optional

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


#: The genetic code each reference transcript declares, from Step 3's
#: sidecar: ``{reference transcript id: table}``, non-standard codes only.
#: Installed once before any thread or fork, read-only afterwards.
_TRANSCRIPT_TABLES = {}


def install_transcript_tables(mapping):
    """Replace the per-transcript declared codes (start of each run: ``{}``)."""
    _TRANSCRIPT_TABLES.clear()
    _TRANSCRIPT_TABLES.update(mapping or {})


def transcript_table(ref_id, default=DEFAULT_TRANSL_TABLE):
    """The code this reference transcript declares, or ``default``."""
    if ref_id is None:
        return default
    return _TRANSCRIPT_TABLES.get(ref_id, default)


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


# ---------------------------------------------------------------------------
# transl_except: codons an annotation declares to translate differently
# ---------------------------------------------------------------------------

#: The GFF3 attribute NCBI writes a codon-level translation exception into:
#: ``transl_except=(pos:complement(30445548..30445550),aa:Sec)`` is SEPHS2's
#: selenocysteine, read through a UGA the genome translates as a stop.
TRANSL_EXCEPT_ATTRIBUTE = 'transl_except'

#: ``aa:TERM`` is a stop completed by polyadenylation (vertebrate mitochondria:
#: a trailing T or TA), not an amino acid.
TRANSL_EXCEPT_TERM = 'TERM'


#: Longest location a transl_except may name (a codon is 3 bp; dog RefSeq has
#: a real 6,910-bp aa:Other range).
MAX_LOCATION_SPAN = 1_000_000


class TranslExcept(NamedTuple):
    """One declared exception, placed on the reference protein."""
    #: 0-based residue index into the reference protein, or None when the
    #: declared codon is not an in-frame codon of the CDS.
    residue: Optional[int]
    #: The amino acid exactly as declared: 'Sec', 'Other', 'Met', 'TERM', ...
    aa: str
    #: Bases the location covers: 3, or 1-2 for a partial TERM.
    length: int
    #: True when the reference codon is a stop the protein continues through
    #: (Sec, Pyl, stop readthrough, an amino acid declared over a stop).
    readthrough: bool
    #: The reference codon's bases in transcript orientation, when known. A
    #: declaration that is neither a read-through nor a start belongs to this
    #: codon, and is carried to a target model only where its codon is the same.
    codon: str = ''


def transl_except_values(attributes):
    """The ``transl_except`` values in one feature's attributes, as strings."""
    if not attributes:
        return []
    value = attributes.get(TRANSL_EXCEPT_ATTRIBUTE)
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(item) for item in value]


def _split_top_level(text):
    parts, depth, current = [], 0, []
    for char in text:
        if char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
        if char == ',' and depth == 0:
            parts.append(''.join(current))
            current = []
        else:
            current.append(char)
    parts.append(''.join(current))
    return [part.strip() for part in parts if part.strip()]


_AA_SEPARATOR = re.compile(r',\s*aa\s*:')


def parse_transl_except(values, errors=None):
    """``[(location, aa), ...]`` for every ``(pos:LOCATION,aa:AA)`` in values.

    gffutils hands back one decoded value per exception
    (``['(pos:25802093..25802095,aa:Sec)', ...]``); a raw GFF3 value keeps
    them in one string with ``%2C`` for the commas. Both are accepted: the
    values are rejoined, decoded, and scanned by parenthesis depth, which also
    keeps the commas inside ``join(...)`` where they belong. Space around
    ``aa:`` is accepted.

    A malformed declaration raises ``ValueError``, unless ``errors`` is a
    list: then it is appended there and the scan carries on with the next one,
    so one bad value no longer costs a transcript all of its others.
    """
    text = ','.join(values).replace('%2C', ',').replace('%2c', ',')
    found = []
    index = 0

    def bad(message):
        if errors is None:
            raise ValueError(message)
        errors.append(message)

    while True:
        start = text.find('(pos:', index)
        if start < 0:
            return found
        depth, end = 0, start
        while end < len(text):
            if text[end] == '(':
                depth += 1
            elif text[end] == ')':
                depth -= 1
                if depth == 0:
                    break
            end += 1
        if depth:
            bad(f'unbalanced parentheses in {text[start:]!r}')
            return found
        body = text[start + len('(pos:'):end]
        separators = list(_AA_SEPARATOR.finditer(body))
        if not separators:
            bad(f'no amino acid in {body!r}')
        else:
            last = separators[-1]
            found.append((body[:last.start()].strip(), body[last.end():].strip()))
        index = end + 1


def parse_location(text):
    """The genomic positions a transl_except location names, as written.

    Accepts ``a..b``, a single base, ``complement(...)``, ``join(...)`` and
    ``order(...)``, nested, with the partial markers ``<``/``>`` ignored.
    Anything else raises ``ValueError``.
    """
    text = text.strip().replace('<', '').replace('>', '')
    if text.startswith('complement(') and text.endswith(')'):
        return parse_location(text[len('complement('):-1])[::-1]
    for prefix in ('join(', 'order('):
        if text.startswith(prefix) and text.endswith(')'):
            positions = []
            for part in _split_top_level(text[len(prefix):-1]):
                positions.extend(parse_location(part))
            return positions
    if '..' in text:
        first, last = (int(value) for value in text.split('..', 1))
        if last < first:
            raise ValueError(f'reversed range {text!r}')
        if last - first + 1 > MAX_LOCATION_SPAN:
            # A declaration names a codon; the longest real one seen (dog
            # RefSeq, aa:Other) spans 6,910 bp. Materialising a malformed
            # chromosome-length range would exhaust memory.
            raise ValueError(f'range {text!r} is longer than {MAX_LOCATION_SPAN:,} bp')
        return list(range(first, last + 1))
    return [int(text)]


def _spliced_offsets(positions, merged_intervals, strand):
    """Offsets of genomic positions within the spliced CDS, in transcript
    orientation, or None if any position lies outside the CDS."""
    total = sum(end - start + 1 for start, end in merged_intervals)
    offsets = []
    for position in positions:
        before, found = 0, None
        for start, end in merged_intervals:
            if start <= position <= end:
                found = before + (position - start)
                break
            before += end - start + 1
        if found is None:
            return None
        offsets.append(total - 1 - found if strand == '-' else found)
    return sorted(offsets)


def transl_excepts_for(values, merged_intervals, strand, phase, protein, context=None,
                       bases_at=None):
    """Place a transcript's declared exceptions on its reference protein.

    ``merged_intervals``, ``strand`` and ``phase`` must be exactly what built
    ``protein`` (``extract_sequence.get_protein_sequence``: the merged
    CDS/start/stop intervals, the initial phase), so that residue ``r`` of the
    protein is the codon at spliced offset ``phase + 3r``. A location is
    placed only when its bases are consecutive in the spliced CDS -- a codon
    split across an exon junction is -- and start in frame; otherwise its
    residue is None and it changes nothing. A malformed value is reported and
    skipped; the others are kept. ``bases_at(positions)``, when given, returns
    the reference bases at those positions in transcript orientation; it
    records each placed codon (:attr:`TranslExcept.codon`).
    """
    placed = []
    errors = []
    declared = parse_transl_except(values, errors)
    if errors:
        from lifton import logger
        for error in errors:
            logger.log_warning(f'{context or "transcript"}: unreadable '
                               f'{TRANSL_EXCEPT_ATTRIBUTE} ignored ({error}).')
    for location, aa in declared:
        try:
            positions = parse_location(location)
        except ValueError as error:
            from lifton import logger
            logger.log_warning(f'{context or "transcript"}: unreadable '
                               f'{TRANSL_EXCEPT_ATTRIBUTE} location {location!r} '
                               f'ignored ({error}).')
            continue
        offsets = _spliced_offsets(positions, merged_intervals, strand)
        residue = None
        if (offsets and 1 <= len(offsets) <= 3
                and offsets == list(range(offsets[0], offsets[0] + len(offsets)))):
            offset = offsets[0] - phase
            if offset >= 0 and offset % 3 == 0:
                residue = offset // 3
        readthrough = (aa != TRANSL_EXCEPT_TERM and len(positions) == 3
                       and residue is not None and residue < len(protein) - 1
                       and protein[residue] == '*')
        codon = ''
        if bases_at is not None and residue is not None and len(positions) == 3:
            try:
                codon = bases_at(positions).upper()
            except Exception:
                codon = ''
        placed.append(TranslExcept(residue, aa, len(positions), readthrough, codon))
    return tuple(placed)
