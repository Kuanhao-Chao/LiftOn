"""Give sparse coding annotations an explicit, deterministic internal hierarchy.

Only parentless CDS and genes with direct CDS children are normalized. Ordinary
transcript hierarchies are returned untouched. Source files/databases are never
modified; generated IDs and original CDS IDs are recorded alongside the copy.
"""
from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

from lifton import coreutils
from lifton.exceptions import LiftOnInputError
from lifton.run_manifest import atomic_write_json


SCHEMA_VERSION = 1


def _feature_key(feature):
    logical_id = feature.attributes.get('ID', [feature.id])[0]
    return (str(logical_id), feature.featuretype, feature.seqid,
            feature.start, feature.end, feature.strand, feature.id)


def _anchors(db):
    # The normal gffutils path answers this with indexed relations lookups,
    # without reconstructing every CDS in a mammalian annotation.
    if type(db).__module__.startswith('gffutils'):
        sql = """SELECT f.id FROM features f WHERE f.featuretype = 'CDS'
                 AND NOT EXISTS (SELECT 1 FROM relations r JOIN features p ON p.id=r.parent
                                 WHERE r.child=f.id AND r.level=1)
                 UNION
                 SELECT DISTINCT p.id FROM relations r
                 JOIN features c ON c.id=r.child JOIN features p ON p.id=r.parent
                 WHERE r.level=1 AND c.featuretype='CDS' AND p.featuretype='gene'"""
        return [db[row[0]] for row in db.execute(sql)]
    result = OrderedDict()
    for cds in db.features_of_type('CDS'):
        parents = list(db.parents(cds, level=1))
        if not parents:
            result[cds.id] = cds
        for parent in parents:
            if parent.featuretype == 'gene':
                result[parent.id] = parent
    return list(result.values())


def normalize_sparse_coding(ref_db, out_dir, strict=False):
    """Return normalization metadata, or None when no sparse model is selected.

    ``strict`` mirrors ``--strict-gff``: a reference whose hierarchy does not
    parse is refused rather than reported.
    """
    db = ref_db.db_connection
    anchors = sorted(_anchors(db), key=_feature_key)
    reserved = set()
    declared_ids = set()
    cds_parent_references = []
    ordinary_transcript_aliases = {}
    for feature in db.all_features():
        reserved.add(feature.id)
        reserved.update(feature.attributes.get('ID', []))
        declared_ids.add(feature.id)
        declared_ids.update(feature.attributes.get('ID', []))
        if feature.featuretype == 'CDS':
            cds_parent_references.extend(
                (feature.id, parent) for parent in feature.attributes.get('Parent', []))
        if feature.featuretype in ('mRNA', 'transcript'):
            for key in ('ID', 'protein_id', 'transcript_id'):
                for alias in feature.attributes.get(key, []):
                    ordinary_transcript_aliases.setdefault(alias, set()).add(feature.id)
            ordinary_transcript_aliases.setdefault(feature.id, set()).add(feature.id)
    dangling = sorted((feature_id, parent) for feature_id, parent in cds_parent_references
                      if parent not in declared_ids)
    if dangling:
        feature_id, parent = dangling[0]
        message = (f'CDS {feature_id!r} names missing Parent {parent!r}; '
                   f'repair the hierarchy first')
        if strict:
            raise LiftOnInputError(message)
        # Fatal by default aborted a 144,415-transcript human lift over 111
        # broken rows in NCBI_RefSeq_no_rRNA.gff -- an artefact of the rRNA
        # filtering that produced it, which left 14 transcripts' CDS behind.
        # v1.0.12 completed that lift. A reference that violates the spec is
        # what --strict-gff is for; the default path reports and continues,
        # exactly as the input validator does a few steps earlier.
        from lifton import logger
        logger.log_warning(
            f'{len(dangling)} CDS name a Parent no row declares (e.g. '
            f'{message}). They are excluded from reference normalization and '
            f'their transcripts cannot be lifted. Pass --strict-gff to make '
            f'this fatal.')
        damaged = {feature_id for feature_id, _ in dangling}
        anchors = [anchor for anchor in anchors if anchor.id not in damaged]
    if not anchors:
        return None

    def allocate(base):
        value, suffix = base, 1
        while value in reserved:
            value = f'{base}_{suffix}'
            suffix += 1
        reserved.add(value)
        return value

    groups = OrderedDict()
    for anchor in anchors:
        if anchor.featuretype == 'CDS':
            if anchor.attributes.get('Parent'):
                # A missing parent is malformed input, not a flat annotation.
                raise LiftOnInputError(f'CDS {anchor.id!r} names a missing Parent; repair the hierarchy first')
            logical = anchor.attributes.get('ID', [anchor.id])[0]
            key = ('CDS', logical)
            groups.setdefault(key, (anchor, []))[1].append(anchor)
        else:
            children = sorted(db.children(anchor, level=1), key=_feature_key)
            unsupported = {c.featuretype for c in children} - {'CDS', 'start_codon', 'stop_codon', 'exon'}
            if unsupported:
                raise LiftOnInputError(
                    f'Gene {anchor.id!r} mixes direct CDS with {sorted(unsupported)} children; '
                    'attach each CDS to its transcript explicitly')
            groups[('gene', anchor.id)] = (anchor, children)

    replacements, consumed, mappings = {}, set(), []
    for (kind, logical), (anchor, children) in groups.items():
        cds = [c for c in children if c.featuretype == 'CDS']
        if not cds:
            continue
        if any(len(c.attributes.get('Parent', [])) > 1 for c in children):
            raise LiftOnInputError(f'Sparse CDS {logical!r} has multiple parents; provide an explicit hierarchy')
        for key in ('protein_id', 'transcript_id', 'transl_table'):
            values = {value for c in cds for value in c.attributes.get(key, [])}
            if len(values) > 1:
                raise LiftOnInputError(f'Sparse CDS {logical!r} has conflicting {key} values')
        ordered_cds = sorted(cds, key=lambda c: (c.start, c.end, c.id))
        if any(a.end >= b.start for a, b in zip(ordered_cds, ordered_cds[1:])):
            raise LiftOnInputError(f'Sparse CDS {logical!r} has overlapping segments; provide an explicit hierarchy')
        if len({(c.seqid, c.strand) for c in children}) != 1 or any(c.strand not in ('+', '-') for c in cds):
            raise LiftOnInputError(f'Sparse CDS {logical!r} has ambiguous sequence/strand; provide an explicit hierarchy')
        start, end = min(c.start for c in children), max(c.end for c in children)
        if kind == 'gene' and (anchor.seqid != cds[0].seqid or anchor.strand != cds[0].strand
                               or anchor.start > start or anchor.end < end):
            raise LiftOnInputError(f'Sparse gene {logical!r} does not contain its children on the same strand/sequence')
        gene = coreutils.clone_feature(anchor)
        gene.featuretype = 'gene'
        gene.start, gene.end, gene.frame = start, end, '.'
        if kind == 'gene':
            gene.start, gene.end = anchor.start, anchor.end
        gene.id = allocate(f'{logical}__lifton_gene') if kind == 'CDS' else anchor.id
        gene.attributes['ID'] = [gene.id]
        gene.attributes['orig_feature_id'] = [logical]
        gene.attributes.setdefault('gene_biotype', ['protein_coding'])
        trans = coreutils.clone_feature(gene)
        trans.id = allocate(f'{logical}__lifton_transcript')
        trans.featuretype = 'mRNA'
        trans.attributes['ID'] = [trans.id]
        trans.attributes['Parent'] = [gene.id]
        block = [gene, trans]
        exons = [c for c in children if c.featuretype == 'exon']
        if exons and any(not any(e.start <= c.start and e.end >= c.end for e in exons) for c in cds):
            raise LiftOnInputError(f'Sparse CDS {logical!r} is not contained in its declared exons')
        if not exons:
            # CDS-only models have no inferred UTR; retain their exact segments.
            for index, child in enumerate(sorted(cds, key=lambda c: (c.start, c.end, c.id)), 1):
                exon = coreutils.clone_feature(child)
                exon.id = allocate(f'{logical}__lifton_exon_{index}')
                exon.featuretype, exon.frame = 'exon', '.'
                exon.attributes['ID'] = [exon.id]
                exon.attributes['Parent'] = [trans.id]
                block.append(exon)
        for child in children:
            copied = coreutils.clone_feature(child)
            copied.attributes['Parent'] = [trans.id]
            block.append(copied)
        consumed.update(c.id for c in children)
        consumed.add(anchor.id)
        replacements[anchor.id] = block
        aliases = {logical}
        for child in cds:
            for key in ('ID', 'protein_id', 'transcript_id'):
                aliases.update(child.attributes.get(key, []))
        mappings.append({'source_id': logical, 'source_type': kind, 'gene_id': gene.id,
                         'transcript_id': trans.id, 'aliases': sorted(aliases),
                         'cds_ids': [c.attributes.get('ID', [c.id])[0] for c in cds]})

    if not replacements:
        return None
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / 'normalized_reference.gff3'
    with output.open('w') as handle:
        directives = getattr(ref_db, 'directives', [])
        handle.write('##gff-version 3\n')
        for directive in directives:
            if not directive.startswith(('##gff-version', '##FASTA')):
                handle.write(directive.rstrip('\n') + '\n')
        for feature in db.all_features():
            if feature.id in replacements:
                for item in replacements[feature.id]:
                    handle.write(str(item) + '\n')
            elif feature.id not in consumed:
                handle.write(str(feature) + '\n')
    mapping_path = out_dir / 'reference_model_mapping.json'
    from lifton.annotation_cache import source_fingerprint
    aliases = ordinary_transcript_aliases
    for model in mappings:
        for alias in [model['transcript_id'], *model['aliases']]:
            aliases.setdefault(alias, set()).add(model['transcript_id'])
    document = {'schema_version': SCHEMA_VERSION, 'source': ref_db.file_name,
                'source_fingerprint': source_fingerprint(ref_db.file_name),
                'normalized_annotation': str(output), 'models': mappings,
                'aliases': {alias: sorted(ids) for alias, ids in sorted(aliases.items())}}
    atomic_write_json(mapping_path, document)
    return {'annotation': str(output), 'mapping': str(mapping_path), 'models': len(mappings),
            'annotation_fingerprint': source_fingerprint(output),
            'mapping_fingerprint': source_fingerprint(mapping_path)}


def alias_fasta(source, mapping, destination):
    """Write supplied sequences using normalized transcript IDs; reject ambiguity."""
    import hashlib
    import json
    import os
    import tempfile
    from Bio import SeqIO
    aliases = json.loads(Path(mapping).read_text())['aliases']
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    seen = {}
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=destination.parent, delete=False) as output:
            temporary = Path(output.name)
            for sequence in SeqIO.parse(str(source), 'fasta'):
                targets = aliases.get(sequence.id, [sequence.id])
                if len(targets) != 1:
                    raise LiftOnInputError(f'Ambiguous supplied FASTA ID {sequence.id!r}: {targets}')
                target = targets[0]
                text = str(sequence.seq)
                digest = hashlib.sha256(text.upper().encode()).hexdigest()
                if target in seen:
                    if seen[target] != digest:
                        raise LiftOnInputError(f'Conflicting supplied FASTA aliases for {target!r}')
                    continue
                seen[target] = digest
                output.write(f'>{target}\n{text}\n')
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return str(destination)
