"""Sparse coding annotations need explicit, unambiguous model hierarchies."""
import json
from pathlib import Path

import pytest

from lifton.annotation import Annotation
from lifton.exceptions import LiftOnInputError


def reference(tmp_path, lines):
    path = tmp_path / 'ref.gff3'
    path.write_text('##gff-version 3\n' + ''.join(lines))
    return Annotation(path, False, False)


def row(kind, start, end, attrs, strand='+', seqid='chr1', phase='.'):
    return f'{seqid}\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n'


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('direct_gene', [False, True])
def test_sparse_segments_keep_cds_semantics_and_form_one_transcript(tmp_path, strand, direct_gene):
    from lifton.reference_models import normalize_sparse_coding
    lines = [row('gene', 10, 50, 'ID=g;gene_biotype=protein_coding', strand)] if direct_gene else []
    parent = ';Parent=g' if direct_gene else ''
    lines += [row('CDS', 10, 21, 'ID=c;protein_id=p;transl_table=11' + parent, strand, phase='1'),
              row('CDS', 40, 50, 'ID=c;protein_id=p;transl_table=11' + parent, strand, phase='1')]
    ann = reference(tmp_path, lines)
    source = Path(ann.file_name).read_bytes()
    result = normalize_sparse_coding(ann, tmp_path / 'normalized')
    normalized = Annotation(result['annotation'], False, False).db_connection
    mapping = json.loads(Path(result['mapping']).read_text())
    assert mapping['schema_version'] == 1 and len(mapping['models']) == 1
    model = mapping['models'][0]
    assert normalized[model['gene_id']].featuretype == 'gene'
    assert normalized[model['transcript_id']].attributes['Parent'] == [model['gene_id']]
    cds = list(normalized.children(model['transcript_id'], featuretype='CDS', order_by='start'))
    assert [(c.start, c.end, c.strand, c.frame) for c in cds] == [(10, 21, strand, '1'), (40, 50, strand, '1')]
    assert all(c.attributes['ID'] == ['c'] and c.attributes['transl_table'] == ['11'] for c in cds)
    assert [(e.start, e.end) for e in normalized.children(model['transcript_id'], featuretype='exon', order_by='start')] == [(10, 21), (40, 50)]
    assert Path(ann.file_name).read_bytes() == source


def test_ordinary_hierarchy_is_a_noop(tmp_path):
    from lifton.reference_models import normalize_sparse_coding
    ann = reference(tmp_path, [row('gene', 1, 30, 'ID=g'), row('mRNA', 1, 30, 'ID=t;Parent=g'),
                               row('exon', 1, 30, 'ID=e;Parent=t'), row('CDS', 1, 30, 'ID=c;Parent=t', phase='0')])
    assert normalize_sparse_coding(ann, tmp_path / 'normalized') is None
    assert not (tmp_path / 'normalized').exists()


def test_partial_cds_phase_and_attributes_are_preserved(tmp_path):
    """Normalizing a partial model must not fill its phase or discard its evidence."""
    from lifton.reference_models import normalize_sparse_coding
    ann = reference(tmp_path, [row('CDS', 10, 24, 'ID=c;partial=true;Note=boundary', phase='.')])
    result = normalize_sparse_coding(ann, tmp_path / 'normalized')
    model = json.loads(Path(result['mapping']).read_text())['models'][0]
    normalized = Annotation(result['annotation'], False, False).db_connection
    cds = list(normalized.children(model['transcript_id'], featuretype='CDS'))
    assert len(cds) == 1
    assert cds[0].frame == '.'
    assert cds[0].attributes['partial'] == ['true']
    assert cds[0].attributes['Note'] == ['boundary']


@pytest.mark.parametrize('damage', ['dangling', 'strand', 'seqid', 'mixed', 'overlap', 'proteins'])
def test_ambiguous_sparse_models_fail_clearly(tmp_path, damage):
    from lifton.reference_models import normalize_sparse_coding
    lines = [row('CDS', 10, 21, 'ID=c', phase='0')]
    if damage == 'dangling':
        lines = [row('CDS', 10, 21, 'ID=c;Parent=missing', phase='0')]
    elif damage == 'mixed':
        lines = [row('gene', 1, 60, 'ID=g'), row('mRNA', 1, 30, 'ID=t;Parent=g'),
                 row('CDS', 40, 51, 'ID=c;Parent=g', phase='0')]
    elif damage == 'proteins':
        lines = [row('gene', 1, 60, 'ID=g'), row('CDS', 10, 21, 'ID=c1;Parent=g;protein_id=p1', phase='0'),
                 row('CDS', 40, 51, 'ID=c2;Parent=g;protein_id=p2', phase='0')]
    else:
        lines += [row('CDS', 20 if damage == 'overlap' else 40, 51, 'ID=c',
                      '-' if damage == 'strand' else '+', 'chr2' if damage == 'seqid' else 'chr1', '0')]
    ann = reference(tmp_path, lines)
    with pytest.raises(LiftOnInputError):
        normalize_sparse_coding(ann, tmp_path / 'normalized')


def test_generated_ids_avoid_real_ids_and_are_deterministic(tmp_path):
    from lifton.reference_models import normalize_sparse_coding
    ann = reference(tmp_path, [row('gene', 100, 120, 'ID=c__lifton_gene'),
                               row('CDS', 10, 21, 'ID=c', phase='0')])
    first = normalize_sparse_coding(ann, tmp_path / 'one')
    second = normalize_sparse_coding(ann, tmp_path / 'two')
    m1 = json.loads(Path(first['mapping']).read_text())['models']
    m2 = json.loads(Path(second['mapping']).read_text())['models']
    assert m1 == m2 and m1[0]['gene_id'] != 'c__lifton_gene'
    assert Path(first['annotation']).read_bytes() == Path(second['annotation']).read_bytes()


def test_supplied_sequences_are_aliased_and_conflicts_rejected(tmp_path):
    from lifton.reference_models import normalize_sparse_coding, alias_fasta
    ann = reference(tmp_path, [row('CDS', 10, 21, 'ID=c;protein_id=p', phase='0')])
    result = normalize_sparse_coding(ann, tmp_path / 'normalized')
    model = json.loads(Path(result['mapping']).read_text())['models'][0]
    supplied = tmp_path / 'proteins.fa'
    supplied.write_text('>p\nMAAA\n>unrelated\nMVVV\n')
    output = alias_fasta(supplied, result['mapping'], tmp_path / 'aliased.fa')
    assert Path(output).read_text() == f'>{model["transcript_id"]}\nMAAA\n>unrelated\nMVVV\n'
    supplied.write_text('>c\nMAAA\n>p\nMVVV\n')
    with pytest.raises(LiftOnInputError, match='[Aa]mbiguous|[Cc]onflict'):
        alias_fasta(supplied, result['mapping'], tmp_path / 'conflict.fa')


def test_dangling_secondary_parent_is_rejected_in_an_ordinary_hierarchy(tmp_path):
    """Removing the dangling-parent validation must make this test fail."""
    from lifton.reference_models import normalize_sparse_coding
    ann = reference(tmp_path, [row('gene', 1, 30, 'ID=g'), row('mRNA', 1, 30, 'ID=t;Parent=g'),
                               row('CDS', 1, 30, 'ID=c;Parent=t,missing', phase='0')])
    with pytest.raises(LiftOnInputError, match='missing'):
        normalize_sparse_coding(ann, tmp_path / 'normalized')


@pytest.mark.parametrize('sparse_first', [False, True])
def test_sparse_alias_collision_with_ordinary_transcript_is_ambiguous(tmp_path, sparse_first):
    """Dropping ordinary transcript IDs from the alias index must make this test fail."""
    from lifton.reference_models import alias_fasta, normalize_sparse_coding
    ordinary = [row('gene', 100, 130, 'ID=g'), row('mRNA', 100, 130, 'ID=c;Parent=g'),
                row('exon', 100, 130, 'ID=e;Parent=c'),
                row('CDS', 100, 129, 'ID=ordinary_cds;Parent=c', phase='0')]
    sparse = [row('CDS', 10, 39, 'ID=c;protein_id=p', phase='0')]
    ann = reference(tmp_path, sparse + ordinary if sparse_first else ordinary + sparse)
    result = normalize_sparse_coding(ann, tmp_path / 'normalized')
    supplied = tmp_path / 'transcripts.fa'
    supplied.write_text('>c\nATGAAATAA\n')
    with pytest.raises(LiftOnInputError, match='Ambiguous'):
        alias_fasta(supplied, result['mapping'], tmp_path / 'aliased.fa')


def test_mixed_ordinary_noncoding_and_sparse_models_preserve_ordinary_rows(tmp_path):
    """Consuming unrelated ordinary/noncoding rows must make this test fail."""
    from lifton.reference_models import normalize_sparse_coding
    lines = [row('gene', 100, 130, 'ID=g;gene_biotype=lncRNA'),
             row('lnc_RNA', 100, 130, 'ID=t;Parent=g'), row('exon', 100, 130, 'ID=e;Parent=t'),
             row('CDS', 10, 39, 'ID=c;protein_id=p', phase='2')]
    ann = reference(tmp_path, lines)
    result = normalize_sparse_coding(ann, tmp_path / 'normalized')
    text = Path(result['annotation']).read_text()
    for original in lines[:3]:
        assert original.rstrip() in text
    model = json.loads(Path(result['mapping']).read_text())['models'][0]
    normalized = Annotation(result['annotation'], False, False).db_connection
    cds = list(normalized.children(model['transcript_id'], featuretype='CDS'))
    assert [(item.id, item.start, item.end, item.frame) for item in cds] == [('c', 10, 39, '2')]


def test_normalization_is_identical_across_database_backends(tmp_path):
    """Backend-dependent anchor iteration must make this test fail."""
    from lifton.reference_models import normalize_sparse_coding
    source = tmp_path / 'ref.gff3'
    source.write_text('##gff-version 3\n' + ''.join([
        row('CDS', 80, 109, 'ID=z;protein_id=pz', phase='1'),
        row('gene', 200, 260, 'ID=g'), row('CDS', 210, 239, 'ID=a;Parent=g', phase='0'),
        row('CDS', 20, 49, 'ID=b;protein_id=pb', phase='2'),
    ]))
    sqlite = normalize_sparse_coding(Annotation(source, False, False, backend='gffutils'), tmp_path / 'sqlite')
    duckdb = normalize_sparse_coding(Annotation(source, False, False, backend='gffbase'), tmp_path / 'duckdb')
    sqlite_map = json.loads(Path(sqlite['mapping']).read_text())
    duckdb_map = json.loads(Path(duckdb['mapping']).read_text())
    assert sqlite_map['models'] == duckdb_map['models']
    assert sqlite_map['aliases'] == duckdb_map['aliases']
    assert Path(sqlite['annotation']).read_bytes() == Path(duckdb['annotation']).read_bytes()


@pytest.mark.parametrize('direct_gene', [False, True])
@pytest.mark.parametrize('options', [[], ['--stream', '--inmemory-liftoff']])
def test_native_sparse_reference_lifts_complete_model(tmp_path, direct_gene, options):
    import os
    import random
    import shutil
    import subprocess
    import sys
    if any(not shutil.which(tool) for tool in ('minimap2', 'miniprot')):
        pytest.skip('requires native minimap2 and miniprot')
    rng = random.Random(14092026)
    coding = 'ATG' + ''.join(rng.choice(['GCT', 'TTC', 'GAA', 'CAA', 'TGG', 'AAC', 'GGT', 'ATT', 'CTG', 'CCA'])
                             for _ in range(199)) + 'TAA'
    sequence = ''.join(rng.choice('ACGT') for _ in range(1000)) + coding
    sequence += ''.join(rng.choice('ACGT') for _ in range(5000 - len(sequence)))
    fasta = tmp_path / 'genome.fa'
    fasta.write_text('>chr1\n' + sequence + '\n')
    lines = [row('gene', 1001, 1603, 'ID=g;gene_biotype=protein_coding')] if direct_gene else []
    lines += [row('CDS', 1001, 1603, 'ID=c;protein_id=p' + (';Parent=g' if direct_gene else ''), phase='0')]
    ref = tmp_path / 'input.gff3'
    ref.write_text('##gff-version 3\n' + ''.join(lines))
    out = tmp_path / 'output.gff3'
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]),
               PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0', LIFTON_MINIPROT_THREADS='1')
    command = [sys.executable, '-c', 'from lifton.lifton import main; main()', str(fasta), str(fasta),
               '-g', str(ref), '-o', str(out), '-t', '1', *options]
    process = subprocess.run(command, env=env, cwd=tmp_path, capture_output=True, text=True, timeout=180)
    assert process.returncode == 0, process.stderr
    manifest = json.loads((tmp_path / 'lifton_output/run_manifest.json').read_text())
    assert manifest['run']['status'] == 'success'
    normal = manifest['input_statistics']['reference_model_normalization']
    assert normal['models'] == 1 and normal['mapping_fingerprint']['sha256']
    mapping = json.loads(Path(normal['mapping']).read_text())['models'][0]
    rows = [line.split('\t') for line in out.read_text().splitlines() if line and not line.startswith('#')]
    assert [r[2] for r in rows] == ['gene', 'mRNA', 'exon', 'CDS']
    cds = rows[-1]
    assert cds[3:5] == ['1001', '1603'] and cds[6:8] == ['+', '0']
    assert 'ID=c;' in cds[8] and f'Parent={mapping["transcript_id"]}' in cds[8]


def _reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def _native_multisegment_case(tmp_path, strand, supplied, extra_options=()):
    import os
    import random
    import shutil
    import subprocess
    import sys
    if any(not shutil.which(tool) for tool in ('minimap2', 'miniprot')):
        pytest.skip('requires native minimap2 and miniprot')
    codons = ('GCT', 'TTC', 'GAA', 'CAA', 'TGG', 'AAC', 'GGT', 'ATT', 'CTG', 'CCA')
    amino_acids = {'GCT': 'A', 'TTC': 'F', 'GAA': 'E', 'CAA': 'Q', 'TGG': 'W',
                   'AAC': 'N', 'GGT': 'G', 'ATT': 'I', 'CTG': 'L', 'CCA': 'P'}
    body = [codons[index % len(codons)] for index in range(199)]
    transcript = 'ATG' + ''.join(body) + 'TAA'
    protein = 'M' + ''.join(amino_acids[codon] for codon in body) + '*'
    rng = random.Random(14092026)
    sequence = list(''.join(rng.choice('ACGT') for _ in range(5000)))
    chunks = (transcript[:300], transcript[300:])
    if strand == '-':
        chunks = (_reverse_complement(transcript[303:]), _reverse_complement(transcript[:303]))
    coordinates = ((1001, 1300), (1401, 1703))
    for (start, end), chunk in zip(coordinates, chunks):
        sequence[start - 1:end] = chunk
    fasta = tmp_path / 'genome.fa'
    fasta.write_text('>chr1\n' + ''.join(sequence) + '\n')
    ref = tmp_path / 'input.gff3'
    ref.write_text('##gff-version 3\n' + ''.join(
        row('CDS', start, end, 'ID=c;protein_id=p', strand, phase='0')
        for start, end in coordinates))
    arguments = []
    if supplied in ('protein', 'both'):
        supplied_proteins = tmp_path / 'supplied-proteins.fa'
        supplied_proteins.write_text(f'>p\n{protein}\n')
        arguments += ['-P', str(supplied_proteins)]
    if supplied in ('transcript', 'both'):
        supplied_transcripts = tmp_path / 'supplied-transcripts.fa'
        supplied_transcripts.write_text(f'>c\n{transcript}\n')
        arguments += ['-T', str(supplied_transcripts)]
    arguments += list(extra_options)
    out = tmp_path / 'output.gff3'
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]),
               PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0',
               LIFTON_MINIPROT_THREADS='1')
    command = [sys.executable, '-c', 'from lifton.lifton import main; main()', str(fasta), str(fasta),
               '-g', str(ref), '-o', str(out), '-t', '1', *arguments]
    process = subprocess.run(command, env=env, cwd=tmp_path, capture_output=True, text=True, timeout=180)
    return process, out, transcript, protein, coordinates


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('supplied', ['protein', 'transcript', 'both'])
def test_native_sparse_aliases_preserve_supplied_and_generated_sequence_truth(tmp_path, strand, supplied):
    """Breaking either -P/-T alias path or strand-aware extraction must make this test fail."""
    process, out, transcript, protein, coordinates = _native_multisegment_case(tmp_path, strand, supplied)
    assert process.returncode == 0, process.stderr
    artifact_dir = tmp_path / 'lifton_output/intermediate_files'
    manifest = json.loads((tmp_path / 'lifton_output/run_manifest.json').read_text())
    normalization = manifest['input_statistics']['reference_model_normalization']
    model = json.loads(Path(normalization['mapping']).read_text())['models'][0]
    transcript_path = (artifact_dir / 'reference_models/supplied_transcripts.fa'
                       if supplied in ('transcript', 'both') else artifact_dir / 'transcripts.fa')
    protein_path = (artifact_dir / 'reference_models/supplied_proteins.fa'
                    if supplied in ('protein', 'both') else artifact_dir / 'proteins.fa')
    assert transcript_path.read_text() == f'>{model["transcript_id"]}\n{transcript}\n'
    assert protein_path.read_text() == f'>{model["transcript_id"]}\n{protein}\n'
    rows = [line.split('\t') for line in out.read_text().splitlines() if line and not line.startswith('#')]
    cds = sorted((int(columns[3]), int(columns[4]), columns[6], columns[7])
                 for columns in rows if columns[2] == 'CDS')
    assert cds == [(start, end, strand, '0') for start, end in coordinates]


@pytest.mark.parametrize('extra_options', [('--gene-only',), ('-f', 'FEATURE_FILE')])
def test_native_sparse_reference_respects_feature_selection_options(tmp_path, extra_options):
    """Skipping flat normalization under --gene-only or breaking explicit -f CDS must make this test fail."""
    if 'FEATURE_FILE' in extra_options:
        feature_file = tmp_path / 'features.txt'
        feature_file.write_text('CDS\n')
        extra_options = tuple(str(feature_file) if value == 'FEATURE_FILE' else value for value in extra_options)
    process, out, _transcript, _protein, coordinates = _native_multisegment_case(
        tmp_path, '+', None, extra_options)
    assert process.returncode == 0, process.stderr
    rows = [line.split('\t') for line in out.read_text().splitlines() if line and not line.startswith('#')]
    cds = sorted((int(columns[3]), int(columns[4])) for columns in rows if columns[2] == 'CDS')
    assert cds == list(coordinates)
