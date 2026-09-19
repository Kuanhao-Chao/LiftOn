"""GTF conversion must preserve biology and keep derived files private."""
from pathlib import Path
import shutil

import pytest

from lifton.annotation import Annotation


@pytest.fixture
def gtf(tmp_path):
    source = tmp_path / 'reference.gtf'
    gene = 'gene_id "g1"; gene_name "coding gene"; gene_biotype "protein_coding";'
    transcript = gene + ' transcript_id "t1"; transcript_biotype "protein_coding";'
    noncoding = 'gene_id "g2"; gene_biotype "lncRNA";'
    source.write_text(
        f'chr1\tEnsembl\tgene\t11\t43\t.\t+\t.\t{gene}\n'
        f'chr1\tEnsembl\ttranscript\t11\t43\t.\t+\t.\t{transcript}\n'
        f'chr1\tEnsembl\texon\t11\t43\t.\t+\t.\t{transcript} exon_id "e1";\n'
        f'chr1\tEnsembl\tCDS\t11\t40\t.\t+\t0\t{transcript} protein_id "p1";\n'
        f'chr1\tEnsembl\tstop_codon\t41\t43\t.\t+\t0\t{transcript}\n'
        f'chr1\tEnsembl\tgene\t71\t100\t.\t-\t.\t{noncoding}\n'
        f'chr1\tEnsembl\ttranscript\t71\t100\t.\t-\t.\t{noncoding} transcript_id "t2";\n'
        f'chr1\tEnsembl\texon\t71\t100\t.\t-\t.\t{noncoding} transcript_id "t2"; exon_id "e2";\n')
    return source


@pytest.mark.skipif(not shutil.which('gffread'), reason='requires native gffread')
@pytest.mark.parametrize('inferred', [False, True])
@pytest.mark.parametrize('backend', ['gffutils', 'gffbase'])
def test_conversion_retains_genes_transcripts_and_child_attributes(gtf, inferred, backend):
    if inferred:
        gtf.write_text(''.join(line for line in gtf.read_text().splitlines(True)
                               if line.split('\t')[2] not in {'gene', 'transcript'}))
    original = gtf.read_bytes()
    ann = Annotation(gtf, False, False, backend=backend)
    db = ann.db_connection
    assert ann.backend == backend
    assert {f.id for f in db.features_of_type('gene')} == {'g1', 'g2'}
    assert db['t1'].attributes['Parent'] == ['g1']
    assert db['t2'].attributes['Parent'] == ['g2']
    assert list(db.children('t1', featuretype='exon'))[0].attributes['exon_id'] == ['e1']
    cds = list(db.children('t1', featuretype='CDS'))
    assert [(c.start, c.end, c.frame) for c in cds] == [(11, 43, '0')]
    assert cds[0].attributes['protein_id'] == ['p1']
    assert not list(db.children('t2', featuretype='CDS'))
    assert gtf.read_bytes() == original


@pytest.mark.skipif(not shutil.which('gffread'), reason='requires native gffread')
def test_conversion_uses_unique_private_paths_and_records_evidence(gtf, tmp_path):
    first = Annotation(gtf, False, False, conversion_dir=tmp_path / 'run')
    second = Annotation(gtf, False, False, conversion_dir=tmp_path / 'run')
    assert first.file_name != second.file_name
    assert Path(first.file_name).is_relative_to(tmp_path / 'run')
    assert not (gtf.parent / 'reference_converted.gff3').exists()
    evidence = first.conversion_provenance
    assert evidence['schema_version'] == 1
    assert evidence['input']['path'] == str(gtf)
    assert evidence['input']['sha256'] and evidence['output']['sha256']
    assert evidence['tool']['sha256'] and evidence['tool']['version']
    assert '--keep-genes' in evidence['argv'] and '--keep-exon-attrs' in evidence['argv']


@pytest.mark.parametrize(('backend', 'from_environment'), [
    ('gffutils', False),
    ('gffbase', False),
    (None, True),
])
def test_direct_gtf_opt_out_keeps_explicit_hierarchy(
        gtf, monkeypatch, backend, from_environment):
    monkeypatch.setattr(Annotation, '_convert_gtf_to_gff3',
                        lambda _: pytest.fail('conversion should be disabled'))
    if from_environment:
        monkeypatch.setenv('LIFTON_USE_GFFBASE', '1')
    ann = Annotation(gtf, False, False, auto_convert_gtf=False, backend=backend)
    assert ann.file_name == str(gtf)
    assert {f.id for f in ann.db_connection.features_of_type('gene')} == {'g1', 'g2'}
    assert [p.id for p in ann.db_connection.parents('t1', level=1)] == ['g1']
    assert ann.db_connection['t1'].attributes['Parent'] == ['g1']
    requested = 'gffbase' if backend == 'gffbase' or from_environment else 'gffutils'
    assert ann.requested_backend == requested
    assert ann.backend == 'gffutils'
    assert ann.backend_fallback_reason == (
        'direct_gtf_hierarchy_requires_gffutils' if requested == 'gffbase' else None)
    assert ann.conversion_provenance is None


@pytest.mark.skipif(not shutil.which('gffread'), reason='requires native gffread')
def test_ensembl_gene_biotype_drives_coding_classification(gtf, tmp_path):
    from types import SimpleNamespace
    from lifton import lifton_utils
    ann = Annotation(gtf, False, False)
    args = SimpleNamespace(annotation_database='Ensembl', evaluation_liftoff_chm13=False)
    features, _, _, _ = lifton_utils.get_ref_liffover_features(['gene'], ann, str(tmp_path), args)
    assert features['g1'].is_protein_coding
    assert features['g2'].is_non_coding


def test_failed_converter_cannot_reuse_stale_success(gtf, tmp_path, monkeypatch):
    import subprocess
    from types import SimpleNamespace
    from lifton.annotation_conversion import convert_gtf
    from lifton.exceptions import LiftOnInputError
    executable = tmp_path / 'fake-gffread'
    executable.write_text('fake executable')
    stale = gtf.with_name('reference_converted.gff3')
    stale.write_text('##gff-version 3\nchr1\ttest\tgene\t1\t99\t.\t+\t.\tID=old\n')
    monkeypatch.setattr(shutil, 'which', lambda _: str(executable))
    def no_output(*args, **kwargs):
        return SimpleNamespace(returncode=0, stdout='', stderr='')
    monkeypatch.setattr(subprocess, 'run', no_output)
    assert convert_gtf(gtf, directory=tmp_path / 'run', gffread=True) == (None, None)
    assert 'ID=old' in stale.read_text()
    def changed_input(command, **kwargs):
        Path(command[-1]).write_text(stale.read_text())
        gtf.write_text(gtf.read_text() + '# modified\n')
        return no_output()
    monkeypatch.setattr(subprocess, 'run', changed_input)
    with pytest.raises(LiftOnInputError, match='changed during conversion'):
        convert_gtf(gtf, directory=tmp_path / 'run', gffread=True)


@pytest.mark.skipif(any(not shutil.which(tool) for tool in ('gffread', 'minimap2', 'miniprot')),
                    reason='requires native gffread, minimap2 and miniprot')
@pytest.mark.parametrize(('options', 'request_gffbase'), [
    ([], False),
    (['--stream'], False),
    (['--no-auto-convert-gtf'], True),
    (['--strict-gff'], False),
])
def test_native_gtf_lift_recovers_coding_hierarchy(tmp_path, options, request_gffbase):
    import json
    import os
    import random
    import subprocess
    import sys
    rng = random.Random(14092026)
    codons = ['GCT', 'TTC', 'GAA', 'CAA', 'TGG', 'AAC', 'GGT', 'ATT', 'CTG', 'CCA']
    coding = 'ATG' + ''.join(rng.choice(codons) for _ in range(199)) + 'TAA'
    sequence = ''.join(rng.choice('ACGT') for _ in range(1000)) + coding
    sequence += ''.join(rng.choice('ACGT') for _ in range(5000 - len(sequence)))
    fasta = tmp_path / 'genome.fa'
    fasta.write_text('>chr1\n' + sequence + '\n')
    gtf = tmp_path / 'reference.gtf'
    attrs = 'gene_id "g1"; gene_biotype "protein_coding";'
    gtf.write_text(''.join(
        f'chr1\tEnsembl\t{kind}\t1001\t1603\t.\t+\t{phase}\t{attrs}{extra}\n'
        for kind, phase, extra in [('gene', '.', ''),
                                   ('transcript', '.', ' transcript_id "t1";'),
                                   ('exon', '.', ' transcript_id "t1"; exon_id "e1";'),
                                   ('CDS', '0', ' transcript_id "t1"; protein_id "p1";')]))
    output = tmp_path / 'out.gff3'
    env = dict(os.environ, PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0',
               PYTHONPATH=str(Path(__file__).resolve().parents[1]), LIFTON_MINIPROT_THREADS='1')
    if request_gffbase:
        env['LIFTON_USE_GFFBASE'] = '1'
    command = [sys.executable, '-c', 'from lifton.lifton import main; main()',
               str(fasta), str(fasta), '-g', str(gtf), '-o', str(output), '-ad', 'Ensembl', '-t', '1', *options]
    result = subprocess.run(command, cwd=tmp_path, env=env, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((tmp_path / 'lifton_output/run_manifest.json').read_text())
    assert manifest['run']['status'] == 'success', manifest.get('failures')
    rows = [line.split('\t') for line in output.read_text().splitlines() if line and not line.startswith('#')]
    assert [r[2] for r in rows].count('gene') == 1
    cds = [r for r in rows if r[2] == 'CDS']
    assert len(cds) == 1 and cds[0][3:5] == ['1001', '1603']
    assert 'Parent=t1' in cds[0][8]
    if request_gffbase:
        assert manifest['run']['backend']['reference_annotation'] == 'gffutils'
        assert manifest['run']['backend']['reference_annotation_requested'] == 'gffbase'
        assert manifest['run']['backend']['reference_annotation_fallback'] == \
            'direct_gtf_hierarchy_requires_gffutils'
    if not options or '--stream' in options:
        assert manifest['input_statistics']['reference_annotation_conversion']['tool']['sha256']


@pytest.mark.parametrize('backend', ['gffutils', 'gffbase'])
def test_direct_gtf_inference_preserves_original_transcript_ids(gtf, backend):
    gtf.write_text(''.join(line for line in gtf.read_text().splitlines(True)
                           if line.split('\t')[2] not in {'gene', 'transcript'}))
    ann = Annotation(gtf, True, True, auto_convert_gtf=False, backend=backend)
    assert {f.id for f in ann.db_connection.features_of_type('gene')} == {'g1', 'g2'}
    assert {f.id for f in ann.db_connection.features_of_type('transcript')} == {'t1', 't2'}
    assert ann.db_connection['t1'].attributes['Parent'] == ['g1']
    assert ann.backend == 'gffutils'
