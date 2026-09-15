"""Independent coding-sequence truth for initial phase and genetic codes."""
import copy

from Bio.Seq import Seq
import gffutils
from pyfaidx import Fasta
import pytest

from lifton import extract_sequence, lifton_class, orf_completion


def model(tmp_path, sequence, phase, strand='+', split=5):
    # Explicit 5'-to-3' transcript; a 10-base intron separates two CDS pieces.
    genomic = sequence[:split] + 'C' * 10 + sequence[split:]
    blocks = [(1, split), (split + 11, len(genomic))]
    if strand == '-':
        genomic = str(Seq(genomic).reverse_complement())
        blocks = [(len(genomic) - end + 1, len(genomic) - start + 1) for start, end in blocks]
    path = tmp_path / 'genome.fa'
    path.write_text('>chr1\n' + genomic + '\n')
    entry = gffutils.Feature(seqid='chr1', source='test', featuretype='mRNA', start=1,
                            end=len(genomic), strand=strand, attributes={'ID': ['t'], 'Parent': ['g']})
    trans = lifton_class.Lifton_TRANS('t', 'g', 'g', 0, entry, copy.deepcopy(entry.attributes))
    cds = []
    for n, (start, end) in enumerate(blocks):
        exon = gffutils.Feature(seqid='chr1', source='test', featuretype='exon', start=start, end=end,
                               strand=strand, attributes={'ID': [f'e{n}'], 'Parent': ['t']})
        trans.add_exon(exon)
        coding = copy.deepcopy(exon)
        coding.featuretype = 'CDS'
        coding.frame = str(phase if n == 0 else (phase - split) % 3)
        coding.attributes['ID'] = [f'c{n}']
        trans.add_cds(coding)
        cds.append(coding)
    return trans, cds, Fasta(str(path))


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('phase', [1, 2])
def test_initial_phase_is_trimmed_once_across_splice_junction(tmp_path, strand, phase):
    sequence = 'C' * phase + 'ATGGCTTAA'
    trans, cds, fasta = model(tmp_path, sequence, phase, strand)
    assert extract_sequence.get_protein_sequence(trans.entry, fasta, cds) == 'MA*'
    coding, transcript = trans.get_coding_trans_seq(fasta)
    assert coding == 'ATGGCTTAA' and transcript == sequence
    assert trans.translate_coding_seq(coding) == 'MA*'
    ordered = sorted((e.cds.entry for e in trans.exons), key=lambda e: e.start, reverse=strand == '-')
    assert [c.frame for c in ordered] == [str(phase), str((phase - 5) % 3)]
    assert trans.get_coding_trans_seq(fasta)[0] == coding  # repeated scoring cannot trim again
    assert trans.get_coding_seq(fasta, include_sequence=False)[2] == [5 - phase, len(sequence) - 5]


def test_phase_trim_does_not_translate_artificial_padding(tmp_path):
    trans, cds, fasta = model(tmp_path, 'CCAT', 2, split=3)
    assert extract_sequence.get_protein_sequence(trans.entry, fasta, cds) == ''
    assert trans.translate_coding_seq(trans.get_coding_trans_seq(fasta)[0]) == ''


def test_phase_aware_terminal_stop_extension(tmp_path):
    from tests.test_orf_completion import _transcript
    path = tmp_path / 'stop.fa'
    path.write_text('>chr1\nCATGGCTTAA\n')
    trans = _transcript('chr1', [(1, 7)], '+')
    trans.exons[0].cds.entry.frame = '1'
    assert orf_completion.complete_terminal_stop(trans, Fasta(str(path)))
    assert trans.exons[0].cds.entry.end == 10


@pytest.mark.parametrize(('phase', 'strand'), [(1, '+'), (2, '-')])
def test_native_partial_reference_keeps_phase_and_exact_protein(tmp_path, phase, strand):
    import json
    import os
    from pathlib import Path
    import random
    import shutil
    import subprocess
    import sys
    if any(not shutil.which(tool) for tool in ('minimap2', 'miniprot')):
        pytest.skip('requires native minimap2 and miniprot')
    rng = random.Random(20260915)
    residues = [('GCT', 'A'), ('TTC', 'F'), ('GAA', 'E'), ('CAA', 'Q'), ('TGG', 'W'), ('AAC', 'N')]
    middle = [rng.choice(residues) for _ in range(199)]
    sequence = 'C' * phase + 'ATG' + ''.join(codon for codon, _ in middle) + 'TAA'
    expected_protein = 'M' + ''.join(amino for _, amino in middle) + '*'
    trans, cds, fasta = model(tmp_path, sequence, phase, strand, split=phase + 301)
    annotation = tmp_path / 'reference.gff3'
    gene = copy.deepcopy(trans.entry)
    gene.featuretype = 'gene'
    gene.attributes = {'ID': ['g'], 'gene_biotype': ['protein_coding']}
    features = [gene, trans.entry] + [e.entry for e in trans.exons] + cds
    annotation.write_text('##gff-version 3\n' + '\n'.join(map(str, features)) + '\n')
    output = tmp_path / 'out.gff3'
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]), PYTHONNOUSERSITE='1',
               PYTHONDONTWRITEBYTECODE='1', LIFTON_MINIPROT_THREADS='1')
    process = subprocess.run([sys.executable, '-c', 'from lifton.lifton import main; main()',
                              str(tmp_path / 'genome.fa'), str(tmp_path / 'genome.fa'),
                              '-g', str(annotation), '-o', str(output), '-t', '1'],
                             cwd=tmp_path, env=env, capture_output=True, text=True, timeout=180)
    assert process.returncode == 0, process.stderr
    manifest = json.loads((tmp_path / 'lifton_output/run_manifest.json').read_text())
    assert manifest['run']['status'] == 'success', process.stderr
    extracted = Fasta(str(tmp_path / 'lifton_output/intermediate_files/proteins.fa'))
    assert str(extracted['t']) == expected_protein
    db = gffutils.create_db(str(output), dbfn=':memory:', force=True, merge_strategy='create_unique',
                           disable_infer_genes=True, disable_infer_transcripts=True)
    output_cds = list(db.features_of_type('CDS', order_by='start', reverse=strand == '-'))
    expected_cds = sorted(cds, key=lambda c: c.start, reverse=strand == '-')
    assert [(c.start, c.end, c.frame) for c in output_cds] == [(c.start, c.end, c.frame) for c in expected_cds]
    assert db['t'].attributes['protein_identity'] == ['1.000']
