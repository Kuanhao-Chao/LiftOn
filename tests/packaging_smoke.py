"""Check an installed artifact outside source; optionally execute fresh lifts.

Run with the installed interpreter from a scratch working directory. This file
is deliberately not a collected unit test or part of the distributed package.
"""
import argparse
from collections import Counter
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(command, directory, environment, expected=0):
    process = subprocess.run(command, cwd=directory, env=environment, capture_output=True,
                             text=True, timeout=180)
    if (expected == 0 and process.returncode != 0) or (expected != 0 and process.returncode == 0):
        raise AssertionError(f'{command}: unexpected exit {process.returncode}\n{process.stderr}')
    return {'command': command, 'returncode': process.returncode,
            'stdout': process.stdout, 'stderr': process.stderr}


def constructed_reference(directory):
    """Two independently specified coding proteins, split across both strands."""
    rng = random.Random(780013)
    alphabet = [('GCT', 'A'), ('TTC', 'F'), ('GAA', 'E'), ('CAA', 'Q'), ('TGG', 'W'), ('AAC', 'N')]
    sequences, lines, expected = [], ['##gff-version 3'], {}
    for number, strand in enumerate(('+', '-'), 1):
        middle = [rng.choice(alphabet) for _ in range(199)]
        coding = 'ATG' + ''.join(codon for codon, _ in middle) + 'TAA'
        expected[f't{number}'] = 'M' + ''.join(residue for _, residue in middle) + '*'
        intron = 'GT' + ''.join(rng.choices('ACGT', k=96)) + 'AG'
        sequence = coding[:301] + intron + coding[301:]
        blocks = [(1, 301), (402, len(sequence))]
        if strand == '-':
            sequence = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            blocks = [(len(sequence) - end + 1, len(sequence) - start + 1) for start, end in blocks]
        sequence = ''.join(rng.choices('ACGT', k=1000)) + sequence + ''.join(rng.choices('ACGT', k=1000))
        blocks = [(start + 1000, end + 1000) for start, end in blocks]
        chrom, gene, trans = f'chr{number}', f'g{number}', f't{number}'
        sequences.append(f'>{chrom}\n{sequence}\n')

        def row(kind, start, end, phase, attrs):
            return f'{chrom}\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}'

        start, end = min(s for s, _ in blocks), max(e for _, e in blocks)
        lines += [row('gene', start, end, '.', f'ID={gene};gene_biotype=protein_coding'),
                  row('mRNA', start, end, '.', f'ID={trans};Parent={gene}')]
        for index, (start, end) in enumerate(blocks):
            lines += [row('exon', start, end, '.', f'ID=e{number}_{index};Parent={trans}'),
                      row('CDS', start, end, str(0 if index == 0 else 2),
                          f'ID=c{number}_{index};Parent={trans}')]
    fasta, annotation = directory / 'reference.fa', directory / 'reference.gff3'
    fasta.write_text(''.join(sequences))
    annotation.write_text('\n'.join(lines) + '\n')
    return fasta, annotation, expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--expected-version', required=True)
    parser.add_argument('--expect-mappy', choices=('present', 'absent'), required=True)
    parser.add_argument('--native-lift', action='store_true')
    parser.add_argument('--compiler-free', action='store_true')
    args = parser.parse_args()
    directory = args.work_dir.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    evidence = {'status': 'running', 'interpreter': sys.executable, 'python': sys.version,
                'checks': [], 'native_runs': []}
    try:
        import lifton
        actual_path = Path(lifton.__file__).resolve()
        assert actual_path.is_relative_to(Path(sys.prefix).resolve()), actual_path
        assert lifton.__version__.lstrip('v') == args.expected_version
        assert importlib.metadata.version('lifton') == args.expected_version
        available = importlib.util.find_spec('mappy') is not None
        assert available == (args.expect_mappy == 'present')
        if available:
            import mappy  # noqa: F401
        evidence.update(package_path=str(actual_path), package_sha256=sha256(actual_path),
                        mappy_available=available,
                        packages={dist.metadata['Name']: dist.version for dist in importlib.metadata.distributions()})
        environment = dict(os.environ, PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1',
                           LIFTON_MINIPROT_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        environment.pop('PYTHONPATH', None)
        if args.compiler_free:
            compilers = {name: shutil.which(name) for name in ('gcc', 'cc', 'clang', 'g++', 'c++')}
            assert not any(compilers.values()), compilers
            evidence['compilers'] = compilers
        cli = [sys.executable, '-c', 'from lifton.lifton import main; main()']
        validator = [sys.executable, '-c', 'from lifton.gff3_validator import _main; raise SystemExit(_main())']
        for command in (cli + ['-V'], cli + ['-h'], validator + ['-h'],
                        [sys.executable, '-m', 'pip', 'check']):
            evidence['checks'].append(run(command, directory, environment))
        fasta, annotation, proteins = constructed_reference(directory)
        evidence['inputs'] = {str(path): sha256(path) for path in (fasta, annotation)}
        # Hide executables without manufacturing a different Python environment.
        absent_tools = dict(environment, PATH=str(directory / 'empty-bin'))
        failed = run(cli + [str(fasta), str(fasta), '-g', str(annotation)],
                     directory, absent_tools, expected=1)
        assert 'miniprot is not installed' in failed['stderr']
        assert 'conda install' in failed['stderr'] and 'pip does not install' in failed['stderr']
        assert not (directory / 'lifton_output').exists()
        evidence['checks'].append(failed)
        if args.native_lift:
            from Bio.Seq import Seq
            import gffutils
            from pyfaidx import Fasta
            for tool in ('minimap2', 'miniprot'):
                binary = shutil.which(tool)
                assert binary, tool
                evidence['checks'].append(run([binary, '--version'], directory, environment))
                evidence.setdefault('tools', {})[tool] = {'path': binary, 'sha256': sha256(binary)}
            for mode, flags in (('standard', []), ('stream-inmemory', ['--stream', '--inmemory-liftoff'])):
                output, intermediates = directory / f'{mode}.gff3', directory / mode
                result = run(cli + [str(fasta), str(fasta), '-g', str(annotation), '-o', str(output),
                                    '-dir', str(intermediates), '-t', '1', '--strict-completeness'] + flags,
                             directory, environment)
                manifest = json.loads((intermediates / 'run_manifest.json').read_text())
                assert manifest['run']['status'] == 'success', result['stderr']
                database = gffutils.create_db(str(output), dbfn=':memory:', force=True,
                                             merge_strategy='create_unique', disable_infer_genes=True,
                                             disable_infer_transcripts=True)
                assert Counter(f.featuretype for f in database.all_features()) == {
                    'gene': 2, 'mRNA': 2, 'exon': 4, 'CDS': 4}
                with Fasta(str(fasta)) as genome:
                    for trans, expected_protein in proteins.items():
                        feature = database[trans]
                        cds = list(database.children(feature, featuretype='CDS', order_by='start'))
                        bases = ''.join(str(genome[c.seqid][c.start - 1:c.end]) for c in cds)
                        if feature.strand == '-':
                            bases = str(Seq(bases).reverse_complement())
                        assert str(Seq(bases).translate()) == expected_protein
                        assert feature.attributes['protein_identity'] == ['1.000']
                evidence['checks'].append(run(validator + [str(output)], directory, environment))
                evidence['native_runs'].append(dict(result, output_sha256=sha256(output), manifest=str(intermediates)))
            assert evidence['native_runs'][0]['output_sha256'] == evidence['native_runs'][1]['output_sha256']
            assert evidence['inputs'] == {str(path): sha256(path) for path in (fasta, annotation)}
        evidence['status'] = 'success'
    except Exception as exc:
        evidence.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        (directory / 'evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')


if __name__ == '__main__':
    main()
