"""Private, recorded GTF conversions used by annotation intake."""
from pathlib import Path
import shutil
import subprocess
import tempfile

from lifton import annotation_cache, logger
from lifton.exceptions import LiftOnInputError


def _fingerprint(path):
    return {'path': str(Path(path).resolve()), **annotation_cache.source_fingerprint(path)}


def convert_gtf(source, *, directory=None, gffread=False, agat=None):
    """Return (GFF3 path, provenance), preserving every conversion attempt.

    Files live under a unique directory, never beside the source annotation.
    The caller retains these files for its database and native aligner inputs.
    """
    source = str(Path(source).resolve())
    original = _fingerprint(source)
    commands = []
    if gffread:
        commands.append(['gffread', '-E', '-F', '--keep-exon-attrs', '--keep-genes', source])
    if agat:
        commands.append([agat, '--gtf' if 'gtf2gff' in agat else '--gff', source])
    if not commands:
        return None, None
    if directory is not None:
        Path(directory).mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='gtf-conversion-', dir=directory))
    for index, command in enumerate(commands, 1):
        output = work / f'converted-{index}.gff3'
        command = command + ['-o', str(output)]
        executable = shutil.which(command[0])
        tool = _fingerprint(executable) if executable else None
        logger.log_info(f'Converting GTF to GFF3: {output}')
        try:
            with (work / f'attempt-{index}.stderr').open('w') as stderr:
                subprocess.run(command, stdout=subprocess.DEVNULL, stderr=stderr, text=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            with (work / f'attempt-{index}.stderr').open('a') as stderr:
                stderr.write(str(getattr(exc, 'stderr', '') or exc))
            logger.log_warning(f'{command[0]} conversion failed: {exc}')
            continue
        if not output.is_file() or not output.stat().st_size:
            logger.log_warning(f'{command[0]} conversion produced an empty file.')
            continue
        from lifton.annotation_validator import scan_annotation
        scan = scan_annotation(output)
        if scan.file_format != 'GFF3' or scan.errors or not dict(scan.feature_type_counts).get('gene'):
            logger.log_warning(f'{command[0]} did not produce GFF3 with a gene hierarchy.')
            continue
        if _fingerprint(source) != original or tool is None or _fingerprint(executable) != tool:
            raise LiftOnInputError('GTF input or conversion executable changed during conversion')
        try:
            version = subprocess.run([executable, '--version'], capture_output=True, text=True, timeout=30)
            tool['version'] = (version.stdout or version.stderr).strip()[:4096]
            tool['version_exit_code'] = version.returncode
        except (OSError, subprocess.TimeoutExpired) as exc:
            tool['version'] = None
            tool['version_error'] = str(exc)
        return str(output), {'schema_version': 1, 'input': original,
                             'output': _fingerprint(output), 'tool': tool,
                             'argv': command, 'stderr': str(work / f'attempt-{index}.stderr')}
    return None, None
