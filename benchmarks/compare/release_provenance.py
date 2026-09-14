"""Content-addressed evidence for release runs; legacy files are never receipts."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from lifton.run_manifest import atomic_write_json


SCHEMA_VERSION = 1


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(path):
    path = Path(path).resolve(strict=True)
    before = path.stat()
    digest = sha256(path)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ino) != (
            after.st_size, after.st_mtime_ns, after.st_ino):
        raise RuntimeError(f'Input changed while hashing: {path}')
    return {'path': str(path), 'size': after.st_size, 'sha256': digest}


def snapshot(root):
    root = Path(root).resolve(strict=True)
    def git(*args):
        return subprocess.run(['git', '-C', str(root), *args], check=True,
                              capture_output=True, text=True).stdout.strip()
    if git('status', '--porcelain', '--untracked-files=normal', '--',
           'lifton', 'setup.py', 'pyproject.toml'):
        raise RuntimeError(f'Release source must be an immutable clean snapshot: {root}')
    return {'root': str(root), 'commit': git('rev-parse', 'HEAD')}


def runtime(python, root, cwd, env):
    code = (
        'import json,sys,lifton; from lifton.run_manifest import collect_dependency_versions; '
        'print(json.dumps(dict(executable=sys.executable,python=sys.version,'
        'module=lifton.__file__,version=lifton.__version__,dependencies=collect_dependency_versions())))'
    )
    result = subprocess.run([python, '-c', code], cwd=cwd, env=env, check=True,
                            capture_output=True, text=True, timeout=60)
    evidence = json.loads(result.stdout)
    if not Path(evidence['module']).resolve().is_relative_to(Path(root).resolve() / 'lifton'):
        raise RuntimeError(f'Wrong imported LiftOn: {evidence["module"]}; expected {root}')
    tools = {}
    for name in ('minimap2', 'miniprot'):
        executable = shutil.which(name, path=env.get('PATH'))
        if executable is None:
            raise RuntimeError(f'Missing required tool: {name}')
        probe = subprocess.run([executable, '--version'], capture_output=True, text=True,
                               check=True, timeout=30, env=env)
        tools[name] = {**fingerprint(executable), 'version': probe.stdout.strip() or probe.stderr.strip()}
    evidence['tools'] = tools
    return evidence


def run_evidence(*, python, root, inputs, argv, cwd, env):
    return {'schema_version': SCHEMA_VERSION, 'source': snapshot(root),
            'runtime': runtime(python, root, cwd, env),
            'inputs': {k: fingerprint(v) for k, v in sorted(inputs.items())},
            'argv': argv, 'environment': dict(sorted(env.items()))}


def read_receipt(path, expected, output, manifest):
    try:
        receipt = json.loads(Path(path).read_text())
        if receipt.get('schema_version') != SCHEMA_VERSION or receipt.get('evidence') != expected:
            return None
        if receipt.get('exit_code') != 0 or receipt.get('manifest_status') != 'success':
            return None
        if receipt.get('output') != fingerprint(output) or receipt.get('manifest') != fingerprint(manifest):
            return None
        if not receipt.get('validation', {}).get('complete'):
            return None
        return receipt
    except (OSError, ValueError, TypeError, KeyError):
        return None


def write_receipt(path, *, expected, output, manifest, profile, validation):
    document = json.loads(Path(manifest).read_text())
    if profile.exit_code != 0 or document.get('run', {}).get('status') != 'success':
        raise RuntimeError('Incomplete/partial run cannot receive a completion receipt')
    if snapshot(expected['source']['root']) != expected['source']:
        raise RuntimeError('Source snapshot changed during the run')
    for value in expected['inputs'].values():
        if fingerprint(value['path']) != value:
            raise RuntimeError(f'Input changed during the run: {value["path"]}')
    if not validation.get('complete'):
        raise RuntimeError('Validator did not finish; refusing completion receipt')
    receipt = {'schema_version': SCHEMA_VERSION, 'evidence': expected,
               'exit_code': profile.exit_code, 'manifest_status': 'success',
               'output': fingerprint(output), 'manifest': fingerprint(manifest),
               'wall_seconds': profile.wall_clock_seconds, 'peak_rss_mb': profile.peak_rss_mb,
               'validation': validation}
    atomic_write_json(path, receipt)
    return receipt


def isolated_env(root):
    # Explicit scientific configuration, independent of the interactive shell's
    # LIFTON_* overrides. A caller can add a documented experimental override.
    return {'PATH': os.environ.get('PATH', os.defpath), 'PYTHONNOUSERSITE': '1',
            'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1',
            'PYTHONPATH': str(Path(root).resolve())}
