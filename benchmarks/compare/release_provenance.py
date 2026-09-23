"""Content-addressed evidence for release runs; legacy files are never receipts."""
from __future__ import annotations

import hashlib
import inspect
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from lifton.run_manifest import atomic_write_json


SCHEMA_VERSION = 2


def install_source_guard(root):
    """Resolve LiftOn imports only from the recorded source, including lazy imports.

    Self-contained because the identical bootstrap runs in released interpreters
    without requiring any changes to their source or installed package.
    """
    import importlib.machinery
    from pathlib import Path
    import sys
    root = Path(root).resolve(strict=True)
    package = root / 'lifton'

    def check_path(value):
        if not Path(value).resolve().is_relative_to(package):
            raise ImportError(f'LiftOn import outside recorded source: {value}; expected {package}')

    for name, module in tuple(sys.modules.items()):
        if name == 'lifton' or name.startswith('lifton.'):
            origin = getattr(module, '__file__', None)
            if origin:
                check_path(origin)
            for location in getattr(module, '__path__', ()):
                check_path(location)
    if any(getattr(finder, '_lifton_source_root', None) == str(root) for finder in sys.meta_path):
        return

    class SourceGuard:
        _lifton_source_root = str(root)

        @staticmethod
        def find_spec(fullname, path=None, target=None):
            if fullname != 'lifton' and not fullname.startswith('lifton.'):
                return None
            search = [str(root)] if fullname == 'lifton' else path
            if fullname != 'lifton':
                for location in search or ():
                    check_path(location)
            spec = importlib.machinery.PathFinder.find_spec(fullname, search, target)
            if spec is None:
                # Returning None would let an editable finder use another tree.
                raise ModuleNotFoundError(f'{fullname} is absent from recorded LiftOn source {root}', name=fullname)
            if spec.origin is not None:
                check_path(spec.origin)
            for location in spec.submodule_search_locations or ():
                check_path(location)
            return spec

    sys.meta_path.insert(0, SourceGuard)


def guarded_code(root, body):
    return inspect.getsource(install_source_guard) + f'\ninstall_source_guard({str(Path(root).resolve())!r})\n' + body


def imported_module_evidence():
    """Fingerprint the actual LiftOn modules loaded by a guarded runtime probe."""
    import hashlib
    from pathlib import Path
    import sys
    modules = {}
    for name, module in sorted(tuple(sys.modules.items())):
        if name != 'lifton' and not name.startswith('lifton.'):
            continue
        origin = getattr(module, '__file__', None)
        if origin:
            path = Path(origin).resolve(strict=True)
            modules[name] = {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
        else:
            modules[name] = {'locations': sorted(str(Path(p).resolve()) for p in getattr(module, '__path__', ()))}
    return modules


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


def verify_fingerprint(record, *, label="Artifact"):
    if not isinstance(record, dict) or not record.get("path") or fingerprint(record["path"]) != record:
        raise ValueError(f"{label} no longer matches recorded input/artifact fingerprint")
    return record


def dependency_evidence():
    # Self-contained so the exact same probe can run in the arm's interpreter.
    import hashlib
    import importlib.metadata
    from pathlib import Path
    from packaging.requirements import Requirement
    from packaging.utils import canonicalize_name
    packages = {}
    names = ("numpy", "biopython", "parasail", "intervaltree", "interlap", "networkx",
             "pyfaidx", "pysam", "gffutils", "ujson", "duckdb", "pyarrow")
    # `mappy` is an OPTIONAL extra (setup.py `extras_require['native']`) as of
    # v1.0.13, and the CI build job asserts it is ABSENT -- so its absence is a
    # supported configuration, not a broken environment. Listing it beside the
    # runtime requirements made every release-evidence test fail on a correctly
    # provisioned machine. Its presence is still worth recording, because it is
    # what makes the native path available, so it is probed and the answer
    # written down either way.
    optional = {canonicalize_name("mappy")}
    pending = [(name, frozenset()) for name in names + tuple(optional)]
    visited = set()
    while pending:
        name, extras = pending.pop()
        name = canonicalize_name(name)
        if (name, extras) in visited:
            continue
        visited.add((name, extras))
        try:
            distribution = importlib.metadata.distribution(name)
        except importlib.metadata.PackageNotFoundError as error:
            if name in optional:
                packages.setdefault(name, {"installed": False})
                continue
            raise RuntimeError(f"Missing required dependency: {name}") from error
        for specification in distribution.requires or []:
            requirement = Requirement(specification)
            if requirement.marker is None or any(
                    requirement.marker.evaluate({"extra": extra}) for extra in {"", *extras}):
                pending.append((requirement.name, frozenset(requirement.extras)))
        if name in packages:
            continue
        if distribution.files is None:
            raise RuntimeError(f"Dependency has no file inventory: {name}")
        # A wheel install lists what it installed (RECORD). A legacy egg-info
        # install can instead yield its SOURCES.txt, which lists the SOURCE
        # tree -- setup.py, LICENSE, the egg-info itself -- most of which were
        # never installed; Python 3.10's importlib does so even when the
        # egg-info also kept installed-files.txt. A legacy interlap install
        # listed nine such files and failed every release-evidence test on an
        # intact 3.10 env. Only when the list IS the source list is an absent
        # file expected: it is skipped and the weaker inventory recorded.
        # Everywhere else a listed file must exist.
        read_text = getattr(distribution, "read_text", None)
        sources = read_text("SOURCES.txt") if callable(read_text) else None
        source_inventory = (
            sources is not None and read_text("RECORD") is None
            and {str(item) for item in distribution.files}
            == {line for line in sources.splitlines() if line.strip()})
        digest = hashlib.sha256()
        count = 0
        for relative in sorted(distribution.files, key=str):
            if str(relative).endswith((".pyc", ".pyo")):
                continue
            path = Path(distribution.locate_file(relative))
            if not path.exists():
                if source_inventory:
                    continue
                raise RuntimeError(
                    f"Dependency {name} lists {relative} as installed, but it "
                    f"is missing: {path}")
            before = path.stat()
            file_hash = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    file_hash.update(chunk)
            after = path.stat()
            if (before.st_size, before.st_mtime_ns, before.st_ino) != (
                    after.st_size, after.st_mtime_ns, after.st_ino):
                raise RuntimeError(f"Dependency changed while hashing: {path}")
            digest.update(str(relative).encode() + b"\0" + file_hash.digest())
            count += 1
        if not count:
            raise RuntimeError(f"Dependency has an empty file inventory: {name}")
        packages[name] = {"version": distribution.version, "files": count, "sha256": digest.hexdigest()}
        if source_inventory:
            packages[name]["inventory"] = "egg-info SOURCES.txt"
    return packages


def evaluation_evidence():
    root = Path(__file__).resolve().parents[2]
    install_source_guard(root)
    import lifton
    package = Path(lifton.__file__).resolve().parent
    compare = Path(__file__).resolve().parent
    sources = {f"lifton/{p.relative_to(package)}": fingerprint(p) for p in sorted(package.rglob("*.py"))}
    sources.update({f"benchmarks/compare/{p.name}": fingerprint(p) for p in sorted(compare.glob("*.py"))})
    return {"sources": sources, "source_root": str(root), "dependencies": dependency_evidence(),
            "python": sys.version, "executable": fingerprint(sys.executable)}


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
        'import json,sys,lifton,lifton.lifton; '
        'print(json.dumps(dict(executable=sys.executable,python=sys.version,'
        'module=lifton.__file__,version=lifton.__version__,modules=imported_module_evidence(),'
        'dependencies=dependency_evidence())))'
    )
    code = guarded_code(root, inspect.getsource(dependency_evidence) + '\n'
                        + inspect.getsource(imported_module_evidence) + '\n' + code)
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
    evidence['interpreter'] = fingerprint(evidence['executable'])
    return evidence


def run_evidence(*, python, root, inputs, argv, cwd, env):
    return {'schema_version': SCHEMA_VERSION, 'source': snapshot(root),
            'runtime': runtime(python, root, cwd, env),
            'inputs': {k: fingerprint(v) for k, v in sorted(inputs.items())},
            'argv': argv, 'cwd': str(Path(cwd).resolve()), 'environment': dict(sorted(env.items()))}


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


def verify_run_evidence(expected):
    """Recheck live evidence before sealing a run or accepting its report."""
    if snapshot(expected['source']['root']) != expected['source']:
        raise RuntimeError('Source snapshot changed during the run')
    observed_runtime = runtime(expected['argv'][0], expected['source']['root'],
                               expected['cwd'], expected['environment'])
    if observed_runtime != expected['runtime']:
        raise RuntimeError('Runtime dependencies or native tools changed during the run')
    for value in expected['inputs'].values():
        if fingerprint(value['path']) != value:
            raise RuntimeError(f'Input changed during the run: {value["path"]}')


def write_receipt(path, *, expected, output, manifest, profile, validation):
    document = json.loads(Path(manifest).read_text())
    if profile.exit_code != 0 or document.get('run', {}).get('status') != 'success':
        raise RuntimeError('Incomplete/partial run cannot receive a completion receipt')
    verify_run_evidence(expected)
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
