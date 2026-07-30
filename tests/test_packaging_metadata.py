"""Release metadata and development-environment consistency checks."""

from __future__ import annotations

import ast
import re
import runpy
import tomllib
from pathlib import Path

from packaging.requirements import Requirement
from packaging.version import Version
import setuptools


ROOT = Path(__file__).resolve().parents[1]


def _package_version() -> tuple[str, str]:
    text = (ROOT / "lifton" / "__init__.py").read_text(encoding="utf-8")
    match = re.search(r"__version__\s*=\s*['\"](?P<version>v?[^'\"]+)['\"]", text)
    assert match is not None
    cli_version = match.group("version")
    return cli_version, cli_version.removeprefix("v")


def _setup_expression(name: str) -> ast.expr:
    tree = ast.parse((ROOT / "setup.py").read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not (
            isinstance(node.func, ast.Attribute)
            and node.func.attr == "setup"
        ):
            continue
        for keyword in node.keywords:
            if keyword.arg == name:
                return keyword.value
    raise AssertionError(f"setup.py does not define {name}")


def _setup_keyword(name: str):
    return ast.literal_eval(_setup_expression(name))


def _pip_environment_requirements() -> dict[str, Requirement]:
    requirements = {}
    in_pip_section = False
    for raw_line in (ROOT / "lifton.yml").read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if stripped == "- pip:":
            in_pip_section = True
            continue
        if not in_pip_section or not stripped.startswith("- "):
            continue
        requirement = Requirement(stripped[2:])
        requirements[requirement.name.lower()] = requirement
    return requirements


def test_current_release_metadata_matches_package_version():
    cli_version, package_version = _package_version()
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    docs_config = runpy.run_path(str(ROOT / "docs" / "source" / "conf.py"))
    badge = f"https://img.shields.io/badge/version-{cli_version}-blue"

    assert re.search(
        rf"(?m)^version:\s*{re.escape(package_version)}\s*$", citation,
    )
    assert re.search(
        r'(?m)^date-released:\s*["\']2026-07-30["\']\s*$', citation,
    )
    assert docs_config["release"] == package_version
    assert docs_config["version"] == package_version
    assert badge in (ROOT / "README.md").read_text(encoding="utf-8")
    assert badge in (ROOT / "docs" / "source" / "index.rst").read_text(
        encoding="utf-8",
    )
    for relative_path in (
        "docs/source/content/function_manual.rst",
        "docs/source/content/installation.rst",
    ):
        assert f"      {cli_version}\n" in (ROOT / relative_path).read_text(
            encoding="utf-8",
        )


def test_duckdb_release_exclusions_match_environment():
    runtime = {
        requirement.name.lower(): requirement
        for requirement in map(Requirement, _setup_keyword("install_requires"))
    }
    environment = _pip_environment_requirements()

    for requirements in (runtime, environment):
        duckdb = requirements["duckdb"]
        assert Version("1.5.2") in duckdb.specifier
        assert Version("1.5.3") not in duckdb.specifier
        assert Version("1.5.4") not in duckdb.specifier
        assert Version("1.5.5") in duckdb.specifier
        assert Version("14") in requirements["pyarrow"].specifier


def test_cigar_is_not_a_dependency():
    """`cigar` 0.1.3 cannot be installed in a clean environment.

    Its setup.py calls `ez_setup.use_setuptools()`, which downloads a setuptools
    egg from a URL that no longer serves one, so `pip install lifton` died with
    `tarfile.ReadError: not a gzip file` for anyone without a cached wheel. It
    was the only dependency that could not build, LiftOn used exactly one call
    from it, and that call now lives in `coreutils.parse_cigar_items`.

    Re-adding it would silently break installation for every new user, so assert
    it stays out of both the runtime requirements and the pinned environment.
    """
    runtime = {
        Requirement(requirement).name.lower()
        for requirement in _setup_keyword("install_requires")
    }
    assert "cigar" not in runtime
    assert "cigar" not in _pip_environment_requirements()

    align_source = (ROOT / "lifton" / "align.py").read_text(encoding="utf-8")
    assert "from cigar import" not in align_source


def test_build_backend_requires_a_modern_setuptools_floor():
    """The former `setuptools<81` ceiling existed only for `cigar`.

    With that dependency gone the cap went with it, but the floor still matters:
    setup.py relies on `find_namespace_packages` and PEP 517 metadata.
    """
    build_system = tomllib.loads(
        (ROOT / "pyproject.toml").read_text(encoding="utf-8"),
    )["build-system"]
    setuptools_requirement = next(
        Requirement(requirement)
        for requirement in build_system["requires"]
        if Requirement(requirement).name.lower() == "setuptools"
    )

    assert Version("60.0.0") not in setuptools_requirement.specifier
    assert Version("61") in setuptools_requirement.specifier
    assert Version("81") in setuptools_requirement.specifier


def test_development_environment_is_relocatable_and_complete():
    environment_text = (ROOT / "lifton.yml").read_text(encoding="utf-8")
    environment = _pip_environment_requirements()

    assert "\nprefix:" not in f"\n{environment_text}"
    assert "lifton" not in environment
    assert {"build", "coverage", "flake8", "hypothesis", "pytest"} <= environment.keys()


def test_test_extra_contains_direct_test_dependencies():
    requirements = {
        Requirement(requirement).name.lower()
        for requirement in _setup_keyword("extras_require")["test"]
    }

    assert {
        "coverage",
        "flake8",
        "hypothesis",
        "packaging",
        "pytest",
    } <= requirements


def test_mappy_is_a_declared_runtime_dependency():
    runtime_names = {
        Requirement(requirement).name.lower()
        for requirement in _setup_keyword("install_requires")
    }
    assert "mappy" in runtime_names


def test_wheel_discovery_excludes_vendored_liftoff_tests():
    package_expression = _setup_expression("packages")
    assert isinstance(package_expression, ast.Call)
    arguments = {
        keyword.arg: ast.literal_eval(keyword.value)
        for keyword in package_expression.keywords
    }
    packages = setuptools.find_namespace_packages(
        where=str(ROOT),
        include=arguments["include"],
        exclude=arguments["exclude"],
    )

    assert "lifton.liftoff" in packages
    assert not any(
        package == "lifton.liftoff.tests"
        or package.startswith("lifton.liftoff.tests.")
        for package in packages
    )

    excluded_data = _setup_keyword("exclude_package_data")
    assert excluded_data["lifton.liftoff"] == ["tests/*", "tests/**/*"]
