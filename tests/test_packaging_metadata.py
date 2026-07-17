"""Release metadata and development-environment consistency checks."""

from __future__ import annotations

import ast
from pathlib import Path

from packaging.requirements import Requirement
from packaging.version import Version
import setuptools


ROOT = Path(__file__).resolve().parents[1]


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
