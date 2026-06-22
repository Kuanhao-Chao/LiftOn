import re
import setuptools
from pathlib import Path

this_directory = Path(__file__).resolve().parent
long_description = (this_directory / "./README.md").read_text()

# Single source of truth for the version: parse lifton/__init__.py and strip a
# leading "v" so the PEP 440 metadata version can never drift from the CLI
# string that `lifton -V` prints (lifton/__init__.py:__version__).
_version_text = (this_directory / "lifton" / "__init__.py").read_text()
_version_match = re.search(
    r"__version__\s*=\s*['\"]v?(?P<v>[^'\"]+)['\"]", _version_text)
if _version_match is None:
    raise RuntimeError("Unable to find __version__ in lifton/__init__.py")
version = _version_match.group("v")

setuptools.setup(
	name="lifton",
	version=version,
	author="Kuan-Hao Chao",
	author_email="kh.chao@cs.jhu.edu",
	description="Combining DNA and protein alignments to improve genome annotation with LiftOn",
	url="https://github.com/Kuanhao-Chao/Lifton",
	license="GPL-3.0-or-later",
	keywords=["genome annotation", "liftover", "lift-over", "bioinformatics",
		"genomics", "miniprot", "minimap2", "gff3", "comparative genomics"],
	project_urls={
		"Documentation": "https://ccb.jhu.edu/lifton/",
		"Source": "https://github.com/Kuanhao-Chao/Lifton",
		"Bug Tracker": "https://github.com/Kuanhao-Chao/Lifton/issues",
		"Publication": "https://doi.org/10.1101/gr.279620.124",
	},
	classifiers=[
		"Development Status :: 5 - Production/Stable",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
		"Operating System :: POSIX :: Linux",
		"Operating System :: MacOS :: MacOS X",
		"Programming Language :: Python :: 3",
		"Programming Language :: Python :: 3.9",
		"Programming Language :: Python :: 3.10",
		"Programming Language :: Python :: 3.11",
		"Programming Language :: Python :: 3.12",
		"Topic :: Scientific/Engineering :: Bio-Informatics",
	],
	install_requires=['numpy>= 1.22.0', "biopython>=1.76", "cigar>=0.1.3", "parasail>=1.2.4", 'intervaltree>=3.1.0', 'interlap>=0.2.6', 'networkx>=3.3', 'pyfaidx>=0.5.8', 'pysam>=0.19.1', 'gffutils>=0.10.1', 'ujson>=3.2.0',
	# Phase 6.1 (vendored gffbase) runtime requirements:
	'duckdb>=1.0', 'pyarrow>=14',
	# Phase 16 Tier 5: mappy unlocks the in-process minimap2 path under
	# `--native` (lifton/liftoff/native_align.py). Without it, --native
	# silently degrades to the subprocess path with a stderr warning
	# per call. Declaring it here makes the advertised flag actually
	# deliver native speed; the runtime still falls back gracefully if
	# the wheel is unavailable on a given platform.
	'mappy'],
	# pytest/hypothesis/coverage are test-only; not pulled by a plain install.
	extras_require={
		'test': ['pytest>=7.0.0', 'hypothesis>=6.0', 'coverage>=6.0'],
	},
	include_package_data=True,
	package_data={
		# gffbase ships the pure-Python fallback parser + the Rust extension
		# SOURCE (no pre-built _native*.so in the tree -> the runtime falls back
		# to _pyfallback). Only declare files that actually exist.
		'lifton.gffbase': ['LICENSE'],
		'lifton.gffbase._pyfallback': ['*.py'],
		'lifton.gffbase._rust': ['Cargo.toml', 'Cargo.lock'],
		'lifton.gffbase._rust.src': ['*.rs'],
	},
	# Realistic floor (was '>=3.6'): the code uses PEP-585 builtin generics
	# (list[...]/dict[...]) which are native in 3.9, and PEP-604 X|None
	# annotations guarded by `from __future__ import annotations`. The dev env
	# and CI run 3.11 (lifton.yml). >=3.9 drops only EOL Pythons 3.6-3.8.
	python_requires='>=3.9',
	packages=setuptools.find_packages(),
	entry_points={'console_scripts': [
            'lifton = lifton.lifton:main',
            'gff3-validate = lifton.gff3_validator:_main',
        ]},
        long_description=long_description,
        long_description_content_type='text/markdown'
)
