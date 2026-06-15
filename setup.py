import setuptools
from pathlib import Path

this_directory = Path(__file__).resolve().parent
long_description = (this_directory / "./README.md").read_text()
setuptools.setup(
	name="lifton",
	version="1.0.8",
	author="Kuan-Hao Chao",
	author_email="kh.chao@cs.jhu.edu",
	description="Combining DNA and protein alignments to improve genome annotation with LiftOn",
	url="https://github.com/Kuanhao-Chao/Lifton",
	install_requires=['numpy>= 1.22.0', "biopython>=1.76", "cigar>=0.1.3", "parasail>=1.2.4", 'intervaltree>=3.1.0', 'interlap>=0.2.6', 'networkx>=3.3', 'pyfaidx>=0.5.8', 'pysam>=0.19.1', 'gffutils>=0.10.1', 'ujson>=3.2.0', 'pytest>=7.0.0',
	# Phase 6.1 (vendored gffbase) runtime requirements:
	'duckdb>=1.0', 'pyarrow>=14',
	# Phase 16 Tier 5: mappy unlocks the in-process minimap2 path under
	# `--native` (lifton/liftoff/native_align.py). Without it, --native
	# silently degrades to the subprocess path with a stderr warning
	# per call. Declaring it here makes the advertised flag actually
	# deliver native speed; the runtime still falls back gracefully if
	# the wheel is unavailable on a given platform.
	'mappy'],
	include_package_data=True,
	package_data={
		# Ship the vendored gffbase Rust extension and its bundled fixtures.
		'lifton.gffbase': ['_native*.so', '_native*.pyd', 'data/*.gff3', 'data/*.gtf', 'LICENSE'],
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