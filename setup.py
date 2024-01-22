import setuptools
from pathlib import Path

this_directory = Path(__file__).resolve().parent
long_description = (this_directory / "./README.md").read_text()
setuptools.setup(
	name="lifton",
	version="0.0.2",
	author="Kuan-Hao Chao",
	author_email="kh.chao@cs.jhu.edu",
	description="A protein-coding gene annotation fixing tool",
	url="https://github.com/Kuanhao-Chao/Lifton",
	install_requires=['numpy', "biopython", "cigar>=0.1.3", "parasail>=1.2.1", 'intervaltree>=3.1.0', 'interlap>=0.2.7', 'pyfaidx', 'pysam', 'gffutils', 'ujson'],
	python_requires='>=3.6',
	packages=setuptools.find_packages(),
	entry_points={'console_scripts': ['lifton = lifton.lifton:main'], },
        long_description=long_description,
        long_description_content_type='text/markdown'
)
