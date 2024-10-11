import setuptools
from pathlib import Path

this_directory = Path(__file__).resolve().parent
long_description = (this_directory / "./README.md").read_text()
setuptools.setup(
	name="lifton",
	version="1.0.5",
	author="Kuan-Hao Chao",
	author_email="kh.chao@cs.jhu.edu",
	description="Combining DNA and protein alignments to improve genome annotation with LiftOn",
	url="https://github.com/Kuanhao-Chao/Lifton",
	install_requires=['numpy>= 1.22.0', "biopython>=1.76", "cigar>=0.1.3", "parasail>=1.2.4", 'intervaltree>=3.1.0', 'interlap>=0.2.6', 'networkx>=3.3', 'pyfaidx>=0.5.8', 'pysam>=0.19.1', 'gffutils>=0.10.1', 'ujson>=3.2.0', 'pytest>=7.0.0'],
	python_requires='>=3.6',
	packages=setuptools.find_packages(),
	entry_points={'console_scripts': ['lifton = lifton.lifton:main'], },
        long_description=long_description,
        long_description_content_type='text/markdown'
)