import setuptools


setuptools.setup(
	name="Lifton",
	version="0.0.1",
	author="Kuan-Hao Chao",
	author_email="kh.chao@cs.jhu.edu",
	description="A protein-coding gene annotation fixing tool",
	url="https://github.com/Kuanhao-Chao/Lifton",
	install_requires=['numpy>=1.22.0', "parasail>=1.2.1"],
	python_requires='>=3.6',
	packages=['lifton'],
	entry_points={'console_scripts': ['lifton = lifton.lifton:main'], },
)
