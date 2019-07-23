import setuptools
from Cython.Build import cythonize

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scDamAndT",
    version="0.0.1",
    author="K. Rooijers, F.J. Rang",
    author_email="j.kind@hubrecht.eu",
    description="Scripts to process raw scDam&T-seq data",
    long_description=long_description,
    url="https://github.com/...",
    packages=setuptools.find_packages(),
    install_requires = [
        "numpy",
        "h5py",
        "pysam",
        "tqdm",
        "sortedcontainers",
        "cython"
    ],
    scripts = [
        "./scripts/add_read_prefix.awk",
        "./scripts/bin_damid_counts.py",
        "./scripts/demultiplex.py",
        "./scripts/fetch_regions.py",
        "./scripts/find_motif_occurrences.py",
        "./scripts/generate_celseq_counts.py",
        "./scripts/generate_damid_counts.py",
        "./scripts/write_posarray.py"
    ],
    py_modules = [
        "support.fastaiter",
        "support.gtfparser",
        "support.stepvector",
        "support.posarray"
    ],
    ext_modules=cythonize(["support/attrsplitter.pyx"]),
)
