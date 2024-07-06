from setuptools import setup, find_packages
import codecs
import os.path
from pathlib import Path


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name="bindcompare",
    version=get_version("bindcompare/__init__.py"),
    description="A Novel DNA/RNA Binding Integration and Analysis Platform.",
    long_description=Path("README.rst").read_text(encoding="utf-8"),
    url="https://github.com/pranavmahabs/bindcompare",
    author="P. Mahableshwarkar, M.Ray, E. Larschan 2024",
    author_email="pranav_mahableshwarkar@brown.edu",
    license="GPL",
    packages=["bindcompare", "bindcompare.bindapp"],
    package_data={"bindcompare": ["reference_files/dm*"]},
    install_requires=[
        "numpy>=1.0",
        "pandas>=1.0",
        "customtkinter>=1.0",
        "matplotlib>=3.0",
        "biopython>=1.0",
        "intervaltree>=1.0",
        "matplotlib_venn>=0.11",
        "scipy>=1.0",
    ],
    entry_points={
        "console_scripts": [
            "bindcompare = bindcompare.bindcompare:main",
            "bindlaunch = bindcompare.bcapp:main",
            "retrievedm6 = bindcompare.retrieve:main",
            "comparexp = bindcompare.comparexp:main",
            "bindexplore = bindcompare.explore:main",
        ],
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
