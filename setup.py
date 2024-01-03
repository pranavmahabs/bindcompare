from setuptools import setup, find_packages
import codecs
import os.path


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
    url="https://github.com/pranavmahabs/bindcompare",
    author="P. Mahableshwarkar, M.Ray, E. Larschan 2024",
    author_email="pranav_mahableshwarkar@brown.edu",
    license="GPL",
    packages=["bindcompare", "bindcompare.bindapp"],
    package_data={"bindcompare": ["reference_files/dm*"]},
    install_requires=[
        "numpy",
        "pandas",
        "customtkinter",
        "matplotlib",
        "biopython",
        "intervaltree",
        "matplotlib_venn",
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
        "Development Status :: 1 - Development",
        "Intended Audience :: Bioinformatics Research",
        "License :: GPL License",
        "Operating System :: Unix :: Linux :: MacOS :: MacOS X",
        "Programming Language :: Python :: 3 :: Only",
    ],
)