import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

exec(open("src/version.py").read())
setuptools.setup(
    name="caspeak",
    version=__version__,
    description="Finding non-reference mobile element insertions (MEIs) based on outer-Cas9 targeted Nanopore sequencing and peak detection",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Rye-lxy/CasPeak",
    author="Xinyi Liu",
    packages=["src"],
    install_requires=["pillow"],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
    ],
    scripts=["caspeak"],
)