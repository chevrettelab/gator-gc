from setuptools import setup, find_packages

VERSION = '0.9.0'
DESCRIPTION = 'GATOR-GC: Genomic Assessment Tool for Orthologous Regions and Gene Clusters'
setup(
    name="gator-gc",
    version=VERSION,
    author="José D. D. Cediel-Becerra",
    author_email="jcedielbecerra@ufl.edu",
    description=DESCRIPTION,
    packages=find_packages(),
    package_data={"": ["modular_domains.hmmdb"]},
    scripts=['bin/gator-gc', 'bin/pre-gator-gc'])

