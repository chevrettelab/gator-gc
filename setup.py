from setuptools import setup

VERSION = '0.9.0'
DESCRIPTION = 'GATOR-GC: Genomic Assessment Tool for Orthologous Regions and Gene Clusters'
setup(
    name="gator-gc",
    version=VERSION,
    author="José D. D. Cediel-Becerra",
    author_email="jcedielbecerra@ufl.edu",
    description=DESCRIPTION,
    packages=['gator-gc'],
    package_data={"flat": ["modular_domains.hmmdb"]},
    scripts=['bin/gator-gc', 'bin/pre-gator-gc'])

