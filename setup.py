import os
from setuptools import setup, find_packages

# version-keeping code based on pybedtools
curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 0
MIN = 1
REV = 0
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'ffqc/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )


setup(
    name='ffqc',
    version=VERSION,
    description='an imitation of FastQC',
    author='Joseph Edmonston',
    author_email='jedmonston@ucsd.edu',
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "ffqc=ffqc.ffqc:main"
        ],
    },
)