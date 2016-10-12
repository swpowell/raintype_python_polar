"""
UW Rain Type
"""

import os
import sys
from setuptools import setup

# Pull the header into a variable
doclines = __doc__.split("\n")

VERSION = '1.0'

# Set variables for setup
PACKAGES = ['uw_raintype_polar']

# Run setup
setup(
      name='uw_raintype_polar',
      version=VERSION,
      url='http://www.atmos.washington.edu/MG/PDFs/JTECH16_Powell-etal_RainCat.pdf',
      author='Scott Powell, Stacy Brodzik',
      author_email='spowell@atmos.colostate.edu',
      description=("Rain-type classification code in native polar coordinates."),
      license='GNU',
      packages=PACKAGES,
      classifiers=["""
          Development Status :: V1.0,
          Programming Language :: Python",
          Topic :: Scientific/Engineering
          Topic :: Scientific/Engineering :: Atmospheric Science
          Operating System :: Unix
          Operating System :: POSIX :: Linux
          Operating System :: MacOS
          """],
      long_description="""
          Python tool for rain type classification of radar data in native polar coordinates.
          """,
      )
