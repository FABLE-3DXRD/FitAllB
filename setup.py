#!/usr/bin/env python
from distutils.core import setup,Extension
import sys

setup(
  name='FitAllB',
  version='0.0.1',
  description='Fitting routines for global parameters (fitgloball), global parameters for each grain (fitglobalgrain) and grain cms, orientations and strain (fitallb)',
  license='GPL', maintainer='Jette Oddershede',
  maintainer_email='jette.oddershede@risoe.dk',
  url='http://fable.wiki.sourceforge.net',
  packages=["FitAllB"],
  package_dir={"FitAllB":"FitAllB"},
  scripts=["scripts/fitallb.py","scripts/fitallb_1D.py","scripts/fitgloball.py","scripts/fitglobalgrain.py","scripts/fitgloball_multidet.py"]
)
