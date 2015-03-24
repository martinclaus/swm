# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 10:17:59 2015

@author: mclaus
"""

from distutils.core import setup

setup(name="KSWM",
      version="0.1dev",
      author="Martin Claus",
      author_email="mclaus@geomar.de",
      packages=["kswm", "kswm.preprocess"],
      description="Utilities for the Kiel Shallow Water Model",
      license="Creative Commons Attribution-Noncommercial-Share Alike license",
      install_requires=[
          "numpy",
          "netCDF4",
          "subprocess32"],
      )
