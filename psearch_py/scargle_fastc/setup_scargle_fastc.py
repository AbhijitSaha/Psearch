#!/usr/bin/env python

# Purpose: Setup.py file for the scargle_fastc module
#  Author: Kenneth Mighell
# Version: 0.3.1  2018MAY05
#
# Build command:
#   python this_file.py build_ext --inplace

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("scargle_fastc",
                             sources=["scargle_fastc.pyx", "scargle_fast_c.c"],
                             include_dirs=[numpy.get_include()])],
)
#EOF
