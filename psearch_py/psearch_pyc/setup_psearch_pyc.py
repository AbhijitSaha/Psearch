#!/usr/bin/env python

# Purpose: setup.py file for the psearch_pyc module
#  Author: Kenneth J. Mighell
# Version: 0.1.1  2018MAY07
#
# Build command:
#   python this_file.py build_ext --inplace

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("psearch_pyc",
                             sources=["psearch_pyc.pyx", "psearch_py_c.c"],
                             include_dirs=[numpy.get_include()])],
)
#EOF
