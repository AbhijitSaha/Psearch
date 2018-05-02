#!/usr/bin/env python

# how to build:
#   python setup_ctheta_slave.py build_ext --inplace

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("ctheta_slavec",
                             sources=["ctheta_slavec.pyx", "ctheta_slave_c.c"],
                             include_dirs=[numpy.get_include()])],
)
#EOF
