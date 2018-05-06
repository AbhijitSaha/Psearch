#!/usr/bin/env python

# Purpose: setup.py file for the ctheta_slavec  module                                   
#  Author: Kenneth J. Mighell                                                           
# Version: 0.3.2  2018MAY06                                                             
#                                                                                       
# Build command:                                                                        
#   python this_file.py build_ext --inplace                                             

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
