# Purpose: Cython wrapper for the scargle_fastc module
#  Author: Kenneth Mighell
# Version: 0.3.1  2018MAY05

import cython
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern void scargle_fast_c(
    double* tAD,
    long sz_tL,
    double* cAD,
    long sz_cL,
    double* omegaAD,
    long sz_omegaL,
    long nfreqL,
    double* pxAD,
    long sz_pxL)

@cython.boundscheck(False)
@cython.wraparound(False)
def scargle_fast(np.ndarray[np.double_t, ndim=1, mode="c"] tAD not None,
                  np.ndarray[np.double_t, ndim=1, mode="c"] cAD not None,
                  np.ndarray[np.double_t, ndim=1, mode="c"] omegaAD not None,
		  long nfreqL
    ):
    cdef long sz_tL = tAD.size
    cdef long sz_cL = cAD.size
    cdef long sz_omegaL = omegaAD.size
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] pxAD = omegaAD.copy(order='K')
    cdef long sz_pxL = pxAD.size
    assert nfreqL == omegaAD.size, "*** ERROR *** nfreqL != omegaAD.size"
    scargle_fast_c(
	&tAD[0], sz_tL,
	&cAD[0], sz_cL,
	&omegaAD[0], sz_omegaL,
        nfreqL,
	&pxAD[0],sz_pxL)
    return pxAD
pass#enddef
#EOF
