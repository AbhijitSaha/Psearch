# Purpose: Cython wrapper for the python_pyc module
#  Author: Kenneth J. Mighell
# Version: 0.1.1  2018MAY07

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

# declare the interface to the C code
cdef extern void ctheta_slave_c(
    double* parrayAD,
    long sz_parrayL,
    double* magAD,
    long sz_magL,
    double* tobsAD,
    long sz_tobsL,
    double* thetaAD,
    long sz_thetaL)

@cython.boundscheck(False)
@cython.wraparound(False)
def ctheta_slave(np.ndarray[np.double_t, ndim=1, mode="c"] parrayAD not None,
    np.ndarray[np.double_t, ndim=1, mode="c"] magAD not None,
    np.ndarray[np.double_t, ndim=1, mode="c"] tobsAD not None
    ):
    cdef long sz_parrayL = parrayAD.size
    cdef long sz_magL = magAD.size
    cdef long sz_tobsL = tobsAD.size
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] thetaAD = parrayAD.copy(order='K')
    cdef long sz_thetaL = thetaAD.size
    ctheta_slave_c(
	&parrayAD[0], sz_parrayL,
	&magAD[0], sz_magL,
	&tobsAD[0], sz_tobsL,
	&thetaAD[0],sz_thetaL)
    return thetaAD

#EOF
