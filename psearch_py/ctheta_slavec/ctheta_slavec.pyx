import cython
import numpy as np
cimport numpy as np

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
pass#enddef
#EOF
