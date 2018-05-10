import numpy as np
import matplotlib.pyplot as plt
import numba
import time as tm
import platform
import os 
import sys

cythonc = True
try:
    import psearch_pyc
except ImportError:
    cythonc = False

# version information:
from collections import namedtuple
version_info = namedtuple('version_info','major minor micro')
version_info = version_info(major=0,minor=19,micro=5) 
__version__ = '%d.%d.%d' % version_info


def reference():
    msg ='pure Python  (*** slow ***)'
    if (cythonc):
        msg = 'Python/Cython/C  (*** fast ***)'
    print ' '
    print 'Saha, A., & Vivas, A. K. 2017, Astronomical Journal, 154, 231;'
    print '    "A Hybrid Algorithm for Period Analysis from Multiband Data with'
    print '    Sparse and Irregular Sampling for Arbitrary Light-curve Shapes"'
    print 'IDL CODE (Abhijit Saha):'
    print '    https://github.com/AbhijitSaha/Psearch'
    print 'PYTHON/CYTHON/C CODE (Kenenth Mighell):'
    print '    https://github.com/AbhijitSaha/Psearch/tree/master/psearch_py'
    print '\nMODULE:'
    print '    %s' % os.path.abspath(__file__)
    print '    [psearch_py (%s)  mode: %s ]' % (__version__,msg)
    print ' '
    return


def psearch_py( hjd, mag, magerr, filts, filtnams, pmin, dphi):
    """
    NAME:
        psearch_py

    INPUTS:
        hjd: time (Heliocentric Julian Day) input data array
        mag: magnitude input data array (co-alligned with hjd)
        magerr: magnitude error input data array (co-alligned with hjd)
        filts: filter input data array (co-aligned with hjd) with
            integer identifier indicating passband
        filtnams =  string array containing character names corresponding to
            coded filts values.  E.g., if you have 5 bands labeled u,g,r,i,z
            with filts values 0,1,2,3,4 respectively, filtnams would be set by:
            filtnams = ['u', 'g', 'r', 'i', 'z']
        pmin: Minimum value of period to be tested.
            E.g., pmin = 0.2
        dphi: Maximum change in relative phase between first and last epoch to
            be permitted when stepping to next test period.
            E.g., dphi = 0.02

    OUTPUTS:
        ptest: 1-D array with N dimensions of periods for which the periodograms
            are computed.  It is the same for ALL bands/channels.
        psi_m: M x N array of the Psi periodogram, where M is the number of
            bands/channels in the input array filtnams
        thresh_m: M x N array containing threshold values of Psi at each period
            and band for assessing significance for psi_m

    ORIGINAL IDL DEFINITION:
        pro Psearch, hjd, mag, magerr, filts, filtnams, pmin, dphi, ptest, $
            psi_m, thresh_m
    """
    assert isinstance(hjd,np.ndarray)
    assert (hjd.dtype == np.float64)
    assert (hjd.ndim == 1)
    hjd_shape = hjd.shape
    assert isinstance(mag,np.ndarray)
    assert (mag.dtype == np.float64)
    assert (mag.shape == hjd_shape)
    assert isinstance(magerr,np.ndarray)
    assert (magerr.dtype == np.float64)
    assert (magerr.shape == hjd_shape)
    assert isinstance(filts,np.ndarray)
    assert (filts.dtype == np.float64)
    assert (filts.shape == hjd_shape)
    print 'psearch: BEGIN ====================================================='
    print '\nREFERENCE:'
    reference()
    nfilts = len(filtnams)
    psiacc = 0.
    confacc = 0.
    for i in xrange(nfilts):
        fwant = i
        print 'psearch: ',filtnams[fwant],' filter'
        x, fy, theta, psi, conf = periodpsi2_py( hjd, mag, magerr, filts, \
            pmin, dphi, fwant )
        if (i == 0):
            # define the output arrays
            # -- needs size of period array from 1st band call to periodpsi2
            psi_m = np.zeros(shape=(nfilts,len(x)))
            thresh_m = np.zeros(shape=(nfilts,len(x)))
        psi_m[i,:] = psi
        thresh_m[i,:]  = conf
        table_psi_kjm_py( xx=x, yy=psi, ee=conf, n=10)
        psiacc = psi + psiacc
        confacc = conf + confacc
    print ' '
    print '========== ALL FILTERS ========== '
    print ' '
    table_psi_kjm_py( xx=x, yy=psiacc, ee=confacc, n=10)
    ptest = x
    print '\nReference:'
    reference()
    print 'psearch: END ======================================================='
    return ptest, psi_m, thresh_m  #KJM


def periodpsi2_py( hjd, mag, magerr, filts, minper, dphi, fwant ):
    """
    NAME:
        periodpsi2_py

    INPUTS:
        hjd: time (Heliocentric Julian Day) input data array
        mag: magnitude input data array (co-alligned with hjd)
        magerr: magnitude error input data array (co-alligned with hjd)
        filts: filter input data array (co-aligned with hjd) with
            integer identifier indicating passband
        minper: minimum period to explore
        dphi: maximum phase change between any two data points resulting from
            one step in frequency or period
        fwant: integer value corresponding to desired passband from among values
            in filts array.

    OUTPUTS: 
        x: period array for which periodograms are computed
        fy: Lomb-Scargle periodogram (co-aligned with x)
        theta: Lafler-Kinman periodogram (co-aligned with x)
        psi: psi periodogram (co-aligned with x)
        conf: simulated PSI periodogram for a non-periodic variable with
            amplitude and noise mimicking real source PLUS of an unvarying
            object with noise mimicking

    ORIGINAL IDL DEFINITION:
        pro periodpsi2, HJD, MAG, MAGERR, FILTS, minper, dphi, fwant, x, fy, $
            theta, psi, conf
    """
    print 'periodpsi2: BEGIN'
    debug1 = False
    #debug1 = True
    debug2 = False
    #debug2 = True
    assert isinstance(hjd,np.ndarray)
    assert (hjd.dtype == np.float64)
    assert (hjd.ndim == 1)
    hjd_shape = hjd.shape
    assert isinstance(mag,np.ndarray)
    assert (mag.dtype == np.float64)
    assert (mag.shape == hjd_shape)
    assert isinstance(magerr,np.ndarray)
    assert (magerr.dtype == np.float64)
    assert (magerr.shape == hjd_shape)
    assert isinstance(filts,np.ndarray)
    assert (filts.dtype == np.float64)
    assert (filts.shape == hjd_shape)
    t0 = np.min(hjd)
    tmax = np.max(hjd)
    tspan = tmax - t0
    ok = (filts == fwant) & (magerr <= 0.2)
    tr = hjd[ok]
    nok = len(tr)
    print 'periodpsi2: ',len(tr),' observations'
    yr = mag[ok]
    er = magerr[ok]*np.random.normal(0.,1.,nok)
    #amp = max(yr)-min(yr)  #UNUSED?  [KJM]
    zr = scramble_py( yr )[0]  # get only the first output
    sss = np.argsort(tr)  #BEWARE of the IDL sort() gotcha!
    tr = tr[sss]
    yr = yr[sss]
    er = er[sss]
    maxfreq = 1./minper
    minfreq = 2./tspan
    deltafreq = dphi/tspan
    nfreq = int( (maxfreq-minfreq)/deltafreq )
    #DEBUG: 
    print 'periodpsi2: number of frequency samples = ', nfreq
    farray = minfreq + np.arange(nfreq)*deltafreq
    x = 1./farray
    #DEBUG: print 'periodpsi2: minimum and maximum periods: ', min(x), max(x)
    omega = farray * 2.0 * np.pi

    time20 = tm.time()
    #om, fy = scargle_py( tr, yr, omega=omega, nfreq=nfreq, old=False )[:2]
    #fy = scargle_fast_py( tr, yr, omega, nfreq )
    #fy = psearch_pyc.scargle_fast( tr, yr, omega, nfreq ) 
    if (cythonc):
        fy = psearch_pyc.scargle_fast( tr, yr, omega, nfreq ) 
    else:
        fy = scargle_fast_py( tr, yr, omega, nfreq )
    time21   = tm.time()
    print 'scargle: DONE  %8.3f seconds' % (time21-time20)
    if (debug1):
        om_, fy_ = scargle_py( tr, yr, omega=omega, nfreq=nfreq, old=False )[:2]
        print np.allclose(fy,fy_),'=np.allclose(fy,fy_)'
        ok1 = np.allclose(fy,fy_)
        if (not ok1):
            print '^--- FY NOT OK!\n'
            plot_absdiff_py( fy_, fy, 'FY' )

    time20   = tm.time()
    #om, fe = scargle_py( tr, er, omega=omega, nfreq=nfreq, old=False )[:2]
    #fe = scargle_fast_py( tr, er, omega, nfreq )
    #fe = psearch_pyc.scargle_fast( tr, er, omega, nfreq )
    if (cythonc):
        fe = psearch_pyc.scargle_fast( tr, er, omega, nfreq ) 
    else:
        fe = scargle_fast_py( tr, er, omega, nfreq )
    time21   = tm.time()
    print 'scargle: DONE  %8.3f seconds' % (time21-time20)
    if (debug1):
        om_, fe_ = scargle_py( tr, er, omega=omega, nfreq=nfreq, old=False )[:2]
        print np.allclose(fe,fe_),'=np.allclose(fe,fe_)'
        ok2 = np.allclose(fe,fe_)
        if (not ok2):
            print '^--- FE NOT OK!\n'
            plot_absdiff_py( fe_, fe, 'FE' )

    time20   = tm.time()
    #om, fz = scargle_py( tr, zr, omega=omega, nfreq=nfreq, old=False )[:2]
    #fz = scargle_fast_py( tr, zr, omega, nfreq )
    #fz = psearch_pyc.scargle_fast( tr, zr, omega, nfreq )
    if (cythonc):
        fz = psearch_pyc.scargle_fast( tr, zr, omega, nfreq ) 
    else:
        fz = scargle_fast_py( tr, zr, omega, nfreq )
    time21   = tm.time()
    print 'scargle: DONE  %8.3f seconds' % (time21-time20)
    if (debug1):
        om_, fz_ = scargle_py( tr, zr, omega=omega, nfreq=nfreq, old=False )[:2]
        print np.allclose(fz,fz_),'=np.allclose(fz,fz_)'
        ok3 = np.allclose(fz,fz_)
        if (not ok3):
            print '^--- FZ NOT OK!\n'
            plot_absdiff_py( fz_, fz, 'FZ' )

    if (debug1):
        ok = ok1 and ok2 and ok3
        if (not ok):
            print '***** ERROR ************ RESULTS ARE NOT CLOSE!  :-(\n\n'
    
    time20 = tm.time()
    #theta = ctheta_slave_py(x, yr, tr, version=1)
    #theta = ctheta_slave_py(x, yr, tr)
    #theta = ctheta_slave_v3_pyjit(x, yr, tr)
    #theta = psearch_pyc.ctheta_slave(x, yr, tr)
    if (cythonc):
        #theta = psearch_pyc.ctheta_slave(x, yr, tr)
        theta = psearch_pyc.ctheta_slave(x, yr, tr)
    else:
        theta = ctheta_slave_v3_pyjit(x, yr, tr)
    time21   = tm.time()
    print 'ctheta_slave: DONE  %8.3f seconds' % (time21-time20)
    if (debug2):
        theta_ = ctheta_slave_py(x, yr, tr, version=1)
        print np.allclose(theta,theta_),'=np.allclose(theta,theta_)'
        ok1 = np.allclose(theta,theta_)
        if (not ok1):
            print '^--- THETA NOT OK!\n'
            plot_absdiff_py( theta_, theta, 'THETA' )

    time20   = tm.time()
    #thetaerr = ctheta_slave_py(x, er, tr, version=1)
    #thetaerr = ctheta_slave_py(x, er, tr)
    #thetaerr = ctheta_slave_v3_pyjit(x, er, tr)
    #thetaerr = psearch_pyc.ctheta_slave(x, er, tr)
    if (cythonc):
        thetaerr = psearch_pyc.ctheta_slave(x, er, tr)
    else:
        thetaerr = ctheta_slave_v3_pyjit(x, er, tr)
    time21   = tm.time()
    print 'ctheta_slave: DONE  %8.3f seconds' % (time21-time20)
    if (debug2):
        thetaerr_ = ctheta_slave_py(x, er, tr, version=1)
        print np.allclose(thetaerr,thetaerr_),'=np.allclose(thetaerr,thetaerr_)'
        ok2 = np.allclose(thetaerr,thetaerr_)
        if (not ok2):
            print '^--- THETAERR NOT OK!\n'
            plot_absdiff_py( thetaerr_, thetaerr, 'THETAERR' )

    time20   = tm.time()
    #thetaz   = ctheta_slave_py(x, zr, tr, version=1)
    #thetaz   = ctheta_slave_py(x, zr, tr)
    #thetaz   = ctheta_slave_v3_pyjit(x, zr, tr)
    #thetaz   = psearch_pyc.ctheta_slave(x, zr, tr)
    if (cythonc):
        thetaz = psearch_pyc.ctheta_slave(x, zr, tr)
    else:
        thetaz = ctheta_slave_v3_pyjit(x, zr, tr)
    time21   = tm.time()
    print 'ctheta_slave: DONE  %8.3f seconds' % (time21-time20)
    if (debug2):
        thetaz_ = ctheta_slave_py(x, zr, tr, version=1)
        print np.allclose(thetaz,thetaz_),'=np.allclose(thetaz,thetaz_)'
        ok3 = np.allclose(thetaz,thetaz_)
        if (not ok3):
            print '^--- THETAZ NOT OK!\n'
            plot_absdiff_py( thetaz_, thetaz, 'THETAZ' )

    if (debug2):
        ok = ok1 and ok2 and ok3
        if (not ok):
            print '***** ERROR ************ RESULTS ARE NOT CLOSE!  :-(\n\n'

    psi = 2.*fy/theta
    conf1 = 2.*fe/thetaerr
    conf1 = conf1*np.sum(psi)/np.sum(conf1)
    conf2 = 2.*fz/thetaz
    conf2 = conf2*np.sum(psi)/np.sum(conf2)
    conf = conf1 + conf2
    print 'periodpsi2: END'
    return x, fy, theta, psi, conf


def scramble_py( inarr ):
    """
    NAME:
        scramble_py

    INPUTS:
        inarr

    OUTPUTS:
        scrambarr
        pickarr

    ORIGINAL IDL DEFINITION:
        pro scramble, inarr, scrambarr, pickarr
    """
    #DEBUG: print 'scramble: BEGIN'
    ns = len(inarr)
    x = np.random.choice( ns, size=ns )
    assert x.dtype == np.int64
    assert (np.min(x) >= 0) & (np.max(x) < ns)
    s = np.argsort(x)  #KJM: BEWARE of the IDL sort() gotcha!
    assert s.dtype == np.int64
    assert (np.min(s) >= 0) & (np.max(s) < ns)
    # outputs:
    scrambarr = inarr[s]
    pickarr = inarr[x]
    #DEBUG: print 'scramble: END'
    return scrambarr, pickarr


def scargle_fast_py( t, c, omega, nfreq ):
    """
    NAME:
        scargle_fast_py

    PURPOSE:
        Compute the Lomb-Scargle periodogram of an unevenly sampled lightcurve

    CATEGORY:
        time series analysis

    INPUTS:
        t: The times at which the time series was measured (e.g. HJD)
        c: counts (corresponding count rates)
        omega: angular frequencies for which the PSD values are desired
            [PSD: Fourier Power Spectral Density]
        nfreq: number of independent frequencies
    
    OUTPUTS:
        px: the psd-values corresponding to omega
    
    DESCRIPTION:
        The Lomb Scargle PSD is computed according to the
        definitions given by Scargle, 1982, ApJ, 263, 835, and Horne
        and Baliunas, 1986, MNRAS, 302, 757. Beware of patterns and
        clustered data points as the Horne results break down in
        this case! Read and understand the papers and this
        code before using it! For the fast algorithm read W.H. Press
        and G.B. Rybicki 1989, ApJ 338, 277.

        The code is still stupid in the sense that it wants normal
        frequencies, but returns angular frequency...

    MODIFICATION HISTORY OF IDL VERSION:
        Version 1.0, 1997, Joern Wilms IAAT
        Version 1.1, 1998.09.23, JW: Do not normalize if variance is 0
            (for computation of LSP of window function...)
        Version 1.2, 1999.01.07, JW: force nfreq to be int
        Version 1.3, 1999.08.05, JW: added omega keyword
        Version 1.4, 1999.08
            KP: significance levels
            JW: pmin,pmax keywords
        Version 1.5, 1999.08.27, JW: compute the significance levels
            from the horne number of independent frequencies, and not from
            nfreq
        Version 1.6, 2000.07.27, SS and SB: added fast algorithm and FAP
            according to white noise lc simulations.
        Version 1.7, 2000.07.28 JW: added debug keyword, sped up
            simulations by factor of four (use /slow to get old
            behavior of the simulations)

    WEBSITE FOR THE IDL VERSION (Version 1.7, 2000.07.28):
        http://astro.uni-tuebingen.de/software/idl/aitlib/timing/scargle.pro

    ORIGINAL IDL DEFINITION:
        PRO scargle,t,c,om,px,fmin=fmin,fmax=fmax,nfreq=nfreq,     $
            nu=nu,period=period,omega=omega,                       $
            fap=fap,signi=signi,simsigni=simsigni,                 $
            pmin=pmin,pmax=pmax,old=old,                           $
            psdpeaksort=psdpeaksort,multiple=multiple,noise=noise, $
            debug=debug,slow=slow
    """
    assert isinstance(t,np.ndarray)
    assert (t.dtype == np.float64)
    assert (t.ndim == 1)
    assert isinstance(c,np.ndarray)
    assert (c.dtype == np.float64)
    assert (c.ndim == 1)
    assert (len(t)==len(c))
    assert isinstance(omega,np.ndarray)
    assert (omega.dtype == np.float64)
    assert (omega.ndim == 1)
    assert (omega.size == nfreq)
    #
    noise = np.sqrt(np.var(c))
    #    
    # make times manageable (Scargle periodogram is time-shift invariant)
    time = t-t[0]
    #
    n0 = len(time)
    assert (n0 > 0)
    #
    om = omega  # alias
    #
    #===== PERIODOGRAM: FAST VERSION ========================================
    # Reference:
    #     Press, W.H., & Rybicki, G.B. 1989, ApJ 338, 277;
    #     FAST ALGORITHM FOR SPECTRAL ANALYSIS OF UNEVENLY SAMPLED DATA"
    #
    # Eq. (6); s2, c2
    s2 = np.zeros(nfreq)
    c2 = np.zeros(nfreq)
    two_time = 2.0*time
    for i in xrange(nfreq):
        s2[i] = np.sum( np.sin(two_time*om[i]) )
        c2[i] = np.sum( np.cos(two_time*om[i]) )
        #s2[i] = np.sum( np.sin(2.*om[i]*time) )
        #c2[i] = np.sum( np.cos(2.*om[i]*time) )
    two_time = 0.0  # clean up
    #
    # Eq. (2): Definition -> tan(2omtau)
    # --- tan(2omtau) =  s2 / c2
    omtau = np.arctan(s2/c2) / (2.)
    # cos(tau), sin(tau)
    cosomtau = np.cos(omtau)
    sinomtau = np.sin(omtau)
    #
    # Eq. (7); total(cos(t-tau)^2) and total(sin(t-tau)^2)
    tmp = c2*np.cos(2.*omtau) + s2*np.sin(2.*omtau)
    tc2 = 0.5*(n0+tmp)
    ts2 = 0.5*(n0-tmp)
    # clean up
    tmp = 0.
    omtau = 0.
    s2 = 0.
    #
    # computing the periodogram for the original light curve
    #
    # subtract mean from data
    cn = c - np.mean(c)
    #   
    # Eq. (5); sh and ch
    sh = np.zeros(nfreq)
    ch = np.zeros(nfreq)
    omi_time = np.zeros(nfreq)
    for i in xrange(nfreq):
        omi_time = om[i] * time
        sh[i] = np.sum( cn*np.sin(omi_time) )
        ch[i] = np.sum( cn*np.cos(omi_time) )
        #sh[i] = np.sum( cn*np.sin(om[i]*time) )
        #ch[i] = np.sum( cn*np.cos(om[i]*time) )
    omi_time = 0.0  # clean up
    #
    # Eq. (3)
    px  = ((ch*cosomtau + sh*sinomtau)**2 / tc2) + \
          ((sh*cosomtau - ch*sinomtau)**2 / ts2)
    # correct normalization
    px = 0.5*px/(noise**2)
    #
    return px


def scargle_py(
    #INPUTS:
    t,c,
    #OPTIONAL INPUTS:
    fmin=None,fmax=None,pmin=None,pmax=None,omega=None,fap=None,noise=None,
    multiple=None,nfreq=None,
    #OPTIONAL BOOLEAN INPUTS (IDL: KEYWORD PARAMETERS):
    old=None,debug=False,slow=None
    ):
    """
   NAME:
        scargle_py

    PURPOSE:
        Compute the Lomb-Scargle periodogram of an unevenly sampled lightcurve

    CATEGORY:
        time series analysis

    INPUTS:
        t: The times at which the time series was measured (e.g. HJD)
        c: counts (corresponding count rates)

    OPTIONAL INPUTS:
        fmin: minimum frequency (NOT ANGULAR FREQ!) to be used
              (has precedence over pmin)
        fmax: maximum frequency (NOT ANGULAR FREQ!) to be used
              (has precedence over pmax)
        pmin: minimum PERIOD to be used
        pmax: maximum PERIOD to be used
        omega: angular frequencies for which the PSD values are desired
               [PSD: Fourier Power Spectral Density]
        fap: false alarm probability desired
             (see Scargle et al., p. 840,a and signi keyword).
             Default equal to 0.01 (99% significance)
        noise: for the normalization of the periodogram and the
               compute (sp?) of the white noise simulations.
               If not set, equal to the variance of the original lc.
        multiple: number of white  noise simulations for the FAP
                  power level. Default equal to 0 (i.e., no simulations).
        nfreq: number of independent frequencies

    OPTIONAL BOOLEAN INPUTS (IDL: KEYWORD PARAMETERS):
        old: if set computing the periodogram according to
            Scargle, J.D. 1982, ApJ 263, 835.
            If not set, compute the periodogram with the fast algorithm of
            Press, W.H., & G.B. Rybicki, G.B. 1989, ApJ 338, 277.
        debug: print out debugging information if set
        slow: if set, a much slower but less memory intensive way to
              perform the white noise simulations is used.

    OUTPUTS:
        om: angular frequency of PSD [PSD: Fourier Power Spectral Density]
        px: the psd-values corresponding to omega
            [KJM: original IDL documentation refers to psd
            --- which did not exist]
        nu: normal frequency  [nu = om/(2.*np.pi)]
        period: period corresponding to each omega  [period = 1./nu]
        signi: power threshold corresponding to the given false alarm
            probabilities fap and according to the desired number of independent
            frequencies
        simsigni: power threshold corresponding to the given false alarm
                  probabilities fap according to white noise simulations
        psdpeaksort: array with the maximum peak pro (sp?) each simulation

    DESCRIPTION:
        The Lomb Scargle PSD is computed according to the
        definitions given by Scargle, 1982, ApJ, 263, 835, and Horne
        and Baliunas, 1986, MNRAS, 302, 757. Beware of patterns and
        clustered data points as the Horne results break down in
        this case! Read and understand the papers and this
        code before using it! For the fast algorithm read W.H. Press
        and G.B. Rybicki 1989, ApJ 338, 277.

        The code is still stupid in the sense that it wants normal
        frequencies, but returns angular frequency...

    MODIFICATION HISTORY OF IDL VERSION:
        Version 1.0, 1997, Joern Wilms IAAT
        Version 1.1, 1998.09.23, JW: Do not normalize if variance is 0
            (for computation of LSP of window function...)
        Version 1.2, 1999.01.07, JW: force nfreq to be int
        Version 1.3, 1999.08.05, JW: added omega keyword
        Version 1.4, 1999.08
            KP: significance levels
            JW: pmin,pmax keywords
        Version 1.5, 1999.08.27, JW: compute the significance levels
            from the horne number of independent frequencies, and not from
            nfreq
        Version 1.6, 2000.07.27, SS and SB: added fast algorithm and FAP
            according to white noise lc simulations.
        Version 1.7, 2000.07.28 JW: added debug keyword, sped up
            simulations by factor of four (use /slow to get old
            behavior of the simulations)

    WEBSITE FOR THE IDL VERSION (Version 1.7, 2000.07.28):
        http://astro.uni-tuebingen.de/software/idl/aitlib/timing/scargle.pro

    ORIGINAL IDL DEFINITION:
        PRO scargle,t,c,om,px,fmin=fmin,fmax=fmax,nfreq=nfreq,     $
            nu=nu,period=period,omega=omega,                       $
            fap=fap,signi=signi,simsigni=simsigni,                 $
            pmin=pmin,pmax=pmax,old=old,                           $
            psdpeaksort=psdpeaksort,multiple=multiple,noise=noise, $
            debug=debug,slow=slow
    """
    #DEBUG: print 'scargle: BEGIN'
    # initial optional output default values
    nu = None
    period = None
    signi = None
    simsigni = None
    psdpeaksort = None
    # defaults
    if noise is None:
        assert c.dtype ==  np.float64
        noise = np.sqrt(np.var(c))
    if multiple is None: multiple = 0
    if fap is None: fap = 0.01
    # make times manageable (Scargle periodogram is time-shift invariant)
    time = t-t[0]
    # number of independent frequencies
    # (Horne and Baliunas, eq. 13)
    n0 = len(time)
    assert n0 > 0
    horne = long(-6.362+1.193*n0+0.00098*(n0**2.))
    if horne < 0: horne = 5
    if nfreq is None:
        nfreq = horne
    else:
        horne = nfreq
    # mininum frequency is 1/T
    if fmin is None:
        #IF (n_elements(pmax) EQ 0) THEN BEGIN
        if pmax is None:
            fmin = 1. / max(time)
        else:
            fmin = 1. / pmax
        pass
    # maximum frequency approximately equal to Nyquist frequency
    if fmax is None:
        #IF (n_elements(pmin) EQ 0) THEN BEGIN
        if pmin is None:
            fmax = n0 / (2.*max(time))
        else:
            fmax = 1. / pmin
        pass
    # if omega is not given, compute it
    #IF (n_elements(omega) EQ 0) THEN BEGIN
    if omega is None:
        om = 2. * np.pi * (fmin+(fmax-fmin)* \
            np.arange(nfreq,dtype=np.float64)/(nfreq-1.))
    else:
        om = omega
    signi = -np.log( 1. - ((1.-fap)**(1./horne)) )
    #
    # Periodogram
    #
    #
    #===== PERIODOGRAM: SLOW VERSION ========================================
    if (old is True):
        # Subtract mean from data
        assert c.dtype == np.float64
        cn = c-np.mean(c)
        # computing the periodogram
        px = np.zeros(nfreq, dtype=np.float64)
        for i in xrange(nfreq):
            tau = np.arctan(np.sum(np.sin(2.*om[i]*time))/\
                np.sum(np.cos(2.0*om[i]*time)))
            tau = tau/(2.*om[i])
            co = np.cos(om[i]*(time-tau))
            si = np.sin(om[i]*(time-tau))
            px[i]=0.5*(np.sum(cn*co)**2/np.sum(c**2)+np.sum(cn*si)**2/\
                np.sum(si**2))
        # correct normalization
        var = np.var(cn)
        if var != 0:
            px = px/var
        else:
            print 'scargle: ***** WARNING ***** Variance is zero (var == 0)'
        # some other nice helpers
        # computed here due to memory usage reasons
        nu = om/(2.*np.pi)
        period = 1./nu
        #DEBUG: print 'scargle: DONE (slow version)'
        #DEBUG: print 'scargle: END (slow version)'
        return om,px,nu,period,signi,simsigni,psdpeaksort
    #
    #
    #===== PERIODOGRAM: FAST VERSION ========================================
    # Reference:
    #     Press, W.H., & Rybicki, G.B. 1989, ApJ 338, 277;
    #     FAST ALGORITHM FOR SPECTRAL ANALYSIS OF UNEVENLY SAMPLED DATA"
    # Eq. (6); s2, c2
    s2 = np.zeros(nfreq)
    c2 = np.zeros(nfreq)
    for i in xrange(nfreq):
        s2[i] = np.sum( np.sin(2.*om[i]*time) )
        c2[i] = np.sum( np.cos(2.*om[i]*time) )
    # Eq. (2): Definition -> tan(2omtau)
    # --- tan(2omtau) =  s2 / c2
    omtau = np.arctan(s2/c2) / (2.)
    # cos(tau), sin(tau)
    cosomtau = np.cos(omtau)
    sinomtau = np.sin(omtau)
    # Eq. (7); total(cos(t-tau)^2) and total(sin(t-tau)^2)
    tmp = c2*np.cos(2.*omtau) + s2*np.sin(2.*omtau)
    tc2 = 0.5*(n0+tmp)
    ts2 = 0.5*(n0-tmp)
    # clean up
    tmp = 0.
    omtau = 0.
    s2 = 0.
    #t2 = 0.  #UNUSED?  [KJM]
    # computing the periodogram for the original lc
    # Subtract mean from data
    cn = c - np.mean(c)
    # Eq. (5); sh and ch
    sh = np.zeros(nfreq)
    ch = np.zeros(nfreq)
    if (multiple > 0) and (slow is None):
        sisi = np.zeros(shape=(n0,nfreq))
        coco = np.zeros(shape=(n0,nfreq))
        for i in xrange(nfreq):
            sisi[:,i] = np.sin(om[i]*time)
            coco[:,i] = np.cos(om[i]*time)
            sh[i] = np.sum(cn*sisi[:,i])
            ch[i] = np.sum(cn*coco[:,i])
    else:
        for i in xrange(nfreq):
            sh[i] = np.sum( cn*np.sin(om[i]*time) )
            ch[i] = np.sum( cn*np.cos(om[i]*time) )
    pass
    # Eq. (3)
    px  = ((ch*cosomtau + sh*sinomtau)**2 / tc2) + \
          ((sh*cosomtau - ch*sinomtau)**2 / ts2)
    # correct normalization
    px = 0.5*px/(noise**2)
    # --- RUN SIMULATIONS for multiple > 0
    if (multiple > 0):
        if (multiple*min(fap) < 10):
            print 'scargle: message: WARNING [/informational]'
            print 'scargle: message: Number of iterations (multiple keyword)'+\
                ' [/informational]'
            print 'scargle: message: not large enough for false alarm '+\
                'probability [/informational]'
            print 'scargle: message: requested (need multiple*FAP > 10 ) '+\
                '[/informational]'
        psdpeak = np.zeros(multiple)
        for m in xrange(multiple):
            if debug is True:
                if ((m+1) % 100 == 0):
                    print 'scargle: working on the', m,'th simulation'
            # white noise simulation
            cn = np.random.normal(0.,1.,n0)*noise
            cn = cn-np.mean(cn)  # .. force OBSERVED count rate to zero
            # Eq. (5); sh and ch
            if slow is not None:
                for i in xrange(nfreq):
                    sh[i] = np.sum(cn*sisi[:,i])
                    ch[i] = np.sum(cn*coco[:,i])
            else:
                for i in xrange(0L, nfreq-1L):
                    sh[i] = np.sum( cn * np.sin(om[i]*time) )
                    ch[i] = np.sum( cn * np.cos(om[i]*time) )
            # Eq. (3) ; computing the periodogram for each simulation
            spud  = ((ch*cosomtau + sh*sinomtau)**2 / tc2) + \
                    ((sh*cosomtau - ch*sinomtau)**2 / ts2)
            psdpeak[m] = max( spud )
        # False Alarm Probability according to simulations
        if len(psdpeak) != 0:
            idx = np.argsort(psdpeak) #BEWARE of the IDL sort() gotcha!
            # correct normalization
            psdpeaksort = 0.5 * psdpeak[idx]/(noise**2)
            simsigni = psdpeaksort[long((1.-fap)*(multiple-1))]
    # some other nice helpers
    # computed here due to memory usage reasons
    nu = om/(2.*np.pi)
    period = 1./nu
    #DEBUG: print 'scargle: DONE (fast version)'
    #DEBUG: print 'scargle: END (fast version)'
    return om,px,nu,period,signi,simsigni,psdpeaksort


def ctheta_slave_py( parray, mag, tobs, version=2 ):
    """
    NAME:
        ctheta_slave_py 

    INPUTS:
        parray
        mag
        tobs
        (version)

    OUTPUTS:
        theta

    DESCRIPTION:
        Computes theta for a pre-specified array of test periods.

    ORIGINAL IDL DEFINITION:
        pro Ctheta_slave, parray, mag, tobs, theta
    """
    assert isinstance(parray,np.ndarray)
    assert (parray.dtype == np.float64)
    assert (parray.ndim == 1)
    assert (parray.size >= 1)
    assert isinstance(mag,np.ndarray)
    assert (mag.dtype == np.float64)
    assert (mag.ndim == 1)
    assert (mag.size >= 1)
    assert isinstance(tobs,np.ndarray)
    assert (tobs.dtype == np.float64)
    assert (tobs.shape == mag.shape)
    assert isinstance(version,int)
    #DEBUG: print 'ctheta_slave: BEGIN'
    #DEBUG: time10 = tm.time()
    t0 = np.min(tobs)
    #tlast = np.max(tobs)  #UNUSED?  [KJM]
    tt = tobs - t0
    theta = 0.*parray
    # Loop over all periods
    if (version == 2):
        # optimized version (about 35% faster than original)
        avm_km = np.sum(mag)/len(mag) 
        denom_km = np.sum( (mag-avm_km)**2 )
        for k in xrange(len(parray)):
            period = parray[k]
            phi = tt/period
            nphi = phi.astype(np.int64)
            phi = phi - nphi
            ss = np.argsort(phi)
            mm  = mag[ss]
            mmplus = np.append(mm[1:], mm[0])
            numer = np.sum( (mmplus - mm)**2 )
            theta[k] = numer/denom_km
    elif (version == 1):
        # original version (line-by-line translation of IDL code)
        for k in xrange(len(parray)):
            period = parray[k]
            phi = tt/period
            nphi = phi.astype(np.int64)
            phi = phi - nphi
            ss = np.argsort(phi)  #KJM: BEWARE of the IDL sort() gotcha!
            phi = phi[ss]  #KJM: Not used -- so why compute?
            mm  = mag[ss]
            avm = np.sum(mm)/len(mm)       #KJM: Suboptimal: move before loop
            denom = np.sum( (mm-avm)**2 )  #KJM: Suboptimal: move before loop
            mmplus = np.append(mm[1:], mm[0])
            numer = np.sum( (mmplus - mm)**2 )
            theta[k] = numer/denom
    else:
        assert( (version == 2) or (version == 1))  #KJM: bad version value
    #DEBUG: time11 = tm.time()
    #DEBUG: print 'ctheta_slave: DONE  [%d] %8.3f seconds' % (version,(time11-time10))
    #DEBUG: print 'ctheta_slave: END'
    return theta


def ctheta_slave_v3_py( parray, mag, tobs ):
    """
    NAME:
        ctheta_slave_v3_py 

    INPUTS:
        parray
        mag
        tobs

    OUTPUTS:
        theta

    DESCRIPTION:
        Computes theta for a pre-specified array of test periods.

    ORIGINAL IDL DEFINITION:
        pro Ctheta_slave, parray, mag, tobs, theta
    """
    t0 = np.min(tobs)
    tt = tobs - t0
    theta = np.zeros_like(parray)
    mmplus_km = np.zeros_like(mag)
    avm_km = np.sum(mag)/len(mag) 
    denom_km = np.sum( (mag-avm_km)**2 )
    m = len(parray)
    for k in xrange(m):
        period = parray[k]
        phi = tt / period
        #nphi = np.fix(phi)          #KJM: literal but slower
        nphi = phi.astype(np.int64)  #KJM: ~25% faster
        phi = phi - nphi
        ss = np.argsort(phi)  #KJM: BEWARE the IDL sort gotcha!
        mm  = mag[ss]
        #mmplus = np.append(mm[1:], mm[0])   #KJM: literal but slower
        #numer = np.sum( (mmplus - mm)**2 )  #KJM: uses mmplus
        mmplus_km[:-1] = mm[1:]  #KJM: Don't use np.append within loops!
        mmplus_km[-1] = mm[0]    #KJM: Don't use np.append within loops!
        #assert np.allclose(mmplus,mmplus_km)  #KJM: NUM? ME VEXO?
        numer = np.sum( (mmplus_km - mm)**2 )  #KJM: uses mmplus_km
        #KJM: using mmplus_km is ~24% faster
        theta[k] = numer/denom_km
    return theta


#KJM: ========== NUMBA JUST-IN-TIME COMPILATION: BEGIN 
    
ctheta_slave_v3_pyjit = numba.jit(nopython=True)(ctheta_slave_v3_py)

#KJM: ========== NUMBA JUST-IN-TIME COMPILATION: END

  
def table_psi_kjm_py( xx=None, yy=None, ee=None, n=None):
    """
    NAME:
        table_psi_kjm_py

    INPUTS:
        xx: periods (e.g., x)
        yy: power   (e.g., psi)
        ee: thresh  (e.g., conf)
        n:  number of ranked periods to show (e.g., 10)
    """
    assert isinstance(xx,np.ndarray)
    assert (xx.ndim == 1)
    xx_shape = xx.shape
    assert isinstance(yy,np.ndarray)
    assert (yy.shape == xx_shape)
    assert isinstance(ee,np.ndarray)
    assert (ee.shape == xx_shape)
    assert isinstance(n,np.int)
    assert (n >= 1)
    sz = len(xx)
    lm_x = np.zeros(sz)
    lm_y = np.zeros(sz)
    lm_k = np.zeros(sz,dtype=np.int_)
    j = 0
    assert (sz >= 3)
    for k in xrange(1,sz-1):
        ym1 = yy[k-1]
        y   = yy[k]
        yp1 = yy[k+1]
        if ((y>ym1) and (y>yp1)):
            lm_y[j] = yy[k]  # local maximum psiacc value
            lm_x[j] = xx[k]  # local maximum period value
            lm_k[j] = k
            j += 1
            lm_n = j
    lm_x = lm_x[:lm_n]
    lm_y = lm_y[:lm_n]
    lm_k = lm_k[:lm_n]
    assert (len(lm_y) >= n)
    idx = (-lm_y).argsort()[:n]  # indices (location) of the n largest values
    print 'TABLE: BEGIN'
    print 'rank   -----Period [days]-----       Psi    index  Frequency  Thresh'
    fmt = '%2d  %12.7f +- %10.7f %9.2f %8d %10.6f %7.2f'
    for j in xrange(n):
        k=idx[j]
        kk = lm_k[k]
        p0 = xx[kk]
        y0 = yy[kk]
        y0err = ee[kk]
        kkp1 = kk + 1
        p1 = xx[kkp1]
        sigma = (p0-p1)/2.  # estimate of error (one standard deviation)
        print fmt % ( j+1, p0, sigma, y0, kk, 1./p0, y0err)
    print 'TABLE: END'
    return


def fig_obs_kjm_py( hjd=None, mag=None, filts=None, filtnams=None, tag=None, \
    plotfile=None ):
    """
    NAME:
        fig_obs_kjm_py

    INPUTS:
        hjd: time (Heliocentric Julian Day) input data array
        mag: magnitude input data array (co-alligned with hjd)
        filts: filter input data array (co-aligned with hjd) with
            integer identifier indicating passband
        filtnams =  string array containing character names corresponding to
            coded filts values.  E.g., if you have 5 bands labeled u,g,r,i,z
            with filts values 0,1,2,3,4 respectively, filtnams would be set by:
            filtnams = []'u', 'g', 'r', 'i', 'z']
        tag: String written in the bottom-left corner (if any)
        plotfile: filename of the output plot (if any)
    """
    assert isinstance(hjd,np.ndarray)
    assert (hjd.ndim == 1)
    hjd_shape = hjd.shape
    assert isinstance(mag,np.ndarray)
    assert (mag.shape == hjd_shape)
    assert isinstance(filts,np.ndarray)
    assert (filts.shape == hjd_shape)
    #DEBUG: print 'fig_obs_kjm: ','OK!  :-)'
    color = ['dodgerblue']  # matplotlib 1.5.0 color names
    nfilts = len(filtnams)
    hjd0 = long(np.min(hjd))  # offset
    x = hjd - hjd0  # days from offset
    dx = max(0.08 * np.max(x),0.25)
    xlim = [-dx, (np.max(x)+dx)]
    xlabel = 'HJD - %d [days]' % hjd0
    dy = 0.5  # delta_mag
    fig, axarr = plt.subplots( nfilts, sharex=True, figsize=(8.5,11) )
    for i in xrange(nfilts):
        fwant = float(i)
        ok = (filts == fwant)
        xx = x[ok]
        yy = mag[ok]
        axarr[i].scatter( xx, yy, color=color[0], alpha=0.5 )
        axarr[i].set_xlim( xlim )  # all subplots have the same X-axis
        # expand and flip Y-axis:
        axarr[i].set_ylim( [(np.max(yy)+dy),(np.min(yy)-dy)] )
        axarr[i].set_ylabel( 'mag', size='x-large' )
        axarr[i].text( 0.97, 0.80, filtnams[i], ha='right', size='x-large', \
            transform=axarr[i].transAxes ) # relative coordinates within subplot
        if (i == (nfilts-1)):  # last subplot needs a label for the X-axis
            axarr[i].set_xlabel( xlabel, size='x-large' )
    if (tag is not None):
        plt.figtext( 0.95, 0.1, tag, ha='right', va='bottom', color='grey', \
            size='large', rotation=90)
    if (plotfile is not None):
        plt.savefig( plotfile, dpi=300 )
    #DEBUG: plt.show()
    plt.close()
    if (plotfile is not None):
        print plotfile, ' <--- plotfile written  :-)'
    return


def fig_psi_kjm_py( freq=None, psi_m=None, thresh_m=None, filtnams=None, \
        tag=None, plotfile=None):
    """
    NAME:
        fig_psi_kjm_py

    INPUTS:
        freq: 1-D array (length of N) frequencies for which the periodograms are
            computed.  It is the same for ALL bands/channels.
        psi_m: M x N array of the Psi periodogram, where M is the number of
            bands/channels in the input array filtnams
        thresh_m: M x N array containing threshold values of Psi at each period
            and band for assessing significance for psi_m
        filtnams =  string array (length of M) containing character names
            corresponding to coded filts values.  E.g., 5 bands labeled
            u,g,r,i,z with filts values:
            filtnams = ['u', 'g', 'r', 'i', 'z']
        tag: String written in the bottom-left corner (if any)
        plotfile: filename of the output plot (if any)
    """
    assert (filtnams is not None)
    nfilts = len(filtnams)
    assert isinstance(freq,np.ndarray)
    assert (freq.ndim == 1)
    ndata = len(freq)
    assert isinstance(psi_m,np.ndarray)
    psi_m_shape = psi_m.shape
    assert (psi_m_shape[0] == nfilts)
    assert (psi_m_shape[1] == ndata)
    assert isinstance(thresh_m,np.ndarray)
    assert (thresh_m.shape == psi_m_shape)
    #DEBUG: print 'fig_psi_kjm: ','OK!  :-)'
    color = ['dodgerblue','salmon']  # matplotlib 1.5.0 color names
    fig, axarr = plt.subplots( nfilts+1, sharex=True, figsize=(8.5,11) )
    for i in xrange(len(filtnams)):
        axarr[i].plot( freq,    psi_m[i], color=color[0], zorder=0 )
        axarr[i].plot( freq, thresh_m[i], color=color[1], zorder=10 )
        axarr[i].set_ylabel( r'${\Psi}$', size=19 )
        axarr[i].text( 0.97, 0.80, filtnams[i], ha='right', size='x-large', \
            transform=axarr[i].transAxes ) # relative coordinates within subplot
    # combine results for all filters
    j = nfilts
    axarr[j].plot( freq,    psi_m.sum(0), color=color[0], zorder=0 )
    axarr[j].plot( freq, thresh_m.sum(0), color=color[1], zorder=10 )
    axarr[j].set_ylabel( r'${\Psi}$', size=19 )
    axarr[j].set_xlabel( r'Frequency [days${^{-1}}$]', size='x-large' )
    axarr[j].text( 0.97, 0.80, 'ALL', ha='right', size='x-large', \
        transform=axarr[j].transAxes ) # relative coordinates within subplot
    if (tag is not None):
        plt.figtext( 0.95, 0.1, tag, ha='right', va='bottom', color='grey', \
            size='large', rotation=90)
    if (plotfile is not None):
        plt.savefig( plotfile, dpi=300 )
    #DEBUG: plt.show()
    plt.close()
    if (plotfile is not None):
        print plotfile, ' <--- plotfile written  :-)'
    return


def fig_phi_kjm_py( hjd=None, mag=None, magerr=None, filts=None, filtnams=None,
        period=None, tag=None, plotfile=None ):
    """
    NAME:
        fig_phi_kjm_py

    INPUTS:
        hjd: time (Heliocentric Julian Day) input data array
        mag: magnitude input data array (co-alligned with hjd)
        magerr: magnitude error input data array (co-alligned with hjd)
        filts: filter input data array (co-aligned with hjd) with
            integer identifier indicating passband
        filtnams =  string array containing character names corresponding to
            coded filts values.  E.g., if you have 5 bands labeled u,g,r,i,z
            with filts values 0,1,2,3,4 respectively, filtnams would be set by:
            filtnams = ['u', 'g', 'r', 'i', 'z']
        period: period to be used to phase up the data [days]
        tag: String written in the bottom-left corner (if any)
        plotfile: filename of the output plot (if any)
    """
    assert isinstance(hjd,np.ndarray)
    assert (hjd.ndim == 1)
    hjd_shape = hjd.shape
    assert isinstance(mag,np.ndarray)
    assert (mag.shape == hjd_shape)
    assert isinstance(magerr,np.ndarray)
    assert (magerr.shape == hjd_shape)
    assert isinstance(filts,np.ndarray)
    assert (filts.shape == hjd_shape)
    assert (period is not None)
    #DEBUG: print 'fig_phi_kjm: ','OK!  :-)'
    color = ['dodgerblue']  # matplotlib 1.5.0 color names
    nfilts = len(filtnams)
    hjd0 = long(np.min(hjd))  # offset
    x = hjd - hjd0  # days from offset
    dx = 0.1
    xlim = [0.0-dx, 2.0+dx]
    xlabel = r'${\phi}$'
    dy = 0.5  # delta_mag
    fig, axarr = plt.subplots( nfilts, sharex=True, figsize=(8.5,11) )
    for i in xrange(nfilts):
        fwant = float(i)
        ok = (filts == fwant)
        xx = x[ok]
        yy = mag[ok]
        ee = magerr[ok]
        phi0 = xx / period
        nphi0 = phi0.astype(np.int64)
        phi = phi0 - nphi0
        axarr[i].errorbar(   phi, yy, yerr=ee, fmt='o', color=color[0], \
            alpha=0.5 )
        axarr[i].errorbar( phi+1, yy, yerr=ee, fmt='o', color=color[0], \
            alpha=0.5 )
        axarr[i].set_xlim( xlim )  # all subplots have the same X-axis
        # expand and flip Y-axis:
        axarr[i].set_ylim( [(np.max(yy+ee)+dy),(np.min(yy-ee)-dy)] )
        axarr[i].set_ylabel( 'mag', size='x-large' )
        axarr[i].text( 0.97, 0.80, filtnams[i], ha='right', size='x-large', \
            transform=axarr[i].transAxes ) # relative coordinates within subplot
        if (i == (nfilts-1)):  # last subplot needs a label for the X-axis
            axarr[i].set_xlabel( xlabel, size=20 )
    plt.figtext( 0.5, 0.93, 'Period: %9.6f days' % period, ha='center', \
         color='black', size='xx-large' )
    if (tag is not None):
        plt.figtext( 0.95, 0.1, tag, ha='right', va='bottom', color='grey', \
            size='large', rotation=90)
    if (plotfile is not None):
        plt.savefig( plotfile, dpi=300 )
    #DEBUG: plt.show()
    plt.close()
    if (plotfile is not None):
        print plotfile, ' <--- plotfile written  :-)'
    return

def do_stats( x, tag=None ):
    if (tag is None):
        tag = ''
    if (len(x)>0):
        print tag, 'median : ', np.median(x)
        print tag, '  mean : ', np.mean(x)
        print tag, '   std : ', np.std(x)
        print tag, '   min : ', np.min(x)
        print tag, '   max : ', np.max(x)
        print tag, '     n : ', len(x)
        print tag, '   NAN : ', np.count_nonzero(np.isnan(x))
    
    
def plot_absdiff_py( fy0_p, fy1_p, title_p ):
    print '\n========== ',title_p,' =========='
    fy0 = fy0_p[:]
    fy1 = fy1_p[:]
    fdy = (fy1 - fy0)
    fady = np.absolute(fdy)
    fn = len(fy0)
    print fn,' : all points  =========='
    do_stats( fady, tag=title_p+'  ' )
    fx = np.arange(fn) + 1
    idx = np.isclose(fy1,fy0)
    ok = idx
    nok = ~idx
    gx = fx[ok]
    gy0 = fy0[ok]
    gy1 = fy1[ok]
    gdy = fdy[ok]
    gady = fady[ok]
    gn = len(gy0)
    hn = 0  # initialize
    if (gn < fn):
        print gn,' : good points  =========='
        do_stats( gady, tag=title_p+'  ' )
        hx = fx[nok]
        hy0 = fy0[nok]
        hy1 = fy1[nok]
        hdy = fdy[nok]
        hady = fady[nok]
        hn = len(hy0)
        if (hn > 0):
            print hn,' : bad points  ==========='
            do_stats( hady, tag=title_p+'  ' )
    fig, ax = plt.subplots()
    if (hn > 0):
        ax.scatter(hx,hady,color='salmon')
    ax.scatter(gx,gady,color='dodgerblue')
    ax.set_ylabel( '|y1-y0|' )
    ax.set_xlabel( 'Frequency bins' )
    ax.set_title( title_p )
    plt.show()
    print '\n===================='


def show_plot_on_mac(plotfile=None):
    assert (plotfile is not None)
    assert (os.path.isfile(plotfile))
    if (platform.system() == 'Darwin'):  # This is a Mac!
        cmd = 'open '+plotfile
        os.system(cmd)
    return


def main():
    import numpy as np
    import time as tm
    import sys

    if (sys.version_info >= (3, 0)):
        sys.stdout.write("Sorry, requires Python 2.X, not Python 3.X\n")
        sys.exit(1)

    # Get some data
    ifile = 'B1392all.tab'
    hjd_, mag_, magerr_, filts_ = np.loadtxt( ifile, unpack=True)[:4]
    
    row= np.arange(len(hjd_))
    bad = (magerr_ > 0.2) | (magerr_ <= 0.)
    nbad = np.count_nonzero(bad)
    if (nbad>0):
        print '\nFound ',nbad,' bad observations:'
        hjd__    = hjd_[bad]
        mag__   = mag_[bad]
        magerr__ = magerr_[bad]
        filts__  = np.fix(filts_[bad])
        row__ = row[bad]+1
        print '***** REJECTED DATA *****: BEGIN'
        print '         HJD      '+'    MAG  '+'   MAGERR'+'  FILTER'+'    row'
        for c1,c2,c3,c4,c5 in zip(hjd__,mag__, magerr__,filts__,row__):
            print '%18.7f %8.3f %8.3f  %6d %6d' % (c1,c2,c3,c4,c5)
        print '***** REJECTED DATA *****: END'
        
    ok = (magerr_ >= 0.0) & (magerr_ <= 0.2 )
    print '\nFound ',np.count_nonzero(ok),' good observations\n'
    hjd0    = hjd_[ok]
    mag0    = mag_[ok]
    magerr0 = magerr_[ok]
    filts0  = filts_[ok]

    prob_cut = 1.000  # 100% --> all the data
    #prob_cut = 0.5    #  50% --> half of the data
    prob = np.random.rand( len(hjd0) )
    idx = (prob <= prob_cut)
    hjd    = hjd0[idx]
    mag    = mag0[idx]
    magerr = magerr0[idx]
    filts  = filts0[idx]

    # Set pmin, dphi, and filtnams
    pmin = .2
    dphi = 0.02
    filtnams = ['u', 'g', 'r', 'i', 'z']
    #filtnams = ['u', 'g']

    # And away we go!
    time00 = tm.time()
    periods, psi_m, thresh_m = \
        psearch_py( hjd, mag, magerr, filts, filtnams, pmin, dphi )
    time01 = tm.time()

    tag = ifile+'     '+'%7.2f%%' % (prob_cut*100.0)  # extra info for plots

    # Plot HJD vs. magnitude for all filters
    plot1 = 'psearch_fig_obs.png'
    fig_obs_kjm_py( hjd, mag, filts, filtnams, tag=tag, plotfile=plot1)

    # Plot Psi vs. Frequency for all filters
    plot2 = 'psearch_fig_psi.png'
    fig_psi_kjm_py( 1/periods, psi_m, thresh_m, filtnams, tag=tag,
        plotfile=plot2 )

    # Period of the strongest peak of the combined Psi distribution
    idx = np.argmax(psi_m.sum(0))
    p_peak = periods[idx]
    
    # Plot phased light curves for all filters
    plot3 = 'psearch_fig_phi.png'
    fig_phi_kjm_py( hjd, mag, magerr, filts, filtnams, period=p_peak, \
        tag=tag, plotfile=plot3 )

    # Show the top 10 peaks of the combined Psi distribution
    #table_psi_kjm_py( xx=periods, yy=psi_m.sum(0), ee=thresh_m.sum(0), n=10 )

    # Show the plotfiles (if using a Mac)
    show_plot_on_mac(plot1)
    show_plot_on_mac(plot2)
    show_plot_on_mac(plot3)

    print '\nPeriod: %9.6f' % p_peak

    print '\nmain: %8.3f seconds [walltime for psearch_py]' % (time01-time00)
    print("\nmain: That's all folks!  :-)")


if __name__ == "__main__":
   main()
