/*
 * Function: scargle_fast_c
 * Language: C 
 *   Author: Kenneth J. Mighell
 *  Version: 0.1.5  2018MAY05
 *
 * This is a C version of the Python function 
 *   def scargle_fast_py( t, c, omega, nfreq ):
 * which is a Python version of the IDL procedure
 *   pro scargle,t,c,om,px,fmin=fmin,fmax=fmax,nfreq=nfreq,nu=nu, $
 *     period=period,omega=omega,fap=fap,signi=signi,simsigni=simsigni, $
 *     pmin=pmin,pmax=pmax,old=old,psdpeaksort=psdpeaksort,multiple=multiple, $
 *     noise=noise,debug=debug,slow=slow
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/*
    NAME:                                                                         
        scargle_fast_c                                                       
                                                                            
    PURPOSE:                                                                      
        Compute the Lomb-Scargle periodogram of an unevenly sampled lightcurve    
                                                                                  
    CATEGORY:                                                                     
        time series analysis                                                      
                                                                                  
    INPUTS:                                                                       
        tAD: times at which the time series was measured (e.g. HJD) (1-d array)
	sz_tL: size of tAD
        cAD: counts (corresponding count rates) (1-d array)
	sz_cL: size of cAD
        omegaAD: angular frequencies for which the PSD values are desired           
            [PSD: Fourier Power Spectral Density] (1-d array)                            
	sz_omegaL: size of omegaAD
        nfreqL: number of independent frequencies                                  
                                                                                  
    OUTPUTS:                                                                      
        pxAD: the psd-values corresponding to omega (1-d array)

    NOTE BENE: 
        (1) tAD and cAD must have the same length (i.e, sz_tL == sz_cL )
        (2) length of omegaAD must be equal to nfreqL (i.e., sz_omegaL == nfreqL )
                                                                                  
    DESCRIPTION:                                                                  
        The Lomb Scargle PSD (Fourier Power Spectral Density) is computed 
        according to the definitions given by

        (1) Scargle, J. D. 1982, ApJ, 263, 835;
            "Studies in astronomical time series analysis. II - Statistical 
            aspects of spectral analysis of unevenly spaced data"

        (2) Horne, J. H., & Baliunas, S. L. 1986, MNRAS, 302, 757;
            "A prescription for period analysis of unevenly sampled time series"

        (3) Press, W.H., & Rybicki, G.B. 1989, ApJ 338, 277;
	    "Fast algorithm for spectral analysis of unevenly sampled data"

    HERITAGE:
        http://astro.uni-tuebingen.de/software/idl/aitlib/timing/scargle.pro      
*/

void scargle_fast_c
(
 double* tAD,
 long sz_tL,
 double* cAD,
 long sz_cL,
 double* omegaAD,
 long sz_omegaL,
 long nfreqL,
 double* pxAD,
 long sz_pxL
)
{
  double* timeAD = NULL;
  double* two_timeAD = NULL;
  double* s2AD = NULL;
  double* c2AD = NULL;
  double* omAD = NULL;
  double* omtauAD = NULL;
  double* cosomtauAD = NULL;
  double* sinomtauAD = NULL;
  double* tmpAD = NULL;
  double* ts2AD = NULL;
  double* tc2AD = NULL;
  double* shAD = NULL;
  double* chAD = NULL;
  double* cnAD = NULL;
  double* omi_timeAD = NULL;

  double chpshD;
  double shmchD;
  double noiseD;
  double tD;
  double varD;
  double meanD;
  double cAD_meanD;
  double sumD;

  long mL =  sz_tL;
  long nL = nfreqL;
  long n0L;
  double omiD;
  long j;
  long i;

  assert( NULL != tAD );
  assert( NULL != cAD );
  assert( NULL != omegaAD );
  assert( NULL != pxAD );
  assert( sz_tL == sz_cL );
  assert( sz_omegaL == nfreqL );
  assert( sz_pxL == nfreqL );

  // zero out the output vector
  for (j=0; j<nL; j++)
    pxAD[j] = 0.0;

  { // compute noise
    long ddofL = 0L; /* Delta Degrees Of Freedom: zero used by numpy */
    //KJM:       ^---- ***** WARNING ***** HACK A FURBALL *****
    sumD = 0.0;
    for (j=0; j<mL; j++)
      sumD += cAD[j];
    meanD = sumD / mL;
    cAD_meanD = meanD;
    varD = 0.0;
    for (j=0; j<mL; j++)
      { 
	tD = cAD[j] - meanD;
	varD += tD * tD;
      }
    varD /= (mL-ddofL);
    noiseD = sqrt(varD);
  }

  // make times manageable (Scargle periodogram is time-shift invariant)
  timeAD = calloc( mL, sizeof(double)); 
  assert( NULL != timeAD );
  tD = tAD[0];
  for (j=0; j<mL; j++)
    timeAD[j] = tAD[j] - tD;

  n0L = mL; /* = sz_tL */;
  assert( n0L > 0L );

  omAD = &omegaAD[0];  // alias

  //
  // ===== PERIODOGRAM: FAST VERSION ========================================
  // Reference:
  //   Press, W.H., & Rybicki, G.B. 1989, ApJ 338, 277;
  //   "Fast algorithm for spectral analysis of unevenly sampled data"
  //     

  // Eq. (6): s2, c2
  s2AD = calloc( nL, sizeof(double)); 
  assert( NULL != s2AD );
  c2AD = calloc( nL, sizeof(double)); 
  assert( NULL != c2AD );
  two_timeAD = calloc( nL, sizeof(double)); 
  assert( NULL != two_timeAD );
  for (j=0; j<mL; j++)
     two_timeAD[j] = 2.0 * timeAD[j];
  for (i=0; i<nL; i++)
    {
      omiD = omAD[i];
      sumD = 0.0;
      for (j=0; j<mL; j++)
	sumD += sin(two_timeAD[j]*omiD);
      s2AD[i] = sumD;
      sumD = 0.0;
      for (j=0; j<mL; j++)
	sumD += cos(two_timeAD[j]*omiD);
      c2AD[i] = sumD;
    }
  free(two_timeAD);

  // Eq. (2): tan(2omtau) =  s2 / c2
  omtauAD = calloc( nL, sizeof(double)); 
  assert( NULL != omtauAD );
  for (j=0; j<nL; j++)
    omtauAD[j] = atan(s2AD[j]/c2AD[j])/2.0;
  cosomtauAD = calloc( nL, sizeof(double)); 
  assert( NULL != cosomtauAD );
  sinomtauAD = calloc( nL, sizeof(double)); 
  assert( NULL != sinomtauAD );
  for (j=0; j<nL; j++)
    {
      cosomtauAD[j] = cos(omtauAD[j]);
      sinomtauAD[j] = sin(omtauAD[j]);
    }

  // Eq. (7)
  tmpAD = calloc( nL, sizeof(double)); 
  assert( NULL != tmpAD );
  for (j=0; j<nL; j++)
    { 
      tD = 2.0*omtauAD[j];
      tmpAD[j] = c2AD[j]*cos(tD) + s2AD[j]*sin(tD);
    }
  tc2AD = calloc( nL, sizeof(double)); 
  assert( NULL != tc2AD );
  for (j=0; j<nL; j++)
    tc2AD[j] = 0.5*(n0L+tmpAD[j]);
  ts2AD = calloc( nL, sizeof(double)); 
  assert( NULL != ts2AD );
  for (j=0; j<nL; j++)
    ts2AD[j] = 0.5*(n0L-tmpAD[j]);

  free(tmpAD);
  free(omtauAD);
  free(s2AD);
  free(c2AD);  

  //
  // computing the periodogram for the original light curve
  //

  // subtract mean from data
  cnAD = calloc( nL, sizeof(double)); 
  assert( NULL != cnAD );
  for (j=0; j<mL; j++)
    cnAD[j] = cAD[j] - cAD_meanD;

  // Eq. (5); sh and ch
  shAD = calloc( nL, sizeof(double)); 
  assert( NULL != shAD );
  chAD = calloc( nL, sizeof(double)); 
  assert( NULL != chAD );
  omi_timeAD = calloc( nL, sizeof(double)); 
  assert( NULL != omi_timeAD );
  for (i=0; i<nL; i++)
    {
      tD = omAD[i];
      for (j=0; j<mL; j++)
	omi_timeAD[j] = tD * timeAD[j];
      sumD = 0.0;
      for (j=0; j<mL; j++)
        sumD += cnAD[j] * sin(omi_timeAD[j]);
      shAD[i] = sumD;
      sumD = 0.0;
      for (j=0; j<mL; j++)
        sumD += cnAD[j] * cos(omi_timeAD[j]);
      chAD[i] = sumD;
    }
  free(cnAD);
  free(omi_timeAD);

  // Eq. (3)
  tD = 0.5/(noiseD*noiseD);
  for (j=0; j<nL; j++)
    {
      chpshD = chAD[j]*cosomtauAD[j] + shAD[j]*sinomtauAD[j];
      shmchD = shAD[j]*cosomtauAD[j] - chAD[j]*sinomtauAD[j];
      pxAD[j] = ((chpshD*chpshD)/tc2AD[j]) +
	        ((shmchD*shmchD)/ts2AD[j]); 
      pxAD[j] *= tD;  /* correct normalization */
    }

  // clean up
  free(chAD);
  free(shAD);
  free(cosomtauAD);
  free(sinomtauAD);
  free(tc2AD);
  free(ts2AD);
  free(timeAD);

  return;
}
//EOF
