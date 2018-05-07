/*
 * Function: ctheta_slave_c
 * Language: C 
 *   Author: Kenneth J. Mighell
 *  Version: 0.1.5  2018MAY07
 *
 * This is a C version of the Python function 
 *   def ctheta_slave_v3_py( parray, mag, tobs ):
 * which is a Python version of the IDL procedure
 *   pro Ctheta_slave, parray, mag, tobs, theta
 * which was written by Abhijit Saha (NOAO).
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define STAND_ALONE
#undef STAND_ALONE
#ifdef STAND_ALONE
extern
void MergeSortIndex( 
  double *phiAD, 
  long *indexAL, 
  long *tmpAL, 
  long leftL, 
  long rightL
  ); 
#endif

#define CYTHON
//#undef CYTHON
#ifdef CYTHON
#include "merge5_c.c"
#endif

/*
    NAME:
        ctheta_slave_v3_py 

    INPUTS:
        parrayAD: periods (1-d array)
        sz_parrayL: size of parrayAD
        magAD: magnitudes (1-d array)
        sz_magL: size of magAD 
        tobsAD: times at which the time series was measured (e.g. HJD) (1-d array)
        sz_tobsL: size of tobsAD

    OUTPUT:
	thetaAD: Lomb-Scargle periodogram (1-d array)
	sz_thetaL: size of thetaAD

    PURPOSE:
        Compute the Lomb-Scargle periodogram of an unevenly sampled lightcurve    

    DESCRIPTION:
        Computes theta for a pre-specified array of test periods.

    NOTE BENE: 
        (1) parrayAD and thetaAD must have the same length (i.e, sz_thetaL == sz_parrayL )
        (2) magaAD and tobsAD must have the same length (i.e., sz_tobsL == sz_magL )
*/

void ctheta_slave_c
(
 double* parrayAD, 
 long sz_parrayL,  
 double* magAD,
 long sz_magL,     
 double* tobsAD,
 long sz_tobsL,
 double* thetaAD,
 long sz_thetaL
)
{
  double *ttAD = NULL;
  double *phiAD = NULL;
  double *mmAD = NULL;
  long   *indexAL = NULL;
  long   *tmpAL = NULL;

  double t0D = 0.0;
  double minD = 0.0;
  double tD = 0.0;
  double avm_kmD = 0.0;
  double sumD = 0.0;
  double denom_kmD = 0.0;
  double periodD = 0.0;
  double numerD = 0.0;

  long tL = 0;
  long j = 0;
  long k = 0;
  long mL = sz_parrayL;
  long nL = sz_magL;

  assert( NULL != parrayAD );
  assert( NULL != magAD );
  assert( NULL != tobsAD );
  assert( NULL != thetaAD );
  assert( sz_thetaL == sz_parrayL );
  assert( sz_tobsL == sz_magL );  
  assert( sz_thetaL >= 1 );
  assert( sz_magL >= 1 );

  ttAD = calloc( nL, sizeof(double)); 
  assert( NULL != ttAD );
  phiAD = calloc( nL, sizeof(double)); 
  assert( NULL != phiAD );
  mmAD = calloc( nL, sizeof(double)); 
  assert( NULL != mmAD );
  indexAL = calloc( nL, sizeof(long)); 
  assert( NULL != indexAL );
  tmpAL = calloc( nL, sizeof(long)); 
  assert( NULL != tmpAL );

  minD = tobsAD[0];
  for (j=1; j<nL; j++)
    {
      tD = tobsAD[j];
      if ( tD < minD )  
	minD = tD;
    }
  t0D = minD;
  for (j=0; j<nL; j++)
    ttAD[j] = tobsAD[j] - t0D;

  for (j=0; j<mL; j++)
    thetaAD[j] = 0.0;
  sumD = 0.0;

  for (j=0; j<nL; j++)
    sumD += magAD[j];
  avm_kmD = sumD / nL;

  sumD = 0.0;
  for (j=0; j<nL; j++)
    {
      tD = magAD[j] - avm_kmD;
      sumD += (tD * tD);
    }
  denom_kmD = sumD;

  assert( mL == sz_parrayL );
  for (k=0; k<mL; k++) 
    {
      periodD = parrayAD[k];
      assert( periodD > 0.0 );  //KJM: NUM ME VEXO?

      for (j=0; j<nL; j++) 
	{
	  tD = ttAD[j] / periodD;
	  tL = (long)tD;
	  phiAD[j] = tD - tL;
	}

      for (j=0; j<nL; j++)
	{ // inititalize arrays
	  indexAL[j] = j;
	  tmpAL[j] = 0;  //KJM: NUM ME VEXO?
	}
      MergeSortIndex(phiAD,indexAL,tmpAL,0,(nL-1)); 
      for (j=0; j<nL; j++)
	mmAD[j] = magAD[indexAL[j]];

      tD = mmAD[0] - mmAD[nL-1];
      numerD = (tD * tD);
      for (j=0; j<(nL-1); j++) 
	{
	  tD = mmAD[j+1] - mmAD[j];
	  numerD += (tD * tD);
	}
      thetaAD[k] = numerD / denom_kmD;
    }

  free(ttAD);
  free(phiAD);
  free(mmAD);
  free(indexAL);
  free(tmpAL);

  return;
}
//EOF
