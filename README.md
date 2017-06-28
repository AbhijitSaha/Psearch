# Psearch
IDL routines for Period search with hybrid algorithm

IDL routines for period finding with multiband data that are sparsely sampled.  Hybridizes Lomb-Scargle and Lafler-Kinman approaches. A paper is nearly ready to submit witha fuller description -- will be cited here when complete.  This repository is intended to include all the needed routines, many of which are already in existence and in the public domain.  


HOW TO USE:

All the necessary is included. There are no dependencies beyond those contained in a regular IDL distribution.
Compatibility with GDL is untested at this time.

There are also 2 data sets as examples:  1) 392work2.sav (which are the data presented in the paper (saha & vivas, 2017) for an RRLyrae star measured in 5 bands (ready for submission at this time)  and 2) simul392data.sav (a simulation of 10 years of observations in the same 5 bands in the LSST survey might look like -- se paper for details).  These are as IDL save files. the IDL RESTORE procedure will ingest the times of observation, measured mags and estimated uncertainties and passband information into IDL variables HJD, MAG, MAGERR and FILTNAMS respectively.

included is a program testrun.pro  (execute @testrun from an IDL prompt) that ingests the data from 392work2.sav, and sets up limits on the period range and period granularity. It then calls Psearch.pro, and receives back ptest, psiall, and threshall. See the explanations at the top of teh PSsearch.pro file for details. Running PSEARCH.PRO with testrun as provided will also produce a plot that should look like Fig. 8 of the paper.  You can chnage testrun to work with the LSST simulation, or provide your own data into the necessary input variables.

The program does not supply and answer. It is up to the user to interprete the periodograms (manipulate the output variables ptest, psiall, threshall as you may wish) to surmise your answer, and fold the data to obtain a phased light curve.
