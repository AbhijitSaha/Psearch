# Psearch
IDL routines for Period search with hybrid algorithm. 

IDL routines for period finding with multiband data that are sparsely sampled.  Hybridizes Lomb-Scargle and Lafler-Kinman approaches. A paper is nearly ready to submit with a fuller description (Saha & Vivas, 2017), henceforth referred to as "the paper".  To be self-complete, this repository includes code provided by others, but within the public domain.


HOW TO USE:

All the necessary code is included. There are no dependencies beyond those contained in a regular IDL distribution.
Compatibility with GDL is untested at this time.

There are also 2 data sets as examples:  1) 392work2.sav (which are the data presented in the paper (saha & vivas, 2017) for an RRLyrae star measured in 5 bands (ready for submission at this time)  and 2) simul392data.sav (a simulation of what 10 years of observations in the same 5 bands for the same object in the LSST survey might look like -- see paper for details).  These data are as IDL save files: the IDL RESTORE procedure will ingest the times of observation, measured mags, estimated uncertainties, passband index,  and passband information into IDL variables HJD, MAG, MAGERR, FILTS and FILTNAMS respectively.

included is a program testrun.pro  (execute @testrun from an IDL prompt) that serves as an example of how to use the code. Running it  ingests the data from 392work2.sav, and sets up limits on the period range and period granularity. It then calls Psearch.pro, and receives back ptest, psiall, and threshall. See the explanations at the top of the Psearch.pro file for details. Running PSEARCH.PRO with testrun as provided will also produce a plot that should look like Fig. 8 of the paper. Running testrun.pro as is should be a good test that the code is successfully installed.  You can change testrun to work with the LSST simulation, or provide your own data into the necessary input variables.

The program does not supply and answer. It is up to the user to interprete the periodograms (manipulate the output variables ptest, psiall and threshall as you may wish) to surmise your answer, and fold the data to obtain a phased light curve.
