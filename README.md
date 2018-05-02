# Psearch
IDL routines and a Python/Cython/C module for Period search with hybrid algorithm. 

IDL routines and the Python/Cython/C module for period finding with multiband data that are sparsely sampled.  Hybridizes Lomb-Scargle and Lafler-Kinman approaches. An article with a fuller description has been published:

    Saha, A., & Vivas, A. K. 2017, Astronomical Journal, 154, 231;
    A HYBRID ALGORITHM FOR PERIOD ANALYSIS FROM MULTI-BAND DATA WITH
    SPARSE AND IRREGULAR SAMPLING FOR ARBITRARY LIGHT CURVE SHAPES

Henceforth, this article will be referred to as "the paper".

If you can not get the article, you can get the preprint at

   https://arxiv.org/abs/1709.10156

To be self-complete, this repository includes code provided by others, but within the public domain.

---

## IDL version: HOW TO USE

**NOTE: All of the IDL code exists in the main directory.**

All the necessary code is included. There are no dependencies beyond those contained in a regular IDL distribution.
For GDL compatibility, the IDL library routines MOMENT.pro and MEAN.pro have also been provided: the code appears to work with GDL Version 0.9.4.

There are also 2 data sets as examples:  1) 392work2.tab (which are the data presented in the paper (saha & vivas, 2017) for an RRLyrae star measured in 5 bands (ready for submission at this time)  and 2) simul392data.tab (a simulation of what 10 years of observations in the same 5 bands for the same object in the LSST survey might look like -- see paper for details).  These data are tabular data: the IDL  procedures testrun and test2run(provided, see below) will ingest the times of observation, measured mags, estimated uncertainties, passband index,  and passband information into IDL variables HJD, MAG, MAGERR, FILTS respectively.

included is a program testrun.pro  (execute @testrun from an IDL prompt) that serves as an example of how to use the code. Running it  ingests the data from 392work2.tab, and sets up limits on the period range and period granularity. It then calls Psearch.pro, and receives back ptest, psiall, and threshall. See the explanations at the top of the Psearch.pro file for details. Running PSEARCH.PRO with testrun as provided will also produce a plot that should look like Fig. 8 of the paper. Running testrun.pro as is should be a good test that the code is successfully installed.  You can change testrun or testrun2 to work with your own data, by replacing the data input specs, so that  it reads your data into the necessary input variables (HJD, MAG, MAGERR and FILTS), and provides a translation from the integer code in FILTS to real filter names, by altering the defintion for FILTNAMS in the testrun (or test2run) code.

The program does not supply an answer. For now is up to the user to interprete the periodograms (manipulate the output variables ptest, psiall and threshall as you may wish) to surmise your answer, and fold the data to obtain a phased light curve. You could add to or modify Psearch.pro, or call a custom routine either from the end of Psearch.pro, or from the IDL prompt when the program returns (the main level retains the output periodograms).


A general tool for helping with the examination and analysis of the periodograms is under construction.

---

## Python/Cython/C version: HOW TO USE

**NOTE: All of the Python/Cython/C code exists in the psearch_py subdirectory and its subdirectories.**

First, go to the psearch_py subdirectory.

The REAME.md file shows how to get the psearch_py.py from GitHub and then how to test the code.

The test input data file (B1392all.tab) is provided.

The main function of the psearch_py module is a built-in demo to test the module. To run the demo, type the following command

    python psearch_py.py

The table presented at the end of the demo indicates that the "best" estimate for the period (based on all of the filters) of the RR Lyrae star B1-392 is 0.5016247 +-  0.0000036 days which has a frequency value of 1.993522.  The phased light curve figure (psearch_fig_phi.png), which is created by running the demo, uses a period value of 0.501625.

If you want to use the much faster Python/Cython/C version, please follow the instructions in the README.md file
in the psearch_py subdirectory.

