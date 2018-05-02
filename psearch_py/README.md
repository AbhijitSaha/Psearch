# psearch_py

This is the Python/Cython/C version of the Psearch algorithm of

    Saha, A., & Vivas, A. K. 2017, Astronomical Journal, 154, 231;
    A Hybrid Algorithm for Period Analysis from Multiband Data with
    Sparse and Irregular Sampling for Abritrary Light-curve Shapes

by Kenneth Mighell (mighell at noao dot edu).

The Python/Cython/C version of the psearch_py module should be faster 
than the original IDL version. Typical speedup factors are about 2.3X.

The Python/Cython/C version of the psearch_py module should be much faster 
than the pure Python version. Typical speedup factors are about 6.4X.

**Note**: psearch_py.py currently works with Python 2.7 .


---

## Fastest way to get psearch_py.py from GitHub

(1) Using your favorite browser, go to the Psearch GitHub repository at

    https://github.com/AbhijitSaha/Psearch

(2) Click on the green "Clone or download" button.

(3) Click on the "Download ZIP" button.

(4) Save the file (Psearch-master.zip).

(5) Go to your download directory.

(6) Move the file to where you want to work (say /tmp).  

        mv Psearch-master.zip /tmp

(7) Unpack the file:

        unzip Psearch-master.zip


---

## Test the codes

(1) Go to the new Psearch-master directory.

        cd Psearch-master

(2) **All of the IDL code exists in the Psearch-master directory.**

(3) To demo the IDL procedure Psearch.pro, type the following command

        @testrun 

from an IDL prompt.

(4) **All of the Python code exists in the Psearch-master/psearch_py 
    directory and its subdirectories.**

(5) To demo the pure Python version of the Psearch algorithm, first 
go to the psearch_py subdirectory:

        cd psearch_py

and then type the following command:

        python psearch_py.py

(6) When done, the directory should have 3 new output figure files:

        psearch_fig_obs.png
        psearch_fig_phi.png
        psearch_fig_psi.png

(7) To demo the much faster Python/Cython/C version of the Psearch algorithm,
    the first thing one must do is ...

(8) **Compile the C/Cython code** with the following command:

        ./build_so_files

If the compilation of the C/Cython code is successful, then 2 new files 
with .so extensions have been created:

        ctheta_slavec.so
        scargle_fastc.so

This build procedure has been tested using the Anaconda Python distribution.

(9) To demo the Python/Cython/C version of psearch_py, retype the command

        python psearch_py.py

The psearch_py module detects the presence of the .so files and then
preferentially uses the compiled C code to give faster (identical) results.

---

**SPEEDUPS**

The Python/Cython/C version should be much faster than the pure Python version.
Typical speedup factors are about 6.4X.

The Python/Cython/C version should be faster than the original IDL version.
Typical speedup factors are about 2.3X.

