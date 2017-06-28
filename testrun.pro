;;;;  Get the test data: 
restore, '392work2.sav' 
;
;;;;  Set pmin,dphi, and filtnams, and define variable names to receive output 
;
pmin = .2
dphi = 0.02
filtnams = ['u', 'g', 'r', 'i', 'z']
psiall = 0.
threshall = 0.
ptest = 0.
; 
;;;;  run the program
;
psearch, hjd,mag,magerr,filts,filtnams, pmin, dphi, ptest, psiall, threshall
;
