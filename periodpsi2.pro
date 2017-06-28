pro periodpsi2, HJD, MAG, MAGERR, FILTS, minper, dphi, fwant, x, fy, theta, psi, conf
;
;;INPUTS
;;;; HJD, MAG, MAGERR, FILTS  — co-aligned input data point arrays:  FILTS carries an 
;;;; integeridentifier indicating passband.
;;;
;;;; minper —>  minimum period to explore
;;;; dphi  —> maximum phase change between any two data points resulting from one step 
;;;;  in frequency or Period 
;;;;  fwant —>  integer value corresponding to desired passband from among values in FILTS 
;;;;    array.
;;
;;
;;OUTPUTS
;;;; x  —> period array for which periodograms are computed  
;;;; fy —>  Lomb-Scargle periodogram (co-aligned with x)
;;;; theta —> Lafler-Kinman periodogram (co-aligned with x)
;;;; psi —>  psi periodogram (co-aligned with x) 
;;;; conf —>  simulated PSI periodogram for a non-periodic variable with amplitude
;;;; and noise mimicking real source PLUS of an unvarying object with noise mimicking
;;;; source to be compared against PSI periodogram of source to estimate confidence 
;;
;; 
HJD = double(HJD)
t0 = min(HJD)
tmax = max(HJD)
tspan = tmax - t0
;
;;;print, ' Enter filter index to use (a coded integer value): '
;;;read, fwant
ok = where(filts eq fwant and magerr le 0.2, nok)
;
tr = HJD(ok)
yr = MAG(ok)
er = MAGERR(ok)*randomn(seed, nok)
amp = max(yr)-min(yr)
scramble, yr, zr, junk
;
sss = sort(tr)
tr = tr(sss)
yr = yr(sss)
er = er(sss) 
;
;
;;;print, 'Enter minimum period (days) desired :'
;;;read, minper
maxfreq = 1.d0/minper
minfreq = 2.d0/tspan
;
;;;print, ' Enter phase accuracy desired: '
;;;read, dphi
deltafreq = dphi/tspan
nfreq = long( (maxfreq-minfreq)/deltafreq )
print, 'No. of frequency samples = ', nfreq
farray = minfreq + dindgen(nfreq)*deltafreq
x = 1./farray
print, ' Min and max periods: ', min(x), max(x)
;
omega = farray*2.d0*!dpi
scargle, tr, yr, om, fy, omega=omega, nfreq=nfreq
scargle, tr, er, om, fe, omega=omega, nfreq=nfreq   ; \rho^{prime}
scargle, tr, zr, om, fz, omega=omega, nfreq=nfreq
pout = 2.d0*!dpi/om
;plot, pout-x, yr=[-10, 10]
;
ctheta_slave, x, yr, tr, theta
ctheta_slave, x, er, tr, thetaerr  ; 
ctheta_slave, x, zr, tr, thetaz    ; \eta 
;
psi = 2.*fy/theta
;
conf1 = 2.*fe/thetaerr
conf1 = conf1*total(psi)/total(conf1)
conf2 = 2.*fz/thetaz
conf2 = conf2*total(psi)/total(conf2)
;
conf = conf1 + conf2
;
;;pause
return
end

