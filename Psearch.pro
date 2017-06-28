
;;;;  INPUTS :   
;;;;;  hjd      =  the time array. Can be in any units. Output period array PTEST (see below) will be in the same units
;;;;;  mag      =  all measured signal values in an array matched to time array hjd. All passbands included.
;;;;;  errmag   =  array of uncertainties associated with mag
;;;;;  filts    =  array of integer values from 0 to n-1 that code the passband or channel, where n = # of bands present 
;;;;;  filtnams =  string array containing character names corresponding to coded filts values.
;;;;;                e.g. if you have 5 bands labeled u,g,r,i,z with filts values 0,1,2,3,4,5 respectively, filtnams would be set by: 
;;;;;                 FILTNAMS = ['u', 'g', 'r', 'i', 'z']
;;;;;  pmin     =  minimum value of period to be tested 
;;;;;  dphi     =  maximum change in relative phase between first and last epoch to be permitted when stepping to next test period
;;;;;
;;;;  OUTPUTS :  
;;;;;   ptest      =  1-D array with N dimensions of periods for which the periodograms are computed. 
;;;;;                 It is the same for ALL bands/channels
;;;;;   psi_m      =  M x N array of the Psi periodogram, where M is the number of bands/channels in the input array FILTNAMS
;;;;;   thresh_m   =  M x N array containing threshold values of Psi at each period and band for assessing significance for psi_m
;;
pro Psearch, hjd, mag, magerr, filts, filtnams, pmin, dphi, ptest, psi_m, thresh_m
;
 if( n_params(0) lt 7) then begin
   print, ' USAGE:  Psearch, hjd, mag, magerr, filts, filtnams, pmin, dphi, [ptest, psi_m, thresh_m] '
   return
 endif
;
;
;
psopen = 0
reply = ' '
print, ' Do you want plots to go to a Postscript file ?  [Y/N] or[y/n] '
read, reply
if(reply eq 'Y' or reply eq 'y') then begin
  set_plot, 'ps'
  device, /portrait, xo=0.5, xs=7.5, yo=0.5, ys=10.0, file='plot.ps', /inches, /color
;  loadct, 13
  tvctedit, 128, 255, 128, 0
  psopen = 1 
endif else begin
  set_plot, 'x'
endelse 
;
nfilts = n_elements(FILTNAMS)
;
;
nppp = max([nfilts+1, 6])
!p.multi = [0, 1, nppp]
;
;
;
psiacc = 0.
confacc = 0.
;
for i = 0, nfilts-1L do begin
  fwant = i
  periodpsi2, hjd, mag, magerr, filts, pmin,dphi,fwant, x, fy, theta, psi, conf
;
  if(i eq 0) then begin
     psi_m = fltarr(nfilts, n_elements(x))   ; define the output arrays -- needs size of period array from 1st band call to periodpsi2
     thresh_m = psi_m
  endif 
;
  psi_m(i,*) = psi
  thresh_m(i,*)  = conf
  ymax = max(psi)*1.2
  yrng = [0., ymax]
  xmax = 1./pmin
  xrng = [0., xmax]
  plot, 1./x, psi, xth=3, yth=3, chars=1.3, charth=2, xr=xrng, yr=yrng, xsty=1, ysty=1
  oplot, 1/x, conf, color=128
  xyouts, 0.95*xmax, 0.8*ymax, filtnams(i), chars=1.5, charth=2
  xyouts, -0.058*xmax, 0.45*ymax, '!7W!X',  orient=90, chars=1.5, charth=2  
;
    psiacc = psi + psiacc
    confacc = conf + confacc
endfor
;
  ymax = max(psiacc)*1.2
  yrng = [0., ymax]
  plot, 1./x, psiacc,  xth=3, yth=3, chars=1.3, charth=2, xr=xrng, yr=yrng, ysty=1, xsty=1
  oplot, 1./x, confacc, color=128
  xyouts, 0.9*xmax, 0.8*ymax, 'ALL', chars=1.5, charth=2
  xyouts, 0.4*xmax, -0.4*ymax, 'Frequency (days!E-1!N)', chars=1.3, charth=2
  xyouts, -0.058*xmax, 0.45*ymax, '!7W!X', orient=90, chars=1.5, charth=2   
;
if(psopen eq 1) then begin
  device, /close
  set_plot, 'x'
endif
;
!p.multi = 0
;
ptest = x
;
return
end

