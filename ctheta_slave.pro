pro Ctheta_slave, parray, mag, tobs, theta
;
;
;;;; this program computes theta for a pre-specified array of test periods
;
t0 = min(tobs)
tlast = max(tobs)
;
tt = tobs - t0
;
theta = 0.*parray
;
;; Loop over all periods…
;
for k = 0L, n_elements(Parray)-1L do begin
     period = parray(k)
     phi = tt/period
     nphi = fix(phi)
     phi = phi - nphi 
     ss = sort(phi)
     phi = phi(ss)
     mm  = mag(ss)
     avm = total(mm)/n_elements(mm)
     denom = total( (mm-avm)^2 ) 
     mmplus = [mm(1:*), mm(0)]
     numer = total( (mmplus - mm)^2 )  
     theta(k)  = numer/denom
endfor
;
;
return
end    