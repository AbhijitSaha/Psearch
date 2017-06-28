PRO scargle,t,c,om,px,fmin=fmin,fmax=fmax,nfreq=nfreq,     $
            nu=nu,period=period,omega=omega,             $
            fap=fap,signi=signi,simsigni=simsigni,       $
            pmin=pmin,pmax=pmax,old=old,                 $
            psdpeaksort=psdpeaksort,multiple=multiple,noise=noise, $
            debug=debug,slow=slow
;+
; NAME:
;         scargle
;
;
; PURPOSE:
;         Compute the lomb-scargle periodogram of an unevenly sampled
;         lightcurve 
;
;
; CATEGORY:
;         time series analysis
;
;
; CALLING SEQUENCE:
;   scargle,t,c,om,px,fmin=fmin,fmax=fmax,nfreq=nfreq,pmin=pmin,pmax=pmax, 
;            nu=nu,period=period,omega=omega,fap=fap,signi=signi
;
; 
; INPUTS:
;         time: The times at which the time series was measured
;         rate: the corresponding count rates
;
;
; OPTIONAL INPUTS:
;         fmin,fmax: minimum and maximum frequency (NOT ANGULAR FREQ!)
;                to be used (has precedence over pmin,pmax)
;         pmin,pmax: minimum and maximum PERIOD to be used
;         omega: angular frequencies for which the PSD values are
;                desired
;         fap : false alarm probability desired
;               (see Scargle et al., p. 840, and signi
;               keyword). Default equal to 0.01 (99% significance)       
;         noise: for the normalization of the periodogram and the
;            compute of the white noise simulations. If not set, equal to
;            the variance of the original lc.   
;         multiple: number of white  noise simulations for the FAP
;            power level. Default equal to 0 (i.e., no simulations).
;         nfreq: number of independent frequencies
;      
;   
; KEYWORD PARAMETERS:
;         old : if set computing the periodogram according to J.D.Scargle
;            1982, ApJ 263, 835. If not set, computing the periodogram
;            with the fast algorithm of W.H. Press and G.B. Rybicki,
;            1989, ApJ 338, 277.
;         debug: print out debugging information if set
;         slow: if set, a much slower but less memory intensive way to
;            perform the white noise simulations is used.
;
; OUTPUTS:
;            om   : angular frequency of PSD
;            psd  : the psd-values corresponding to omega
;
;
; OPTIONAL OUTPUTS:
;            nu    : normal frequency  (nu=omega/(2*!DPI))
;            period: period corresponding to each omega
;            signi : power threshold corresponding to the given 
;                    false alarm probabilities fap and according to the
;                    desired number of independent frequencies
;            simsigni : power threshold corresponding to the given 
;                    false alarm probabilities fap according to white
;                    noise simulations
;            psdpeaksort : array with the maximum peak pro each simulation    
;
; PROCEDURE:
;         The Lomb Scargle PSD is computed according to the
;         definitions given by Scargle, 1982, ApJ, 263, 835, and Horne
;         and Baliunas, 1986, MNRAS, 302, 757. Beware of patterns and
;         clustered data points as the Horne results break down in
;         this case! Read and understand the papers and this
;         code before using it! For the fast algorithm read W.H. Press
;         and G.B. Rybicki 1989, ApJ 338, 277.
;
;         The code is still stupid in the sense that it wants normal
;         frequencies, but returns angular frequency...   
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;          Version 1.0, 1997, Joern Wilms IAAT
;          Version 1.1, 1998.09.23, JW: Do not normalize if variance is 0
;             (for computation of LSP of window function...)
;          Version 1.2, 1999.01.07, JW: force nfreq to be int
;          Version 1.3, 1999.08.05, JW: added omega keyword   
;          Version 1.4, 1999.08
;              KP: significance levels   
;              JW: pmin,pmax keywords
;          Version 1.5, 1999.08.27, JW: compute the significance levels
;               from the horne number of independent frequencies, and not from
;               nfreq
;          Version 1.6, 2000.07.27, SS and SB: added fast algorithm and FAP
;               according to white noise lc simulations.    
;          Version 1.7, 2000.07.28 JW: added debug keyword, sped up 
;               simulations by factor of four (use /slow to get old
;               behavior of the simulations)
;-
   
   ;; defaults
   IF n_elements(noise) EQ 0 THEN noise = double(sqrt((moment(c))[1]))
   IF n_elements(multiple) EQ 0 THEN multiple = 0
   IF n_elements(fap) EQ 0 THEN fap = 0.01
   
   
   ;; make times manageable (Scargle periodogram is time-shift invariant)
   time = t-t[0]

   
   ;; number of independent frequencies
   ;;  (Horne and Baliunas, eq. 13)
   n0    = n_elements(time)
   horne = long(-6.362+1.193*n0+0.00098*n0^2.)
   IF (horne LT 0) THEN horne=5
   
   IF (n_elements(nfreq) EQ 0) THEN nfreq = horne ELSE horne = nfreq

   ;; min.freq is 1/T
   IF (n_elements(fmin) EQ 0) THEN BEGIN 
       IF (n_elements(pmax) EQ 0) THEN BEGIN 
           fmin = 1.D0 / max(time)
       END ELSE BEGIN 
           fmin = 1.D0 / pmax
       END
   ENDIF 

   ;; max. freq: approx. to Nyquist frequency
   IF (n_elements(fmax) EQ 0) THEN BEGIN 
       IF (n_elements(pmin) EQ 0) THEN BEGIN 
           fmax = n0 / (2.D0*max(time))
       END ELSE BEGIN 
           fmax = 1.D0 / pmin
       END
   ENDIF 
   
   ;; if omega is not given, compute it
   IF (n_elements(omega) EQ 0) THEN BEGIN 
       om = 2.D0 * !DPI* (fmin+(fmax-fmin)*findgen(nfreq)/(nfreq-1.D0))
   END ELSE BEGIN 
       om = omega
   END
   
   ;; False Alarm Probability according to Numf
   signi = -alog(  1.D0 - ((1.D0-fap)^(1./horne))  )
      
   ;; Periodogram
   IF keyword_set(old) THEN BEGIN ;; slow version
       ;; Subtract mean from data
       cn   = c-mean(c)
       
       ;; computing the periodogram
       px = fltarr(nfreq)
       FOR i=0L,nfreq-1L DO BEGIN 
           tau = atan(total(sin(2.D0*om[i]*time))/total(cos(2.D0*om[i]*time)))
           tau = tau/(2.*om[i])
           
           co = cos(om[i]*(time-tau))
           si = sin(om[i]*(time-tau))
           
           px[i]=0.5D0*(total(cn*co)^2/total(co^2)+total(cn*si)^2/total(si^2))
       ENDFOR 
       
       ;; correct normalization
       var = (moment(cn))[1]
       IF var NE 0. THEN BEGIN 
           px = px/var
       END ELSE BEGIN
           print,'Scargle Warning: Variance is zero'
       END
       
       ;; some other nice helpers
       ;; computed here due to memory usage reasons
       nu     = om/(2.D0*!dpi)
       period = 1.D0/nu
       
       return
       
   ENDIF 
   
   ;; Ref.: W.H. Press and G.B. Rybicki, 1989, ApJ 338, 277
   
   ;; Eq. (6); s2, c2
   s2 = dblarr(nfreq) 
   c2 = dblarr(nfreq) 
   FOR i=0L,nfreq-1L DO BEGIN 
       s2[i] = total( sin(2.D0*om[i]*time) )
       c2[i] = total( cos(2.D0*om[i]*time) )
   ENDFOR
       
   ;; Eq. (2): Definition -> tan(2omtau)
   ;; --- tan(2omtau)  =  s2 / c2
   omtau = atan(s2/c2) / (2.D0)       
   
   ;; cos(tau), sin(tau)
   cosomtau= cos(omtau)  
   sinomtau= sin(omtau)
   
   ;; Eq. (7); total(cos(t-tau)^2)  and total(sin(t-tau)^2) 
   tmp = c2*cos(2.D0*omtau) + s2*sin(2.D0*omtau)
   tc2 = 0.5D0*(n0+tmp)         ; total(cos(t-tau)^2)       
   ts2 = 0.5D0*(n0-tmp)         ; total(sin(t-tau)^2) 
   
   ;; clean up
   tmp = 0. & omtau= 0.
   s2  = 0. & t2  = 0.
   
   ;; computing the periodogram for the original lc
   
   ;; Subtract mean from data
   cn = c - mean(c)

   ;; Eq. (5); sh and ch
   sh = dblarr(nfreq) 
   ch = dblarr(nfreq) 
   
   IF (multiple GT 0 AND NOT keyword_set(slow)) THEN BEGIN
       sisi=dblarr(n0,nfreq)
       coco=dblarr(n0,nfreq)
       FOR i=0,nfreq-1L DO BEGIN 
           sisi[*,i]=sin(om[i]*time)
           coco[*,i]=cos(om[i]*time)
           
           sh[i]=total(cn*sisi[*,i])
           ch[i]=total(cn*coco[*,i])
       END 
   END ELSE BEGIN 
       FOR i=0L,nfreq-1L DO BEGIN 
           sh[i] = total( cn*sin(om[i]*time) )
           ch[i] = total( cn*cos(om[i]*time) )
       ENDFOR 
   END 
   
   ;; Eq. (3)
   px = (ch*cosomtau + sh*sinomtau)^2  / tc2 + $
     (sh*cosomtau - ch*sinomtau)^2  / ts2     
   
   ;; correct normalization 
   px = 0.5D0*px/(noise^2)
   
   ;; --- RUN SIMULATIONS for multiple > 0
   IF multiple GT 0 THEN BEGIN
       IF (multiple*min(fap) LT 10) THEN BEGIN 
           message,'WARNING',/informational
           message,'Number of iterations (multiple keyword)',/informational
           message,'not large enough for false alarm probability',/informational
           message,'requested (need multiple*FAP > 10 )',/informational
       ENDIF 
       
       
       IF (keyword_set(debug)) THEN BEGIN 
           t0=systime(1)
       END 
       
       psdpeak = dblarr(multiple)
       FOR m=0L,multiple-1L DO BEGIN
           IF (keyword_set(debug)) THEN BEGIN 
               IF ((m+1) MOD 100 EQ 0) THEN BEGIN 
                   message,'...working on '+strtrim(m+1,2)+'th simul. ('+$
                     strtrim(string(format='(F8.3)',100.*m/multiple),2)+ $
                     '% done)',/informational
                   exec=systime(1)-t0
                   message,'   time expired : '+strtrim(exec,2)+'s',/info
                   message,'   time per step: '+ $
                     string(format='(F10.6)',exec/(m+1))+'s',/info
                   message,'   time remaining: '+ $
                     strtrim(exec/(m+1)*(multiple-m),2)+'s',/info
               END 
           END 

           ;; white noise simulation
           cn = randomn(anyseed,n0)*noise
           cn = cn-mean(cn) ;; .. force OBSERVED count rate to zero
           
           ;; Eq. (5); sh and ch
           IF (NOT keyword_set(slow)) THEN BEGIN 
               FOR i=0L,nfreq-1L DO BEGIN 
                   sh[i]=total(cn*sisi[*,i])
                   ch[i]=total(cn*coco[*,i])
               END 
           END ELSE BEGIN 
               FOR i=0L,nfreq-1L DO BEGIN 
                   sh[i] = total( cn * sin(om[i]*time) )
                   ch[i] = total( cn * cos(om[i]*time) )
               ENDFOR
           END 
           
           ;; Eq. (3) ; computing the periodogram for each simulation
           psdpeak[m] = max ( (ch*cosomtau + sh*sinomtau)^2 / tc2 + $
             (sh*cosomtau - ch*sinomtau)^2 / ts2 )
       ENDFOR                
           
       ;; False Alarm Probability according to simulations
       IF n_elements(psdpeak) NE 0  THEN BEGIN 
           idx = sort(psdpeak)
           ;; correct normalization 
           psdpeaksort = 0.5D0 * psdpeak[idx]/(noise^2)
           simsigni = psdpeaksort[long((1.-fap)*(multiple-1))]
       ENDIF
   END        
   
   ;; some other nice helpers
   ;; computed here due to memory usage reasons
   nu     = om/(2.D0*!dpi)
   period = 1.D0/nu
   
END 



