pro tvctedit,level,r,g,b
; replace level of color table with r,g,b values
; for graphical overlay

	tvlct,rtab,gtab,btab,/get
	rtab(level)=r
	gtab(level)=g
	btab(level)=b
	tvlct,rtab,gtab,btab
	return
end
