pro rcols, fname, c1, c2, c3, c4
;
c1 = fltarr(10000L)
c2 = c1
c3 = c1
c4 = c1
;
fff = '(2x, F16.7, 3x, F8.3, 3x, F8.3, 3x, F8.3)'
;
openr, u201, fname, /get_lun
;
 i = 0
 while ~EOF(u201) do begin 
     readf, u201, cc1, cc2, cc3, cc4, format=fff
     c1(i) = cc1 & c2(i) = cc2 & c3(i) = cc3 & c4(i) = cc4
     i = i + 1
 endwhile
 close, u201
 free_lun, u201
;
 nrows = i
 c1 = c1(0:nrows-1L)
 c2 = c2(0:nrows-1L)
 c3 = c3(0:nrows-1L)
 c4 = c4(0:nrows-1L)
;
;
return
end
   