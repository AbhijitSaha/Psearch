pro wcols, fname, c1, c2, c3, c4
;
nrows = n_elements(c1)
;
fff = '(2x, F16.7, 3x, F8.3, 3x, F8.3, 3x, F8.3)'
;
openw, u101, fname, /get_lun
;
 for i = 0, nrows-1L do begin 
     printf, u101, c1(i), c2(i), c3(i), c4(i), format=fff
 endfor
 close, u101
 free_lun, u101
;
return
end
   