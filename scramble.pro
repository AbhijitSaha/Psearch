pro scramble, inarr, scrambarr, pickarr
;
ns = n_elements(inarr)
x = randomu(seed, ns)*ns
x = long(x+0.5)
s = sort(x)
scrambarr = inarr(s)
pickarr = inarr(x)
;
return
end


