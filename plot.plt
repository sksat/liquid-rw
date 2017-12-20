if(exist("n")==0 || n<0) n = n0
fname = sprintf("out/output_%05d.prof", n)
plot fname every ::4 using 2:3

n = n + dn
if(n < n1) reread
undefine n
