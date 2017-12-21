if(exist("n")==0 || n<0) n = n0
fname = sprintf("out/output_%05d.prof", n)

rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b)

unset colorbox
plot fname every ::5 u 2:3:(int($1)>1 ? rgb(255,0,0) : rgb(0,0,255)) w p pt 7 ps 0.5 fc rgb variable

n = n + dn
if(n < n1) reread
undefine n
