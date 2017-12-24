set nokey
set term gif animate
set output "out/out.gif"
set xrange[-20:20]
set yrange[-20:20]
n0 = 0
n1 = 319
dn = 1
load "plot.plt"
