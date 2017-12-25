set nokey
set term gif animate
set output "out/out.gif"
set xrange[-0.2:0.2]
set yrange[-0.2:0.2]
n0 = 0
n1 = int(system("find out -type f -name \"[^.]*\" | wc -l"))-1
dn = 1
load "plot.plt"
