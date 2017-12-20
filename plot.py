#!/usr/bin/python3

import PyGnuplot as gp

gp.c("set nokey")
gp.c("set term gif animate")
gp.c("set output \"out/out.gif\"")
gp.c("n0 = 0")
gp.c("n1 = 100")
gp.c("dn = 1")
gp.c("load \"plot.plt\"")
