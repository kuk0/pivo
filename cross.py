#!/usr/bin/python
import sys
import os
from copy import *
from types import *
from random import *
from Tree import Tree
import History
from Genome import *
from DCJ import DCJ
from MultiLinDCJ import MultiLinDCJ
from LinRev import LinRev
from CircRev import CircRev

if len(sys.argv) != 4:
  print "Usage: python cross.py working-dir n"
  sys.exit(1)

History.Genome = DCJ #CircRev #MultiLinDCJ

wdir = sys.argv[1]
n = int(sys.argv[2])
if n < 1: n = 500
m = int(sys.argv[3])
if m < 1: m = 500

T = Tree(wdir)
h = History.History(T, wdir)
l = sorted(os.listdir('data/'+wdir+'/'))
l = [s for s in l if s[0] != '.' and s != 'T' and s != 'G']
n = min(n, len(l))
l = l[0:n]

for i in xrange(0,n):
  h.read_cand(l[i])

h.neigh2()

for i in xrange(m,n):
  h.read_cand(l[i])

h.opt_neigh()
h.write()
print h
print h.score()
