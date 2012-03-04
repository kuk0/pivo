#!/usr/bin/python
import sys
import os
from Tree import Tree
import History
from Genome import *
from DCJ import DCJ
from MultiLinDCJ import MultiLinDCJ
from TestDCJ2 import TestDCJ2
from LinRev import LinRev
from CircRev import CircRev
from copy import copy

if len(sys.argv) != 4:
    print "Usage: python genm.py working-dir n p"
    sys.exit(1)

History.Genome = MultiLinDCJ #TestDCJ2
#TestDCJ2.nMedianCalls = 0

print len(sys.argv)
print sys.argv

wdir = sys.argv[1]
n = int(sys.argv[2])
if n < 1: n = 500
p = int(sys.argv[3])
if p < 0: p = 10

T = Tree(wdir)
h = History.History(T, wdir)
h.rand_hist(2)
suf = "" #'-'+str(n)+'-'+str(p)+'.txt'
init = 0

#h.rand_hist(2)
#h.local_opt = h.local_opt_desc
#h.local_opt()
h.local_opt = h.local_opt_medians
for i in xrange(0,n):
    if p==0:
        h.rand_hist(2)
    elif i%p == 0:
        l = os.listdir('data/'+wdir+'/')
        l = [s for s in l if '0' <= s[0] <= '9']
        print l
        if len(l) > 0:
            best = min(l)
            h.read (best)
            print best
            print h.score()
        else: h.rand_hist(2)
        h.mutate()
    h.local_opt()
    h.write()
