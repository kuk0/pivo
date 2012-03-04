#!/usr/bin/python
import sys
import os
from Tree import Tree
import History
from Genome import *
from DCJ import DCJ
from MultiLinDCJ import MultiLinDCJ
from LinRev import LinRev
from CircRev import CircRev
from copy import copy

if len(sys.argv) != 3:
    print "Usage: python prepare.py working-dir n"
    sys.exit(1)

History.Genome = MultiLinDCJ

print len(sys.argv)
print sys.argv

wdir = sys.argv[1]
n = int(sys.argv[2])
if n < 1: n = 500

T = Tree(wdir)
h = History.History(T, wdir)

for i in xrange(0,n):
    h.rand_hist(2)
    h.write("init"+str(i))
