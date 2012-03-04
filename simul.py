#!/usr/bin/python
import sys
from copy import *
from types import *
from random import *
from Genome import *
from Tree import Tree
import History
from LinRev import LinRev

if len(sys.argv) != 4:
  print "Usage: python simul.py working-dir  mean  #genes"
  sys.exit(1)

wdir = sys.argv[1]
mean = int(sys.argv[2])
n = int(sys.argv[3])

History.Genome = LinRev
T = Tree(wdir)
h = History.History(T, wdir, True)
h.simul_data(mean,n)
