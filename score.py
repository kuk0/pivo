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
from CircDCJ import CircDCJ
from LinRev import LinRev
from CircRev import CircRev

if len(sys.argv) != 3:
  print "Usage: python score.py working-dir file"
  sys.exit(1)

History.Genome = CircRev #DCJ #CircDCJ #Rev #MultiLinDCJ

wdir = sys.argv[1]
f = sys.argv[2]

T = Tree(wdir)
h = History.History(T, wdir)
h.read (f)
print h.score()
#h.draw_tree(5, 0.5)
