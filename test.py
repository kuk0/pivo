#!/usr/bin/python
import sys
import os
#import argparse
from Tree import Tree
import History
from Genome import *
from DCJ import DCJ
from MultiLinDCJ import MultiLinDCJ
from CircDCJ import CircDCJ
from LinRev import LinRev
from CircRev import CircRev
from copy import copy
from random import randint

wdir = 't100n20m'
History.Genome = CircDCJ

n, mean = 30, 8
suf = str(randint(0,999))
F = open('tst3-'+suf, 'a')
F.write("""
# 100 vrcholov; 30 markerov, priemerne 8 reverzov na hranu (2 hrany z korena su 1 hrana)
# 1 historiu vytvorime nahodne a na kazdej spustime 10x steinerizaciu;
# zaciname z historie, ktoru vyrobime tak, ze z nahodnych 3 potomkov spocitame median
# nasledne po kazdej stein. skusime better_neighbours a all_neighbours, ci to nevieme zlepsit
""")
stat = []
for i in xrange(0,10):
    T = Tree(wdir)
    h = History.History(T, wdir, True)
    h.x = -1
    h.chunks = 20
    h.simul_data(mean,n)
    a = h.score()
    #print h.score() 
    #print h
    

    Genome.median = Genome.single_median

    F.write(str(a)+'   ')
    for k in xrange(0,10):
    	h.rand_hist(1)
    	h.cand()
    	#h.fill_in_desc()
    	h.init_hist()
    	h.opt_neigh()
    
    	h.local_opt = h.local_opt_steinerization
    	h.local_opt()
    	print h.score(), '------------------------------------------------'
    	b = h.score()
    	#h.local_opt = h.local_opt_neighbours 
    	h.local_opt = h.local_opt_better_neighbours
    	h.local_opt()
    	print h.score(), '================================================'
    	c = h.score()

    	h.local_opt = h.local_opt_neighbours 
    	h.local_opt()
    	print h.score(), '================================================'
    	d = h.score()
    	#F.write(str(b)+' '+str(c)+' '+'   ') #str(d)+'   ')
    	F.write(str(b)+' '+str(c)+' '+str(d)+'   ')
    F.write('\n')
#   stat.append([a,b,c])
F.close()
#aa, bb, cc = zip(*stat)
#nn = len(stat)
#print sum(aa)/float(nn), sum(bb)/float(nn), sum(cc)/float(nn)

