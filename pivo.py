#!/usr/bin/python
import sys
import os
import argparse
from Tree import Tree
import History
from Genome import *
from DCJ import DCJ
from MultiLinDCJ import MultiLinDCJ
from CircDCJ import CircDCJ
from LinRev import LinRev
from CircRev import CircRev
from copy import copy


def chain (g2, h2): # we assume that the last
    g, h = copy(g2), copy(h2)
    d = g[-1] - 1
    h = [x+d for x in h]
    g.extend(h[1:])
    return g

def concat (g2, h2):
    g, h = copy(g2), copy(h2)
    d = g[-1]
    h = [x+d for x in h]
    g.extend(h)
    return g

def insert (f, g, h): # ZLEEE
    return concat(concat(f,g), h)

def LR (g):
    return LinRev(" ".join([str(x) for x in g])+" $")

#H, SH = [1, 3, 2, 4], [1, 7, 3, 5, 4, 6, 2, 8]
#H1 = [1, 3]; H2 = [2, 4]

#g = CircRev (id_cgenome(23))
#h = CircRev (random_circular_genome(23))
#print g.dist (h)
#print g.dist2 (h)

#g1 = MultiLinDCJ ("-1 -7 4 -2 -6 9 -3 5 -8 $") 
#g2 = MultiLinDCJ ("1 7 8 9 2 3 4 5 6 $")
#g3 = MultiLinDCJ ("2 9 4 6 -3 8 1 5 -7 $")
#sys.exit();

parser = argparse.ArgumentParser(prog='pivo')
#parser.add_argument('integers', metavar='N', type=int, nargs='+', help='an integer for the accumulator')
#parser.add_argument('--sum', dest='accumulate', action='store_const',
#                   const=sum, default=max,
#                   help='sum the integers (default: find the max)')
parser.add_argument('wdir',
                    help = 'the working directory')
parser.add_argument('n',
                    type=int,
                    help = 'the number of reconstructed histories')
parser.add_argument('p',
                    type=int,
                    help = 'how often we start from the best history')
parser.add_argument('-m', '--model', 
                    choices = ['dcj', 'lin-dcj', 'circ-dcj', 'rev', 'circ-rev'],
                    default = 'dcj',
                    #metavar = 'M',
                    help = 'the genome model used; either dcj, multilinear dcj,'
                          +'circular dcj, reversal, or circular reversal '
                          +'(dcj is the default)')
parser.add_argument('-s', '--strategy',
                    type = int, 
                    choices = [0, 1, 2, 3, 4, 5, 6, 7],
                    default = 2,
                    #metavar = 'M',
                    help = '0 - steinerization\r'
                          +'1 - medians\n'
                          +'2 - medians with neighbours\n'
                          +'3 - all neighbours\n'
                          +'4 - better neighbours\n'
                          +'5 - median sample\n')
parser.add_argument('-x',
                    type=int,
                    default = -1,
                    help = 'start [when partitioning] ')
parser.add_argument('-y',
                    type=int,
                    default = 20,
                    help = 'chunks size [when partitioning] ')
parser.add_argument('-b',
                    type=int,
                    default = 0,
                    help = 'start [when partitioning] ')

args = parser.parse_args()

print args

wdir = args.wdir
n = args.n
if n < 1: n = 500
p = args.p
if p < 0: p = 10

print 'model:', args.model
if args.model == 'dcj':
    History.Genome = DCJ
elif args.model == 'lin-dcj':
    History.Genome = MultiLinDCJ
elif args.model == 'circ-dcj':
    History.Genome = CircDCJ
elif args.model == 'rev':
    History.Genome = LinRev
elif args.model == 'circ-rev':
    History.Genome = CircRev

# 0 = po 1000
History.Genome.nMedianCalls = -1
    
T = Tree(wdir)
h = History.History(T, wdir)
#h.rand_hist(2)

s = args.strategy
if s == 0:
    h.local_opt = h.local_opt_steinerization
    Genome.median = Genome.single_median
    s = 'steinerization'
elif s == 1:
    h.local_opt = h.local_opt_medians
    Genome.median = Genome.single_median
    s = 'medians'
elif s == 2:
    h.local_opt = h.local_opt_medians
    Genome.median = Genome.median_cloud
    s = 'medians with neighbours'
elif s == 3:
    h.local_opt = h.local_opt_neighbours
    s = 'all neighbours'
elif s == 4:
    h.local_opt = h.local_opt_better_neighbours
    s = 'better neighbours'
elif s == 5:
    h.local_opt = h.local_opt_medians
    Genome.median = Genome.median_sample
    s = 'random sample of medians'
elif s == 6:
    h.local_opt = h.local_opt_neighbours2
    s = 'all neighbours'
elif s == 7:
    h.local_opt = h.local_opt_circ
    s = 'get rid of extra chromosomes'
print 'strategy:', s

h.x = args.x
h.chunks = args.y

#l = os.listdir('data/'+wdir+'/')
#l = sorted ([s for s in l if '0' <= s[0] <= '9'])
#l = l[args.b : args.b+20]
#for hh in l:
    #h.read(hh)
    #s = h.score()
    #h.local_opt()
    #s2 = h.score()
    #if s2 < s:
       #h.write()
       #print "YAY"
       #print '-------------------------------------------------------------------'
#sys.exit()

for i in xrange(0,n):
    if p==0:
        h.rand_hist(2)
    elif i%p == 0:
        l = os.listdir('data/'+wdir+'/')
        l = [s for s in l if '0' <= s[0] <= '9']
        #print l
        if len(l) > 0:
            best = min(l)
            h.read (best)
            print best
            #h.mutate()
            print h.score()
        else: h.rand_hist(2)
    h.local_opt()
    #h.write()
    #print "TABU:", [len(x) for x in h.TABU[T.n:T.N]]
    
#h.local_opt_tabu()
#h.write()