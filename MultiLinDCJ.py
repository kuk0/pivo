import os
import sys
from DCJ import DCJ
from random import choice, randint, shuffle
from copy import deepcopy

class MultiLinDCJ(DCJ):
    def median_solver (self, g1, g2, g3, errors=0):
        """
        Use ASMedian solver to find median of g1, g2, and g3.
        """
        if not (g1.multilin() and g2.multilin() and g3.multilin()):
            print g1
            print g2
            print g3
            raise TypeError("The input genomes are not multilinear")
        if MultiLinDCJ.nMedianCalls > 0:
            MultiLinDCJ.nMedianCalls -= 1
            if MultiLinDCJ.nMedianCalls == 0:
                sys.exit()
        
        r = str(randint(0,99999))
        inf = "input" + r
        outf = "output" + r
        f = open (inf, 'w')
        f.write('> Genome 1\n')
        f.write(g1.__str__().replace("$", "$\n"))
        f.write('> Genome 2\n')
        f.write(g2.__str__().replace("$", "$\n"))
        f.write('> Genome 3\n')
        f.write(g3.__str__().replace("$", "$\n"))
        f.close()
        os.system('java -cp ms/ASM ASMedian '+inf+' >> '+outf)
        f = open (outf, 'r')
        s = ''
        zac = False
        for line in f:
            if not zac:
                if line[0] == '>': zac = True
                continue
            pos = line.find(':')
            if pos == -1: break
            if line[pos - 1] == '-': s += line[pos + 1:] + ' $ '
            else:                    s += line[pos + 1:] + ' @ '
        f.close()
        os.system('rm '+inf)
        os.system('rm '+outf)
        print '.'
        M = DCJ(s)
        if self.n != M.n:
            if errors == 3:
                #TODO: najlepsie odchytit exception a vypisat historiu
                #sys.exit()
                return [g1, g2, g3]
            else:
                print "ASM ERROR"
                return self.median_solver(g1, g2, g3, errors+1)
        r = self.linearize(M)
        return r
    
    def linearize (self, g):
        """
        Returns a set of multilinear genomes close to the genome.
        """
        L, C = g.numch()
        if C == 0: return [MultiLinDCJ(g._g), ]
        while C > 1: g = choice ([h for h in g.neigh() if h.numcirc() < C]); C = g.numcirc()
        return [MultiLinDCJ(h._g) for h in g.neigh() if h.multilin()]
    
    def rand_neigh(self):
        u = deepcopy (self)
        i, j = randint(0, self.n-1), randint(0, self.n-1)
        #i, j = min (i, self._g[i]), min (j, self._g[j])
        u.swap (i, self._g[i], j, self._g[j])
	while not u.multilin():
	    u = deepcopy (self)
            i, j = randint(0, self.n-1), randint(0, self.n-1)
            #i, j = min (i, self._g[i]), min (j, self._g[j])
            u.swap (i, self._g[i], j, self._g[j])
        return u #MultiLinDCJ(self.rand_neigh2(3)._g)
    
    def neigh(self):
        c = []
        for i in xrange(0, self.n):
            if self._g[i] >= i:
                for j in xrange(i + 1, self.n):
                    if self._g[j] >= j and j != self._g[i]:
                        u = deepcopy(self)
                        u.swap (i, self._g[i], j, self._g[j])
                        if u.multilin():
                            c.append(u)
                if self._g[i] != i:
                    j = self._g[i]
                    u = deepcopy(self)
                    u.swap (i, j, i, j)
                    if u.multilin():
                        c.append(u)
        return c

    def permute (self, g, pi):
        g = g.__str__().split()
        for i in xrange(0,len(g)):
            if g[i] != '$' and g[i] != '@':
                x = int(g[i])
                if x > 0:
                    g[i] = str(pi[x])
                else:
                    g[i] = str(-pi[-x])
        return MultiLinDCJ (" ".join(g))
