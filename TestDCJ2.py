import os
import sys
from DCJ import DCJ
from random import choice, randint
from copy import deepcopy

class TestDCJ2(DCJ):    
    def median (self, g1, g2, g3):
        """
        Use ASMedian solver to find median of g1, g2, and g3.
        """
        TestDCJ2.nMedianCalls += 1
        if TestDCJ2.nMedianCalls > 100000:
            sys.exit()

        if not (g1.multilin() and g2.multilin() and g3.multilin()):
            print g1
            print g2
            print g3
            raise TypeError("The input genomes are not multilinear")
        
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
        if self.n != M.n: print 'MS: ' + s; sys.exit()
        r = self.linearize(M)
        return r
    
    def linearize (self, g):
        """
        Returns a set of multilinear genomes close to the genome.
        """
        L, C = g.numch()
        if C == 0: return [TestDCJ2(g._g), ]
        while C > 1: g = choice ([h for h in g.neigh() if h.numcirc() < C]); C = g.numcirc()
        return [TestDCJ2(h._g) for h in g.neigh() if h.multilin()]
    
    def rand_neigh(self):
        return TestDCJ2(self.rand_neigh2(3)._g)
    
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
