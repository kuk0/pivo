import os
import sys
from DCJ import DCJ
from random import choice, randint
from copy import deepcopy

class TestDCJ(DCJ):    
    def median2 (self, g1, g2, g3):
        """
        Use ASMedian solver to find median of g1, g2, and g3.
        """
        TestDCJ.nMedianCalls += 1
        if TestDCJ.nMedianCalls > 10000:
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
        if C == 0: return [TestDCJ(g._g), ]
        while C > 1: g = choice ([h for h in g.neigh() if h.numcirc() < C]); C = g.numcirc()
        return [TestDCJ(h._g) for h in g.neigh() if h.multilin()]
    
    def median (self, g1, g2, g3):
        g = self.median2(g1, g2, g3)[0]
        t = g.median_score(g1, g2, g3)
        S = [set(), set(), set(), set()] 
        a = [[g], [], [], []]
        i = 0
        while i < 10000 and len(S[0]) <= 300: # len(a) > 0 and len(S) <= 300:
            if len(a[0]) > 0: g = a[0].pop()
            elif len(a[1]) > 0: g = a[1].pop()
            elif len(a[2]) > 0: g = a[2].pop()
            elif len(a[3]) > 0: g = a[3].pop()
            else: break
            for h in g.neigh():
                k = h.median_score(g1, g2, g3) - t
                i += 1
                if i >= 10000: break
                if k > 3: continue
                if k < 0: k = 0
                if not h in S[k]:
                    S[k].add(h)
                    a[k].append(h)
                    if len(S[0]) > 300: break
        #print len(S[0]), len(S[1]), len(S[2]), len(S[3]), i
        l = list(S[0])
        if len(l) < 300: l.extend(list(S[1])) 
        if len(l) < 300: l.extend(list(S[2])) 
        if len(l) < 300: l.extend(list(S[3]))
        if len(l) > 300: l = l[0:300]
        return l

    def rand_neigh(self):
        return TestDCJ(self.rand_neigh2(3)._g)
    
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
