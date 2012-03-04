import os
import sys
from DCJ import DCJ
from random import choice, randint
from copy import deepcopy

class CircDCJ(DCJ):
    #def __init__ (self, s):
        #g = DCJ(s)
        #if not g.circular():
            #print "Warning:", g, "not circular"
        #self = self.circularize(g)

    def median (self, g1, g2, g3, errors=0): # pomocou median solvera
        """
        Use BIOMedian solver to find median of g1, g2, and g3.

        BIOMedian solves the median problem for multicircular DCJ model, 
        however, the input genomes have to be circular.
        """
        if not (g1.circular() and g2.circular() and g3.circular()):
            print g1
            print g2
            print g3
            raise TypeError("The input genomes are not circular")
        r = str(randint(0,99999))
        inf = "input" + r
        outf = "input" + r + ".rst"
        f = open (inf, 'w')
        f.write('> Genome 1\n')
        f.write('C: ' + g1.__str__()[:-2] + '\n')
        f.write('> Genome 2\n')
        f.write('C: ' + g2.__str__()[:-2] + '\n')
        f.write('> Genome 3\n')
        f.write('C: ' + g3.__str__()[:-2] + '\n')
        f.close()
        #print g1
        #print g2
        #print g3
        os.system('java -cp ms/BIO BIOMedian '+inf+' >> /dev/null')
        f = open (outf, 'r')
        s = ''
        for line in f:
            if line[0] == '>': continue
            if line[0] == '#': break
            pos = line.find(':')
            s += line[pos + 1:] + ' @ '
        f.close()
        os.system('rm '+inf)
        os.system('rm '+outf)
        print '.',
        M = DCJ(s)
        if self.n != M.n:
            if errors == 3:
                #TODO: najlepsie odchytit exception a vypisat historiu
                #sys.exit()
                return [g1, g2, g3]
            else:
                print "BIO ERROR"
                return self.median_solver(g1, g2, g3, errors+1)
	r = self.circularize(M)
        return [r,]

    def circularize (self, g):
        """
        Returns a set of circular genomes close to the genome.
        """
        L, C = g.numch()
        if L == 0 and C == 1: return CircDCJ(g._g)
        while L > 0:
            h = g.rand_neigh()
            while h.numlin() >= L:
                h = g.rand_neigh()
            g = h
            L = g.numlin()
        while C > 1:
            h = g.rand_neigh()
            while h.numcirc() >= C:
                h = g.rand_neigh()
            g = h
            C = g.numcirc()
        return CircDCJ(h._g)

    def circularize2 (self):
        """
        Returns a set of circular genomes close to the genome.
        """
        print self
        L, C = self.numch()
        if L == 0 and C == 1: return # [self, ]
        while L > 0:
            h = self.rand_neigh()
            while h.numlin() >= L:
                h = self.rand_neigh()
            self = h
            L = self.numlin()
            print self
        while C > 1:
            h = self.rand_neigh()
            while h.numcirc() >= C:
                h = self.rand_neigh()
            self = h
            C = self.numcirc()
            print self
        print "vysledok:", self
        #return [CircDCJ(h._self)]

    def rand_neigh(self):
        g = DCJ(self._g)
        h = g.rand_neigh()
        #while not h.circular():
            #h = g.rand_neigh()
        return CircDCJ(h._g)
    
    def neigh(self):
        c = []
        for i in xrange(0, self.n):
            if self._g[i] > i:
                for j in xrange(i + 1, self.n):
                    if self._g[j] > j:
                        u = deepcopy(self)
                        u.swap (i, self._g[i], j, self._g[j])
                        if u.circular():
                            c.append(u)
        return c
    
    def betterer_neigh (self, z):
        g = DCJ(self._g)
        z = [DCJ(g._g) for g in z]
        neigh = g.betterer_neigh(z)
        if len(neigh) == 0:
            print '-----'
            print z[0]
            print z[1]
            print z[2]
        return [self.circularize(h) for h in neigh]
