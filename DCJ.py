from Genome import Genome
from copy import deepcopy, copy
from random import randint, choice
import os
import sys

class DCJ(Genome):
    def __init__ (self, s):
        if isinstance(s, str):
            self.set (s)
        else:
            self.n = len(s)
            self._g = s

    def set (self, s):
        """
        Set the genome to the genome described in string s.
        """
        s = s.split()
        self.n = 2 * (len(s) - s.count("$") - s.count("@"))
        self._g = self.n * [0]
        first, last = -1, -1
        emptychr = True
        for x in s:
            if x == '$':
                if not emptychr: self._g[last], first, last, emptychr = last, -1, -1, True
            elif x == '@':
                if not emptychr: self._g[last], self._g[first], first, last, emptychr = first, last, -1, -1, True
            else:
                x = int(x)
                if x > 0:
                    x = 2 * (x - 1)
                    y = x + 1
                else:
                    x = -2 * x - 1
                    y = x - 1
                if first == -1: first, last = x, x
                self._g[last], self._g[x] = x, last
                last, emptychr = y, False
#    print self._g

    def __len__ (self):
        return self.n/2
        
    def _compstr (self, v, i):
        """
        Auxilliary function creating the string description of the genome.
        """
        j = i
        s = ""
        while not v[j]:
            v[j] = True
            if j % 2 == 0:
                s = s + str(j / 2 + 1) + " "
                v[j + 1] = True
                j = self._g[j + 1]
            else:
                s = s + str(-(j / 2) - 1) + " "
                v[j - 1] = True
                j = self._g[j - 1]
        return s

    def __str__ (self):
        v = self.n * [False]
        s = ""
        for i in xrange(0, self.n):
            if not v[i] and self._g[i] == i:
                s = s + self._compstr (v, i) + "$ "
        for i in xrange(0, self.n):
            if not v[i]:
                s = s + self._compstr (v, i) + "@ "
        return s

    def _mark (self, v, i):
        j = i
        while not v[j]:
            v[j] = True
            if j % 2 == 0: j += 1
            else: j -= 1
            v[j] = True
            j = self._g[j]

    def numch (self):
        """
        Compute the number of chromosomes.
        
        returns (L, C) where L is the number of linear and C the number of circular chromosomes.
        """
        v = self.n * [False]
        s = ""
        L, C = 0, 0
        for i in xrange(0, self.n):
            if not v[i] and self._g[i] == i:
                self._mark (v, i)
                L += 1
        for i in xrange(0, self.n):
            if not v[i]:
                self._mark (v, i)
                C += 1
        return L, C

    def numlin (self):
        """
        Returns the number of linear chromosomes
        """
        L, C = self.numch()
        return L

    def numcirc (self):
        """
        Returns the number of circular chromosomes
        """
        L, C = self.numch()
        return C

    def _comp_size (self, g, v, i, p):
        j = i
        while not v[j]:
            v[j] = 1
    #    print j,
            if p == 0: j = g[j]; p = 1
            else: j = self._g[j]; p = 0
        return p

    def dist (self, B):
        """
        Compute the DCJ distance from genome B.
        """
        if self.n != B.n:
            raise TypeError ("Computing distance of genomes with different number of markers.")
            #print 'BUUG, rozne dlzky', self.n, B.n; print self; print B
        op, c = 0, 0
        v = self.n * [False]
        for i in xrange(0, self.n):
            if not v[i] and self._g[i] == i:
    #     print "("+str(comp_size (A, B, v, i, 0))+")"
                op = op + self._comp_size (B._g, v, i, 0)
            elif not v[i] and B._g[i] == i:
    #     print "("+str(1-comp_size (A, B, v, i, 1))+")"
                op = op + 1 - self._comp_size (B._g, v, i, 1)
        for i in xrange(0, self.n):
            if not v[i]:
    #      print "("+str(comp_size (A, B, v, i, 1))+")"
                c = c + self._comp_size (B._g, v, i, 1)
    #  print n/2, c, op
        return self.n / 2 - c - op / 2

    def swap (self, p, q, p2, q2):
        if p == q and p2 == q2: # {p} {p2} --> {p,p2}
            self._g[p], self._g[p2] = p2, p
        elif p == q: # {p} {p2,q2} --> {p,p2}, {q2}
            self._g[p], self._g[p2], self._g[q2] = p2, p, q2
        elif p2 == q2: # {p,q} {p2} --> {p,p2}, {q}
            self._g[p], self._g[p2], self._g[q] = p2, p, q
        else: # {p,q} {p2,q2} --> {p,p2} {q,q2}
            self._g[p], self._g[p2], self._g[q], self._g[q2] = p2, p, q2, q

    def unswap (self, p, q, p2, q2):
        # nefunguje pre {p,q} {pp}
        if p == q and p2 == q2: # {p} {p2} <-- {p,p2}
            self._g[p], self._g[p2] == p, p2
        elif p == q: # {p} {p2,q2} <-- {p,p2}, {q2}
            self._g[p], self._g[p2], self._g[q2] = p, q2, p2
        elif p2 == q2: # {p,q} {p2} <-- {p,p2}, {q}
            self._g[p], self._g[p2], self._g[q] = q, p2, p
        else: # {p,q} {p2,q2} --> {p,p2} {q,q2}
            self._g[p], self._g[p2], self._g[q], self._g[q2] = q, q2, p, p2

    def linear (self):
        """
        Is the genome linear?
        """
        L, C = self.numch()
        return (C == 0) and (L == 1)

    def multilin (self):
        """
        Is the genome multilinear?
        """
        L, C = self.numch()
        return (C == 0)

    def circular (self):
        """
        Is the genome circular?
        """
        L, C = self.numch()
        return (L == 0) and (C == 1)

    def multicirc (self):
        """
        Is the genome multicircular?
        """
        L, C = self.numch()
        return (L == 0)

    def tempcirc (self):
        """
        Does the genome contain only one (temporary) circular chromosome?
        """
        L, C = self.numch()
        return (C <= 1)

    def unitempcirc (self):
        """
        Is the genome linear, with possibly one (temporary) circular chromosome?
        """
        L, C = self.numch()
        return (L == 1) and (C <= 1)

    def linearize2 (self):
        """
        Returns a set of multilinear genomes close to the genome.
        """
        L, C = self.numch()
        if C == 0: return [self, ]
        g = DCJ(self._g)
        while C > 1: g = choice ([h for h in g.neigh() if h.numcirc() < C]); C = g.numcirc()
        for h in g.neigh():
            if h.multilin(): print h 
        return [h for h in g.neigh() if h.multilin()]

    def neigh (self):
        """
        List all genomes within one DCJ operation from the genome.
        
        The list has size O(n^2).
        """
        c = []
        for i in xrange(0, self.n):
            if self._g[i] >= i:
                for j in xrange(i + 1, self.n):
                    if self._g[j] >= j and j != self._g[i]:
                        u = deepcopy(self)
                        u.swap (i, self._g[i], j, self._g[j])
                        c.append(u)
                if self._g[i] != i:
                    j = self._g[i]
                    u = deepcopy(self)
                    u.swap (i, j, i, j)
                    c.append(u)
        return c

    def neigh_part (self, i1, i2, j1, j2):
        """
        List all genomes within one DCJ operation from the genome.
        
        The list has size O(n^2).
        """
        c = []
        i2 = min(i2, self.n)
        j2 = min(j2, self.n)
        for i in xrange(i1, i2):
            if self._g[i] >= i:
                for j in xrange(j1, j2):
                    if j > i and self._g[j] >= j and j != self._g[i]:
                        u = deepcopy(self)
                        u.swap (i, self._g[i], j, self._g[j])
                        c.append(u)
                if self._g[i] != i:
                    j = self._g[i]
                    u = deepcopy(self)
                    u.swap (i, j, i, j)
                    c.append(u)
        return c

    def neigh2 (self):
        """
        List only multilinear neighbours

        The list has size O(n^2).
        """
        return [g for g in self.neigh() if self.multilin(g)]

    def better_neigh (self, z):
        """
        List only neighbours which improve the overall distance from the neighbouring vertices (in z).
        
        z is a list of genomes
        The returned list has size O(n^2) in the worst case, in practice, it should be shorter. 
        """
        c = self.neigh()
        D = sum([self.dist(g) for g in z])
        return [g for g in c if sum([g.dist(g2) for g2 in z]) < D]

    def even_better_neigh (self, g2):
        """
        List neighbours that are closer to genome g2.
        
        It goes through all adjacencies in g2 and if they are breakpoints, it creates
        a new genome with healed bp in the list of neighbours. The list has size O(n).     
        """
        c = []
        for i in xrange(0, self.n):
            if g2._g[i] >= i:
                j = g2._g[i]
                if self._g[i] != j:
                    u = deepcopy (self)
                    u.swap(i, self._g[i], j, self._g[j])
                    c.append(u)
        return c

    def betterer_neigh (self, z):
        """
        List neigbours that are closer to one of the neighbouring vertices.
        """
        c = []
        c.extend (self.even_better_neigh(z[0]))
        c.extend (self.even_better_neigh(z[1]))
        c.extend (self.even_better_neigh(z[2]))
        return c
        D = sum([self.dist(g) for g in z])
        return [g for g in c if sum([g.dist(g2) for g2 in z]) < D]

    def median_score(self, g1, g2, g3):
        return self.dist(g1) + self.dist(g2) + self.dist(g3)

    # dake pokusy o heuristicke mediany
    def median2 (self, g1, g2, g3):
        z = [g1, g2, g3]
        G = self
        while True:
            c = G.neigh();
            D = G.median_score(g1, g2, g3)
            D2, G = min([(g.median_score(g1, g2, g3), g) for g in c])
            if D2 >= D and randint(0, 3) == 0: break
        return [G, ]

    def _next (self, j):
        if j % 2 == 0: j += 1
        else: j -= 1
        return self._g[j]

    def circularize (self, i1, i2, j1, j2):
        # zatial predpokladam, ze to ma len 2C
        # a spravim vsetky moznosti ako vrazit jeden do druheho
        v = self.n * [False]
        c = []
        self._mark (v, 0)        
        for i in xrange(1, self.n):
            if not v[i]:
                self._mark (v, i)
                k = 0
                while True:
                    j = i
                    if i1 <= k < i2:
                        while True:
                            if j1 <= j < j2:
                                g = deepcopy (self)
                                g.swap(k, self._g[k], j, self._g[j])
                                c.append(g)
                                g = deepcopy (self)
                                g.swap(k, self._g[k], self._g[j], j)
                                c.append(g)
                            j = self._next(j)
                            if j == i: break 
                    k = self._next(k)
                    if k == 0: break
        return c
        
    def median3 (self, g1, g2, g3):
        g = [0, 0, 0]; s = [set(), set(), set()]
        n, g[0], g[1], g[2] = g1.n, copy(g1._g), copy(g2._g), copy(g3._g)
        oth2 = [(0, 1, 2), (1, 0, 2), (2, 0, 1)]
        perm3 = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
        for i in xrange(0, n):
            for j in xrange(0, 3):
                if g[j][i] <= i: s[j].add((i, g[j][i]))
        ii = s[0] & s[1] & s[2]
        s[0] -= ii; s[1] -= ii; s[2] -= ii;
        while True:
            z = False
            for (i, j, k) in oth2:
                for a in xrange(0, n):
                    if g[i][a] < a:
                        b = g[i][a]; c, d = g[j][a], g[j][b]
                        if g[i][c] == d and g[j][a] == g[k][a] and g[j][b] == g[k][b] and g[j][a] != g[i][a] and g[j][b] != g[i][b]:
                            swap2(g[i], a, b, c, d)
                            z = True
            if not z: break
        while True:
            ss = set()
            z = False
            for (i, j, k) in oth2:
                for a in xrange(0, n):
                    if g[i][a] < a and g[j][a] == g[k][a] != g[i][a]:
                        b, c = g[i][a], g[j][a]; d = g[i][c]
                        ss.add((i, a, b, c, d))
            if len(ss) > 0:
                (i, a, b, c, d) = choice(list(ss))
                swap2(g[i], a, b, c, d)
                continue
            for (i, j, k) in perm3:
                for a in xrange(0, n):
                    if g[i][a] < a:
                        b = g[i][a]
                        x, w, y, z = g[j][a], g[k][a], g[j][b], g[k][b]
                        if len(set([a, b, x, y, w, z])) == 6:
                            ss.add((i, j, k, a, b, x, w, y, z))
            if len(ss) > 0:
                (i, j, k, a, b, x, w, y, z) = choice(list(ss))
                swap2(g[j], a, x, b, y)
                swap2(g[k], a, w, b, z)
                continue
            else: break
        for (i, j, k) in oth2:
            for a in xrange(0, n):
                if g[i][a] < a and g[j][a] == g[k][a] == a:
                    swap2(g[i], a, g[i][a], a, g[i][a])
#    print Genome(g[0])
#    print Genome(g[1])
#    print Genome(g[2])
        return DCJ(g[0]).neigh()

    def rand_neigh (self):
        c = deepcopy(self)
        i = randint(0, c.n - 1)
        j = randint(0, c.n - 1)
        if i != j and i != c._g[j]:
            c.swap(i, c._g[i], j, c._g[j])
        return c

    def rand_neigh2 (self, restr):
#    print self, '-->',
        c = self.neigh()
        if restr == 1 or restr == 4:
            c = [g for g in c if g.multicirc()]
        elif restr == 2:
            c = [g for g in c if g.unitempcirc()]
        elif restr == 3:
            c = [g for g in c if g.tempcirc()]
        g = choice(c)
#    print g, '-->',
        L, C = g.numch()
        if restr == 1:
            if C >= 2: g = choice ([h for h in g.neigh() if h.circular()])
        elif restr == 2:
            if C >= 1: g = choice ([h for h in g.neigh() if h.linear()])
        elif restr == 3:
            if C >= 1: g = choice ([h for h in g.neigh() if h.multilin()])
#    print g
        return g

    def path (self, B):
        P = []
        C = deepcopy (self)
        s = set()
        for p in xrange(0, B.n):
            q = B._g[p]
            if q < p: continue
            if p == q:
                u = C._g[p]
                if u != p:
                    s.add((p, p))
            else:
                u, v = min(p, C._g[p]), min(q, C._g[q])
                if u != v:
                    s.add((p, q))
        k = C.dist(B) - 1
        i = 0
        while i < k:
            p, q = choice (list(s))
            u, v = C._g[p], C._g[q]
            if p == q:
                if B._g[u] == u:
                    C.swap(p, u, p, u)
                    s.remove((u, u))
                    s.remove((p, p))
                else: continue
            else:
                s.remove((p, q))
                if u == p and q == v:
                    C.swap(p, q, q, p)
                elif u == p:
                    C.swap(q, v, p, v)
                    if B._g[v] == v: s.remove((v, v))
                elif v == q:
                    C.swap(p, u, q, u)
                    if B._g[u] == u: s.remove((u, u))
                else:
                    C.swap(p, u, q, v)
                    if B._g[u] == v: s.remove((min(u, v), max(u, v)))
            i += 1
            P.append(deepcopy(C))
        return P

    def restrict (self, n):
        s = self.__str__().split()
        t = []
        for x in s:
            if x == '$' or x == '@': t.append(x)
            elif abs(int(x)) <= n: t.append(x)
        self.set(" ".join(t))
        if not self.OK(): print 'BUUUG restrict'

    def add (self):
        self._g.append(self.n)
        self._g.append(self.n + 1)
        self.n += 2
        if not self.OK(): print 'BUUUG add'

    def __hash__ (self):
        return hash(self.__str__())

    def __cmp__ (self, B):
        return cmp (self._g, B._g)

    def OK (self):
        for i in xrange(0, self.n):
            if i != self._g[self._g[i]]: return False
        return True


def swap2 (g, p, q, p2, q2):
    if p == q and p2 == q2: # {p} {p2} --> {p,p2}
        g[p], g[p2] = p2, p
    elif p == q: # {p} {p2,q2} --> {p,p2}, {q2}
        g[p], g[p2], g[q2] = p2, p, q2
    elif p2 == q2: # {p,q} {p2} --> {p,p2}, {q}
        g[p], g[p2], g[q] = p2, p, q
    else: # {p,q} {p2,q2} --> {p,p2} {q,q2}
        g[p], g[p2], g[q], g[q2] = p2, p, q2, q
