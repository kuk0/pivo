from Genome import Genome
from DCJ import DCJ
from copy import deepcopy, copy
from random import randint
import os

class LinRev(Genome):
    def __init__ (self, s):
        if isinstance(s, str):
            self.set (s)
        else:
            self.n = s.n
            self.g = copy(s.g)
            self.s = copy(s.s)

    def set (self, s):
        s = s.split()
        if s[-1] != '$':
            raise TypeError("The input genome is not linear") 
        s.pop()
        self.n = len(s)
        self.g, self.s, check = [0], [True], (self.n + 1) * [False]
        for x in s:
            x = int(x) # try
            self.s.append(x > 0)
            x = abs(x)
            if not (0 < x <= self.n): print 'mimo'
            if check[x]: print 'uz bolo...'
            check[x] = True
            self.g.append(x)
        self.g.append (self.n + 1)
        self.s.append (True)

    def __len__ (self): return self.n

    def __str__ (self):
        return self.toString() + ' $'
    
    def toString(self):
        sg = [('' if self.s[i] else '-') + str(self.g[i]) for i in xrange(1, self.n + 1)]
        return " ".join(sg)

    def numch (self):
        return (1, 0)

    def _find_components (self, pi, si):
        n, badc = len(pi) - 2, 0
        M1, M2, S1, S2, M, m, comp = [n + 1], [0], [0], [0], [n + 1], [0], []
        Max, Min, mark1, mark2 = (n + 2) * [0], (n + 2) * [n + 1], (n + 2) * [False], (n + 2) * [False]
        for i in xrange(1, n + 2):
            # compute M[i]
            if pi[i - 1] > pi[i]:
                M1.append (pi[i - 1])
            else:
                while M1[-1] < pi[i]: M1.pop()
            M.append (M1[-1])
            # find direct components
            #print S1[-1], pi, M
            s = S1[-1]
            while s > 0 and (pi[s] > pi[i] or M[s] < pi[i]):
                S1.pop()
                s2 = S1[-1]
                Max[s2] = max(Max[s], Max[s2])
                mark1[s2] |= mark1[s]
                s = s2
            pii, pis = pi[i], pi[s]            
            if si[i] and M[i] == M[s] and i - s == pii - pis:
                if pii < pis: pii, pis = pis, pii
                good = mark1[s] or pii==pis+1
                comp.append ((pis, pii, good))
                if not good:
                    badc += 1
                mark1[s] = False
            # compute m[i]
            if pi[i - 1] < pi[i]:
                M2.append (pi[i - 1])
            else:
                while M2[-1] > pi[i]: M2.pop()
            m.append (M2[-1])
            # find reversed components
            s = S2[-1]
            while s > 0 and (pi[s] < pi[i] or m[s] > pi[i]):
                S2.pop()
                s2 = S2[-1]
                Min[s2] = min(Min[s], Min[s2])
                mark2[s2] |= mark2[s]
                s = s2
            pii, pis = pi[i], pi[s]            
            if not si[i] and m[i] == m[s] and i - s == pis - pii:
                if pii < pis: pii, pis = pis, pii
                good = mark2[s] or pii==pis+1
                comp.append ((pis, pii, good))
                if not good:
                    badc += 1
                mark2[s] = False
            # update stacks and marks
            if si[i]: S1.append (i)
            else:     S2.append (i)
            Max.append(pi[i])
            Min[S2[-1]] = min(Min[S2[-1]], pi[i])
            Min.append(pi[i])
            Max[S1[-1]] = max(Max[S2[-1]], pi[i])
            mark1[S1[-1]] = not si[i];
            mark2[S2[-1]] = si[i];
        return badc, comp

    def _cover(self, comp):
        n = m = len(comp) # n = #components; m = #all nodes (including square nodes - these have numbers n...m-1)
        parent = 2*n*[-1]
        node, first, last = [], [], [-1]
        i = 0
        for (x, y, good) in comp:
            if x < last[-1]: # component containing 
                while len(first) > 0 and first[-1] > x:
                    parent[node[-1]] = i
                    node.pop()
                    first.pop()
                    last.pop()
            if x > last[-1]: # new circle node
                node.append(i)
                first.append(x)
                last.append(y)
            if x == last[-1]: # chain
                last_node = node[-1]
                if last_node < n: # last_node is a circle node - we create a new square node
                    node[-1] = parent[last_node] = parent[i] = m
                    m += 1
                else: # square node
                    parent[i] = last_node
            last[-1] = y
            i += 1
        #print 'parents:', parent[0:m]
        
        last_deg3 = -1
        leaves = 0
        j = n
        deg = m*[0]  # degree after cutting good branches
        branch = m*[0] # number of bad components on the branch
        for i in xrange(0,n):
            if parent[i] > j: # process square node j first
                branch[parent[j]] += branch[j]
                if deg[j] > 0:
                    deg[parent[j]] += 1
                if deg[j] > 2:
                    last_deg3 = j                    
                j += 1
            good = comp[i][2]
            branch[parent[i]] += branch[i]
            if not good or deg[i] > 0: # if it is bad or it has bad descendant 
                deg[parent[i]] += 1
            if not good:
                branch[i] += 1
                branch[parent[i]] += branch[i]
                if deg[i] == 0: # bad leaf
                    leaves += 1                
            if deg[i] > 2:
                last_deg3 = i
        while j<m:
            if deg[j] > 2:
                last_deg3 = j
            if parent[j] != -1:
                branch[parent[j]] += branch[j]
                if deg[j] > 0:
                    deg[parent[j]] += 1
            j += 1
        #print 'degs:', deg[0:m]
        #print 'branches:', branch[0:m]

        v = last_deg3
        if last_deg3 == -1:
            #print 'no deg3 vertex'
            return 2 # we assume that badc > 2
        else:
            badc = 0
            v = parent[v]
            while v != -1:
                if v < n and not comp[v][2]:
                    badc += 1
                    #print 'toto:', v
                v = parent[v]
            if badc > 0:
                leaves += 1
                if badc == 1:
                    return leaves # short branch
        
        #print  'leaves: ', leaves
        if leaves%2==0:
            return leaves

        short_branch = False
        for i in xrange(0,m-1):
            if deg[parent[i]] >= 3 and branch[i] == 1:
                short_branch = True
                break
        if short_branch:
            return leaves
        else:
            return leaves+1
#        print deg                

    def cappedString (self):
        sg = ["1"]
        sg.extend([('' if self.s[i] else '-') + str(self.g[i]+1) for i in xrange(1, self.n + 1)])
        sg.append(str(self.n+2))
        return " ".join(sg) + " $"
        
    def dist (self, B):
        if self.n != B.n:
            raise TypeError ("Computing distance of genomes with different number of markers.")
#             print 'BUUG, rozne dlzky', self.n, B.n; print self; print B
        inv, invs = (self.n + 2) * [0], (self.n + 2) * [True]
        for i in xrange(1, self.n + 1):
            inv[B.g[i]] = i
            invs[B.g[i]] = B.s[i]
        inv[self.n + 1] = self.n + 1
        pi, s = [], []
        for i in xrange(0, self.n + 2):
            pi.append(inv[self.g[i]])
            s.append (invs[self.g[i]] == self.s[i])
        # print pi, s
        d = DCJ(self.cappedString()).dist(DCJ(B.cappedString()))
        #print 'dcj: ', d
        badc, comp = self._find_components (pi, s)
        #print 'bad components: ', badc
        #print 'comp: ', comp
        if badc <= 2:
            return d + badc
        #print 'cover:'
        return d + self._cover(comp)
    
    def dist2 (self, B):
        f = open ('input', 'w')
        print >>f, '> g'
        print >>f, self
        print >>f, '> h'
        print >>f, B
        f.close()
        os.system('rm output')
        os.system('./grimm -f input -o output -Ld')
        f = open ('output', 'r')
        for l in f:
            if l[:18] == 'Reversal Distance:':
                f.close()
                return int(l[19:])
        raise Error("grimm")

    def rev (self, i, j): # reverse interval from i to j (inclusive)
        while i < j:
            self.g[i], self.s[i], self.g[j], self.s[j] = self.g[j], not self.s[j], self.g[i], not self.s[i]
            i, j = i + 1, j - 1
        if i == j: self.s[i] = not self.s[i]

    def neigh (self):
        c = []
        for i in xrange (1, self.n + 1):
            for j in xrange (i, self.n + 1):
                t = deepcopy (self)
                t.rev(i, j)
                c.append(t)
        return c
    
    def neigh_part (self, i1, i2, j1, j2):
        """
        List all genomes within one DCJ operation from the genome.
        
        The list has size O(n^2).
        """
        c = []
	i1 = max(i1, 1)
	i2 = max(i2, 1)
        i2 = min(i2, self.n+1)
        j2 = min(j2, self.n+1)
        for i in xrange (i1, i2):
            for j in xrange (j1, j2):
                t = deepcopy (self)
                t.rev(i, j)
                c.append(t)
        return c

    def betterer_neigh (self, z):
        g = DCJ(self.g)
        z = [DCJ(g.__str__()) for g in z]
        return [LinRev(h.__str__()) for h in g.betterer_neigh(z) if h.linear()]
        
    def rand_neigh (self):
        i = randint (1, self.n)
        j = randint (i, self.n)
        g = deepcopy (self)
        g.rev(i, j)
        return g

    def median (self, g1, g2, g3):
        f = open ('input', 'w')
        print >>f, '> g1'
        print >>f, g1.toString()
        print >>f, '> g2'
        print >>f, g2.toString()
        print >>f, '> g3'
        print >>f, g3.toString()
        f.close()
        print os.system("ms/siepel/inv_medians -f input -n " + str(len(g1)) + " -a -o output")
        ret = []
        f = open ('input', 'r')
        for l in f:
            if l[0] != '>':
                ret.append(LinRev(l+' $'))
        f.close()
        return ret
