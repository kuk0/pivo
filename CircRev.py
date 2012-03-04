from LinRev import LinRev
from DCJ import DCJ
import os

class CircRev(LinRev):
    def __init__ (self, s):
        if isinstance(s, str):
            self.set (s)
        else:
            self.n = len(s)
            self._g = s

    def set (self, s):
        s = s.split()
        if s[-1] != '@':
            raise TypeError("The input genome is not circular")
        s.pop()
        s = [int(x) for x in s]
        n = len(s)
        for i in xrange(0,n):
            if s[i] == n or s[i] == -n:
                break
        if s[i] == n:
            s = s[i+1:n] + s[0:i]
        else:
            l1, l2 = s[0:i], s[i+1:n] 
            l1.reverse()
            l2.reverse()
            s = l1 + l2
            s = [-x for x in s]
        #print s
        self.n = n-1
        self.g, self.s, check = [0], [True], (self.n + 1) * [False]
        for x in s:
            self.s.append(x > 0)
            x = abs(x)
            if not (0 < x < n): print 'mimo'
            if check[x]: print 'uz bolo...'
            check[x] = True
            self.g.append(x)
        self.g.append (self.n + 1)
        self.s.append (True)

    def __len__ (self): return self.n+1

    def __str__ (self):
        return self.toString() + ' @'

    def toString (self):
        return " ".join([('' if self.s[i] else '-') + str(self.g[i]) for i in xrange(1, self.n+1)]) + ' ' + str(self.n+1)
   
    def numch (self):
        return (1, 0)

    def dist2 (self, B):
        f = open ('input', 'w')
        print >>f, '> g'
        print >>f, self.toString()
        print >>f, '> h'
        print >>f, B.toString()
        f.close()
        os.system('rm output')
        os.system('./grimm -f input -o output -Cd')
        f = open ('output', 'r')
        for l in f:
            if l[:18] == 'Reversal Distance:':
                f.close()
                return int(l[19:])
        raise Error("grimm")

    def better_neigh (self, z):
        g = DCJ(self.g)
        z = [DCJ(g.__str__()) for g in z]
        return [CircRev(h.__str__()) for h in g.betterer_neigh(z) if h.circular()]
