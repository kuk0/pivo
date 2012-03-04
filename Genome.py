# TODO: median_sample by nemal byt len pre MultiLinDCJ 
#from MultiLinDCJ import MultiLinDCJ
from random import choice, shuffle, sample, randint

class Genome:
    def __init__ (self, s):
        raise NotImplementedError( "Should have implemented this" )
    def set (self, s):
        raise NotImplementedError( "Should have implemented this" )
    def __len__ (self):
        raise NotImplementedError( "Should have implemented this" )
    def __str__ (self):
        raise NotImplementedError( "Should have implemented this" )
    def numch (self):
        raise NotImplementedError( "Should have implemented this" )
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

    def dist (self, B):
        raise NotImplementedError( "Should have implemented this" )
    def neigh (self):
        raise NotImplementedError( "Should have implemented this" )
    def median_solver (self, g1, g2, g3, errors=0):
        raise NotImplementedError( "Should have implemented this" )

    def median_cloud (self, g1, g2, g3, MAX=300):
        g = self.median_solver(g1, g2, g3)[0]
        t = g.median_score(g1, g2, g3)
        S = [set(), set(), set(), set()] 
        a = [[g], [], [], []]
        i = 0
        while i < 10000 and len(S[0]) <= MAX: # len(a) > 0 and len(S) <= MAX:
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
                    if len(S[0]) > MAX: break
        #print len(S[0]), len(S[1]), len(S[2]), len(S[3]), i
        l = list(S[0])
        if len(l) < MAX: l.extend(list(S[1])) 
        if len(l) < MAX: l.extend(list(S[2])) 
        if len(l) < MAX: l.extend(list(S[3]))
        if len(l) > MAX: l = l[0:MAX]
        return l

    def single_median (self, g1, g2, g3, MAX=300):
	r = self.median_solver(g1, g2, g3)
	if len(r) <= MAX: return r
        else: return sample(r, MAX)
        
    def random_median (self, g1, g2, g3):
        pi = range(1,len(g1)+1)
        shuffle(pi)
        pi = [0] + pi
        ro = range(0,len(g1)+1)
        for i in xrange(1, len(pi)):
            ro[pi[i]] = i
            if randint(0,1) == 0:
                ro[pi[i]] = -i
                pi[i] = -pi[i]
        h1 = self.permute(g1, pi)
        h2 = self.permute(g2, pi)
        h3 = self.permute(g3, pi)
        m = self.median_solver(h1, h2, h3)
        for i in xrange(0, len(m)):
            m[i] = self.permute(m[i], ro) 
        return m
    
    def median_sample (self, g1, g2, g3, MAX=5):
        s = set()
        for i in xrange(0, 2*MAX):
            for m in self.random_median(g1, g2, g3):
                s.add(m)
            if len(s) >= MAX: break
        return list(s)

def id_genome (n):
    """
    Linear genome '1 2 3 ... n $'
    """
    return " ".join([str(x) for x in xrange(1, n + 1)]) + " $"

def id_cgenome (n):
    return " ".join([str(x) for x in xrange(1, n + 1)]) + " @"

def random_linear_genome (n):
    g = [str(x * choice([-1, 1])) for x in xrange(1, n+1)]
    shuffle(g)
    return " ".join(g) + " $"

def random_circular_genome (n):
    g = [str(x * choice([-1, 1])) for x in xrange(1, n+1)]
    shuffle(g)
    return " ".join(g) + " @"
