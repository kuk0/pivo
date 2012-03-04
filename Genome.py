from random import choice, shuffle, sample

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
    def dist (self, B):
        raise NotImplementedError( "Should have implemented this" )
    def neigh (self):
        raise NotImplementedError( "Should have implemented this" )
    def median_solver (self, g1, g2, g3):
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
