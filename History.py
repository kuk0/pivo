import os
from Tree import Tree
from Genome import *
from copy import deepcopy, copy
from random import choice, randint, shuffle, sample
from numpy.random import poisson
import time
#from pyxfigs import *
from CircDCJ import CircDCJ

chngf = open ('changes', 'a')
heurf = open ('heuristics', 'a')

def penalty (p):
#    return 0
    L, C = p
    if C > 1: return 1
    return 0
#    if C > 0: return 10
#    else: return 0
#  if C > 1 or (L > 0 and C > 0): return 10
#  else: return max (0, L-2)

def _c(g, i):
    return "C%d,%d" % (g, i)
def _g(g, i, j):
    assert i != j
    if i > j: i, j = j, i
    return "G%d,%d,%d" % (g, i, j)
def _p(g, i, j):
    assert i != j
    if i > j: i, j = j, i
    return "P%d,%d,%d" % (g, i, j)
def _u(g, i):
    return "U%d,%d" % (g, i)

class History:
    def __init__ (self, T, wdir, simul=False):
        "Zo suboru T vo wdir nacita genomy v listoch (v liste moze byt aj viac moznosti)"
        self.T = T
        N, n = T.N, T.n
        self.wdir = 'data' + os.sep + wdir
        if simul:
            self.g = [0] * N
            return
        f = open(self.wdir + os.sep + 'G', 'r')
        self.g = []
        self.input = []
        for i in xrange(0, n):
            self.input.append([])
        i = 0
        for line in f:
            line = line.strip()
            if line == '' or line[0] == '#': continue
            pos = line.find(' ')
            self.input[self.T.leaf[line[0:pos]]].append(Genome(line[pos + 1:]))
        f.close()
        #for i in xrange(0, n):
            #print self.T.name[i], self.input[i]
        #print n, N
        for i in xrange(0, n):
            self.g.append(deepcopy(choice(self.input[i])))
        for i in xrange(n, N):
            if randint(0, 1): self.g.append (self.g[self.T.right[i]])
            else: self.g.append (self.g[self.T.left[i]])
        self.chng = N * [0]
        self.cand()
        self.TABU = []
        for i in xrange(0,N):
            self.TABU.append({})
        self.tmpTABU = []
        for i in xrange(0,N):
            self.tmpTABU.append(set())

    def __str__ (self):
        s = ""
        for i in xrange(0, self.T.N):
            s = s + self.T.name[i] + ' ' + str(self.g[i].dist(self.g[self.T.parent[i]])) + ' ' + self.g[i].__str__() + '\n' #' ' + str(self.g[i].numch()) +'\n'
        return s #+str(self.score())+'\n'

    def read (self, fn):
        """
        Zo suboru fn nacita historiu v standardnom formate
        (meno, vzdialenost od otca, genom). Pozor, mena sa
        ignoruju - vrcholy musia ist v takom poradi,
        ako su v strome.
        """
        f = open(self.wdir + os.sep + fn, 'r')
        #g = []
        i = 0
        for line in f:
            line = line.strip()
            if line == '' or line[0] == '#': continue
            #print line
            self.g[i] = Genome(line.split(None, 2)[2])
            i += 1
        f.close()

    def read_cand (self, fn):
        """
        Zo suboru fn nacita historiu ako dalsieho kandidata.
        """
        f = open(self.wdir + os.sep + fn, 'r')
        #g = []
        i = 0
        tmp = CircDCJ("1 2 3 @")
        for line in f:
            line = line.strip()
            if line == '' or line[0] == '#': continue
            if i >= self.T.n:
                #self.c[i].append(Genome(line.split(None, 2)[2]))
                self.c[i].append(tmp.circularize(Genome(line.split(None, 2)[2])))
                ##self.heur[i].append(-1)
            i += 1
        f.close()

    def cand (self):
        """
        Vytvori nove zoznamy kandidatov, inicializuje ich genomami v historii
        """
        self.c = []
        ##self.heur = []
        for i in xrange(0, self.T.n):
            self.c.append(deepcopy(self.input[i]))
            ##self.heur.append([0]*len(self.input[i]))
        for i in xrange(self.T.n, self.T.N):
#            if self.g[i].circular():
                self.c.append([self.g[i], ])
#            else:
#                self.c.append([])
            ##self.heur.append([0,])
            ##if len(self.heur[i]) != len(self.c[i]): print 'cand'

    def cand2 (self, m): # dana dlzka genomu
        """
        Vytvori nove zoznamy kandidatov, inicializuje ich genomami v historii
        """
        self.c = []
        ##self.heur = []
        for g in self.g:
            if g.n > 2 * m: g.restrict(m)
            if 2 * m > g.n: g.add()
        for i in xrange(0, self.T.n):
            self.c.append(deepcopy(self.input[i]))
            for g in self.c[i]: g.restrict(m)
            ##self.heur.append([0]*len(self.input[i]))
        for i in xrange(self.T.n, self.T.N):
            self.c.append([self.g[i], ])
            ##self.heur.append([0,])
            ##if len(self.heur[i]) != len(self.c[i]): print 'cand'

    def fill_in_desc (self):
        "Do kazdeho zoznamu kandidatov prida tych, ktori su v zozname pod nim"
        for i in xrange(self.T.n, self.T.N):
            ##if len(self.heur[i]) != len(self.c[i]): print 'pred desc'
            self.c[i].extend (self.c[self.T.left[i]])
            ##self.heur[i].extend ([1]*len(self.c[self.T.left[i]]))
            self.c[i].extend (self.c[self.T.right[i]])
            if len(self.c[i]) > 300:
                self.c[i] = sample(self.c[i], 300)
            ##self.heur[i].extend ([1]*len(self.c[right[i]]))
            ##if len(self.heur[i]) != len(self.c[i]): print 'desc'

    def fill_in_middle (self):
        "Do zoznamu kandidatov prida genomy na nahodnej ceste medzi lavym a pravym synom"
        for i in xrange(self.T.n, self.T.N):
            r, l = self.g[self.T.right[i]], self.g[self.T.left[i]]
            p = l.path(r)
            self.c[i].extend (p)
            ##self.heur[i].extend ([2]*len(p))
            ##if len(self.heur[i]) != len(self.c[i]): print 'middle'

    def fill_paths (self):
        "Do zoznamu kandidatov prida nahodne cesty medzi vsetkymi dvojicami listov v podstromoch"
        desc = []
        for i in xrange(0, self.T.n):
            desc.append([i, ])
        for i in xrange(self.T.n, self.T.N):
            l, r = self.T.left[i], self.T.right[i]
            for j in xrange(0, len(desc[l])):
                for k in xrange(0, len(desc[r])):
                    e = self.g[desc[l][j]].path(self.g[desc[r][k]])
                    self.c[i].extend(e)
                    ##self.heur[i].extend([3]*len(e))
                    ##if len(self.heur[i]) != len(self.c[i]): print 'paths'
            desc.append (desc[l] + desc[r])

    def score (self):
        """
        Score of the evolutionary history
        """
        s = 0
        for i in xrange(0, self.T.N):
            s = s + self.g[i].dist(self.g[self.T.parent[i]])
        for i in xrange(self.T.n, self.T.N):
            s = s + penalty (self.g[i].numch())
        return s

    def write (self, fn = None): #, suffix):
        """
        Write the history into a file in the working directory.
        
        The filename is of the form score-time.
        """
        #now = time.localtime(time.time())
        if fn == None:
            fn = ("%04d" % self.score()) + '-' + str(time.time())
            #+time.strftime("%d-%m-%y_%H:%M:%S", now)
        f = open (self.wdir + os.sep + fn, 'w')
        f.write (self.__str__())
        f.close()

    def neigh (self):
        "Do zoznamu kandidatov prida okolia aktualnych genomov"
        for i in xrange(self.T.n, self.T.N - 1):
            e = []
            e.extend (self.g[i].neigh())
            #e.extend (self.g[i].even_better_neigh(self.g[self.T.left[i]]))
            #e.extend (self.g[i].even_better_neigh(self.g[self.T.right[i]]))
            if i == self.T.left[self.T.N - 1]: p = self.T.right[self.T.N - 1]
            elif i == self.T.right[self.T.N - 1]: p = self.T.left[self.T.N - 1]
            else: p = self.T.parent[i]
            #e.extend (self.g[i].even_better_neigh(self.g[p]))
            #e = self.g[i].betterer_neigh([self.g[self.T.left[i]], self.g[self.T.right[i]], self.g[p]])
            self.c[i].extend(e)
            #print len(self.c[i])
            ##self.heur[i].extend([4]*len(e))
#   e = self.g[N-1].better_neigh([self.g[self.T.left[N-1]], self.g[self.T.right[N-1]]])
#   self.c[N-1].extend(e)
        ##self.heur[N-1].extend([4]*len(e))

    def neigh_part (self, i1, i2, j1, j2):
        "Do zoznamu kandidatov prida okolia aktualnych genomov"
        for i in xrange(self.T.n, self.T.N - 1):
            e = []
            if self.g[i].circular():
                e.extend (self.g[i].neigh_part(i1, i2, j1, j2))
            else:
                e.extend (self.g[i].circularize(i1, i2, j1, j2))
            self.c[i].extend(e)

    def neigh2 (self):
        "Do zoznamu kandidatov prida okolia kandidatov"
        for i in xrange(self.T.n, self.T.N):
            k = len(self.c[i])
            for j in xrange(0, k):
                e = self.c[i][j].betterer_neigh([self.g[self.T.left[i]], self.g[self.T.right[i]], self.g[self.T.parent[i]]])
                self.c[i].extend(e)
                ##self.heur[i].extend([5]*len(ebetter_neigh))

    def neigh3 (self):
        "Do zoznamu kandidatov prida okolia kandidatov"
        for i in xrange(self.T.n, self.T.N):
            e = self.g[i].better_neigh([self.g[self.T.left[i]], self.g[self.T.right[i]], self.g[self.T.parent[i]]])
            self.c[i].extend(e)

    def median (self):
        for i in xrange(self.T.n, self.T.N - 1):
            if i == self.T.left[self.T.N - 1]:
                self.c[i].extend(self.g[i].median(self.g[self.T.left[i]], self.g[self.T.right[i]], self.g[self.T.right[self.T.N - 1]]))
            elif i == self.T.right[self.T.N - 1]:
                self.c[i].extend(self.g[i].median(self.g[self.T.left[i]], self.g[self.T.right[i]], self.g[self.T.left[self.T.N - 1]]))
            else:
                self.c[i].extend(self.g[i].median(self.g[self.T.left[i]], self.g[self.T.right[i]], self.g[self.T.parent[i]]))

    def _replace (self, M, v, i):
        global chng
        self.chng[v] = self.g[v].dist(self.c[v][i])
        self.g[v] = self.c[v][i]
        if v < self.T.n: return
        self._replace (M, self.T.left[v], M[v][i][0])
        self._replace (M, self.T.right[v], M[v][i][1])
        
    def _purge (self):
        N, n = self.T.N, self.T.n
        for i in xrange(n,N):
            self.c[i] = [g for g in self.c[i] if g not in self.TABU[i] and g not in self.tmpTABU[i]]
            #+[self.c[i][0]]
        
    def _opt_sols (self, B, M, u, i, m, vis, S):
        if u < self.T.n: return
        if (u,i) in vis: return
        vis.add((u,i))
        self.TABU[u][deepcopy(self.c[u][i])] = S
        v, w = self.T.left[u], self.T.right[u]
        lv, lw = len(self.c[v]), len(self.c[w])
        bl = B[v][M[u][i][0]]
        dl = self.c[v][M[u][i][0]].dist(self.c[u][i])
        br = B[w][M[u][i][1]]
        dr = self.c[w][M[u][i][1]].dist(self.c[u][i])
        for j in xrange(0, lv):
            dd = self.c[v][j].dist(self.c[u][i])
            if B[v][j] + dd == bl + dl:
                self._opt_sols (B, M, v, j, bl, vis, S)
        for j in xrange(0, lw):
            dd = self.c[w][j].dist(self.c[u][i])
            if B[w][j] + dd == br + dr:
                self._opt_sols (B, M, w, j, br, vis, S)

    def opt_neigh (self, tabu=True, newtabu=False):
        "Najde historiu v okoli s najmensim skore"
        global chngf#, heurf
        N, n = self.T.N, self.T.n
        left, right = self.T.left, self.T.right
        if not tabu:
            self._purge()
        c = self.c
        if len(c[left[N - 1]]) < len(c[right[N - 1]]):
            c[N - 1] = c[left[N - 1]]
##      self.heur[self.T.N-1] = self.heur[self.T.left[self.T.N-1]]
        else:
            c[N - 1] = c[right[N - 1]]
##      self.heur[self.T.N-1] = self.heur[self.T.right[self.T.N-1]]
        B, M = [], []
        print [len(c[i]) for i in xrange(n, N)]
##    ns = [0]*self.T.N
##    for u in xrange0,self.T.N):
##      ns[u] = [1]*len(self.c[u])
##      if len(self.heur[u]) != len(self.c[u]): print u, '!!!!!'
        for u in xrange(0, n):
            B.append(len(c[u]) * [0])
            M.append([])
            for i in xrange(0, len(c[u])):
                M[u].append([0, 0])
        for u in xrange(n, N):
            v, w = left[u], right[u]
            ul, vl, wl = len(c[u]), len(c[v]), len(c[w])
            B.append([]); M.append([])
            for i in xrange(0, ul):
                B[u].append(penalty(c[u][i].numch()))
                M[u].append([0, 0])
            for i in xrange(0, ul):
                m, mi, q = 999999, 0, 1
                for j in xrange(0, vl):
                    t = B[v][j] + c[u][i].dist(c[v][j])
##          if t == m: q += ns[v][j]
                    if t < m: m, mi = t, j # , q = t, j, ns[v][j]
                B[u][i], M[u][i][0] = B[u][i] + m, mi #, ns[u][i] = B[u][i]+m, mi, q
                m, mi, q = 999999, 0, 1
                for j in xrange(0, wl):
                    t = B[w][j] + c[u][i].dist(c[w][j])
##          if t == m: q += ns[w][j]
                    if t < m: m, mi = t, j # , q = t, j, ns[w][j]
                B[u][i], M[u][i][1] = B[u][i] + m, mi
##        ns[u][i] *= q
        m, mi, q = 999999, 0, 1
        for i in xrange(0, len(c[N - 1])):
##      if B[self.T.N-1][i] == m: q += ns[self.T.N-1][i]
            if B[N - 1][i] < m: m, mi = B[N - 1][i], i #, q = B[self.T.N-1][i], i, ns[self.T.N-1][i]

        if newtabu:
            for i in xrange(0,len(c[N-1])):
                if B[N-1][i] == m:
                    ## ns[N-1][i] = -ns[N-1][i]
                    self._opt_sols(B, M, N-1, i, m, set(), m)
##    hsug = [0]*9
##    hopt = [0]*9
##    for v in xrange(self.T.n,self.T.N):
##      for i in xrange(0,len(self.c[v])):
##        hsug[self.heur[v][i]] += 1
##        if ns[v][i] < 0: hopt[self.heur[v][i]] += 1
##    print >> heurf, q, "optimalnych rieseni"
##    for i in xrange(0,9): print >> heurf, "%d/%d" % (hopt[i],hsug[i]),
##    print >> heurf
        self.chng = [0] * self.T.N
        self._replace (M, self.T.N - 1, mi)
        print >> chngf, '\n', self.chng[self.T.n:self.T.N], m,
        if m != self.score(): print "BUGGGGGGGGGGG", m, self.score()
        return m

    def local_opt (self):
        raise NotImplementedError( "Should have implemented this" )
#        "Najde lokalne minimum historii"
#        self.cand()
#        self.fill_in_desc()
#        self.rand_neigh()
#    self.fill_in_middle()
#    self.fill_fitch()
        #s, s2 = 9999, self.opt_neigh()
#        s, s2 = 9999, 9998
#    print >> chngf, ' *',
#    self.cand()
#    s, s2 = s2, self.opt_neigh()
#    print s2
        #while s2 < s:
            #print "pred"
            #print self
        #    self.cand()
#      self.fill_in_middle()
#      self.fill_fitch()
#            self.median()
         #   self.median()
            #self.neigh()
            #print "po"
            #print self
          #  s, s2 = s2, self.opt_neigh()
            #print s2
        #print self

    def local_opt_steinerization (self):
        "Najde lokalne minimum historii"
        s, s2 = 9999, 9998
        while s2 < s:
            s = s2
            for i in xrange(self.T.n, self.T.N - 1):
                self.cand()
                if i == self.T.left[self.T.N - 1]:
                    self.c[i].extend(self.g[i].median(self.g[self.T.left[i]], self.g[self.T.right[i]], self.g[self.T.right[self.T.N - 1]]))
                elif i == self.T.right[self.T.N - 1]:
                    self.c[i].extend(self.g[i].median(self.g[self.T.left[i]], self.g[self.T.right[i]], self.g[self.T.left[self.T.N - 1]]))
                else:
                    self.c[i].extend(self.g[i].median(self.g[self.T.left[i]], self.g[self.T.right[i]], self.g[self.T.parent[i]]))
                s2 = self.opt_neigh()
            print s2
        print self
        
    def local_opt_medians (self):
        "Najde lokalne minimum historii"
        s, s2 = 9999, 9998
        while s2 < s:
            self.cand()
            self.median()
            s, s2 = s2, self.opt_neigh()
            print s2
        print selfobject
        
    def local_opt_better_neighbours (self):
        "Najde lokalne minimum historii"
        s, s2 = 9999, 9998
        while s2 < s:
            self.cand()
            self.neigh2()
            s, s2 = s2, self.opt_neigh()
            print s2
        print self
        
    def local_opt_neighbours2 (self):
        "Najde lokalne minimum historii"
        self.cand()
        self.neigh()
        s, s2 = 9999, self.opt_neigh()
        print s2
        while s2 < s:
            self.cand()
            self.neigh2()
            s, s2 = s2, self.opt_neigh()
            print s2
        self.cand()
        self.neigh()
        self.opt_neigh(newtabu=True)
        print self
        
    def local_opt_neighbours (self):
        "Najde lokalne minimum historii"
        m = len(self.g[0])
        #self.cand()
        #self.init_hist()
        s, s2 = 9999, 9998 #self.opt_neigh()
        while s2 < s:
            s = s2
#            for i in xrange(0,m,20):
            y = self.chunks
            if self.x == -1:
                rr = xrange(0,m,y)
            else:
                rr = [self.x]
            for i in rr:
                for j in xrange(i,m,y):
                    self.cand()
                    self.neigh_part(i, i+y, j, j+y)
                    s2 = self.opt_neigh()
                    #print i, j, '>', s2, self.noncirc()            
        print self
    
    def local_opt_tabu (self):
        "Najde lokalne minimum historii"
        self.cand()
        for i in xrange(self.T.n,self.T.N):
            self.c[i] = self.TABU[i].keys()
        s, s2 = 9999, self.opt_neigh(True)
        while s2 < s:
            self.cand()
            self.neigh2()
            s, s2 = s2, self.opt_neigh(True)
            print s2
        print self
        
    def noncirc (self):
        num = 0
        for i in xrange(self.T.n, self.T.N):
            if not self.g[i].circular():
                num += 1
        return num

    def local_opt_desc (self):
        "Najde lokalne minimum historii"
        self.cand()
        self.fill_in_desc()
        s, s2 = 9999, self.opt_neigh()
        while s2 < s:
            self.cand()
            self.median()
            s, s2 = s2, self.opt_neigh()
            print s2
        print self

    def local_opt_incr_lengths (self):
        "Najde lokalne minimum historii"
        gl = self.input[0][0].n / 2
        for i in xrange(5, gl):
            print 'LEN:', i
            self.cand2(i);
            self.fill_in_middle(); self.fill_in_desc()
            s, s2 = 9999, self.opt_neigh()
            while s2 < s:
                self.cand2(i); self.fill_in_middle(); self.fill_in_desc(); self.neigh()
                s, s2 = s2, self.opt_neigh()
                print s2
        print "---"
        
    def local_opt_circ (self):
        N, n = self.T.N, self.T.n
        g = self.g
        for i in xrange(n,N):
            if g[i].circular(): continue
            gg = g[i]
            a, b = gg._circ_chromosomes()
            if len(a) < len(b):
                a, b = b, a
            print i-n, ". genom ma", gg.numch()
            for p in a:
                self.cand();
                q = gg._g[p]
                for s in b:
                    t = gg._g[s]
                    h = deepcopy(g[i])
                    h.swap(p,q,s,t)
                    self.c[i].append(h)
                    for j in xrange(n,N):
                        h = deepcopy(g[j])
                        qq, pp, tt, ss = h._g[p], h._g[q], h._g[s], h._g[t]
                        if q!=qq or s!=ss:
                            h.swap(p,qq,ss,t)
                            if h.circular(): self.c[j].append(h);
                        if p!=pp or s!=ss: 
                            h = deepcopy(g[j]); h.swap(pp,q,ss,t);
                            if h.circular(): self.c[j].append(h);
                        if q!=qq or t!=tt:
                            h = deepcopy(g[j]); h.swap(p,qq,s,tt);
                            if h.circular(): self.c[j].append(h);
                        if p!=pp or t!=tt:
                            h = deepcopy(g[j]); h.swap(pp,q,s,tt);
                            if h.circular(): self.c[j].append(h);
                self.opt_neigh()
                if g[i].circular():
                    print "HURA"
                    break
        print self
    
    def _rand_leaf(self, desc): # z int 0,...,a-1,b+1,...,n-1
        a, b = desc
        x = randint (0,a+self.T.n-b-2)
        if x < a: return x
        else: return (x-a) + (b+1) 
         
    def init_hist (self):
        n, N = self.T.n, self.T.N
        left, right, parent = self.T.left, self.T.right, self.T.parent
        desc = []
        for i in xrange(0, n):
            desc.append((i,i))
        for i in xrange(n, N):
            desc.append((desc[left[i]][0],desc[right[i]][1]))
        for i in xrange(n, N-1):
            for k in xrange(0,1): # TODO
                z0 = self.g[randint(desc[left[i]][0],desc[left[i]][1])]
                z1 = self.g[randint(desc[right[i]][0],desc[right[i]][1])]
                z2 = self.g[self._rand_leaf(desc[i])]
#                print i, ':', choice(desc[left[i]]), choice(desc[right[i]]), choice(desc2[i])
                self.c[i].extend (self.g[i].median(z0, z1, z2)) #, MAX=20
            
#    def init_hist (self):
#        n, N = self.T.n, self.T.N
#        left, right, parent = self.T.left, self.T.right, self.T.parent
#        desc = []
#        for i in xrange(0, n):
#            desc.append([i])
#        for i in xrange(n, N):
#            desc.append(desc[left[i]]+desc[right[i]])
##	    print 'desc', i, desc[i]
#    	desc2 = N*[0]
#        ll, rr = left[N-1], right[N-1]
#        desc2[ll] = copy(desc[rr])
#        desc2[rr] = copy(desc[ll])
##	print 'desc2', ll, desc2[ll]
##	print 'desc2', rr, desc2[rr]
#        for i in reversed(xrange(n, N-1)):
#            if i == ll or i == rr: continue
#            desc2[i] = copy(desc2[parent[i]])
#            if left[parent[i]] == i: # lavy syn
#                desc2[i].extend(desc[right[parent[i]]])
#            else:
#                desc2[i].extend(desc[left[parent[i]]])
##	    print 'desc2', i, desc2[i]
#        for i in xrange(n, N-1):
#            for k in xrange(0,10):
#                z0 = self.g[choice(desc[left[i]])]
#                z1 = self.g[choice(desc[right[i]])]
#                z2 = self.g[choice(desc2[i])]
##                print i, ':', choice(desc[left[i]]), choice(desc[right[i]]), choice(desc2[i])
#                self.c[i].extend (self.g[i].median(z0, z1, z2)) #, MAX=20

    def paths_opt (self):
        self.cand()
        self.fill_paths()
        self.opt_neigh()
        print self.score()
        print self

    def median_opt (self):
        s = 10000
        s2 = s - 1
        while s2 < s:
            self.fill_in_oe(0)
            s, s2 = s2, self.opt_neigh()
            self.fill_in_oe(1)
            s2 = self.opt_neigh()
            print s2
        print "---"

    def swap_subtree (self, h, v):
        if v < self.T.n: return
        self.g[v], h.g[v] = h.g[v], self.g[v]
        self.swap_subtree(h, self.T.left[v])
        self.swap_subtree(h, self.T.right[v])

    def mutate (self):
        # i = randint (self.T.n,self.T.N-1)
        for i in xrange(0,self.T.N):
            self.tmpTABU[i] = set()
        for i in xrange(0, self.T.N):
            self.g[i] = self.g[i].rand_neigh()
            self.tmpTABU[i].add(deepcopy(self.g[i]))
#            self.g[i] = self.g[i].rand_neigh()
#            self.tmpTABU[i].add(deepcopy(self.g[i]))
            self.g[i] = self.g[i].rand_neigh()
            #self.g[i] = self.g[i].rand_neigh().rand_neigh() #.rand_neigh() #.rand_neigh().rand_neigh()

    def rand_hist (self, restr=None):
        g = range(1, len(self.input[0][0]) + 1)
        for i in xrange(self.T.n, self.T.N):
            for j in xrange(0, len(g)):
                g[j] = g[j] * choice([-1, 1])
            shuffle(g)
            if restr == 1 or restr == 4:
                s = " ".join([str(x) for x in g]) + " @"
            elif restr == 2 or restr == 3:
                s = " ".join([str(x) for x in g]) + " $"
            else: s = " ".join([str(x) for x in g]) + " " + choice('$@')
            self.g[i] = Genome(s)

    def rand_neigh (self):
        for i in xrange(self.T.n, self.T.N):
            self.c.append([self.g[i], ])
##      self.heur.append([0,]) #TODO ???
            for j in xrange(0, 100):
                self.c[i].append(self.g[i].rand_neigh())
            for j in xrange(0, 200):
                self.c[i].append(self.g[i].rand_neigh().rand_neigh())
##      self.heur[i].extend([6]*100)
##      self.heur[i].extend([7]*200)
            #for j in xrange(0,300):
            #  self.c[i].append(self.g[i].rand_neigh().rand_neigh().rand_neigh())

    def rand_opt (self):
        s = self.score()
        s2 = s - 1
        print s
        while s2 < s:
            self.cand()
            self.fill_in_desc()
            self.rand_neigh()
            s, s2 = s2, self.opt_neigh()
            print s2

    def fill_in_rand (self):
        for i in xrange(self.T.n, self.T.N):
            for k in xrange(0, 300):
                g = range(1, 41)
                for j in xrange(0, len(g)):
                    g[j] = g[j] * choice([-1, 1])
                shuffle(g)
                s = " ".join([str(x) for x in g]) + " " + choice('$@')
                self.c[i].append (Genome(s))
##        self.heur[i].append (8)

    def fill_fitch (self):
        gl = self.g[0].n
        S = []
        for i in xrange(0, self.T.n):
            S.append([])
            for j in xrange(0, gl):
                S[i].append(set([self.g[i]._g[j], ]))
        for i in xrange(self.T.n, self.T.N):
            S.append([])
            u, v = self.T.left[i], self.T.right[i]
            for j in xrange(0, gl):
                s = S[u][j] & S[v][j]
                if len(s) == 0: S[i].append(S[u][j] | S[v][j])
                else: S[i].append(s)
        for i in reversed(xrange(self.T.n, self.T.N)):
            u, v = self.T.left[i], self.T.right[i]
            for j in xrange(0, gl):
                if not((S[u][j] | S[v][j]) <= S[i][j]): S[u][j], S[v][j] = S[i][j], S[i][j]
        for i in xrange(self.T.n, self.T.N):
            for ii in xrange(0, 20):
                f = set(range(0, gl))
                cand = Genome(copy(self.g[i]._g))
                for j in xrange(0, gl):
                    if j not in f: continue
                    g = S[i][j] & f
                    if len(g) == 0: k = choice(list(f)); #print "u",
                    else: k = choice(list(g)); #print "O",
                    cand._g[j], cand._g[k] = k, j
                    f.discard(k)
                    f.discard(j)
                self.c[i].append(cand)

###    self.c[i].append (Genome(s))
##        self.heur[i].append (8)

    def fill_in_oe (self, p):
        self.cand()
        for i in xrange(self.T.n, self.T.N):
            if level[i] % 2 == p:
                s = []
#        s = set()
                neigh = []
                for k in xrange(0, len(self.c[i])):
                    neigh.extend (self.c[i][k].neigh())
                self.c[i].extend (neigh)
                print i, len(neigh)
                for k in xrange(0, len(neigh)):
                    s.extend (neigh[k].neigh())
#           s.union(set(neigh[k].neigh()))  # preco to tam nic nehadze?
#           if k % 50 == 0:
#             print '.', len(s)
                self.c[i].extend (list(s))

    def _simul (self, v, mean):
        T = self.T
        if T.parent[v] == T.N-1:
            d = poisson (mean*0.5, 1)[0]
        else:
            d = poisson (mean, 1)[0]
        G = self.g[T.parent[v]] # Genome()
        for i in xrange(0, d):
            G = G.rand_neigh()
        self.g[v] = G
        if T.left[v] < v: self._simul(T.left[v], mean)
        if T.right[v] < v: self._simul(T.right[v], mean)

    def simul_data (self, mean, m):
        T = self.T
        self.g[T.N - 1] = Genome (id_cgenome (m))             # TODO id_genome alebo id_cgenome?
        self._simul (T.left[T.N - 1], mean)
        self._simul (T.right[T.N - 1], mean)
        self.input = []
        for i in xrange(0, T.n):
            self.input.append([self.g[i]])
        f = open (self.wdir + os.sep + 'true-hist.txt', 'w')
        f.write ('# ' + str(self.score()) + '\n')
        f.write (self.__str__())
        f.close ()
        f = open (self.wdir + os.sep + 'G', 'w')
        for i in xrange(0, self.T.n):
            f.write(self.T.name[i]+' ' + self.g[i].__str__() + '\n')
        f.close()

    def dist_tab (self):
        T = []
        for i in xrange(0, self.T.n):
            self.T.append([])
            for j in xrange(0, i):
                T[i].append(self.g[i].dist(self.g[j]))
        return T

    def _path (self, i, j):
        s = []
        while level[i] > level[j]:
            s.append('h' + str(i))
            i = self.T.parent[i]
        while level[j] > level[i]:
            s.append('h' + str(j))
            j = self.T.parent[j]
        while i != j:
            s.append('h' + str(i))
            s.append('h' + str(j))
            i = self.T.parent[i]
            j = self.T.parent[j]
        return s

    def _path2 (self, i, j):
        return '+'.join(self._path(i, j))

    def lower_bound (self, fn):
        f = open (fn, 'w')
        print >> f, "minimize"
        print >> f, "obj:"
        all = ['h' + str(i) for i in xrange(0, self.T.N - 1)]
        print >> f, '+'.join(all)
        print >> f, "subject to"
        T = self.dist_tab()
        for i in xrange(0, self.T.n):
            for j in xrange(0, i):
                print >> f, 'c_' + str(i) + ',' + str(j) + ':', self._path2(i, j), '>=', T[i][j]
        print >> f, "bounds"
        for h in all:
            print >> f, h, '>= 0'
        print >> f, 'generals'
        print >> f, ' '.join(all)
        print >> f, "end"
        f.close()

    def _path3 (self, i, j, k):
        t = self._path(i, j) + self._path(i, k) + self._path(j, k)
        print t
        t.sort()
        print t
        last = t[0]
        lasti = i = 1
        self.T.n = len(t)
        while i < self.T.n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        t = t[:lasti]
        return '+'.join(t)

    def lower_bound2 (self, fn):
        f = open (fn, 'w')
        print >> f, "minimize"
        print >> f, "obj:"
        all = ['h' + str(i) for i in xrange(0, self.T.N - 1)]
        print >> f, '+'.join(all)
        print >> f, "subject to"
        T = self.dist_tab()
        for i in xrange(0, self.T.n):
            for j in xrange(0, i):
                print >> f, 'c_' + str(i) + ',' + str(j) + ':', self._path2(i, j), '>=', T[i][j]
        for i in xrange(0, self.T.n - 2):
            for j in xrange(i + 1, self.T.n - 1):
                for k in xrange(j + 1, self.T.n):
                    fm = open ('med/' + str(i) + '-' + str(j) + '-' + str(k) + '.md', 'r')
                    print >> f, 'm_' + str(i) + ',' + str(j) + ',' + str(k) + ':', self._path3(i, j, k), '>=', int(fm.read())
                    fm.close()
        print >> f, "bounds"
        for h in all:
            print >> f, h, '>= 0'
        print >> f, 'generals'
        print >> f, ' '.join(all)
        print >> f, "end"
        f.close()

    def LP (self, fn):
        f = open (fn, 'w')
        print >> f, "maximize"
        print >> f, "obj:"

        L = self.g[0].n
        cycles = [_c(g, i) for g in xrange(0, self.T.N - 1) for i in xrange(0, L)]
        print >> f, '+'.join(cycles)

        print >> f, "subject to"
        for g in xrange(0, self.T.n):
            for i in xrange(0, L):
                j = self.g[g]._g[i]
                if i < j:
                    print g, i, j
                    print >> f, _g(g, i, j) + " = 1"
        for g in xrange(0, self.T.N):
            for i in xrange(0, L - 1):
                print >> f, '+'.join([_g(g, i, j) for j in xrange(0, L) if i != j]) + " = 1"
        for g in xrange(0, self.T.N - 1):
            h = self.T.parent[g]
            for i in xrange(0, L - 1):
                for j in xrange(i + 1, L):
                    print >> f, _p(g, i, j), '-', _g(g, i, j), ">= 0"
                    print >> f, _p(g, i, j), '-', _g(h, i, j), ">= 0"
                    for k in xrange(0, L):
                        if k != i and k != j:
                            print >> f, _p(g, i, j), '-', _p(g, i, k), '-', _g(g, k, j), ">= -1"
                            print >> f, _p(g, i, j), '-', _p(g, i, k), '-', _g(h, k, j), ">= -1"
                    print >> f, _c(g, j), '+', _p(g, i, j) + "<= 1"

        print >> f, "binary"
        print >> f, '\n'.join(cycles)
        print >> f, '\n'.join([_g(g, i, j) for g in xrange(0, self.T.N - 1) for i in xrange(0, L - 1) for j in xrange(i + 1, L)])
        print >> f, '\n'.join([_p(g, i, j) for g in xrange(0, self.T.N - 1) for i in xrange(0, L - 1) for j in xrange(i + 1, L)])
        print >> f, "end"
        f.close()


    def LPnove (self, fn):
        f = open (fn, 'w')
        print >> f, "maximize"
        print >> f, "obj:"

        L = self.g[0].n
        cycles = [_c(g, i) for g in xrange(0, self.T.N - 1) for i in xrange(0, L)]
        print >> f, '+'.join(cycles)

        print >> f, "subject to"
        for g in xrange(0, self.T.n):
            for i in xrange(0, L):
                j = self.g[g]._g[i]
                if i < j:
                    print g, i, j
                    print >> f, _g(g, i, j) + " = 1"
        for g in xrange(0, self.T.N):
            for i in xrange(0, L - 1):
                print >> f, '+'.join([_g(g, i, j) for j in xrange(0, L) if i != j]) + " = 1"
        for g in xrange(0, self.T.N - 1):
            h = self.T.parent[g]
            for i in xrange(0, L):
                for j in xrange(0, L):
                    if i != j:
                        print >> f, _u(g, i), '-', _u(g, j), '+', L, _g(g, i, j), '-', L, _c(g, i), "<=", L - 1
                        print >> f, _u(g, j), '-', _u(g, i), '+', L, _g(g, i, j), '-', L, _c(g, i), "<=", L - 1
                        print >> f, _u(g, i), '-', _u(g, j), '+', L, _g(h, i, j), '-', L, _c(g, i), "<=", L - 1
                        print >> f, _u(g, j), '-', _u(g, i), '+', L, _g(h, i, j), '-', L, _c(g, i), "<=", L - 1
        for g in xrange(0, self.T.N - 1):
            for i in xrange(0, L):
                print >> f, _u(g, i), '+', _c(g, i), ' >= 1'
#          print >> f, _p(g,i,j), '-', _g(g,i,j), ">= 0"
#          print >> f, _p(g,i,j), '-', _g(h,i,j), ">= 0"
#          for k in xrange(0,L):
#            if k != i and k != j:
#              print >> f, _p(g,i,j), '-', _p(g,i,k), '-', _g(g,k,j), ">= -1"
#              print >> f, _p(g,i,j), '-', _p(g,i,k), '-', _g(h,k,j), ">= -1"
#      for i in xrange(0,L-1):
#        for j in xrange(i+1,L):
#          print >> f, _c(g,j), '+', _p(g,i,j) + "<= 1"

        print >> f, "binary"
        print >> f, '\n'.join(cycles)
        print >> f, '\n'.join([_g(g, i, j) for g in xrange(0, self.T.N - 1) for i in xrange(0, L - 1) for j in xrange(i + 1, L)])
#    print >> f, '\n'.join([_p(g, i, j) for g in xrange(0,self.T.N-1) for i in xrange(0,L-1) for j in xrange(i+1,L)])
        print >> f, "bounds"
        for g in xrange(0, self.T.N - 1):
            for i in xrange(0, L):
                print >> f, _u(g, i), '<=', L - 1
        print >> f, "end"
        f.close()

    def draw_tree (self, treewidth, dy, lengths=False):
        n, N = self.T.n, self.T.N
        height = N*[0]
        for i in xrange(n,N):
            height[i] = max(height[self.T.left[i]], height[self.T.right[i]]) + 1
        new()
        posx, posy = N*[0], N*[0]
        for i in xrange(0,n):
            posy[i] = i*dy*1j
            labelr(posy[i], "\\it "+self.T.name[i])
            #print i, height[i], self.T.level[i]
            #labelr(treewidth*height[i]+posy[i], self.T.name[i])
        for i in xrange(n,N):
            l, r, h, lev = self.T.left[i], self.T.right[i], height[i], self.T.level[i]
            #print i, h, lev
            posx[i] = -float(treewidth) / (lev+h) * h
            posy[i] = (posy[l]+posy[r]) / 2.0
            line (posx[i]+posy[l], posx[l]+posy[l])
            line (posx[i]+posy[r], posx[r]+posy[r])
            line (posx[i]+posy[l], posx[i]+posy[r])
            if lengths:
                setcolor (color.grey(1))
                if i == N-1:
                    p = posx[i]+(posy[l]+posy[r])/2
                    node (p, 0.15, str(self.g[i].dist(self.g[self.T.left[i]])
                                     + self.g[i].dist(self.g[self.T.right[i]])))
                else:
                    node ((posx[i]+posx[l])/2 +posy[l], 0.15, str(self.g[i].dist(self.g[self.T.left[i]])))
                    node ((posx[i]+posx[r])/2 +posy[r], 0.15, str(self.g[i].dist(self.g[self.T.right[i]])))
                setcolor (color.grey(0))
            #line (treewidth*height[i]+posy[l], treewidth*height[l]+posy[l])
            #line (treewidth*height[i]+posy[r], treewidth*height[r]+posy[r])
            #line (treewidth*height[i]+posy[l], treewidth*height[i]+posy[r])
        save("T.pdf")

    # verzia pre veroniku 
    def opt_neigh2 (self):
        "Najde historiu v okoli s najmensim skore"
        if len(self.c[self.T.left[self.T.N - 1]]) < len(self.c[self.T.right[self.T.N - 1]]):
            self.c[self.T.N - 1] = self.c[self.T.left[self.T.N - 1]]
        else:
            self.c[self.T.N - 1] = self.c[self.T.right[self.T.N - 1]]
        B, M = [], []
        print [len(self.c[i]) for i in xrange(self.T.n, self.T.N)]
        for u in xrange(0, self.T.n):
            B.append(len(self.c[u]) * [0])
            M.append([])
            for i in xrange(0, len(self.c[u])):
                M[u].append([0, 0])
        for u in xrange(self.T.n, self.T.N):
            v, w = self.T.left[u], self.T.right[u]
            ul, vl, wl = len(self.c[u]), len(self.c[v]), len(self.c[w])
            B.append([]); M.append([])
            for i in xrange(0, ul):
                B[u].append(penalty(self.c[u][i].numch()))
                M[u].append([0, 0])
            chl = []
            chl.extend(self.c[v])
            chl.extend(self.c[w])
            DM = dist_matrix(self.c[u], chl) 
            for i in xrange(0, ul):
                m, mi, q = 999999, 0, 1
                for j in xrange(0, vl):
                    t = B[v][j] + DM[i][j]
                    if t < m: m, mi = t, j # , q = t, j, ns[v][j]
                B[u][i], M[u][i][0] = B[u][i] + m, mi #, ns[u][i] = B[u][i]+m, mi, q
                m, mi, q = 999999, 0, 1
                for j in xrange(0, wl):
                    t = B[w][j] + DM[i][vl + j]
                    if t < m: m, mi = t, j # , q = t, j, ns[w][j]
                B[u][i], M[u][i][1] = B[u][i] + m, mi
        m, mi, q = 999999, 0, 1
        for i in xrange(0, len(self.c[self.T.N - 1])):
            if B[self.T.N - 1][i] < m: m, mi = B[self.T.N - 1][i], i #, q = B[self.T.N-1][i], i, ns[self.T.N-1][i]
        self._replace (M, self.T.N - 1, mi)
        if m != self.score(): print "BUGGGGGGGGGGG", m, self.score()
        return m



#
# def lower_bound (self):
#   global self.T.n, self.T.N, self.T.left, self.T.right
#   desc,  = []
#   for i in xrange(0,self.T.n):
#     desc.append([i,])
#   for i in xrange(self.T.n,self.T.N):
#     l, r = self.T.left[i], self.T.right[i]
#     for j in xrange(0,len(desc[l])):
#       for k in xrange(0,len(desc[r])):
#
#         self.c[i].extend(self.g[desc[l][j]].path(self.g[desc[r][k]]))
#     desc.append (desc[l]+desc[r])



#  def fill_in_rand (self):
#    global self.T.n, self.T.N
#    g = range(1,27)
#    for i in xrange(self.T.n,self.T.N):
#      for j in xrange(0,len(g)):
#        g[j] = g[j]*choice([-1,1])
#      shuffle(g)
#      s = " ".join([str(x) for x in g]) + " " + choice('$@')
#      self.g[i] = Genome(s)
