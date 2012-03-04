import re

class Tree:
    # 0...n-1 su listy, n...N-1 su vnutorne vrcholy
    #n = N = 0
    #name, left, right, parent, level = [], [], [], [], []
    #leaf = {}
    
    def __init__(self, wdir):
        f = open ('data/' + wdir + '/T', 'r')
        s = f.readline ()
        f.close ()
        self._set_tree(s)
        self._levels (self.N - 1, 0)
        #print self.n, self.N

    def _set_tree(self, s):
        #print s
        s = re.sub(r':[0-9.]*', '', s)
        s = s.translate (None, "(.;\n").replace(",", " ").replace(")", " * ").split()
        #print s
        self.n = 0
        self.name = []
        self.left, self.right, self.parent = [], [], []
        self.leaf = {}
        namel, namer = [], []
        for i in xrange(0, len(s)):
            if s[i] == '*':
                s[i] = -1
            else:
                self.name.append(s[i])
                namel.append(s[i])
                namer.append(s[i])
                self.left.append(self.n)
                self.right.append(self.n)
                self.parent.append(self.n)
                self.leaf[s[i]] = self.n
                s[i] = self.n
                self.n += 1
        self.N = self.n
        z = []
        for i in xrange(0, len(s)):
            if s[i] >= 0:
                z.append(s[i])
            else:
                r, l = z.pop(), z.pop()
                namel.append(namel[l])
                namer.append(namer[r])
                self.name.append(namel[self.N] + '-' + namer[self.N])
                self.left.append(l)
                self.right.append(r)
                self.parent.append(self.N)
                self.parent[l] = self.parent[r] = self.N
                z.append(self.N)
                self.N += 1
        self.level = self.N * [0]
        #  parent[left[N-1]] = right[N-1]
        #  parent[right[N-1]] = left[N-1]

    def _levels (self, v, l):
        self.level[v] = l
        if v < self.n: return
        self._levels (self.left[v], l + 1)
        self._levels (self.right[v], l + 1)
