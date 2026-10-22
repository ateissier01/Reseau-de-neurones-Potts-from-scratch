import math
import pickle

class pot_neurone :
    def __init__(self,activate) :
        self.activate=activate
        self.state=0
    
    def output(self,field_matrix) :
        h=[]
        for L in field_matrix :
            h.append(0)
            for W in L :
                h[-1]+=W
        self.state=self.activate(h)
        return self.state
        
    

class Combination:
    def __init__(self, *values):
        self.values = tuple(values)  # ordre conservé

    def __eq__(self, other):
        if not isinstance(other, Combination):
            return False
        return self.values == other.values

    def __hash__(self):
        return hash(self.values)  # hashable comme un tuple

    def __repr__(self):
        return f"Comb{self.values}"
    
    def orga(self) :
        l=[]
        for v in self.values :
            l.append(v)
        return l



class PairMap:
    """
    Associe une valeur à une paire de combinaisons.
    L’ordre de la paire n’a pas d’importance :
    (a,b) == (b,a)
    """
    def __init__(self):
        self._store = {}

    def _key(self, a, b):
        # clé = frozenset de 2 Combination
        return frozenset([a, b])

    def associer(self, a, b, valeur):
        """Associer une valeur à la paire (a,b)."""
        self._store[self._key(a, b)] = valeur

    def retrouver(self, a, b, default=None):
        """Retrouver la valeur associée à (a,b), ou default si absente."""
        return self._store.get(self._key(a, b), default)

    def __repr__(self):
        return repr(self._store)

    def getmatrix(self) :
        matrix=[]
        Max=0
        
        c=0
        Map={}
        for k in list(self._store.keys()) :
            for n in k :
                if not n in list(Map.keys()) :
                    Map[n]=c
                    c+=1
                for i in n.values :
                    if i>Max :
                        Max=i
        
        Max=Max+1
        Max=(Max**2)
        
        
        for i in range(Max) :
            matrix.append([])
            for j in range(Max) :
                matrix[i].append(0)
                
                
        for k in list(self._store.keys()) :
            val=self._store[k]
            key=[]
            for n in k :
                key.append(n)
            matrix[Map[key[0]]][Map[key[1]]]=val
            matrix[Map[key[1]]][Map[key[0]]]=val
        
        return matrix

class UndirectedGraph:
    def __init__(self):
        self.adj = {}       # dictionnaire d’adjacence

    def add_edge(self, a, b):
        self.adj.setdefault(a, set()).add(b)
        self.adj.setdefault(b, set()).add(a)

    def voisins(self, a):
        return self.adj.get(a, set())


class pott_network :
    def __init__(self, size,activate,Tokens=2,maxcodist=float('inf')) :
        self.neurone=[]
        self.weight=PairMap()
        self.nei=UndirectedGraph()
        self.limit=size-1
        self.size=(self.limit+1)**2
        self.learn_count = 0          # nombre de motifs appris
        self.number_token=Tokens
        self.maxdist=maxcodist
        for i in range(size) :
            self.neurone.append([])
            for j in range(size) :
                n=Combination(i,j)
                for x in range(size) :
                    for y in range(size) :
                        m=Combination(x,y)
                        if not n.__eq__(m) and math.sqrt((i-x)**2+(j-y)**2)<self.maxdist :
                            comatrix=[]
                            for k in range(self.number_token):
                                mat=[]
                                for q in range(self.number_token):
                                    mat.append(0)
                                comatrix.append(mat)
                            self.weight.associer(n, m,comatrix)
                            self.nei.add_edge(n,m)
                self.neurone[i].append(pot_neurone(activate))

    def add_hebb(self,motif) :
        self.learn_count+=1
        for i in range(len(motif)) :
            for j in range(len(motif[i])) :
                n=Combination(i,j)
                for x in range(len(motif)) :
                    for y in range(len(motif[i])) :
                        m=Combination(x,y)
                        if not n.__eq__(m) and math.sqrt((i-x)**2+(j-y)**2)<self.maxdist :
                            self.weight.retrouver(n,m)[motif[i][j]][motif[x][y]]+=1
    

    def retrieve(self,motif,maxiter,stoppingfactor) :
        state = [row[:] for row in motif]
        size = self.limit + 1

        for iteration in range(maxiter):
            prev_state = [row[:] for row in state]
            new_state = [[0 for _ in range(size)] for _ in range(size)]

            # --- mise à jour synchrone ---
            for i in range(size):
                for j in range(size):
                    comb = Combination(i, j)
                    voisins = self.nei.voisins(comb)
                    field_matrix=[]
                    for T in range(self.number_token) :
                        field_matrix.append([])
    
                    for v in voisins :
                        for x,y in [v.values] :
                            W=self.weight.retrouver(comb, v)
                            S=state[x][y]
                            for T in range(self.number_token) :
                                field_matrix[T].append(W[T][S])

                    
                    new_state[i][j] = self.neurone[i][j].output(field_matrix)

            # --- appliquer simultanément ---
            state = new_state

            # --- vérifier convergence ---
            diff = sum((state[i][j] - prev_state[i][j])**2 for i in range(size) for j in range(size))
            if diff < stoppingfactor:
                print("dif break : diff = {diff}")
                break
        
        return state
    
    def save(self,directory) :
        with open(directory, "wb") as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    @staticmethod
    def load(directory) :
        with open(directory, "rb") as f:
            return pickle.load(f)

def activation_higherscore(h) :
    higher=0
    index=0
    for i in range(len(h)) :
        if h[i]>higher :
            higher=h[i]
            index=i
    return index
