#!/usr/bin/env python
# coding: utf-8

# In[26]:


from itertools import combinations, product

class Operator(object):
    '''
    Models a simplicial operator of the form:

    s ... s d ... d

    represented in the canonical form

    s > ... > s d < ... < d

    Parameters
    ----------
    deg_maps  : tuple or list
                An ordered collection of integers representing the degeneracy maps of the operator

    face_maps : tuple or list
                An ordered collection of integers representing the face maps of the operator

    '''

    def __init__(self, deg_maps = [], face_maps = []):

        self.deg_maps  = Operator.deg_maps_sort(deg_maps)
        self.face_maps = Operator.face_maps_sort(face_maps)

    def __repr__(self):

        s, d = repr(self.deg_maps), repr(self.face_maps)

        return f'Operator(deg_maps={s}, face_maps={d})'

    def __str__(self):

        if not self.deg_maps and not self.face_maps:

            return 'id'

        if not self.deg_maps:
            s = ''
        else:
            s = f's_{"s_".join(str(i) for i in self.deg_maps)}'

        if not self.face_maps:
            d = ''
        else:
            d = f'd_{"d_".join(str(i) for i in self.face_maps)}'

        return(s+d)

    def __call__(self, simplex):
        '''applies the operator to a simplex represented by a list or a tuple'''

        simplex = list(simplex)

        if self.face_maps and len(simplex)-1 < self.face_maps[-1]:
            raise ValueError('simplex not in the domain of the operator')

        # face maps
        simplex = ([v for i, v in enumerate(simplex)
                    if i not in self.face_maps])

        # degeneracy maps
        for i in sorted(self.deg_maps):
            try:
                simplex.insert(i+1, simplex[i])
            except IndexError:
                raise ValueError('simplex not in the domain of the operator')

        return simplex

    @property
    def degree(self):
        '''returns the degree of the operator as an int'''

        return len(self.deg_maps) - len(self.face_maps)

    @property
    def is_degenerate(self):
        '''returns True if the operator is degenerate and False if not'''

        return bool(self.deg_maps)

    def compose(self, other):
        '''returns the operator self other'''

        s_left, d_left  = list(self.deg_maps), list(self.face_maps)
        s_right, d_right = list(other.deg_maps), list(other.face_maps)

        # getting s_right passed d_left
        d = d_left
        s = []
        for j in s_right:
            if j+1 in d or j in d:
                p = max(filter(lambda k: k == j or k == j+1, d))
                d = d[:d.index(p)] + [k-1 for k in d[d.index(p)+1:]]

            else:
                lt = list(filter(lambda k: k < j, d))
                gt = list(filter(lambda k: k > j+1, d))
                d = lt + [k-1 for k in gt]
                s.append(j - len(lt))

        # One can maybe optimize the following using that
        # s_left, s, d, and d_right are totally ordered

        new_deg_maps  = Operator.deg_maps_sort(s_left + s)
        new_face_maps = Operator.face_maps_sort(d + d_right)

        return Operator(new_deg_maps, new_face_maps)

    def prime(self):
        '''adds 1 to all the numbers defining the face and degeneracy maps 
        of the operator'''
        
        return Operator((v+1 for v in self.deg_maps),
                        (w+1 for w in self.face_maps))
    
    @staticmethod
    def deg_maps_sort(deg_maps):
        '''puts the degeneracy maps in canonical order s > ... > s using
        the simplicial identity s_i s_j = s_{j+1} s_i if i <= j'''

        deg_maps = list(deg_maps)
        for index in reversed(range(len(deg_maps)-1)):

            currentvalue = deg_maps[index]
            position = index

            while (position<len(deg_maps)-1 and
                   currentvalue <= deg_maps[position+1]):

                deg_maps[position] = deg_maps[position+1]+1
                position = position+1

            deg_maps[position] = currentvalue

        return tuple(deg_maps)

    @staticmethod
    def face_maps_sort(face_maps):
        '''puts the face maps in canonical order d < ... < d using the
        simplicial identity d_i d_j = d_j d_{i+1} if i >= j'''

        face_maps = list(face_maps)
        for index in range(1, len(face_maps)):

            currentvalue = face_maps[index]
            position = index

            while (position > 0 and
                   face_maps[position-1] >= currentvalue):

                face_maps[position] = face_maps[position-1]+1
                position = position-1

            face_maps[position] = currentvalue

        return tuple(face_maps)
    

class Bioperator:
    '''
    Models a simplicial bioperator of form:

    s ... s d ... d (x) s ... s d ... d

    represented in the canonical form

    s > ... > s d < ... < d (x) s > ... > s d < ... < d
    '''

    def __init__(self, deg_maps1=[], face_maps1=[], 
                       deg_maps2=[], face_maps2=[]):
        
        if isinstance(deg_maps1, Operator) and isinstance(face_maps1, Operator):
            self.op1 = deg_maps1
            self.op2 = face_maps1
        else:
            self.op1 = Operator(deg_maps1, face_maps1)
            self.op2 = Operator(deg_maps2, face_maps2)

    def __repr__(self):

        s1, d1 = self.op1.deg_maps, self.op1.face_maps
        s2, d2 = self.op2.deg_maps, self.op2.face_maps

        return (f'Bioperator(deg_maps1={s1}, face_maps1={d1}, '
                         + f'deg_maps2={s2}, face_maps2={d2})')

    def __str__(self):

        return(str(self.op1) + ' x ' + str(self.op2))

    def __call__(self, spx1, spx2=None):
        
        if spx2 is None:
            return self.op1(spx1), self.op2(spx1)
        else:
            return self.op1(spx1), self.op2(spx2)
    
    @property
    def bidegree(self):
        '''returns the bidegree of the operator as a pair of int'''

        return self.op1.degree(), self.op2.degree()
    
    @property
    def is_degenerate(self):
        '''returns True if the operator is degenerate and False if not'''

        s1 = set(self.op1.deg_maps)
        s2 = set(self.op2.deg_maps)

        return bool(s1.intersection(s2))

    def compose(self, other):
        '''returns the operator self other'''

        op1 = self.op1.compose(other.op1)
        op2 = self.op2.compose(other.op2)

        return Bioperator(op1.deg_maps, op1.face_maps,
                          op2.deg_maps, op2.face_maps)
    
    def prime(self):
        '''adds 1 to all the numbers defining the face and degeneracy maps 
        of the operator'''
        
        return Bioperator(self.op1.prime(), self.op2.prime())


####### EZ-AW #########

def alexander_whitney(n, q=None):
    '''if an integer n is passed it returns the linear combination of bioperators defining the 
    restriction of AW to degree n. If two integers (p, q) are passed it provides the single 
    bioperad defining the projection to bidegree (p, q) of the AW map'''
    
    if not q is None:
        p = n
        return Bioperator(face_maps1 = range(p+1,p+q+1), 
                          face_maps2 = range(0,p))
    
    # dictionary of bioperators indexed by (-n,0), (-n+1, -1), ... , (0,-n)
    if q == None:
        answer = {}
        for i in range(n+1):
            answer[(-n+i, -i)] = Bioperator(face_maps1 = range(i+1,n+1), 
                                            face_maps2 = range(0,i))

        return answer
    
def eilenberg_zilber(n, q=None, all_bidegrees=False):
    '''if a single integer n is passed, it returns a dictionary whose keys are bidegrees adding to n
    and values the bioperators defining the restriction of EZ to that bidegree. If all_bidegrees 
    is True it gives the same but the condition is that bidegrees add up to less than or equal to n.
    
    If a pair of integers p,q is passed it returns the bioperators defining EZ in bidegree (p,q). 
    If all_bidegrees is true it acts as if a single integer p+q was passed.'''
    
    # operators in bidegree (q,p) acting on elements of bidegree (p,q)
    if q != None:
        p = n
        if (p,q) == (0,0):
            return {Bioperator()}

        answer = set()
        if p > 0:
            west = eilenberg_zilber(p-1, q)
            for biop in west:
                answer ^= {Bioperator(
                                deg_maps1 = biop.op1.deg_maps,
                                deg_maps2 = [p+q-1] + list(biop.op2.deg_maps))}
        if q > 0:
            south = eilenberg_zilber(p, q-1)
            for biop in south:
                answer ^= {Bioperator(
                                deg_maps1 = [p+q-1] + list(biop.op1.deg_maps),
                                deg_maps2 = biop.op2.deg_maps)}
        return answer
    
    # dictionary of all bioperators indexed by bidegrees
    # (memory intensive)
    if all_bidegrees:
        if q:
            n += q
            
        answer = {(0,0): {Bioperator()}}
        for m in range(1, n+1):
            answer.update({(p,m-p): set() for p in range(m+1)})
            for i in range(m):
                j = m-i-1
                for biop in bioperators[(i,j)]:
                    answer[(i,j+1)] ^= {Bioperator(
                                            deg_maps1 = biop.op1.deg_maps,
                                            deg_maps2 = [i+j] + list(biop.op2.deg_maps))}

                    answer[(i+1,j)] ^= {Bioperator(
                                            deg_maps1 = [i+j] + list(biop.op1.deg_maps),
                                            deg_maps2 = biop.op2.deg_maps)}
            
        return bioperators
    
    # dictionary of bioperators indexed by bidegrees (0,n), (1,n-1), ... , (n,0)
    if q == None:
        if n == 0:
            return {(0,0) : {Bioperator()}}

        if n > 0:
            answer = {(i,n-i): set() for i in range(n+1)}
            for bidegree, biops in eilenberg_zilber(n-1).items():
                i, j = bidegree
                answer[(i,j+1)] ^= {Bioperator(
                                        deg_maps1 = biop.op1.deg_maps,
                                        deg_maps2 = [i+j] + list(biop.op2.deg_maps))
                                            for biop in biops}

                answer[(i+1,j)] ^= {Bioperator(
                                        deg_maps1 = [i+j] + list(biop.op1.deg_maps),
                                        deg_maps2 = biop.op2.deg_maps)
                                            for biop in biops}
            return answer

def shih(n):
    '''returns all bioperators defining the chain homotopy between 
    EZAW and the identity. Some of them are degenerate.'''
    
    if n == 0:
        return set()
    
    s_0s_0 = Bioperator(deg_maps1=[0], deg_maps2=[0])
    
    ez = eilenberg_zilber(n)
    aw = alexander_whitney(n)
    ezaw = set()
    for i in range(n+1):
        ezaw ^= {biop.compose(aw[(-i,-n+i)]) for biop in ez[(i,n-i)]}
        
    answer = {biop.prime().compose(s_0s_0) for biop in ezaw}
    
    if n == 1:
        return answer

    if n > 1:
        return answer^{biop.prime() for biop in shih(n-1)}
        

######## STEENROD ########

def steenrod_diagonal(n,i):
    '''...'''
    if type(n) != int or n < 0 :
        raise ValueError('The first entry (dimension) must be a non-negative integer')
    
    elif type(i) != int :
        raise ValueError('The second entry (i) must be an integer')

    elif n < i or i < 0 :
        return set() # no operators model for 0 operator
    
    answer = set()
    for U in combinations(range(n+1), n-i):
        
        U_minus, U_plus = [], []
        for u in U :
            if (U.index(u) + u) % 2 == 1:
                U_minus.append(u)
            else :
                U_plus.append(u)
        
        answer ^= {Bioperator(face_maps1 = U_minus, 
                              face_maps2 = U_plus)}
    
    return answer

