from itertools import product
from .base_class import Operator

def alexander_whitney(n, q=None):
    '''if an integer n is passed it returns the linear combination of bioperators defining the 
    restriction of AW to degree n. If two integers (p, q) are passed it provides the single 
    bioperad defining the projection to bidegree (p, q) of the AW map'''
    
    if not q is None:
        p = n
        return (Operator(face_maps = range(p+1,p+q+1)), 
                Operator(face_maps = range(0,p)))
    
    # dictionary of bioperators indexed by (-n,0), (-n+1, -1), ... , (0,-n)
    if q == None:
        answer = {}
        for i in range(n+1):
            answer[(-n+i, -i)] = (Operator(face_maps = range(i+1,n+1)), 
                                  Operator(face_maps = range(0,i)))

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
            return {(Operator(), Operator())}

        answer = set()
        if p > 0:
            west = eilenberg_zilber(p-1, q)
            for biop in west:
                answer ^= {(Operator(deg_maps = biop[0].deg_maps),
                            Operator(deg_maps = [p+q-1] + list(biop[1].deg_maps)))}
        if q > 0:
            south = eilenberg_zilber(p, q-1)
            for biop in south:
                answer ^= {(Operator(deg_maps = [p+q-1] + list(biop[0].deg_maps)),
                            Operator(deg_maps = biop[1].deg_maps))}
        return answer
    
    # dictionary of all bioperators indexed by their bidegrees
    if all_bidegrees:
        if q:
            n += q
            
        answer = {(0,0): {(Operator(), Operator())}}
        for m in range(1, n+1):
            answer.update({(p,m-p): set() for p in range(m+1)})
            for i in range(m):
                j = m-i-1
                for biop in answer[(i,j)]:
                    answer[(i,j+1)] ^= {(Operator(deg_maps = biop[0].deg_maps),
                                         Operator(deg_maps = [i+j] + list(biop[1].deg_maps)))}
                    
                    answer[(i+1,j)] ^= {(Operator(deg_maps = [i+j] + list(biop[0].deg_maps)),
                                         Operator(deg_maps = biop[1].deg_maps))}
        return answer
    
    # dictionary of bioperators indexed by bidegrees (0,n), (1,n-1), ... , (n,0)
    if q == None and not all_bidegrees:
        if n == 0:
            return {(0,0): {(Operator(), Operator())}}

        if n > 0:
            answer = {(i,n-i): set() for i in range(n+1)}
            for bidegree, biops in eilenberg_zilber(n-1).items():
                i, j = bidegree
                answer[(i,j+1)] ^= {(Operator(deg_maps = biop[0].deg_maps),
                                     Operator(deg_maps = [i+j] + list(biop[1].deg_maps)))
                                         for biop in biops}

                answer[(i+1,j)] ^= {(Operator(deg_maps = [i+j] + list(biop[0].deg_maps)),
                                     Operator(deg_maps = biop[1].deg_maps))
                                         for biop in biops}
            return answer

def shih(n):
    '''returns all bioperators defining the chain homotopy between 
    EZAW and the identity. Some of them are degenerate.'''
    
    if n == 0:
        return set()
    
    ez = eilenberg_zilber(n)
    aw = alexander_whitney(n)
    ezaw = set()
    for i in range(n+1):
        a0, a1 = aw[(-i,-n+i)]
        ezaw ^= {(op0.compose(a0), op1.compose(a1)) for op0, op1 in ez[(i,n-i)]}
        
    s_0 = Operator(deg_maps=[0])
    answer = {(op0.prime().compose(s_0), 
               op1.prime().compose(s_0)) for op0, op1 in ezaw}
    
    if n == 1:
        return answer

    if n > 1:
        return answer^{(op0.prime(), op1.prime()) for op0, op1 in shih(n-1)}