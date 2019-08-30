from itertools import combinations
from .base_class import Operator

def steenrod_diagonal(i, n):
    '''...'''
    if type(n) != int or n < 0 :
        raise ValueError('The first entry (dimension) must be a non-negative integer')
    
    elif type(i) != int :
        raise ValueError('The second entry (i) must be an integer')

    elif n < i or i < 0 :
        return set()
    
    answer = set()
    for U in combinations(range(n+1), n-i):
        
        U_minus, U_plus = [], []
        for u in U :
            if (U.index(u) + u) % 2 == 1:
                U_minus.append(u)
            else :
                U_plus.append(u)
        
        answer ^= {(Operator(face_maps = U_minus), 
                    Operator(face_maps = U_plus))}
    
    return answer