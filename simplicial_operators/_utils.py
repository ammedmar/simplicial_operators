from itertools import combinations, product, chain, permutations
from math import floor

def partitions(n, k, smallest_value=1, largest_value=None, ordered=False):
    '''n is the integer to partition and k is the length of partitions.
    It returns all k tuples of integers greater or equal to smallest_value 
    and less than or equal to largest_value that add to to n. 
    If ordered == True it returns all tuples if False it returns those 
    in non-decreassing order '''
    if largest_value == None:
        largest_value = n
    
    def unordered_partitions(n,k,l=smallest_value, m=largest_value):
        if k == 1:
            if l <= n <= m:
                yield (n,)
            return
        for i in range(l,m+1):
            for result in unordered_partitions(n-i,k-1,i,m):
                yield (i,)+result
                
    if ordered:
        return chain.from_iterable(set(permutations(p)) for p 
                                   in unordered_partitions(n,k))
    if not ordered:
        return unordered_partitions(n,k)
    
def harmonic_partitions(n,k,smallest_value=1):
    '''returns tuple (a_1,...,a_k) of non-negative integer greater than 
    or equal to smallest_value such that a_1 + 2a_2 + ... + ka_k = n'''
    if k == 1:
        if n >= smallest_value:
            yield (n,)
        return
    for i in range(smallest_value, floor(n/k)+1):
        for result in harmonic_partitions(n-(k*i), k-1, smallest_value):
            yield (i,)+result