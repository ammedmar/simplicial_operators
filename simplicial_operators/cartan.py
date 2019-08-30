from .base_class import Operator
from .aw_ez_shih import shih
from .barratt_eccles import barratt_eccles_operator

def cartan_first_homotopy(n):
    '''applies the first homotopy to the element 
    \tilde x_n = (e, (12), ..., (12)^n)

    \sum_{i=0}^n ( (23)f(e), \dots,  (23)f((12)^i), g((12)^i), \dots, g((12)^n))

    (23)f(12) = (23)(13)(24) = (2,4,1,3)
        g(12) = (12)(34)     = (2,1,4,3)

    '''

    # permutations (23), (12)(34)(23) and e, (12)(34) 
    a = {0: (1,3,2,4), 1: (2,4,1,3)}
    b = {0: (1,2,3,4), 1: (2,1,4,3)}

    x = tuple(a[i%2] for i in range(n+1))
    y = tuple(b[i%2] for i in range(n+1))

    answer = set()
    for i in range(n+1):
        answer ^= {x[:i+1]+y[i:]}

    return answer

def cartan_second_homotopy(n):
    '''
    applies the second homotopy to the element 
    \tilde x_n = (e, (12), ..., (12)^n)

    \circ_{\mathcal E}( (e, \dots, e) \otimes SHI (x \otimes x) )
    '''

    sigma_two = {0: (1,2), 1: (2,1)}
    x = [sigma_two[i%2] for i in range(n+1)]

    if n == 0:
        return set()

    # composition (e, x, y) with x,y in sigma_two
    composition = {((1,2), (2,1)) : (1,2,4,3),
                   ((2,1), (1,2)) : (2,1,3,4),
                   ((2,1), (2,1)) : (2,1,4,3),
                   ((1,2), (1,2)) : (1,2,3,4)}

    values = set()
    for biop in shih(n):
        if not Operator.is_degenerate(biop):
            values ^= {(biop[0](x), biop[1](x))}

    answer = set()
    for value in values:

        table = tuple()
        for i in range(n+2):
            table += tuple((composition[(value[0][i], value[1][i])],))

        answer ^= {table}

    return answer

def cartan_operator(i, n):
    '''
    it returns the multioperators defining the i-th cartan coboundary in degree n when 
    applied to homogeneous cocycles
    '''
    if i >= n:
        return set()
    
    first  = cartan_first_homotopy(i)
    second = cartan_second_homotopy(i)
    all_operators = barratt_eccles_operator(first^second, n)
    
    filtered_by_degree = {multiop for multiop in all_operators if 
                          multiop[0].degree == multiop[1].degree and 
                          multiop[2].degree == multiop[3].degree}
    
    return filtered_by_degree