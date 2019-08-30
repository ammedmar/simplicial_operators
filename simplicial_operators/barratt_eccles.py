from ._utils import partitions
from .surjections import surjection_operator

def table_reduction(bar_ecc_elements):
    '''given a set of basis element in the Barratt-Eccles operad, it returns the set of 
    surjections in its image via the table reduction morphism'''
    
    if not isinstance(bar_ecc_elements, set):
        bar_ecc_elements = {bar_ecc_elements}
    
    answer = set()
    for bar_ecc_element in bar_ecc_elements:
        
        d, a = len(bar_ecc_element)-1, max(bar_ecc_element[0]) #dimension and arity

        for pi in partitions(d+a, d+1, ordered=True):
            
            surjection, removed = [], []
            degenerate = False
            for idx, i in enumerate(pi):
                filtered =  [i for i in bar_ecc_element[idx] if i not in removed]

                if idx > 0 and surjection[-1] == filtered[0]:
                    degenerate = True
                    break    

                if i > 1:
                    removed += filtered[:i-1]

                surjection += filtered[:i]

            if not degenerate:
                answer ^= {tuple(surjection)}
 
    return answer

def barratt_eccles_operator(bar_ecc_elements, n):
    '''returns the multioperator defining the action of a barratt-eccles 
    element on simplices of dimension d'''
    
    if not isinstance(bar_ecc_elements, set):
        bar_ecc_elements = {bar_ecc_elements}
        
    answer = set()
    for bar_ecc_element in bar_ecc_elements:
        surjections = table_reduction(bar_ecc_element)
        for surjection in surjections:
            answer ^= surjection_operator(surjection, n)
        
    return answer