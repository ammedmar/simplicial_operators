from .base_class import Operator

def surjection_operator(surjections, n):
    '''returns the set of multioperators representing the action of the passed 
    set of surjection on a n-simplex'''
    
    def _new_term(term, num_to_append, pos_to_append, k):
            tuple_to_replace = term['seq'][pos_to_append]+(num_to_append,)

            return {'pos': term['pos']+k, 
                    'seq': term['seq'][:pos_to_append] 
                            + (tuple_to_replace,)
                            + term['seq'][pos_to_append+1:]}
        
    if not isinstance(surjections, set):
        surjections = {surjections}
     
    for surj in surjections:
        after = [{'pos': 0, 
                  'seq': ((),)*(surj[0]-1) + ((0,),) + ((),)*(max(surj)-surj[0])}]

        for i in range(n+len(surj)-1):
            before = after[:]
            after = []
            for term in before:
                if term['pos'] < len(surj)-1:
                    num_to_append = term['seq'][surj[term['pos']]-1][-1]
                    pos_to_append = surj[term['pos']+1]-1

                    empty = bool(term['seq'][pos_to_append])
                    if not empty or num_to_append != term['seq'][pos_to_append][-1]:
                        after.append(_new_term(term, num_to_append, pos_to_append,1))

                if term['seq'][surj[term['pos']]-1][-1] < n:
                    num_to_append = term['seq'][surj[term['pos']]-1][-1] + 1
                    pos_to_append = surj[term['pos']]-1

                    after.append(_new_term(term, num_to_append, pos_to_append,0))

        answer = set()
        for term in after:
            operators = tuple()
            for seq in term['seq']:
                operators += (Operator(face_maps=[i for i in range(n+1) if i not in seq]),)
            answer ^= {operators}

    return answer