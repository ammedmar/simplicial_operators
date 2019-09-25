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

        return tuple(simplex)

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
    
    @staticmethod
    def is_degenerate(multiop):
        '''returns True if a multioperator is degenerate'''
        if isinstance(multiop, Operator):
            return multiop.is_degenerate
        if isinstance(multiop, tuple):
            deg = set(multiop[0].deg_maps)
            for op in multiop:
                deg = deg.intersection(set(op.deg_maps))
            return bool(deg)
        else:
            raise TypeError('expected either Operator or tuple of Operators')
    
    @staticmethod
    def is_nondegenerate(multiop):
        return not Operator.is_degenerate(multiop)
    
    @staticmethod
    def action(multiop, lincomb):
        '''modeling the action of operators'''
        # single Operator and single simplex
        if isinstance(multiop, Operator) and not isinstance(lincomb, set):
            return multiop(lincomb)
        
        # single Operator and set of simplices
        if isinstance(multiop, Operator) and isinstance(lincomb, set):
            answer = set()
            for spx in lincomb:
                answer ^= {multiop(spx)}
            return answer
        
        # single multioperator and single multisimplex
        if isinstance(multiop, tuple) and not isinstance(lincomb, set):
            if len(multiop) == len(lincomb): 
                return tuple(op(spx) for op, spx in zip(multiop, lincomb))
            else:
                raise TypeError('arities do not match')
        
        # single multioperator and set of multisimplices
        if isinstance(multiop, tuple) and isinstance(lincomb, set):
            answer = set()
            for spx in lincomb:
                answer ^= {Operator.action(multiop, spx)}
            return answer
        
        # set of Operators or multioperators and single spx or multisimplex
        if isinstance(multiop, set) and not isinstance(lincomb, set):
            try:
                return Operator.action(multiop, {tuple(lincomb)})
            except TypeError:
                raise(TypeError('Cannot be acted on'))
        
        # set of Operators or multioperators and set of simplices or multisimplices
        if isinstance(multiop, set) and isinstance(lincomb, set):
            answer = set()
            for mop in multiop:
                answer ^= Operator.action(mop, lincomb)
            return answer
        
        else:
            raise TypeError('cannot act: must be Operator, '+
                            'tuple of them, or a set of such tuples')
            
    @staticmethod
    def display(multiop):
        '''tool to aid visualization of operators'''
        
        if isinstance(multiop, Operator):
            return str(multiop)
        
        if isinstance(multiop, tuple):
            string = ''
            for op in multiop:
                string += str(op) + ' x '
            return string[:-3]
        
        if isinstance(multiop, set):
            string = '  '
            for op in multiop:
                string += Operator.display(op) + '\n+ '
            return string[:-3]
        
        else:
            raise TypeError('Expected types: Operator, tuple of '+
                            'Operator, or set of tuple of Operator')
            
    @staticmethod
    def display_action(multiop, lincomb):
        '''modeling the action of operators'''
        string = str(Operator.action(multiop, lincomb))
        
        if '{((' in string:
            string = string.replace(')), ((', ')\n+ (')
            string = string.replace('{(', '  ')
            string = string.replace(')}', '')
            string = string.replace('),', ') x')
            return string
        
        elif '((' in string:
            string = string.replace('), (', ') x (')
            string = string.replace('((', '  (')
            string = string.replace('))', ')')
            return string  
                
        elif '{(' in string:
            string = string.replace('), (', ')\n+ (')
            string = string.replace('{', '  ')
            string = string.replace('}', '')
            return string
        
        elif '(' in string:
            return '  ' + string
        
        else:
            return '0'
        
    # TODO general composition