

class SimplicialOperator(object):
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

        self.deg_maps  = SimplicialOperator.deg_maps_sort(deg_maps)
        self.face_maps = SimplicialOperator.face_maps_sort(face_maps)

    def __repr__(self):

        s, d = repr(self.deg_maps), repr(self.face_maps)

        return f'SimplicialOperator(deg_maps={s}, face_maps={d})'

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

        new_deg_maps  = SimplicialOperator.deg_maps_sort(s_left + s)
        new_face_maps = SimplicialOperator.face_maps_sort(d + d_right)

        return SimplicialOperator(new_deg_maps, new_face_maps)

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


class SimplicialBioperator:
    '''
    Models a simplicial bioperator of form:

    s ... s d ... d (x) s ... s d ... d

    represented in the canonical form

    s > ... > s d < ... < d (x) s > ... > s d < ... < d
    '''

    def __init__(self, deg_maps1=[], face_maps1=[], deg_maps2=[], face_maps2=[]):

        self.op1 = SimplicialOperator(deg_maps1, face_maps1)
        self.op2 = SimplicialOperator(deg_maps2, face_maps2)

    def __repr__(self):

        s1, d1 = self.op1.deg_maps, self.op1.face_maps
        s2, d2 = self.op2.deg_maps, self.op2.face_maps

        return (f'Bioperator(deg_maps1={s1}, face_maps1={d1}, '
                         + f'deg_maps2={s2}, face_maps2={d2}')

    def __str__(self):

        return(str(self.op1) + ' x ' + str(self.op2))

    def __call__(self, spx1, spx2):

        return self.op1(spx1), self.op2(spx2)

    def bidegree(self):
        '''returns the bidegree of the operator as a pair of int'''

        return self.op1.degree(), self.op2.degree()

    def is_degenerate(self):
        '''returns True if the operator is degenerate and False if not'''

        s1 = set(self.op1.deg_maps)
        s2 = set(self.op2.deg_maps)

        return bool(s1.intersection(s2))

    def compose(self, other):
        '''returns the operator self other'''

        op1 = self.op1.compose(other.op1)
        op2 = self.op2.compose(other.op2)

        return SimplicialBioperator(op1.deg_maps, op1.face_maps,
                                    op2.deg_maps, op2.face_maps)
