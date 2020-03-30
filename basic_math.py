#
# basic_math.py
#

# helper functions for printing

def intify(a):
    if a.imag == 0:
        a = a.real
        if int(a) == a:
            return int(a)
    return a

def num_print(a):
    x = intify(a.real); y = intify(a.imag)
    if y == 0: 
        return str(x)
    if x == 0: 
        return str(y) + 'i'
    if y < 0: return str(x) + ' - ' + str(abs(y)) + 'i'
    else: return str(x) + ' + ' + str(y) + 'i'

# permutations

def sign(perm):
    copy = perm[:]
    copy.sort()
    assert copy == range(len(perm))
    out = 1
    for i in range(len(perm)):
        for j in range(i):
            if perm[i] < perm[j]:
                out *= -1
    return out

# linear algebra

class vector(tuple):

    def __add__(self, other):
        return self.__class__(map(lambda x,y: x+y, self, other)) # self.__class__ is vector, unless i am a polynomial or something!

    def __neg__(self):
        return self.__class__(map(lambda x: -x, self))

    def __sub__(self, other):
        return self + (-other)
#        return self.__class__(map(lambda x,y: x-y, self, other))

    def __mul__(self, other): # left mul == scalar mul on left
        if isinstance(other, Number):
            return self.__class__(map(lambda x: x*other, self))
        else: raise TypeError
        
    def __rmul__(self, other):
        if isinstance(other, Number): # multiplication by scalars is commutative... so
            return self * other
        elif isinstance(other, matrix):
            p, q = self
            a, b, c, d = other
            # [a b][p] = [ap + bq]
            # [c d][q]   [cp + dq]
            return self.__class__(a*p + b*q, c*p + d*q)
        else: raise TypeError

#    def __repr__(self): 
#        pretty = [num_print(x) for x in self]
#        return '(' + ', '.join(printable) + ')'

def matrix(tuple):

    def __repr__(self): 
        A, B, C, D = (num_print(x) for x in self)
        return ''.join(['[[', A, ', ', B, '], [', C, ', ', D,']]'])

    def cast(self, other):
        if isinstance(other, Number):
            return matrix(other, 0, 0, other)
        else: return other

    def __mul__(self, other):
        other = self.cast(other)
        if isinstance(other, matrix):
            a, b, c, d = self
            w, x, y, z = other
            # [a b][w x] = [aw + by ax + bz]
            # [c d][y z]   [cw + dy cx + dz]
            return matrix(a*w + b*y, a*x + b*z, 
                          c*w + d*y, c*x + d*z)
        else: return other.__rmul__(self)

    def trace(self):
        a, b, c, d = self
        return a + d

    def det(self):
        a, b, c, d = self
        return a*d - b*c

    def unit(self):
        D = self.det()
        return self * intify(sqrt(D)**(-1))

    def inverse(self):
        D = self.det()
        a, b, c, d = self
        # [a b]*[ d -b] = [D 0]
        # [c d] [-c  a]   [0 D]. 
        # One way to remember the formula is to think about how you
        # invert elliptic, parabolic, and hyperbolic elements.  
        return intify(D**(-1)) * matrix(d, -b, -c, a) 
        # NB - I intified the inverse of the determinant because one
        # use case is det == 1.  In that case using D**-1 (ie a float)
        # will lose accuracy even for very modestly sized matrices.

    def transpose(self):
        a, b, c, d = self
        return matrix(a, c, b, d)

    def conjugate_transpose(self):
        return (self.conjugate()).transpose()

    def __call__(self, z): # z is a number
        v = vector(z, 1)
        w = self * v
        return w.frac()

    
#   projective geometry   

class CP1(tuple):

    def complex(self):
        if abs(a[1]) < 0.000001 * abs(a[0]):
            return complex(10,10) ### hack, useful for debugging
        else:
            return a[0]/a[1]

    def preferred_rep(a):
      if abs(a[1]) < 0.00001:
        return a
      else:
        return [a[0]/a[1], complex(1.0,0.0)]   ### saul: perhaps should divide by cmath.sqrt(a[0]*a[1])

  ## g[p_,q_,r_,s_,t_,u_] is a matrix that takes p,q,r to s,t,u

  ## h[z1_, z2_, z3_] := {{z2 - z3, -z1 (z2 - z3)}, {z2 - z1, -z3 (z2 - z1)}};
  ## g[p_, q_, r_, u_, v_, w_] := Inverse[h[u, v, w]].h[p, q, r];


def inf_zero_one_to_triple(p,q,r):
    p1, p2 = p
    q1, q2 = q
    r1, r2 = r
    M = matrix(p1, q1, p2, q2)
    Minv = M.inverse() # matrix2_inv(M)
    mu, lam = Minv * r
    # [mu, lam] = matrix_mult_vector(matrix2_inv([[p1,q1],[p2,q2]]), [r1,r2])
    return matrix(mu*p1, lam*q1, mu*p2, lam*q2)


def two_triples_to_PSL(a1,b1,c1,a2,b2,c2):
    return matrix_mult( inf_zero_one_to_triple(a2,b2,c2), matrix2_inv(inf_zero_one_to_triple(a1,b1,c1) ) ) 


if __name__ == '__main__':
    print sign([0,1,2,3])
    print sign([0,2,1,3])
    print sign([0,1,2,2])
