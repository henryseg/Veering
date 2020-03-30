#
# basic_math.py
#

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

#   linear algebra

class vector(list):
  def __add__(self, other):
    return self.__class__(map(lambda x,y: x+y, self, other)) #self.__class__ is vector, unless i am a polynomial or something!

  def __neg__(self):
    return self.__class__(map(lambda x: -x, self))

  def __sub__(self, other):
    return self.__class__(map(lambda x,y: x-y, self, other))

  def __mul__(self, other): #mult by scalar
    return self.__class__(map(lambda x: x*other, self))

  def __rmul__(self, other):
    return (self*other)

def zero_matrix(m,n):
  """create zero matrix"""
  return [[0 for row_entry in range(n)] for row in range(m)]

def matrix_mult(matrix1,matrix2):
  """matrix multiplication"""
  m1, m2 = matrix1, matrix2
  assert len(m1[0]) == len(m2), 'Matrices must be m*n and n*p to multiply!'
  new_matrix = zero_matrix(len(m1),len(m2[0]))
  for i in range(len(m1)):
    for j in range(len(m2[0])):
      for k in range(len(m2)):
        new_matrix[i][j] += m1[i][k]*m2[k][j]
  return new_matrix 

def matrix_mult_vector(M,v): #treats v as if it were a vertical vector, when it is actually just a list
  u = [[entry] for entry in v]
  out = matrix_mult(M, u)
  return [entry for [entry] in out]

def matrix2_det(M):
  """determinant of 2x2 matrix"""
  return M[0][0]*M[1][1] - M[0][1]*M[1][0]
  
def matrix2_inv(M):
  inv_det = 1.0/matrix2_det(M)
  return [[inv_det*M[1][1],-inv_det*M[0][1]],[-inv_det*M[1][0],inv_det*M[0][0]]]

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
  p1,p2=p
  q1,q2=q
  r1,r2=r
  M = [[p1,q1],[p2,q2]]
  Minv = matrix2_inv(M)
  [mu,lam] = matrix_mult_vector(matrix2_inv([[p1,q1],[p2,q2]]), [r1,r2])
  return [[mu*p1, lam*q1],[mu*p2, lam*q2]]

def two_triples_to_PSL(a1,b1,c1,a2,b2,c2):
  return matrix_mult( inf_zero_one_to_triple(a2,b2,c2), matrix2_inv(inf_zero_one_to_triple(a1,b1,c1) ) ) 

if __name__ == '__main__':
    print sign([0,1,2,3])
    print sign([0,2,1,3])
    print sign([0,1,2,2])
