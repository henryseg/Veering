
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

def times(a,b): ### 2d coordinatewise multiplication
  return [a[0]*b[0], a[1]*b[1]]

def add(a,b):
  return [a[0]*b[1]+a[1]*b[0], a[1]*b[1]]

def neg(b):
  return [-b[0], b[1]]

def invert(a):
  return [a[1],a[0]]

def preferred_rep(a):
  if abs(a[1]) < 0.00001:
    return a
  else:
    return [a[0]/a[1], complex(1.0,0.0)]   ### saul: perhaps should divide by cmath.sqrt(a[0]*a[1])

def developed_position(A1,A2,A3,z): #use Feng's "solving Thurston's equations in a commutative ring"
  # print A1, A2, A3, z
  a1,b1 = A1 ## eg 1,0
  a2,b2 = A2 ## eg 0,1
  a3,b3 = A3 ## eg 1,1
  assert (a1*b2-b1*a2) != complex(0,0), ",".join([str(a1), str(a2), str(b1), str(b2)]) 
  invdet = 1.0/(a1*b2-b1*a2) ## eg 1
  X = [[b2*invdet, -a2*invdet],[-b1*invdet, a1*invdet]] ## eg [[1,0],[0,1]]
  [a3p,b3p] = matrix_mult_vector(X,[a3,b3]) ## eg [1,1]
  Xinv = matrix2_inv(X)
  assert abs(b3p) > 0.00000001 and abs(a3p) > 0.00000001, str(b3p) + ' ' + str(a3p) + ' ' + str([A1,A2,A3])
  out = matrix_mult_vector(Xinv, [-z/b3p,-1/a3p]) ## eg [-z,-1] = [z,1]
  return preferred_rep(out) 

  ## g[p_,q_,r_,s_,t_,u_] is a matrix that takes p,q,r to s,t,u

  ## h[z1_, z2_, z3_] := {{z2 - z3, -z1 (z2 - z3)}, {z2 - z1, -z3 (z2 - z1)}};
  ## g[p_, q_, r_, u_, v_, w_] := Inverse[h[u, v, w]].h[p, q, r];

def inf_zero_one_to_triple(p,q,r):
  #print 'pqr', p, q, r
  p1,p2=p
  q1,q2=q
  r1,r2=r
  M = [[p1,q1],[p2,q2]]
  Minv = matrix2_inv(M)
  #print 'M, matrix2_det(M), Minv',M, matrix2_det(M), Minv
  [mu,lam] = matrix_mult_vector(matrix2_inv([[p1,q1],[p2,q2]]), [r1,r2])
  # print 'mu,lam', mu, lam
  #print [[mu*p1, lam*q1],[mu*p2, lam*q2]]
  return [[mu*p1, lam*q1],[mu*p2, lam*q2]]

def two_triples_to_PSL(a1,b1,c1,a2,b2,c2):
  return matrix_mult( inf_zero_one_to_triple(a2,b2,c2), matrix2_inv(inf_zero_one_to_triple(a1,b1,c1) ) ) 
  
infty = complex(10,10)  ### hack

def convert_to_complex(a):
  if abs(a[1]) < 0.000001 * abs(a[0]):
    return infty
  else:
    return a[0]/a[1]

unknown_vert_to_known_verts_ordering = {0:(3,2,1), 1:(2,3,0), 2:(1,0,3), 3:(0,1,2)}
### if we don't know the position of vert i, this gives the order of the vertices to put in for inf, 0, 1
### Note that the complex angles we have are associated to the edges 01|23
### This is consistent with the orientation given by veering 

def develop_vert_posns(vertex_posns, gluing, face_vertex, tet_shape):
  """Get vert posns for a new tetrahedron given vert posns for an existing tetrahedron and gluing info"""
  next_tet_vertex_posns = []
  current_vertex_nums = [0,1,2,3]
  current_vertex_nums.remove(face_vertex)
  next_tet_vertex_posns = [None, None, None, None]
  for i in current_vertex_nums:
    next_tet_vertex_posns[gluing[i]] = vertex_posns[i]
  #get shape of next tet
  next_tet_new_vertex_num = gluing[face_vertex] #the number of the vertex we don't have the position of
  ordering = unknown_vert_to_known_verts_ordering[next_tet_new_vertex_num]
  new_position = developed_position(next_tet_vertex_posns[ordering[0]], next_tet_vertex_posns[ordering[1]], next_tet_vertex_posns[ordering[2]], tet_shape)
  next_tet_vertex_posns[next_tet_new_vertex_num] = new_position
  return next_tet_vertex_posns



  