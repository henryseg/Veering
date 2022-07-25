#
# Finding a path of mobius transformations
#

from sage.symbolic.all import I
from sage.functions.other import sqrt
from sage.functions.log import log
from sage.functions.log import exp

from basic_math import matrix

# M is the original matrix.

# [1/2*I*sqrt(3) + 1/2      -sqrt(3) - 2*I]
# [1/2*sqrt(3) + 1/2*I      -I*sqrt(3) - 2]

# M.eigenvalues() gives

e = -1/4*I*sqrt(3) - 1/4*sqrt(6*I*sqrt(3) - 10) - 3/4                                                                
f = -1/4*I*sqrt(3) + 1/4*sqrt(6*I*sqrt(3) - 10) - 3/4

e = e.n()
f = f.n()

le = log(e)
lf = log(f)

# Diagonalise using
# var("p, q, r, s")                                                                                                    
# eqns = [p*s - q*r == 1,
#         e*p*s - f*q*r == 1/2*I*sqrt(3) + 1/2,
#         (e - f)*q*s == -sqrt(3) - 2*I,
#         (-e + f)*p*r == 1/2*sqrt(3) + 1/2*I,
#         -e*q*r + f*p*s == -I*sqrt(3) - 2,
#         r == 1]
# solve(eqns, p, q, r, s) 

p = -1/52*(sqrt(3) + 7*I)*sqrt(6*I*sqrt(3) - 10)                                                                     
q = 1/104*((15*I*sqrt(3) - 1)*sqrt(6*I*sqrt(3) - 10) - 52)                                                           
r = 1                                                                                                                
s = 1/8*(sqrt(3) - I)*sqrt(6*I*sqrt(3) - 10) - 1/8*(8*sqrt(3) + 4*I)

C = matrix((p, q, r, s))
C = C.n()

def many_matrices(n): 
    matrices = [] 
    for k in range(n): 
        D = matrix((exp(k*le/n), 0, 0, exp(k*lf/n)))
        D = D.n() 
        matrices.append(C.inverse()*D*C) 
    return matrices 
