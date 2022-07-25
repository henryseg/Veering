#
# Finding a path of mobius transformations
#

from sage.symbolic.all import I
from sage.symbolic.constants import pi
from sage.functions.other import sqrt
from sage.functions.log import log
from sage.functions.log import exp

from basic_math import matrix

# M is the original matrix.

# [1/2*I*sqrt(3) + 1/2      -sqrt(3) - 2*I]
# [1/2*sqrt(3) + 1/2*I      -I*sqrt(3) - 2]

# M.eigenvalues() gives

e = -1/4*I*sqrt(3) - 1/4*sqrt(6*I*sqrt(3) - 10) - 3/4
# f = -1/4*I*sqrt(3) + 1/4*sqrt(6*I*sqrt(3) - 10) - 3/4

e = e.n()
# f = f.n()

# Due to the branch of the logarithm, and because e^2 = e/f is the
# stretch factor, we need to fix the argument.
le = complex(log(e) + pi*I)
lf = -le

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

p = complex(p.n())
q = complex(q.n())
r = complex(r)
s = complex(s.n())

C = matrix((p, q, r, s))

def many_matrices(n): 
    matrices = [] 
    for k in range(n+1): 
        D = matrix((exp(k*le/n), 0, 0, exp(k*lf/n)))
        matrices.append(C.inverse()*D*C) 
    return matrices 

# More generally, this should be "path of matrices" and take as
# arguments the matrix M and the integer n.  Then we do all of the
# work (eigenvalues, diagonalisation, etc, inside of the function.
# I'll fix this if we ever need it.
