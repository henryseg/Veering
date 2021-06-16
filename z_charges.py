#
# z_charges
#

# Goal - generate z_charges, following [BB]

import regina

from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector
from sage.numerical.mip import MIPSolverException, MixedIntegerLinearProgram
from sage.misc.misc import powerset

from taut import charges_to_taut_struct
from taut_polytope import dot_prod, extract_solution
from z2_taut import is_trivial_in_cohomology

def tet_vector(i, num_tet):
    """
    Gives the tet equation for the i^th tetrahedron.
    """
    out = []
    for j in range(num_tet): 
        if i == j: out.extend([1]*3)
        else: out.extend([0]*3)
    return out

def solution_vector(num_tet, num_cusps):
    """
    The desired sums of the tet, edge, and holonomy equations.
    """
    out = []
    out.extend([1]*num_tet) # pi
    out.extend([2]*num_tet) # 2 pi
    out.extend([0]*2*num_cusps) # zero
    return vector(out)

def sol_and_kernel(M):
    """
    Given a snappy manifold M, returns a solution to the tet, edge,
    and holonomy equations, as well as a basis for the kernel.
    """
    G = M.gluing_equations()
    T = Matrix(tet_vector(i, M.num_tetrahedra()) for i in range(M.num_tetrahedra()))
    E = T.transpose().augment(G.transpose()).transpose()
    S = solution_vector(M.num_tetrahedra(), M.num_cusps())
    return E.solve_right(S), E.right_kernel().basis()

def real_sol_and_kernel(M):
    nt = M.num_tetrahedra()
    nc = M.num_cusps()
    G = M.gluing_equations() # edge and cusp coeffs
    T = Matrix(tet_vector(i, nt) for i in range(nt)) # tet coeffs
    E = T.transpose().augment(G.transpose()).transpose() 
    S = solution_vector(nt, nc) # the desired sums

    q = MixedIntegerLinearProgram( maximization = False, solver = "GLPK" )
    w = q.new_variable(real = True, nonnegative = True)

    for i, v in enumerate(E.rows()):
        q.add_constraint( dot_prod(v, w) == S[i] )
    s = sum( w[i] for i in range(len(E.columns())) )
    q.add_constraint( s == nt )
    return q.polyhedron(), E.right_kernel().basis()

def int_sol_and_kernel(M):
    nt = M.num_tetrahedra()
    nc = M.num_cusps()
    G = M.gluing_equations() # edge and cusp coeffs
    T = Matrix(tet_vector(i, nt) for i in range(nt)) # tet coeffs
    E = T.transpose().augment(G.transpose()).transpose() 
    S = solution_vector(nt, nc) # the desired sums

    q = MixedIntegerLinearProgram( maximization = False, solver = "GLPK" )
    w = q.new_variable(integer = True, nonnegative = False)

    for i, v in enumerate(E.rows()):
        q.add_constraint( dot_prod(v, w) == S[i] )
    s = sum( w[i] for i in range(len(E.columns())) )
    q.add_constraint( s == nt )
    q.set_objective(None)
    q.solve()
    sol = extract_solution(q, w)
    assert all(a.is_integer() for a in sol)
    sol = vector(int(a) for a in sol)
    return sol, E.right_kernel().basis()

def reduce(u):
    return vector(a % 2 for a in u)

def reduced_charges(M):
    """
    Given a snappy manifold M, find all reduced charges so that:
    (1) no tetrahedron has three pi's and 
    (2) no loop in the triangulation passes an odd number of pi's.
    """
    t, A = int_sol_and_kernel(M)
    nt = M.num_tetrahedra()
    out = [reduce(t + sum(B)) for B in powerset(A)]
    out = [v for v in out if sum(v) == nt] # reject if there are three pi's in any tet.

    tri = regina.Triangulation3.fromSnapPea(M._to_string())
    out = [charges_to_taut_struct(v) for v in out]
    out_new = [angle for angle in out if is_trivial_in_cohomology(tri, angle)]
    if len(out_new) < len(out): print("yay")
    return out_new
