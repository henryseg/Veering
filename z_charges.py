#
# z_charges
#

# Goal - generate z_charges, following [BB]

import regina

from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector
from sage.numerical.mip import MIPSolverException, MixedIntegerLinearProgram
from sage.misc.misc import powerset
from sage.rings.integer_ring import ZZ

from taut import charges_to_taut_struct, lex_smallest_angle_structure, is_taut
from taut_polytope import dot_prod, extract_solution, is_layered
from z2_taut import is_trivial_in_cohomology
from veering import is_veering

def tet_vector(i, num_tet):
    """
    Gives the tet equation for the i^th tetrahedron.
    """
    out = []
    for j in range(num_tet): 
        if i == j: out.extend([1]*3)
        else: out.extend([0]*3)
    return out

def solution_vector(M):
    """
    Given a snappy manifold, returns the desired sums of the tet,
    edge, and holonomy equations.
    """
    num_tet = M.num_tetrahedra()
    num_cusps = M.num_cusps()
    sol = []
    sol.extend([1]*num_tet) # pi
    sol.extend([2]*num_tet) # 2 pi
    sol.extend([0]*2*num_cusps) # zero
    return vector(ZZ, sol)

def angle_equations(M):
    """
    Given a snappy manifold M, returns the matrix of angle equations
    (that is, the tet, edge, and cusp equations).
    """
    num_tet = M.num_tetrahedra()
    G = M.gluing_equations()
    T = Matrix(ZZ, [tet_vector(i, num_tet) for i in range(num_tet)])
    return T.transpose().augment(G.transpose()).transpose() # sigh

def sol_and_kernel(M):
    """
    Given a snappy manifold M, returns a solution to the angle
    equations, as well as a basis for the kernel.
    """
    A = angle_equations(M)
    b = solution_vector(M)
    return A.solve_right(b), A.right_kernel().basis()

def last_row_with_non_zero_ith_entry(A, i):
    out = None
    for row in A:
        if row[i] != 0:
            out = row
    return out

def intify_sol(t, A):
    """
    t is rational, A is a list of integer vectors of same length as t. 
    Use multiples of rows of A to intify t
    """
    for i in range(len(t)):
        if not t[i].is_integer():
            row = last_row_with_non_zero_ith_entry(A, i)
            t = t - (t[i]/row[i]) * row  ### make ith entry zero, which is an integer!
    for x in t:
        if not x.is_integer():
            return None
    return t

def better_int_sol_and_kernel(M):
    t, A = sol_and_kernel(M)
    t = intify_sol(t, A)
    if t == None:
        return None
    else:
        return t, A

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
    Given a snappy manifold M, we find all reduced charges so that:
    (1) no tetrahedron has three pi's and 
    (2) no loop in the triangulation passes an odd number of pi's.
    We return these after converting them to "angle structures".  
    The quotes are there because the edge equations are only satisfied 
    modulo two. 
    """
    out = final_int_sol_and_kernel(M)
    # out = better_int_sol_and_kernel(M)
    if out == None:
        return None
    t, A = out
    nt = M.num_tetrahedra()
    charges = [reduce(t + sum(B)) for B in powerset(A)] 
    charges = [c for c in charges if sum(c) == nt] # reject if there are three pi's in any tet.

    tri = regina.Triangulation3(M)
    angles = [charges_to_taut_struct(c) for c in charges] 

    lex_angles = [lex_smallest_angle_structure(tri, angle) for angle in angles]
    angles = []
    for angle in lex_angles:  ## remove symmetries
        if angle not in angles:
            angles.append(angle) 

    return [angle for angle in angles if is_trivial_in_cohomology(tri, angle)]

def can_deal_with_reduced_charges(M):
    """
    Returns True or False, answering the question of whether our techniques 
    recognise each of the reduced charges we find.
    """
    angles = reduced_charges(M)
    if angles == None:
        return False
    tri = regina.Triangulation3(M)
    for angle in angles:
        if not is_taut(tri, angle):
            return False
        elif not (is_veering(tri, angle) or is_layered(tri, angle)):
            return False
    return True

def final_int_sol_and_kernel(M):
    A = angle_equations(M)
    b = solution_vector(M)
    D, U, V = A.smith_form()
    min_dim, max_dim = A.dimensions()
    assert min_dim <= max_dim
    c = U*b
    if not all(D[i][i].divides(c[i]) for i in range(min_dim)):
        return None
    c = [c[i] / D[i][i] if D[i][i] != 0 else c[i] for i in range(min_dim)]
    padding = vector(ZZ, [0]*(max_dim - min_dim))
    c.extend(padding)
    c = vector(ZZ, c)
    x = V * c
    assert A*x == b
    return x, A.right_kernel().basis()
