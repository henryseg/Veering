#
# z_charge.py
#

# Goal - generate z_charges, following [BB] and use them to prove that
# some manifolds have finitely many veering structures. 


from sage.rings.integer_ring import ZZ
from sage.misc.misc import powerset
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector
from sage.modules.free_module import VectorSpace
from sage.numerical.mip import MIPSolverException, MixedIntegerLinearProgram

import regina
import snappy
import flipper

from snappy.snap import t3mlite as t3m

from .taut import is_taut, charge_to_angle, angle_to_charge, lex_smallest_angle_structure, unsorted_vert_pair_to_edge_pair
from .taut_polytope import dot_prod, extract_solution, is_layered
from .veering_tri import is_veering
from .z2_taut import is_trivial_in_cohomology

ZZ2 = ZZ.quotient(ZZ(2))


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
    Given a snappy manifold M, returns the matrix of left-hand-sides
    of angle equations (that is, the tet, edge, and cusp equations).
    """
    num_tet = M.num_tetrahedra()
    G = M.gluing_equations()
    T = Matrix(ZZ, [tet_vector(i, num_tet) for i in range(num_tet)])
    return T.transpose().augment(G.transpose()).transpose() # sigh

def sol_and_kernel(M):
    """
    Given a snappy manifold M, returns an integer solution to the
    angle equations, as well as an integer basis for the kernel.
    """
    A = angle_equations(M)
    b = solution_vector(M)
    D, U, V = A.smith_form() 
    # UAV = D so Uinv D Vinv = A
    # want to solve Ax = b
    # that is Uinv D Vinv x = b
    # that is D Vinv x = Ub
    c = U * b
    
    min_dim, max_dim = A.dimensions()
    assert min_dim <= max_dim
    # D is diagonal, so Dinv means "divide".
    # by [L] + [BB] there is always a solution:
    assert all(D[i][i].divides(c[i]) for i in range(min_dim))
    # so divide, and set 0 / 0 equal to 0
    c = [c[i] / D[i][i] if D[i][i] != 0 else 0 for i in range(min_dim)]

    # need to have the correct dimension
    padding = vector(ZZ, [0]*(max_dim - min_dim))
    c.extend(padding)
    c = vector(ZZ, c)
    # thus Vinv x = Dinv U b
    # thus x = V Dinv U b
    x = V * c
    assert A*x == b
    return x, A.right_kernel().basis()

def leading_trailing_deformations(M):
    tri = regina.Triangulation3(M)
    num_tet = tri.countTetrahedra()
    out = []
    for e in tri.edges():
        defm = [0] * (3 * num_tet)
        for i in range(e.degree()):
            emb = e.embedding(i)
            tet_num = emb.simplex().index()
            v0, v1, v2 = emb.vertices()[0], emb.vertices()[1], emb.vertices()[2]
            a = unsorted_vert_pair_to_edge_pair[(v0, v2)]
            b = unsorted_vert_pair_to_edge_pair[(v1, v2)]
            defm[3*tet_num + a] += 1
            defm[3*tet_num + b] -= 1
        out.append(vector(defm))
    return out

def has_pi_triple(u):
    n = len(u)
    assert n % 3 == 0
    one = vector([1, 1, 1])
    for i in range(int(n/3)):
        if u[3*i:3*i + 3] == one: 
            return True
    return False

def reduced_charges(M):
    """
    Given a snappy manifold M, we find all reduced charges so that:
    (1) no loop in the triangulation passes an odd number of pi's.
    (2) no tetrahedron has three pi's and 
    """
    out = sol_and_kernel(M)
    x, A = out
    nt = M.num_tetrahedra()
    dim = 3*nt
    V = VectorSpace(ZZ2, dim)
    AA = V.subspace(A) # the reduced kernel
    xx = V(x) # the reduced solution
    
    charges = [xx + sum(B) for B in powerset(AA.basis())]    
    charges = [c for c in charges if not has_pi_triple(c)] # reject if there are three pi's in any tet.
    return charges

def reduced_angles(M):
    """
    Given a snappy manifold M, compute the reduced charges, convert to
    angle structures, remove repeated structures (using symmetries of
    the triangulation), remove non-trivial structures (in cohomology),
    and return what remains.
    """
    charges = reduced_charges(M)    
    tri = regina.Triangulation3(M)
    angles = [charge_to_angle(c) for c in charges] 

    # remove symmetries
    lex_angles = [lex_smallest_angle_structure(tri, angle) for angle in angles]
    angles = []
    for angle in lex_angles:  
        if angle not in angles:
            angles.append(angle) 

    # remove the angles that flip a triangle over 
    angles = [angle for angle in angles if is_trivial_in_cohomology(tri, angle)] 
    return angles

def can_deal_with_reduced_angle(tri, angle):
    """
    Returns True or False, as our techniques can recognise the given
    reduced angle structure.
    """
    if not is_taut(tri, angle):
        return False  
    if is_veering(tri, angle):
        return True
    if is_layered(tri, angle):  # has_internal_singularities is not needed here.
        return True
    return False

def can_deal_with_reduced_angles(M, report = False):
    """
    Returns True if we can deal with all of the reduced angles. 
    """
    angles = reduced_angles(M)
    tri = regina.Triangulation3(M)
    if report:
        nv = num_veering_structs(M, angles = angles, use_flipper = False)
        return all(can_deal_with_reduced_angle(tri, angle) for angle in angles), len(angles), nv
    else: 
        return all(can_deal_with_reduced_angle(tri, angle) for angle in angles) 

def num_veering_structs(M, angles = None, use_flipper = True):
    """
    Tries to count them (in a very naive way). 
    Tries to eliminate overcounting due to symmetries (in reduced_angles)
    but will fail if one of the symmetries is hidden by retriangulation. 
    If use_flipper = False then we get a (true) lower bound, since then it only
    counts veering structures on the given triangulation.
    """
    if angles == None:
        angles = reduced_angles(M)
    tri = regina.Triangulation3(M)
    for angle in angles:
        if not is_taut(tri, angle):
            return None
    for angle in angles:
        if not (is_veering(tri, angle) or is_layered(tri, angle)):
            return None
    total = 0
    for angle in angles:
        if is_veering(tri, angle):
            total = total + 1
        elif use_flipper:
            assert is_layered(tri, angle)
            print(M.name(), angle, "needs flipper")
            try: 
                if not has_internal_singularities(tri, angle):
                    total = total + 1
            except:
                print("flipper failed")
    return total

def has_internal_singularities(tri, angle):
    """
    Given a regina manifold tri and an angle structure (assumed to be
    layered), convert tri to a t3m triangulation, convert angle to the
    correct flipper format, use them both to get a flipper
    TautStructure, find the monodromy, and then check the stratum.  If
    there are internal singularities, return True.
    """
    # See
    # https://github.com/MarkCBell/flipper/blob/master/flipper/kernel/taut.py
    # for the relevant code in flipper
    T = t3m.Mcomplex(snappy.Manifold(tri))
    angle_vector = angle_to_charge(angle, flipper_format = True)
    taut_struct = flipper.kernel.taut.TautStructure(T, angle_vector)
    strat = taut_struct.monodromy().stratum()
    return any(punc.filled for punc in strat.keys())
