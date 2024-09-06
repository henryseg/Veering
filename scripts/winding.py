#
# winding.py
#

# Goal - generate windings (called \ZZ--charges by [BB]) and use them
# to prove that some manifolds have finitely many veering structures.


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

from veering.taut import pi_edgepair_dict, is_taut, lex_smallest_angle_structure, unsorted_vert_pair_to_edge_pair, isosig_to_tri_angle
from veering.taut_polytope import dot_prod, extract_solution, is_layered, is_torus_bundle
from veering.veering_tri import is_veering
from veering.z2_taut import is_trivial_in_z2_cohomology

ZZ2 = ZZ.quotient(ZZ(2))


def winding_to_preangle(w):
    """
    Given a list of 3*n integers with each triple of the form (1,0,0),
    (0,1,0), or (0,0,1), convert it to our angle structure format.
    Note that we may obtain a "preangle" structure as we may not have
    the edge condition.
    """
    assert len(w) % 3 == 0
    n = int(round(len(w)/3))
    out = []
    for i in range(n):
        tet = w[ 3*i : 3*i+3 ]
        out.append( pi_edgepair_dict[tuple(tet)] )
    return out


def preangle_to_winding(angle, flipper_format = False):
    """
    Given a list of n integers in [0,2], convert to winding format.
    """
    if flipper_format:
        # The veering code uses "vertex with 0 (minus one)". On the
        # other hand, flipper and t3m use "vertex with 3".
        # See line 25 of
        # https://github.com/MarkCBell/flipper/blob/master/flipper/kernel/taut.py
        angle = [2 - a for a in angle]
    out = [0] * (3*len(angle))
    for i, a in enumerate(angle):
        out[3*i + a] = 1
    if flipper_format:
        # flipper adds a variable to homogenise, so we do the same.
        out.append(-1)
    return out


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


def tet_edge_cusp_equations(M):
    """
    Given a snappy manifold M, returns the matrix of left-hand-sides
    of angle equations (that is, the tet, edge, and cusp equations).
    """
    num_tet = M.num_tetrahedra()
    G = M.gluing_equations()
    T = Matrix(ZZ, [tet_vector(i, num_tet) for i in range(num_tet)])
    return T.transpose().augment(G.transpose()).transpose() # sigh


def is_trivial_in_bdy_cohomology(M, w):
    """
    Given a snappy manifold M and a winding w, decides if w induces
    the trivial cohomology class on the boundary.  (In fact we also
    check that w is a winding over ZZ.)
    """
    A = tet_edge_cusp_equations(M)
    b = solution_vector(M)
    return A*x == b


def windings_vanishing_on_bdy(M):
    """
    Given a snappy manifold M, returns an integer solution to the
    tet_edge_cusp_equations, as well as an integer basis for the
    kernel.  The resulting lattice is the set of windings with trivial
    boundary cohomology.
    """
    A = tet_edge_cusp_equations(M)
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


def is_flat(u):
    """
    Given a (possibly reduced) winding u, decide if it is flat.
    Note that this is not an affine condition.
    """
    n = len(u)
    assert n % 3 == 0
    v = [a % 2 for a in u]
    for i in range(int(n/3)):
        a, b, c = v[3*i:3*i + 3]
        if a == b == c:  # at least one is odd...
            return False
    return True


def flat_reduced_windings(M):
    """
    Given a snappy manifold M, we reduce (mod 2) the windings given by
    windings_vanishing_on_bdy and eliminate the non-flat ones.
    """
    x, A = windings_vanishing_on_bdy(M)
    dim = 3*M.num_tetrahedra()
    V = VectorSpace(ZZ2, dim)
    AA = V.subspace(A)  # the reduced kernel
    xx = V(x)  # the reduced solution

    diffs = [sum(B) for B in powerset(AA.basis())]
    windings = [xx + v for v in diffs]
    windings = [c for c in windings if is_flat(c)]
    return windings


def preangles_from_frws(M):
    """
    Given a snappy manifold M, compute the flat reduced windings,
    convert to pre-angle structures (as the edge equations may be
    violated), remove repeated structures (using symmetries of the
    triangulation), remove non-trivial structures (in ZZ2 cohomology),
    and return what remains.
    """
    windings = flat_reduced_windings(M)
    tri = regina.Triangulation3(M)
    angles = [winding_to_preangle(w) for w in windings] 

    # remove symmetries
    lex_angles = [lex_smallest_angle_structure(tri, angle) for angle in angles]
    angles = []
    for angle in lex_angles:  
        if angle not in angles:
            angles.append(angle) 

    # remove the angles that flip a triangle over
    angles = [angle for angle in angles if is_trivial_in_z2_cohomology(tri, angle)] 
    return angles


def can_deal_with_reduced_angle(tri, angle):
    """
    Returns True or False, as our techniques can recognise the given
    reduced angle structure.
    """
    if not is_taut(tri, angle):
        return (False, "not taut")
    if is_torus_bundle(tri, angle):
        return (True, "torus bundle")
    if is_veering(tri, angle) and is_layered(tri, angle):
        return (True, "veering+layered")
    if is_veering(tri, angle):
        return (True, "veering")
    if is_layered(tri, angle):  # has_internal_singularities is not needed here.
        return (True, "layered")
    return (False, None)


def can_deal_with_reduced_angles(M):
    """
    Returns True if we can deal with all of the reduced angles. 
    """
    angles = preangles_from_frws(M)
    tri = regina.Triangulation3(M)
    results = [can_deal_with_reduced_angle(tri, angle) for angle in angles]
    if any(r[1] == "torus bundle" for r in results):
        return (True, "torus bundle")
    return (all(r[0] for r in results), None)


def get_some_sigs(M, tries = 20):
    sigs = set()
    for i in range(tries):
        sigs.add(M.triangulation_isosig())
        M.randomize()
    return sigs


def check_some_sigs(sigs):
    """
    Given a list of sigs for a particular manifold, check if, for any
    one of them, we can deal with all of its angle structures.
    """
    for sig in sigs:
        M = snappy.Manifold(sig)
        b, r = can_deal_with_reduced_angles(M)
        if b:
            return (True, sig, r)
    return (False, None, None)


def check_snappy_manifold(M):
    ids = M.identify()
    s = M.triangulation_isosig()
    b, r = can_deal_with_reduced_angles(M)
    if b:
        print("can deal with", ids, s, r)
    else:
        sigs = get_some_sigs(M)
        b, s, r = check_some_sigs(sigs)
        if b:
            print("can deal with", ids, s, r)

            
def check_veering_sig(veering_sig):
    tri, angle = isosig_to_tri_angle(veering_sig)
    M = snappy.Manifold(tri)
    check_snappy_manifold(M)


def num_veering_structs(M, angles = None, use_flipper = True):
    """
    Tries to count them (in a very naive way). 
    Tries to eliminate overcounting due to symmetries (in preangles_from_frws)
    but will fail if one of the symmetries is hidden by retriangulation. 
    If use_flipper = False then we get a (true) lower bound, since then it only
    counts veering structures on the given triangulation.
    """
    if angles == None:
        angles = preangles_from_frws(M)
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
