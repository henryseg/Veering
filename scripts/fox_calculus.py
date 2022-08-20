
#
#  fox_calculus.py
#

from sage.arith.misc import gcd
from sage.groups.free_group import FreeGroup

from veering.taut import liberal, vert_pair_to_edge_num
from veering.veering_tri import veering_triangulation
from veering.transverse_taut import edge_side_face_collections, top_bottom_embeddings_of_faces, get_tet_top_vert_nums
from veering.taut_homology import group_ring, faces_in_laurent, matrix_laurent_to_poly, normalise_poly
from veering.fundamental_domain import spanning_dual_tree, non_tree_edge_loops_oriented


def fundamental_group(tri, angle, simplified = True):
    """Returns the finite presentation of the fundamental group arising
    from the canonical spanning tree.
    
    Input 

    * "tri" -- a regina triangulation.

    * "angle" -- an veering angle structure

    * "simplified" -- returns the Tietze simpification [Havas, 1994]
      of the fundamental group.

    Examples:

    sage: sig = "cPcbbbdxm_10"  
    sage: tri, angle = isosig_to_tri_angle(sig)
    sage: fundamental_group(sig, False)
    Finitely presented group < x0, x1, x2, x3 | x1^-1*x0^-1*x1^-1*x3*x2*x3, x2^-1*x1^-1*x2^-1*x0*x3*x0, x0 >

    This is the unsimplified presentation of the fundamental group of
    the figure-eight sibling.

    sage: sig = "cPcbbbiht_12"  # figure-eight
    sage: tri, angle = isosig_to_tri_angle(sig)
    sage: fundamental_group(sig)
    Finitely presented group < x1, x2 | x1^-1*x2^-2*x1^-1*x2*x1*x2^-1*x1*x2 >
    
    This is the simplified presentation of the fundamental group of
    the figure-eight knot complement.
    """
    tree = spanning_dual_tree(tri)
    tree_edges = tree[0]
    F = FreeGroup(tri.countTriangles())
    # We add all of the faces as generators.  This is useful, for
    # example, when defining homomorphisms from \pi_1(M).
    gens = F.gens()

    esfc = edge_side_face_collections(tri, angle)
    relators = []
    for edge in esfc:  # relations for the edges
        rel = F.one()
        left, right = edge
        left.reverse()
        for term in left:
            rel = rel * gens[term[0]].inverse()
        for term in right:
            rel = rel * gens[term[0]]
        relators.append(rel)

    gens = F.gens()
    for f in tree_edges:  # killing the tree faces
        relators.append(gens[f])
    G = F/relators
    if simplified:
        G = G.simplified()
    return G


def is_AB_turn(vt, top_bottom_embeddings, face0, face1, face0_dir, face1_dir):
    """
    For the "turn" (face0, face1, face0_dir, face1_dir) in the veering
    triangulation vt, we decide if it is an "anti-branching" (AB) turn
    as defined at the top of page 16 of arxiv:2008.04836.  In more
    detail: we travel along face0 in direction face0_dir (+1 if with
    the coorientation) into a tet t. We then leave through face1 in
    direction face1_dir.  We return True if this turn is an AB turn:
    the triangles are adjacent along an equatorial edge of t of the
    same colour as the top diagonal of the edge.
    """
    top, bottom = top_bottom_embeddings
    if face0_dir == 1:
        embed0 = bottom[face0]
    else:
        embed0 = top[face0]
    if face1_dir == -1:
        embed1 = bottom[face1]
    else:
        embed1 = top[face1]
    t0 = embed0.tetrahedron()
    t1 = embed1.tetrahedron()
    # print(t0.index(), t1.index())
    assert(t0 == t1)
    t = t0
    if face0_dir != face1_dir:
        return False
    f0 = embed0.face()
    f1 = embed1.face()
    equatorial_nums = [0,1,2,3]
    equatorial_nums.remove(f0)
    equatorial_nums.remove(f1)
    equatorial_colour = vt.get_edge_between_verts_colour(t.index(), equatorial_nums)
    top_vert_nums = get_tet_top_vert_nums(vt.coorientations, t.index())
    top_colour = vt.get_edge_between_verts_colour(t.index(), top_vert_nums)
    return top_colour == equatorial_colour


def loop_twistednesses(tri, angle):
    """
    Returns the (list of the) images of the face generators under the
    edge-orientation homomorphism.  For each face of the triangulation
    we take its canonical loop (using the canonical spanning tree).
    We then compute the image (either +1 or -1) by counting the parity
    of number of AB turns along the loop.  We return these as a list.
    (Note that the homomorphism is trivial on tree faces.)  See
    Proposition 5.7 of arxiv:2101.12162.
    """
    vt = veering_triangulation(tri, angle)
    top_bottom_embeddings = top_bottom_embeddings_of_faces(tri, angle)
    twistednesses_dict = {}
    oriented_loops, all_signs = non_tree_edge_loops_oriented(tri, angle)
    for i in range(len(oriented_loops)):
        loop = oriented_loops[i]
        signs = all_signs[i]
        count = 0
        for j in range(len(loop)):
            f0, f1 = loop[j], loop[(j+1)%len(loop)]
            f0d, f1d = signs[j], signs[(j+1)%len(loop)]
            if is_AB_turn(vt, top_bottom_embeddings, f0, f1, f0d, f1d):
                count += 1
        twistednesses_dict[loop[0]] = (-1)**(count % 2)  # first in loop is the non-tree face
    for i in range(tri.countTriangles()):
        if i not in twistednesses_dict:
            twistednesses_dict[i] = 1
    
    return [twistednesses_dict[i] for i in range(tri.countTriangles())]


@liberal
def taut_polynomial_via_fox_calculus(tri, angle, simplified = True):
    """
    Computes the taut polynomial.  Under the hood we are relying on
    (1) sage's Fox calculus to find the Alexander matrix of the
    fundamental group and (2) the result that the taut polynomial is
    the Alexander polynomial twisted by the edge-orientation
    homomorphism.  See Proposition 5.7 of arxiv:2101.12162.
    
    Input 

    * "tri" -- a regina triangulation.

    * "angle" -- an veering angle structure

    * "simplified" -- uses the Tietze simpification [Havas, 1994]
      of the fundamental group. 

    Examples:

    sage: sig = "cPcbbbdxm_10"
    sage: taut_polynomial_via_fox_calculus(sig)  
    a^2 - 3*a + 1

    The above works because of the "liberal" decorator.

    sage: sig = "eLMkbcddddedde_2100"
    sage: tri, angle = isosig_to_tri_angle(sig)
    sage: taut_polynomial_via_fox_calculus(tri, angle)
    a^2*b - a^2 - a*b - b^2 + b

    """
    lt = loop_twistednesses(tri, angle)  # images under the edge-orientation homo
    ZH = group_ring(tri, angle, [], alpha = True)
    # Note that alpha = True is "ok" here because in the census we
    # have b_1 \leq 4.
    P = ZH.polynomial_ring() 
    fl = faces_in_laurent(tri, angle, [], ZH)  # images in ZZ[H_1/torsion]
    flt = [a * x for a, x in zip(fl, lt)]  # twisting

    G = fundamental_group(tri, angle)
    if simplified:
        G = G.simplified()
        indices = [int(str(x)[1:]) for x in G.gens()]  # the hackest of hacks
        flt = [flt[i] for i in indices]
    M = G.alexander_matrix(flt)
    N = matrix_laurent_to_poly(M, ZH, P)
    n = len(G.gens()) - 1
    poly = gcd(N.minors(n))
    poly = normalise_poly(poly, ZH, P)
    return poly
