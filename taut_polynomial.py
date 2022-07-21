#
# taut_polynomial.py
#

# compute taut polynomial as defined by Sam Taylor et al.  We also
# compute the alexander polynomial here, to ensure that we are using
# the same spanning tree.

import regina

from sage.arith.misc import gcd
from sage.rings.rational_field import QQ
from sage.matrix.constructor import Matrix

from taut import liberal, vert_pair_to_edge_num
from transverse_taut import is_transverse_taut
from fundamental_domain import spanning_dual_tree
from taut_homology import (edge_equation_matrix_taut, group_ring,
                           faces_in_laurent, matrix_laurent_to_poly,
                           normalise_poly, epimorphism_in_laurent)


verbose = 0


# upper and lower diagonals of the tet

def tet_lower_upper_edges(tetrahedron, coorientations):
    tet_coor = coorientations[tetrahedron.index()]
    lower_edge_endpoints = [i for i in range(4) if tet_coor[i] == +1]
    lower_edge_endpoints.sort()
    lower_edge_num = vert_pair_to_edge_num[ tuple(lower_edge_endpoints) ]
    upper_edge_num = 5 - lower_edge_num
    return ( tetrahedron.face(1, lower_edge_num), tetrahedron.face(1, upper_edge_num) )


# Polynomials that come with snappy (in sage) - delete these or move
# them to the correct place.

def alex_is_monic(M):
    p = M.alexander_polynomial()
    return p.is_monic()


def hyper_is_monic(M): # add a way to dial up the precision
    p = M.hyperbolic_torsion()
    lead = p.coefficients(sparse=False)[-1]
    return abs(1 - lead) < 0.000001  # worry about lead = -1


# computing the taut polynomial


@liberal
def edges_to_triangles_matrix(triangulation, angle_structure, cycles, ZH, P, mode = "taut"):
    # In mode alexander, we are computing the transpose of the
    # boundary operator from relative 1-chains to relative 0-chains of
    # the dual 2-complex (relative to its vertex set).  Note that we
    # technically should only be taking a single dual vertex
    # downstairs and lifting that... but then we would have to compute
    # a kernel.  By using all of the dual vertices we are splitting
    # off a free summand of the correct rank.

    # In mode "taut" we are computing the matrix which assigns to a
    # triangle the switch condition on its dual lower track
    coorientations = is_transverse_taut(triangulation, angle_structure, return_type = "tet_vert_coorientations")
    if verbose > 0:
        print(("coorientations", coorientations))
    face_laurents = faces_in_laurent(triangulation, angle_structure, cycles, ZH)
    if verbose > 0:
        print(("face_laurents", face_laurents))

    ET_matrix = [] # now to find the face coefficients relative to each edge
    for tet in triangulation.tetrahedra():
        # get its upper edge - we iterate over upper edges of tetrahedra
        if verbose > 0:
            print(("tet_index", tet.index()))
        edge = tet_lower_upper_edges(tet, coorientations)[1]
        if verbose > 0:
            print(("edge_index", edge.index()))
        embeddings = list(edge.embeddings())

        # find index of tet in the list of embeddings of edge
        for i, embed in enumerate(embeddings):
            tet = embed.tetrahedron()
            if verbose > 0:
                print(("current_tet", tet.index()))
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]

            if (coorientations[tet.index()][trailing_vert_num] == +1 and
                coorientations[tet.index()][leading_vert_num]  == +1):
                bottom_index = i
                break

        # rotate
        embeddings = embeddings[bottom_index:] + embeddings[:bottom_index]

        face_coeffs = [ ZH(0) ] * 2 * triangulation.countTetrahedra()
        sign = -1 # we are going up the left side of the edge
        current_coeff = ZH(1)

        embed = embeddings[0]
        assert tet.index() == embed.tetrahedron().index() # sanity check

        if mode == "taut":
            # because of a sign change, we install the first and last by hand.
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            leading_face = tet.triangle(trailing_vert_num)
            face_coeffs[leading_face.index()] = face_coeffs[leading_face.index()] + current_coeff
            trailing_face = tet.triangle(leading_vert_num)
            face_coeffs[trailing_face.index()] = face_coeffs[trailing_face.index()] + current_coeff
            if verbose > 0:
                print(("face_coeffs", face_coeffs))
            embeddings = embeddings[1:-1] # and remove them

        for embed in embeddings:
            tet = embed.tetrahedron()
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            leading_face = tet.triangle(trailing_vert_num)
            trailing_face = tet.triangle(leading_vert_num)

            if sign == -1:
                use_face = trailing_face
            else: # sign == -1
                use_face = leading_face
            current_coeff = current_coeff * (face_laurents[use_face.index()])**sign
            if verbose > 0:
                print(("current_coeff", current_coeff))

            if (coorientations[tet.index()][trailing_vert_num] == -1 and 
                coorientations[tet.index()][leading_vert_num]  == -1): ## we are the top embed
                sign = 1 # so now are going down the right side of the edge
                current_coeff = current_coeff * (face_laurents[leading_face.index()])**sign
                if verbose > 0:
                    print(("top current_coeff", current_coeff))
            if mode == "taut":
                face_coeffs[leading_face.index()] = face_coeffs[leading_face.index()] - current_coeff
            elif mode == "alexander":
                face_coeffs[leading_face.index()] = face_coeffs[leading_face.index()] + sign*current_coeff
            if verbose > 0:
                print(("face_coeffs", face_coeffs))

        ET_matrix.append(face_coeffs)

    # convert and return
    return matrix_laurent_to_poly(ET_matrix, ZH, P)


@liberal
def edges_to_triangles_matrix_wrapper(tri, angle, cycles):
    ZH = group_ring(tri, angle, cycles, alpha = True)
    P = ZH.polynomial_ring()
    return edges_to_triangles_matrix(tri, angle, cycles, ZH, P, mode = "taut")


@liberal
def taut_polynomial(tri, angle, cycles = [], alpha = True, mode = "taut"):
    # set up
    ZH = group_ring(tri, angle, cycles, alpha = alpha)
    P = ZH.polynomial_ring()

    ET = edges_to_triangles_matrix(tri, angle, cycles, ZH, P, mode = mode)

    # compute via minors
    minors = ET.minors(tri.countTetrahedra())
    return normalise_poly(gcd(minors), ZH, P)


@liberal
def taut_polynomial_via_tree(tri, angle, cycles = [], alpha = True, mode = "taut"):
    """
    If cycles = [] then this is the taut polynomial.  If cycles != []
    then this is the zero-Fitting invariant of the module obtained by
    'extension of scalars' - that is, form the epimorphism obtained by
    killing the given boundary cycles, apply it to the presentation
    matrix, and only then compute the gcd of the (correct minors).       
    """

    # set up 
    ZH = group_ring(tri, angle, cycles, alpha = alpha)
    P = ZH.polynomial_ring()

    ET = edges_to_triangles_matrix(tri, angle, cycles, ZH, P, mode = mode)
    _, non_tree_faces, _ = spanning_dual_tree(tri)

    ET = ET.transpose()
    ET = Matrix([row for i, row in enumerate(ET) if i in non_tree_faces]).transpose()

    # compute via minors
    minors = ET.minors(tri.countTetrahedra())
    return normalise_poly(gcd(minors), ZH, P)


@liberal
def taut_polynomial_image(tri, angle, cycles = [], alpha = True, mode = "taut"):
    """
    If cycles = [] then this is the taut polynomial.  If cycles != []
    then this is the image of the taut polynomial under the
    epimorphism obtained by killing the given boundary cycles.
    """
    taut_poly = taut_polynomial_via_tree(tri, angle, alpha = alpha, mode = mode)
    ZG = group_ring(tri, angle, [], alpha = alpha)
    ZH = group_ring(tri, angle, cycles, alpha = alpha)
    P = ZH.polynomial_ring()

    image_in_laurent = epimorphism_in_laurent(tri, angle, cycles, ZH)
    epi = ZG.Hom(ZH)(image_in_laurent)
    return normalise_poly( epi(taut_poly), ZH, P ) 


@liberal
def taut_polynomial_via_smith(tri, angle, cycles = [], alpha = True, mode = "taut"):
    # set up
    assert tri.homology().rank() == 1 # need the polynomial ring to be a PID
    ZH = group_ring(tri, angle, cycles, alpha = alpha, ring = QQ) # ditto
    P = ZH.polynomial_ring()

    ET = edges_to_triangles_matrix(tri, angle, cycles, ZH, P, mode = mode)

    # compute via smith normal form
    ETs = ET.smith_form()[0]
    a = tri.countEdges()
    ETs_reduced = Matrix([row[:a] for row in ETs])
    return normalise_poly(ETs_reduced.determinant(), ZH, P)


@liberal
def taut_polynomial_via_tree_and_smith(tri, angle, cycles = [], alpha = True, mode = "taut"):
    # set up
    assert tri.homology().rank() == 1 # need the polynomial ring to be a PID
    ZH = group_ring(tri, angle, cycles, alpha = alpha, ring = QQ) # ditto
    P = ZH.polynomial_ring()

    ET = edges_to_triangles_matrix(tri, angle, cycles, ZH, P, mode = mode)
    _, non_tree_faces, _ = spanning_dual_tree(tri)

    ET = ET.transpose()
    ET = Matrix([row for i, row in enumerate(ET) if i in non_tree_faces]).transpose()

    # compute via smith normal form
    ETs = ET.smith_form()[0]
    a = tri.countEdges()
    ETs_reduced = Matrix([row[:a] for row in ETs])
    return normalise_poly(ETs_reduced.determinant(), ZH, P)


# Remarks on taut_polynomial_via_smith and
# taut_polynomial_via_tree_and_smith

# 1. Running on "gLLAQacdefefjkaaqks_200210" appears to involve very
# large intermediate computations - I've never waited long enough for
# the function to actually return.  So there is a serious dichotomy:
# the version using smith normal form is exponentially faster most of
# the time, but is exponentially slower every once in a while.

# 2. Also, smith normal form wants the polynomial ring to be a PID, so
# we are restricted to using QQ[x].  So we need b_1 = 1 or we need to
# specialise the matrix ET before computing.  Both of these argue
# strongly against using smith normal form...


@liberal
def taut_polynomial_via_interpolate(tri, angle, cycles = [], alpha = True):
    # In the one-variable case: (1) get an upper bound on the
    # degree of the answer (2) plug integers into the matrix and
    # compute smith normal form (over ZZ).  (3) Write down the
    # function you get. (4) use difference equations to compute the
    # derivatives, back compute, and win.
    # This might also work in the multivariable case.
    return None
