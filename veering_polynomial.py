#
# veering_polynomial.py
#

# Goal - compute the "big" and "small" veering polynomials as defined
# by Sam Taylor et al.

import regina

from sage.arith.misc import gcd
from sage.rings.rational_field import QQ
from sage.matrix.constructor import Matrix

from taut import liberal
from transverse_taut import is_transverse_taut
from taut_homology import (build_spanning_dual_tree, edge_equation_matrix_taut,
                           group_ring, faces_in_laurent, matrix_laurent_to_poly,
                           normalise_poly)
from veering import is_veering

verbose = 0

# Polynomials that come with snappy (in sage) - delete these or move
# them to the correct place.


def alex_is_monic(M):
    p = M.alexander_polynomial()
    return p.is_monic()


def hyper_is_monic(M): # add a way to dial up the precision
    p = M.hyperbolic_torsion()
    lead = p.coefficients(sparse=False)[-1]
    return abs(1 - lead) < 0.000001  # worry about lead = -1


# veering structure


edge_vert_index_map = {(0, 1):0, (0, 2):1, (0, 3):2, (1, 2): 3, (1, 3):4, (2, 3):5 }


def tet_lower_upper_edges(tetrahedron, coorientations):
    tet_coor = coorientations[tetrahedron.index()]
    lower_edge_endpoints = [i for i in range(4) if tet_coor[i] == +1]
    lower_edge_endpoints.sort()
    lower_edge_num = edge_vert_index_map[ tuple(lower_edge_endpoints) ]
    upper_edge_num = 5 - lower_edge_num
    return ( tetrahedron.face(1, lower_edge_num), tetrahedron.face(1, upper_edge_num) )


### computing the big veering polynomial


def has_red_lower_edge(tetrahedron, coorientations, edge_colours):
    lower_edge = tet_lower_upper_edges(tetrahedron, coorientations)[0]
    return edge_colours[lower_edge.index()] == 'R'


@liberal
def edges_to_tetrahedra_matrix(triangulation, angle_structure, ZH, P):
    coorientations = is_transverse_taut(triangulation, angle_structure, return_type = 'tet_vert_coorientations')
    if verbose > 0: print 'coorientations', coorientations
    edge_colours = is_veering(triangulation, angle_structure, return_type = 'veering_colours')
    if verbose > 0: print 'edge_colours', edge_colours
    red_tetrahedra = []
    blue_tetrahedra = []
    for tet in triangulation.tetrahedra():
        if has_red_lower_edge(tet, coorientations, edge_colours):
            red_tetrahedra.append(tet)
        else:
            blue_tetrahedra.append(tet)
    if verbose > 0: print 'how many reds and blues', len(red_tetrahedra), len(blue_tetrahedra)

    face_laurents = faces_in_laurent(triangulation, angle_structure, ZH)
    if verbose > 0: print 'face_laurents', face_laurents

    ET_matrix = [] # now to find the tet coefficients relative to each edge
    for tet in triangulation.tetrahedra():
        if verbose > 0: print 'tet_index', tet.index()
        edge = tet_lower_upper_edges(tet, coorientations)[1]
        if verbose > 0: print 'edge_index', edge.index()
        edge_colour = edge_colours[edge.index()]
        if verbose > 0: print 'edge_colour', edge_colour
        embeddings = edge.embeddings()
        tet_coeffs = [ZH(0)] * triangulation.countTetrahedra()
        tet_coeffs[tet.index()] = 1  # bottom tet around the edge gets a 1
        if verbose > 0: print 'initial tet_coeffs', tet_coeffs
        current_coeff = ZH(1)
        # find index of bottom embedding in the list of embedding
        for i, embed in enumerate(embeddings):
            tet = embed.tetrahedron()
            if verbose > 0: print 'current_tet', tet.index()
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            if coorientations[tet.index()][trailing_vert_num] == +1 and coorientations[tet.index()][leading_vert_num] == +1:
                bottom_index = i
                break
        embeddings = embeddings[bottom_index:] + embeddings[:bottom_index]
        sign = 1 # we are going up the left side of the edge
        for embed in embeddings[1:]:  # skipping the first
            tet = embed.tetrahedron()
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            current_coeff = current_coeff * face_laurents[tet.face(2,leading_vert_num).index()]**sign
            if coorientations[tet.index()][trailing_vert_num] == -1 and coorientations[tet.index()][leading_vert_num] == -1:
                # we are the top embed so:
                tet_coeffs[tet.index()] = tet_coeffs[tet.index()] - current_coeff
                sign = -1 # now we go down the right side
            elif edge_colour == 'L' and tet in red_tetrahedra or edge_colour == 'R' and tet in blue_tetrahedra:
                tet_coeffs[tet.index()] = tet_coeffs[tet.index()] - current_coeff
            if verbose > 0: print 'current tet_coeffs', tet_coeffs

        ET_matrix.append(tet_coeffs)

    # convert and return
    return matrix_laurent_to_poly(ET_matrix, ZH, P)


@liberal
def big_polynomial(tri, angle, alpha = True):
    # set up
    ZH = group_ring(tri, angle, alpha = alpha)
    P = ZH.polynomial_ring()
    if verbose > 0: print 'angle', angle

    ET = edges_to_tetrahedra_matrix(tri, angle, ZH, P)
    return normalise_poly(ET.determinant(), ZH, P)


# computing the small veering polynomial


@liberal
def edges_to_triangles_matrix(triangulation, angle_structure, ZH, P, mode = 'veering'):
    # In mode alexander, we are computing the transpose of the
    # boundary operator from relative 1-chains to relative 0-chains of
    # the dual 2-complex (relative to its vertex set).  Note that we
    # technically should only be taking a single dual vertex
    # downstairs and lifting that... but then we would have to compute
    # a kernel.  By using all of the dual vertices we are splitting
    # off a free summand of the correct rank.

    # In mode veering we are computing the matrix which assigns to a triangle the switch condition on its dual lower track
    coorientations = is_transverse_taut(triangulation, angle_structure, return_type = 'tet_vert_coorientations')
    if verbose > 0: print 'coorientations', coorientations
    face_laurents = faces_in_laurent(triangulation, angle_structure, ZH)
    if verbose > 0: print 'face_laurents', face_laurents

    ET_matrix = [] # now to find the face coefficients relative to each edge
    for tet in triangulation.tetrahedra():
        # get its upper edge - we iterate over upper edges of tetrahedra
        if verbose > 0: print 'tet_index', tet.index()
        edge = tet_lower_upper_edges(tet, coorientations)[1]
        if verbose > 0: print 'edge_index', edge.index()
        embeddings = edge.embeddings()

        # find index of tet in the list of embeddings of edge
        for i, embed in enumerate(embeddings):
            tet = embed.tetrahedron()
            if verbose > 0: print 'current_tet', tet.index()
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            if coorientations[tet.index()][trailing_vert_num] == +1 and coorientations[tet.index()][leading_vert_num] == +1:
                bottom_index = i
                break
        # rotate
        embeddings = embeddings[bottom_index:] + embeddings[:bottom_index]

        face_coeffs = [ ZH(0) ] * 2 * triangulation.countTetrahedra()
        sign = 1 # we are going up the left side of the edge
        current_coeff = ZH(1)

        # because of a sign change, we install the first and last by hand.
        embed = embeddings[0]
        assert tet.index() == embed.tetrahedron().index() # sanity check

        if mode == 'veering':
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            leading_face = tet.triangle(trailing_vert_num)
            face_coeffs[leading_face.index()] = face_coeffs[leading_face.index()] + current_coeff
            trailing_face = tet.triangle(leading_vert_num)
            face_coeffs[trailing_face.index()] = face_coeffs[trailing_face.index()] + current_coeff
            if verbose > 0: print 'face_coeffs', face_coeffs
            embeddings = embeddings[1:-1] # so we can skip the first and last

        for embed in embeddings:
            tet = embed.tetrahedron()
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            leading_face = tet.triangle(trailing_vert_num)
            trailing_face = tet.triangle(leading_vert_num)
            if sign == 1:
                use_face = trailing_face
            else: # sign == -1
                use_face = leading_face
            current_coeff = current_coeff * (face_laurents[use_face.index()])**sign
            if verbose > 0: print 'current_coeff', current_coeff
            # have we reached the top?
            if coorientations[tet.index()][trailing_vert_num] == -1 and coorientations[tet.index()][leading_vert_num] == -1: ## we are the top embed
                sign = -1 # we are going down the right side of the edge
                current_coeff = current_coeff * (face_laurents[leading_face.index()])**sign
                if verbose > 0: print 'top current_coeff', current_coeff
            if mode == 'veering':
                face_coeffs[leading_face.index()] = face_coeffs[leading_face.index()] - current_coeff
            elif mode == 'alexander':
                face_coeffs[leading_face.index()] = face_coeffs[leading_face.index()] + sign*current_coeff
            if verbose > 0: print 'face_coeffs', face_coeffs

        ET_matrix.append(face_coeffs)

    # convert and return
    return matrix_laurent_to_poly(ET_matrix, ZH, P)


@liberal
def small_polynomial(tri, angle, alpha = True, mode = 'veering'):
    # set up
    ZH = group_ring(tri, angle, alpha = alpha)
    P = ZH.polynomial_ring()

    ET = edges_to_triangles_matrix(tri, angle, ZH, P, mode = mode)

    # compute via minors
    minors = ET.minors(tri.countTetrahedra())
    return normalise_poly(gcd(minors), ZH, P)


@liberal
def small_polynomial_via_tree(tri, angle, alpha = True, mode = 'veering'):
    # set up
    ZH = group_ring(tri, angle, alpha = alpha)
    P = ZH.polynomial_ring()

    ET = edges_to_triangles_matrix(tri, angle, ZH, P, mode = mode)
    tree_faces, non_tree_faces = build_spanning_dual_tree(tri)

    ET = ET.transpose()
    ET = Matrix([row for i, row in enumerate(ET) if i in non_tree_faces]).transpose()

    # compute via minors
    minors = ET.minors(tri.countTetrahedra())
    return normalise_poly(gcd(minors), ZH, P)


@liberal
def small_polynomial_via_smith(tri, angle, alpha = True, mode = 'veering'):
    # set up
    assert tri.homology().rank() == 1 # need the polynomial ring to be a PID
    ZH = group_ring(tri, angle, alpha = alpha, ring = QQ) # ditto
    P = ZH.polynomial_ring()

    ET = edges_to_triangles_matrix(tri, angle, ZH, P, mode = mode)

    # compute via smith normal form
    ETs = ET.smith_form()[0]
    a = tri.countEdges()
    ETs_reduced = Matrix([row[:a] for row in ETs])
    return normalise_poly(ETs_reduced.determinant(), ZH, P)


@liberal
def small_polynomial_via_tree_and_smith(tri, angle, alpha = True, mode = 'veering'):
    # set up
    assert tri.homology().rank() == 1 # need the polynomial ring to be a PID
    ZH = group_ring(tri, angle, alpha = alpha, ring = QQ) # ditto
    P = ZH.polynomial_ring()

    ET = edges_to_triangles_matrix(tri, angle, ZH, P, mode = mode)
    tree_faces, non_tree_faces = build_spanning_dual_tree(tri)

    ET = ET.transpose()
    ET = Matrix([row for i, row in enumerate(ET) if i in non_tree_faces]).transpose()

    # compute via smith normal form
    ETs = ET.smith_form()[0]
    a = tri.countEdges()
    ETs_reduced = Matrix([row[:a] for row in ETs])
    return normalise_poly(ETs_reduced.determinant(), ZH, P)


# Remarks on small_polynomial_via_smith and
# small_polynomial_via_tree_and_smith

# 1. Running on 'gLLAQacdefefjkaaqks_200210' appears to involve very
# large intermediate computations - I've never waited long enough for
# the function to actually return.  So there is a serious dichotomy:
# the version using smith normal form is exponentially faster most of
# the time, but is exponentially slower every once in a while.

# 2. Also, smith normal form wants the polynomial ring to be a PID, so
# we are restricted to using QQ[x].  So we need b_1 = 1 or we need to
# specialise the matrix ET before computing.  Both of these argue
# strongly against using smith normal form...


@liberal
def small_polynomial_via_interpolate(tri, angle, alpha = True):
    # In the one-variable case: (1) get an upper bound on the
    # degree of the answer (2) plug integers into the matrix and
    # compute smith normal form (over ZZ).  (3) Write down the
    # function you get. (4) use difference equations to compute the
    # derivatives, back compute, and win.
    # This might also work in the multivariable case.
    return None
