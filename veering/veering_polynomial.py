#
# veering_polynomial.py
#

# compute the lower (and upper) veering polynomials as defined by Sam
# Taylor et al.

from sage.arith.misc import gcd
from sage.rings.rational_field import QQ
from sage.matrix.constructor import Matrix

from .basic_math import sign
from .sage_tools import matrix_laurent_to_poly, normalise_poly
from .taut import liberal
from .transverse_taut import is_transverse_taut
from .taut_homology import edge_equation_matrix_taut, group_ring, faces_in_laurent
from .taut_polynomial import tet_lower_upper_edges
from .veering_tri import is_veering

verbose = 0

# computing the veering polynomial

def has_red_lower_edge(tetrahedron, coorientations, edge_colours):
    lower_edge = tet_lower_upper_edges(tetrahedron, coorientations)[0]
    return edge_colours[lower_edge.index()] == "red"


@liberal
def edges_to_tetrahedra_matrix(triangulation, angle_structure, coorientations, normalisation, ZH, P):

    edge_colours = is_veering(triangulation, angle_structure, return_type = "veering_colours")
    if verbose > 0:
        print(("edge_colours", edge_colours))
    red_tetrahedra = []
    blue_tetrahedra = []

    for tet in triangulation.tetrahedra():
        if has_red_lower_edge(tet, coorientations, edge_colours):
            red_tetrahedra.append(tet)
        else:
            blue_tetrahedra.append(tet)
    if verbose > 0:
        print(("how many reds and blues", len(red_tetrahedra), len(blue_tetrahedra)))

    face_laurents = faces_in_laurent(triangulation, angle_structure, [], ZH)  # empty list of cycles.
    if verbose > 0:
        print(("face_laurents", face_laurents))

    ET_matrix = []  # now to find the tet coefficients relative to each edge
    for tet in triangulation.tetrahedra():
        if verbose > 0:
            print(("tet_index", tet.index()))
        edge = tet_lower_upper_edges(tet, coorientations)[1]
        if verbose > 0:
            print(("edge_index", edge.index()))
        edge_colour = edge_colours[edge.index()]
        if verbose > 0:
            print(("edge_colour", edge_colour))
        embeddings = list(edge.embeddings())
        tet_coeffs = [ZH(0)] * triangulation.countTetrahedra()
        tet_coeffs[tet.index()] = 1  # bottom tet around the edge gets a 1
        if verbose > 0:
            print(("initial tet_coeffs", tet_coeffs))
        current_coeff = ZH(1)

        # find index of bottom embedding in the list of embedding
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

        embeddings = embeddings[bottom_index:] + embeddings[:bottom_index]
        sign = -1  # we are going up the left side of the edge
        for embed in embeddings[1:]:  # skipping the first
            tet = embed.tetrahedron()
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            current_coeff = current_coeff * face_laurents[tet.face(2,leading_vert_num).index()]**sign

            if (coorientations[tet.index()][trailing_vert_num] == -1 and
                coorientations[tet.index()][leading_vert_num]  == -1):
                # we are the top embed so:
                tet_coeffs[tet.index()] = tet_coeffs[tet.index()] - current_coeff
                sign = 1  # now we go down the right side
            elif ((edge_colour == "blue" and tet in red_tetrahedra) or
                  (edge_colour == "red" and tet in blue_tetrahedra)): 
                tet_coeffs[tet.index()] = tet_coeffs[tet.index()] - current_coeff
            if verbose > 0:
                print(("current tet_coeffs", tet_coeffs))

        ET_matrix.append(tet_coeffs)        
        
    if normalisation == True:    
    # convert and return
        return matrix_laurent_to_poly(ET_matrix, ZH, P)
    else: 
        return Matrix(ET_matrix)

def permutation_tet_top_diagonal(tri, angle, coorientations):
    # for the lower veering polynomial we identify a tetrahedron with its top diagonal (we always compute the lower - the upper is computed by switching coorientations)
    # to compute the unnormalised veering polynomial not only up to a sign, we need to find the parity of the permutation tet_index -> top_edge_index
    perm = []
    for tet in tri.tetrahedra():
        top_edge = tet_lower_upper_edges(tet, coorientations)[1]
        index_of_top = top_edge.index()
        perm.append(index_of_top)
    return perm
        

@liberal
def veering_polynomial(tri, angle, alpha = True, normalisation = True, mode = "lower"):
    # set up
    ZH = group_ring(tri, angle, [], alpha = alpha)
    P = ZH.polynomial_ring()
    if verbose > 0:
        print(("angle", angle))
           
    coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    
    if mode == "upper":
        coorientations = [[-x for x in coor] for coor in coorientations]
    
    ET = edges_to_tetrahedra_matrix(tri, angle, coorientations, normalisation, ZH, P)
    #print(ET)
    if normalisation == True:
        return normalise_poly(ET.determinant(), ZH, P)
    else:
        perm = permutation_tet_top_diagonal(tri, angle, coorientations)
        #print(perm, "sign", sign(perm))
        return sign(perm)*ET.determinant()


