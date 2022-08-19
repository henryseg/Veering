#
# taut_homology.py
#

# Various computations in homology, given a regina triangulation
# and an angle_structure.

from string import ascii_lowercase as ascii

from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.matrix.constructor import Matrix
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing

from .fundamental_domain import spanning_dual_tree
from .transverse_taut import is_transverse_taut, edge_side_face_collections

verbose = 0


def edge_equation_matrix_taut(triangulation, angle_struct):
    """
    For each edge, find the face numbers on either side of the pis,
    put +1 for one side and -1 for the other.
    """
    matrix = []

    edge_sides = edge_side_face_collections(triangulation, angle_struct)
    for left_faces, right_faces in edge_sides:
        row = [0] * triangulation.countTriangles()
        for (face_num, vert) in left_faces:
            row[face_num] = row[face_num] + 1
        for (face_num, vert) in right_faces:
            row[face_num] = row[face_num] - 1
        matrix.append(row)
    return matrix


def edge_equation_matrix_taut_reduced(triangulation, angle_structure, cycles):
    _, non_tree_faces, _ = spanning_dual_tree(triangulation)
    N = edge_equation_matrix_taut(triangulation, angle_structure)
    N = N + cycles
    if verbose > 0:
        print(("edge equation matrix", N))
    # N is a 2n by n matrix - the rows are longer than the columns are
    # high.  So we delete the columns corresponding to the tree edges.
    reduced = [[a for i, a in enumerate(row) if i in non_tree_faces] for row in N]
    if verbose > 0:
        print(("edge equation matrix, reduced", reduced))
    return reduced


def elem_vector(i, dim):
    vec = [0] * dim
    vec[i] = 1
    return vector(vec)


def faces_in_smith(triangulation, angle_structure, cycles):
    N = edge_equation_matrix_taut_reduced(triangulation, angle_structure, cycles)
    N = Matrix(N)
    N = N.transpose()
    return N.smith_form()


def rank_of_quotient(S):
    image_dim = sum(1 for a in S.diagonal() if a != 0)
    ambient_dim = S.dimensions()[0]
    return ambient_dim - image_dim, ambient_dim, image_dim


def faces_in_homology(triangulation, angle_structure, cycles):
    S, U, V = faces_in_smith(triangulation, angle_structure, cycles)
    rank, ambient_dim, image_dim = rank_of_quotient(S)

    if len(cycles) == 0:
        assert rank == triangulation.homology().rank()

    zero_vec = vector([0] * rank)

    tree_faces, non_tree_faces, _ = spanning_dual_tree(triangulation)
    # We recompute the tree :( but we don't have to pass it around :)

    n = len(tree_faces) + len(non_tree_faces)
    face_vecs = []
    for i in range(n):
        if i in tree_faces:
            face_vecs.append(zero_vec)
        else:
            j = non_tree_faces.index(i)
            face_vecs.append( (U * elem_vector(j, ambient_dim))[image_dim:] )
    if verbose > 0:
        print(("face_vecs", face_vecs))
    return [tuple(vec) for vec in face_vecs]


def group_ring(triangulation, angle_structure, cycles, alpha = False, ring = ZZ):
    S, U, V = faces_in_smith(triangulation, angle_structure, cycles)
    rank, _, _ = rank_of_quotient(S)
    if alpha:
        assert rank < 26
        return LaurentPolynomialRing(ring, list(ascii)[:rank], rank)
    else:
        return LaurentPolynomialRing(ring, "x", rank)


def faces_in_laurent(triangulation, angle_structure, cycles, ZH):
    face_vecs = faces_in_homology(triangulation, angle_structure, cycles)
    return [ ZH( {vec:1} ) for vec in face_vecs]


def epimorphism_in_laurent(tri, angle, cycles, ZH):
    """
    The argument cycles specifies a group epimorphism from the
    manifold to the filled manifold.  This function returns the image
    of the generators of the group ring under the induced epimorphism.
    """
    n = tri.countTetrahedra()
    S,U,V = faces_in_smith(tri, angle, []) # basis before filling, so no cycles 
    r = rank_of_quotient(S)[0]
    S2, U2, V2 = faces_in_smith(tri, angle, cycles) # basis after filling
    r2 = rank_of_quotient(S2)[0]

    A = U.inverse().delete_columns(range(n+1-r)) 
    B = U2.delete_rows(range(n+1-r2))

    image_on_gens = (B*A).columns()
    image_on_gens = [tuple(col) for col in image_on_gens]

    image_in_laurent = [ZH( { image_on_gens[i]:1 } ) for i in range(r)]
    return image_in_laurent


# The code below is copied and modified (with permission) from
# https://github.com/3-manifolds/SnapPy/blob/master/python/snap/nsagetools.py


def join_lists(list_of_lists):
    out = []
    for l in list_of_lists:
        out = out + l
    return out


def uniform_exponents(poly):
    return [list(e) if hasattr(e, "__getitem__") else (e,) for e in poly.exponents()]


def monomial_multiplier(elts, ZH):
    if all(elt == 0 for elt in elts):
        # Zero (unlike other constants) has valuation -\infty.  This
        # can show up as an issue when computing the big polynomial.
        # For an example, compute ET for
        # "kLLLMPPkcdgfehijjijhshassqhdqr_1222011022"
        return ZH(1)
    elts = [ZH(elt) for elt in elts]
    A =  Matrix(ZZ, join_lists([ uniform_exponents(p) for p in elts]))
    min_exp = tuple( [min(row) for row in A.transpose()] )
    return ZH( {min_exp:1} )


def laurent_to_poly(elt, P):
    if type(elt) is int:
        return P(elt)
    return P( elt.dict() )


def matrix_laurent_to_poly(M, ZH, P):
    # convert to polynomials after shifting rows
    if verbose > 0:
        print("matrix")
        print(M)
    muls = [ monomial_multiplier(row, ZH) for row in M ]
    if verbose > 0:
        print(("muls", muls))
    return Matrix( [ [ laurent_to_poly(p / mul, P) for p in row ] for row, mul in zip(M, muls) ] )


def normalise_poly(poly, ZH, P):
    if verbose > 0:
        print(("poly", poly))
    if poly == 0:
        return poly
    mul = monomial_multiplier([poly], ZH)
    if verbose > 0:
        print(("mul", mul))
    poly = laurent_to_poly(poly / mul, P)
    # if poly.coefficients()[-1] < 0:  # I don't trust this.
    # if str(poly)[0] == "-":  # This is a hack, of course.
    if poly.lt().coefficients()[0] < 0:  # lt = leading term
        poly = -poly
    return poly

### end of copied/modified code
