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

verbose = 0


def build_spanning_dual_tree(triangulation):
    """
    Returns three lists - (dual) edges in the spanning tree, (dual)
    edges not in the spanning tree, and the distance of each tet (in
    the tree) to the root. We use the regina numbering to determine
    the tree.
    """
    explored_tetrahedra = [0]
    distances_to_root = [0] + [None]*(triangulation.countTetrahedra() - 1)
    frontier_tet_faces = [(0, 0), (0, 1), (0, 2), (0, 3)]
    tree_faces = []
    non_tree_faces = []
    
    while len(frontier_tet_faces) > 0:
        my_tet_num, my_face_num = frontier_tet_faces.pop()
        my_tet = triangulation.tetrahedron(my_tet_num)
        my_face = my_tet.face(2, my_face_num)
        neighbour_tet = my_tet.adjacentSimplex(my_face_num)
        neighbour_face_num = my_tet.adjacentGluing(my_face_num)[my_face_num]
        if neighbour_tet.index() in explored_tetrahedra:
            non_tree_faces.append(my_face.index())
            frontier_tet_faces.remove((neighbour_tet.index(), neighbour_face_num))
        else:
            tree_faces.append(my_face.index())
            frontier_tet_faces.extend( [(neighbour_tet.index(), i) for i in range(4) if i != neighbour_face_num] )
            explored_tetrahedra.append(neighbour_tet.index())
            distances_to_root[neighbour_tet.index()] = distances_to_root[my_tet.index()] + 1
    tree_faces.sort()
    non_tree_faces.sort()
    if verbose > 0:
        print(("tree and non-tree faces", (tree_faces, non_tree_faces)))
    return (tree_faces, non_tree_faces, distances_to_root)


def walk_towards_root(tet, tree_faces, distances_to_root):
    """
    returns the tet closer to the root, and the face we crossed
    """
    my_dist_to_root = distances_to_root[tet.index()]
    assert my_dist_to_root > 0 # if we are at the root, blow up
    for i in range(4):
        if tet.triangle(i).index() in tree_faces:
            neighbour = tet.adjacentSimplex(i)
            if distances_to_root[neighbour.index()] < my_dist_to_root:
                return (neighbour, tet.triangle(i).index())
    assert False

    
def build_non_tree_edge_cycles(triangulation):
    """
    For every non-tree (dual) edge, find the corresponding cycle.
    """
    tree_faces, non_tree_faces, distances_to_root = build_spanning_dual_tree(triangulation)
    cycles = []
    for orig_face_ind in non_tree_faces:
        embeddings = triangulation.triangle(orig_face_ind).embeddings()
        tet0, tet1 = (embed.simplex() for embed in embeddings)
        zero_faces = []
        one_faces = []
        while tet0.index() != tet1.index():
            dist0 = distances_to_root[tet0.index()]
            dist1 = distances_to_root[tet1.index()]
            if dist0 < dist1:
                tet1, face1_ind = walk_towards_root(tet1, tree_faces, distances_to_root)
                one_faces.append(face1_ind)
            else:
                tet0, face0_ind = walk_towards_root(tet0, tree_faces, distances_to_root)
                zero_faces.append(face0_ind)
        zero_faces.reverse()
        cycles.append([orig_face_ind] + one_faces + zero_faces)
    return cycles
                

def there_is_a_pi_here(angle_struct, embed):
    """
    Given an embedding of an edge in a tetrahedron, tells us if there
    is a pi at that edge.
    """
    tet = embed.tetrahedron()
    vert_perm = embed.vertices()
    vert_nums = [vert_perm[0], vert_perm[1]]
    vert_nums.sort()
    in_tet_edge_num = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)].index(tuple(vert_nums))
    if in_tet_edge_num <= 2:
        edgepair = in_tet_edge_num
    else:
        edgepair = 5 - in_tet_edge_num
    return angle_struct[tet.index()] == edgepair

def edge_side_face_collections(triangulation, angle_struct, tet_vert_coorientations = None):
    """
    For each edge e, find the two collections of (face numbers, edge_index in that face corresponding to e)
    on either side of the pis, ordered from bottom to top if we have coorientations
    """

    out = []
    for edge_num in range(triangulation.countEdges()):
        edge = triangulation.edge(edge_num)
        embeddings = edge.embeddings()

        triangle_num_sets = [[],[],[]]
        which_set_to_add_to = 0
        for embed in embeddings:
            tet = embed.tetrahedron()
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            # as we walk around the edge, leading is in front of us, trailing is behind us
            # see http://regina.sourceforge.net/engine-docs/classregina_1_1NTetrahedron.html#a54d99721b2ab2a0a0a72b6216b436440
            f = tet.triangle(leading_vert_num)  ### this is the face behind the tetrahedron as we walk around
            fmapping = tet.faceMapping(2, leading_vert_num)
            index_of_opposite_vert_in_f = fmapping.inverse()[trailing_vert_num]
            assert f.edge(index_of_opposite_vert_in_f) == edge
            triangle_num_sets[which_set_to_add_to].append((f.index(), index_of_opposite_vert_in_f))
            if there_is_a_pi_here(angle_struct, embed):
                which_set_to_add_to += 1
        triangle_num_sets = [triangle_num_sets[2] + triangle_num_sets[0], triangle_num_sets[1]]
        ## we wrap around, have to combine the two lists on the side we started and ended on

        triangle_num_sets[1].reverse() ## flip one so they are going the same way.
        ## Now if we have coorientations, make them both go up
        if tet_vert_coorientations != None:
            embed = embeddings[0]
            tet = embed.tetrahedron()
            vert_perm = embed.vertices()
            leading_vert_num = vert_perm[3]
            if tet_vert_coorientations[tet.index()][leading_vert_num] == +1: 
            ### then coorientation points out of this tetrahedron through this face,
            ### which is behind the tetrahedron as we walk around
            ### so we are the wrong way round
                triangle_num_sets[0].reverse()
                triangle_num_sets[1].reverse()
        out.append(triangle_num_sets)
    return out


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
    tree_faces, non_tree_faces, distances_to_root = build_spanning_dual_tree(triangulation)
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

    tree_faces, non_tree_faces, distances_to_root = build_spanning_dual_tree(triangulation)
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
        return LaurentPolynomialRing(ring, list(ascii)[:rank])
    else:
        return LaurentPolynomialRing(ring, "x", rank)


def faces_in_laurent(triangulation, angle_structure, cycles, ZH):
    face_vecs = faces_in_homology(triangulation, angle_structure, cycles)
    if len(face_vecs[0]) == 1:
        face_vecs = [vec[0] for vec in face_vecs]
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

    if len(image_on_gens[0]) == 1: # flatten if neccesary
        image_on_gens= [vec[0] for vec in image_on_gens]
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
    if len(min_exp) == 1:
        min_exp = min_exp[0]
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
    if poly.coefficients()[-1] < 0:
        poly = -poly
    return poly


### end of copied/modified code
