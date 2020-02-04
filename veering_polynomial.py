#
# veering_polynomial.py
#

# Goal - compute the "big" veering polynomial as defined by Sam Taylor et al.

import regina

from string import ascii_lowercase as ascii

from sage.arith.misc import gcd
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.matrix.constructor import Matrix
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from taut import isosig_to_tri_angle
from transverse_taut import is_transverse_taut
from veering import is_veering
from edge_equation_matrix_taut import edge_equation_matrix_taut

verbose = 0

def build_spanning_dual_tree(triangulation):
    ### return list of dual edges in the tree
    explored_tetrahedra = [0]
    frontier_tet_faces = [(0,0), (0,1), (0,2), (0,3)]
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
    tree_faces.sort()
    non_tree_faces.sort()
    if verbose > 0: print 'tree and non-tree faces', (tree_faces, non_tree_faces)
    return (tree_faces, non_tree_faces)

def reduced_edge_equation_matrix_taut(triangulation, angle_structure):
    tree_faces, non_tree_faces = build_spanning_dual_tree(triangulation)
    N = edge_equation_matrix_taut(triangulation, angle_structure)
    if verbose > 0: print 'edge_equation_matrix', N
    # N is a 2n by n matrix - the rows are longer than the columns
    # are high.  So we have to zero out the _columns_ for the tree
    # edges.
    if verbose > 0: print 'reduced_edge_equation_matrix', [[a for i, a in enumerate(row) if i in non_tree_faces] for row in N]
    return [[a for i, a in enumerate(row) if i in non_tree_faces] for row in N]

def elem_vector(i, dim): 
    vec = [0]*dim
    vec[i] = 1
    return vector(vec)

def faces_in_homology(triangulation, angle_structure):
    tree_faces, non_tree_faces = build_spanning_dual_tree(triangulation)
    N = reduced_edge_equation_matrix_taut(triangulation, angle_structure)
    N = Matrix(N)
    N = N.transpose()
    S, U, V = N.smith_form()
    image_dim = sum(1 for a in S.diagonal() if a != 0)
    ambient_dim = S.dimensions()[0]
    betti = ambient_dim - image_dim
    zero_vec = vector([0]*betti)
    n = len(tree_faces) + len(non_tree_faces)
    face_vecs = []
    for i in range(n):
        if i in tree_faces:
            face_vecs.append(zero_vec)
        else:
            j = non_tree_faces.index(i)
            face_vecs.append( (U*elem_vector(j, ambient_dim))[image_dim:] )
    if verbose > 0: print 'face_vecs', face_vecs
    return face_vecs

def monomial(vector, L):
    mon = 1
    for i,e in enumerate(vector):
        mon = mon * L.gen(i)**e
    return mon

def faces_in_laurent(triangulation, angle_structure, alpha = False):
    face_vecs = faces_in_homology(triangulation, angle_structure)
    b = len(face_vecs[0]) # betti number
    if alpha:
        assert b < 26
        L = LaurentPolynomialRing(ZZ, list(ascii)[:b])
    else:
        L = LaurentPolynomialRing(ZZ, 'x', b)
    return [monomial(vec, L) for vec in face_vecs]

edge_vert_index_map = {(0,1):0, (0,2):1, (0,3):2, (1,2): 3, (1,3):4, (2,3):5 }

def tet_lower_upper_edges(tetrahedron, coorientations):
    tet_coor = coorientations[tetrahedron.index()]
    lower_edge_endpoints = [i for i in range(4) if tet_coor[i] == +1]
    lower_edge_endpoints.sort()
    lower_edge_num = edge_vert_index_map[ tuple(lower_edge_endpoints) ]
    upper_edge_num = 5 - lower_edge_num
    return ( tetrahedron.face(1, lower_edge_num), tetrahedron.face(1, upper_edge_num) )

def has_red_lower_edge(tetrahedron, coorientations, edge_colours):
    lower_edge = tet_lower_upper_edges(tetrahedron, coorientations)[0]
    return edge_colours[lower_edge.index()] == 'R'

def edges_to_tetrahedra_matrix(triangulation, angle_structure, alpha = False):
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
            
    face_laurents = faces_in_laurent(triangulation, angle_structure, alpha = alpha)
    if verbose > 0: print 'face_laurents', face_laurents

    matrix = []
    for tet in triangulation.tetrahedra():
        if verbose > 0: print 'tet_index', tet.index()
        edge = tet_lower_upper_edges(tet, coorientations)[1]
        if verbose > 0: print 'edge_index', edge.index()
        edge_colour = edge_colours[edge.index()]
        if verbose > 0: print 'edge_colour', edge_colour
        embeddings = edge.embeddings()
        tet_coeffs = [0] * triangulation.countTetrahedra()
        tet_coeffs[tet.index()] = 1  # bottom tet around the edge gets a 1
        if verbose > 0: print 'initial tet_coeffs', tet_coeffs
        current_coeff = 1
        ### find index of bottom embedding in the list of embedding
        for i, embed in enumerate(embeddings):
            tet = embed.tetrahedron()
            if verbose > 0: print 'current_tet', tet.index()
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            if coorientations[tet.index()][trailing_vert_num] == +1 and coorientations[tet.index()][leading_vert_num] == +1:
                bottom_index = i
                break
        embeddings = embeddings[bottom_index:] + embeddings[:bottom_index]
        sign = 1 ### are we going up or coming down the other side of the edge
        for embed in embeddings[1:]:  ## first entry already dealt with
            tet = embed.tetrahedron()
            vert_perm = embed.vertices() 
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            current_coeff = current_coeff * face_laurents[tet.face(2,leading_vert_num).index()]**sign
            if coorientations[tet.index()][trailing_vert_num] == -1 and coorientations[tet.index()][leading_vert_num] == -1: ## we are the top embed
                tet_coeffs[tet.index()] = tet_coeffs[tet.index()] - current_coeff
                sign = -1
            elif edge_colour == 'L' and tet in red_tetrahedra or edge_colour == 'R' and tet in blue_tetrahedra:
                tet_coeffs[tet.index()] = tet_coeffs[tet.index()] - current_coeff
            if verbose > 0: print 'current tet_coeffs', tet_coeffs
                
        matrix.append(tet_coeffs)
    return matrix

def big_polynomial(veer_sig):
    tri, angle = isosig_to_tri_angle(veer_sig)
    if verbose > 0: print 'angle', angle
    ET = edges_to_tetrahedra_matrix(tri, angle, alpha=True)
    ET = Matrix(ET)
    if verbose > 0: print 'edges to tetrahedra'
    if verbose > 0: print ET
    return ET.determinant()

def edges_to_triangles_matrix(triangulation, angle_structure, alpha = False):
    coorientations = is_transverse_taut(triangulation, angle_structure, return_type = 'tet_vert_coorientations')
    if verbose > 0: print 'coorientations', coorientations             
    face_laurents = faces_in_laurent(triangulation, angle_structure, alpha = alpha)
    if verbose > 0: print 'face_laurents', face_laurents

    matrix = [] # this will get the face coefficients relative to each edge
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
        
        face_coeffs = [0] * 2 * triangulation.countTetrahedra()
        sign = 1 ### are we going up or coming down the other side of the edge
        current_coeff = 1

        # because of a sign change, we install the first and last by hand.
        embed = embeddings[0]
        assert tet.index() == embed.tetrahedron().index() # sanity
        vert_perm = embed.vertices()
        trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
        leading_face = tet.triangle(trailing_vert_num)
        face_coeffs[leading_face.index()] = face_coeffs[leading_face.index()] + current_coeff
        trailing_face = tet.triangle(leading_vert_num)
        face_coeffs[trailing_face.index()] = face_coeffs[trailing_face.index()] + current_coeff
        if verbose > 0: print 'face_coeffs', face_coeffs
        
        for embed in embeddings[1:-1]:
            tet = embed.tetrahedron()
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            leading_face = tet.triangle(trailing_vert_num)
            trailing_face = tet.triangle(leading_vert_num)
            if sign == 1:
                use_face = trailing_face
            else:
                use_face = leading_face
            current_coeff = current_coeff * (face_laurents[use_face.index()])**sign
            if verbose > 0: print 'current_coeff', current_coeff
            # have we reached the top?
            if coorientations[tet.index()][trailing_vert_num] == -1 and coorientations[tet.index()][leading_vert_num] == -1: ## we are the top embed
                sign = -1
                current_coeff = current_coeff * (face_laurents[leading_face.index()])**sign
                if verbose > 0: print 'top current_coeff', current_coeff
            face_coeffs[leading_face.index()] = face_coeffs[leading_face.index()] - current_coeff
            if verbose > 0: print 'face_coeffs', face_coeffs
            
        # one last sanity check
        embed = embeddings[0]
        tet = embed.tetrahedron()        
        vert_perm = embed.vertices()
        trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
        leading_face = tet.triangle(trailing_vert_num)
        trailing_face = tet.triangle(leading_vert_num)
        assert current_coeff == face_laurents[trailing_face.index()]
            
        matrix.append(face_coeffs)
    return matrix

def laurent_to_poly(p):
    if verbose > 0: print p
    L = p.parent()
    gens = L.gens()
    if verbose > 0: print L, gens
    exps = p.exponents()
    if verbose > 0: print exps
    if len(gens) == 1:
        n = -min(exps)
    else: # num gens > 1
        exps = [list(e) for e in exps]
        flat = []
        for e in exps:
            flat = flat + e
        n = -min(flat)
    n = max(n, 0)
    mul = 1
    for gen in L.gens():
        mul =  mul * gen**n
    if verbose > 0: print mul
    q = p * mul
    L = q.parent()
    P = L.polynomial_ring()
    if verbose > 0: print 'q', q
    return P(q)

def small_polynomial(veer_sig):
    tri, angle = isosig_to_tri_angle(veer_sig)
    if verbose > 0: print 'angle', angle
    ET = edges_to_triangles_matrix(tri, angle, alpha=True)
    ET = Matrix(ET)
    if verbose > 0: print 'edges to triangles'
    if verbose > 0: print ET
    minors = ET.minors(tri.countTetrahedra())
    if verbose > 0: print minors
    minors = [laurent_to_poly(minor) for minor in minors if minor != 0]
    if verbose > 0: print minors
    return gcd(minors)

# veering_polynomial.big_polynomial('cPcbbbiht_12')
# veering_polynomial.big_polynomial('gLLPQcdfefefuoaaauo_022110')

