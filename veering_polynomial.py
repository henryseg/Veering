#
# veering_polynomial.py
#

# Goal - compute the "big" and "small" veering polynomials as defined
# by Sam Taylor et al.

import regina

from string import ascii_lowercase as ascii

from sage.arith.misc import gcd
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector
from sage.matrix.constructor import Matrix
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing

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
    assert betti == triangulation.homology().rank() 

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
    return [tuple(vec) for vec in face_vecs]

def group_ring(triangulation, angle_structure, alpha = False, ring = ZZ):
    b = triangulation.homology().rank() # betti number
    if alpha:
        assert b < 26
        return LaurentPolynomialRing(ring, list(ascii)[:b])
    else:
        return LaurentPolynomialRing(ring, 'x', b)

def faces_in_laurent(triangulation, angle_structure, ZH):
    face_vecs = faces_in_homology(triangulation, angle_structure)
    if len(face_vecs[0]) == 1:
        face_vecs = [vec[0] for vec in face_vecs]
    return [ ZH( {vec:1} ) for vec in face_vecs]

edge_vert_index_map = {(0,1):0, (0,2):1, (0,3):2, (1,2): 3, (1,3):4, (2,3):5 }

def tet_lower_upper_edges(tetrahedron, coorientations):
    tet_coor = coorientations[tetrahedron.index()]
    lower_edge_endpoints = [i for i in range(4) if tet_coor[i] == +1]
    lower_edge_endpoints.sort()
    lower_edge_num = edge_vert_index_map[ tuple(lower_edge_endpoints) ]
    upper_edge_num = 5 - lower_edge_num
    return ( tetrahedron.face(1, lower_edge_num), tetrahedron.face(1, upper_edge_num) )

### Code copied and modified from
### https://github.com/3-manifolds/SnapPy/blob/master/python/snap/nsagetools.py

def join_lists(list_of_lists):
    out = []
    for l in list_of_lists:
        out = out + l
    return out

def uniform_exponents(poly):
    return [list(e) if hasattr(e, "__getitem__") else (e,) for e in poly.exponents()]
    
def monomial_multiplier(elts, ZH):
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
    muls = [ monomial_multiplier(row, ZH) for row in M ]
    return Matrix( [ [ laurent_to_poly(p / mul, P) for p in row ] for row, mul in zip(M, muls) ] )

# computing the big veering polynomial

def has_red_lower_edge(tetrahedron, coorientations, edge_colours):
    lower_edge = tet_lower_upper_edges(tetrahedron, coorientations)[0]
    return edge_colours[lower_edge.index()] == 'R'

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
        for embed in embeddings[1:]:  # first entry already dealt with
            tet = embed.tetrahedron()
            vert_perm = embed.vertices() 
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            current_coeff = current_coeff * face_laurents[tet.face(2,leading_vert_num).index()]**sign
            if coorientations[tet.index()][trailing_vert_num] == -1 and coorientations[tet.index()][leading_vert_num] == -1:
                # we are the top embed so:
                tet_coeffs[tet.index()] = tet_coeffs[tet.index()] - current_coeff
                sign = -1 # we are going down the right side
            elif edge_colour == 'L' and tet in red_tetrahedra or edge_colour == 'R' and tet in blue_tetrahedra:
                tet_coeffs[tet.index()] = tet_coeffs[tet.index()] - current_coeff
            if verbose > 0: print 'current tet_coeffs', tet_coeffs
                
        ET_matrix.append(tet_coeffs)

    # convert and return
    return matrix_laurent_to_poly(ET_matrix, ZH, P)

def big_polynomial(veer_sig, alpha = True):
    # set up
    tri, angle = isosig_to_tri_angle(veer_sig)
    ZH = group_ring(tri, angle, alpha = alpha)
    P = ZH.polynomial_ring()
    if verbose > 0: print 'angle', angle

    ET = edges_to_tetrahedra_matrix(tri, angle, ZH, P)

    out_poly = ET.determinant()
    if out_poly == 0:
        return out_poly
    
    # normalize
    mul = monomial_multiplier([out_poly], ZH)
    out_poly = laurent_to_poly(out_poly / mul, P)
    if out_poly.coefficients()[-1] < 0:
        out_poly = -out_poly
    return out_poly
    
# computing the small veering polynomial

def edges_to_triangles_matrix(triangulation, angle_structure, ZH, P):
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
        sign = 1 ### are we going up or coming down the other side of the edge
        current_coeff = ZH(1)

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
            
        ET_matrix.append(face_coeffs)

    # convert and return
    return matrix_laurent_to_poly(ET_matrix, ZH, P)

def small_polynomial(veer_sig, alpha = True):
    # set up
    tri, angle = isosig_to_tri_angle(veer_sig)
    ZH = group_ring(tri, angle, alpha = alpha)
    P = ZH.polynomial_ring()

    ET = edges_to_triangles_matrix(tri, angle, ZH, P)

    # compute via minors
    minors = ET.minors(tri.countTetrahedra())
    out_poly = gcd(minors)
    if out_poly == 0:
        return out_poly

    # normalize
    mul = monomial_multiplier([out_poly], ZH)
    out_poly = laurent_to_poly(out_poly / mul, P)
    if out_poly.coefficients()[-1] < 0:
        out_poly = -out_poly
    return out_poly

### End of copied/modified code

def small_polynomial_via_smith(veer_sig, alpha = True):
    # set up
    tri, angle = isosig_to_tri_angle(veer_sig)
    assert tri.homology().rank() == 1 # need the polynomial ring to be a PID
    ZH = group_ring(tri, angle, alpha = alpha, ring = QQ) # ditto
    P = ZH.polynomial_ring()

    ET = edges_to_triangles_matrix(tri, angle, ZH, P)

    # compute via smith normal form
    ETs = ET.smith_form()[0]
    a = tri.countEdges()
    out_poly = Matrix([row[:a] for row in ETs]).determinant()
    if out_poly == 0:
        return out_poly

    # normalize
    mul = monomial_multiplier([out_poly], ZH)
    out_poly = laurent_to_poly(out_poly / mul, P)
    if out_poly.coefficients()[-1] < 0:
        out_poly = -out_poly
    return out_poly

# Remarks on the smith normal form version'gLLAQacdefefjkaaqks_200210'

def small_polynomial_via_interpolate(veer_sig, alpha = True):
    # In the one-variable case: (1) get an upper bound on the
    # degree of the answer (2) plug integers into the matrix and
    # compute smith normal form (over ZZ).  (3) Write down the
    # function you get. (4) use difference equations to compute the
    # derivatives, back compute, and win.
    # This might also work in the multivariable case... 
    return None
