#
# transverse_taut.py
#

# functions for working with transverse taut ideal triangulations.

import regina # needed inside of imported files

from .taut import liberal, is_taut, apply_isom_to_angle_struct_list, vert_pair_to_edge_num, there_is_a_pi_here


vertexSplit = [[0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]]  

def add_coors_with_check(triangulation, coorientations, tet_num, edgepair, direction):
    """
    Adds coorientation directions to coorientations for this tet_num,
    checks if they match neighbours' coorientations
    """
    # edgepair is the taut angle structure pi angle direction is +1 if
    # from the 0,1 faces to the 2,3 faces if we are in vertexSplit[0],
    # and similarly for the other cases
    for i in range(0, 2):
        if coorientations[tet_num][vertexSplit[edgepair][i]] == -direction:
            return False
        else:
            coorientations[tet_num][vertexSplit[edgepair][i]] = direction
    for i in range(2, 4):
        if coorientations[tet_num][vertexSplit[edgepair][i]] == direction:
            return False
        else:
            coorientations[tet_num][vertexSplit[edgepair][i]] = -direction
    # now check whether these new coorientations match the tetrahedra
    # these are glued to:
    tet = triangulation.tetrahedron(tet_num)
    for face in range(4):
        adjtet, adjgluing = tet.adjacentTetrahedron(face), tet.adjacentGluing(face)
        adjtet_num = adjtet.index()
        if coorientations[tet_num][face] * coorientations[adjtet_num][adjgluing[face]] == 1:
            return False # failed neighbours check
    return True # passed neighbours checks

@liberal
def is_transverse_taut(triangulation, taut_angle_struct, return_type = "boolean"):
    assert is_taut(triangulation, taut_angle_struct)
    assert triangulation.isOriented()
    coorientations = []
    for i in range(triangulation.countTetrahedra()):
        coorientations.append([0,0,0,0]) # +1 is out of tet, -1 is into tet
    # first decide coorientation for tet 0
    edgepair = taut_angle_struct[0]
    if not add_coors_with_check(triangulation, coorientations, 0, edgepair, +1): # choose arbitrary direction for tet 0
        return False # tet 0 could have self gluings incompatible with the transverse taut structure, so we could fail here

    explored_tetrahedra = [0]
    frontier_tet_faces = [(0, 0), (0, 1), (0, 2), (0, 3)]
    while len(frontier_tet_faces) > 0:
        my_tet_num, my_face_num = frontier_tet_faces.pop()
        my_tet = triangulation.tetrahedron(my_tet_num)
        neighbour_tet = my_tet.adjacentSimplex(my_face_num)
        neighbour_face_num = my_tet.adjacentGluing(my_face_num)[my_face_num]
        if neighbour_tet.index() in explored_tetrahedra:
            frontier_tet_faces.remove((neighbour_tet.index(), neighbour_face_num))
        else:
            newcoor = -coorientations[my_tet_num][my_face_num] # opposite of the adjacent face's coorientation
            edgepair = taut_angle_struct[neighbour_tet.index()]
            # now calculate the direction that the flow points through neighbour_tet
            if neighbour_face_num in [vertexSplit[edgepair][i] for i in range(2)]:
                direction = newcoor # agrees with newcoor
            else:
                direction = -newcoor
            if not add_coors_with_check(triangulation, coorientations, neighbour_tet.index(), taut_angle_struct[neighbour_tet.index()], direction):
                return False

            frontier_tet_faces.extend( [(neighbour_tet.index(), i) for i in range(4) if i != neighbour_face_num] )
            explored_tetrahedra.append(neighbour_tet.index())
    if return_type == "tet_vert_coorientations":
        return coorientations
    elif return_type == "face_coorientations":
        return convert_tetrahedron_coorientations_to_faces(triangulation, coorientations)
    else:
        return True

def convert_tetrahedron_coorientations_to_faces(triangulation, coorientations):
    """
    returns a list of +-1 for each face. +1 if face numbering
    orientation, converted into a coorientation via the right hand
    rule, agrees with the coorientation coming from the transverse
    taut structure. In our veering tetrahedron pictures, +1 means that
    the face numbering orientation is anticlockwise.
    """
    out = []
    for face in triangulation.faces(2):
        does_face_agree_with_coorientation = []
        for i in range(2):
            face_embedding = face.embedding(i)
            tet = face_embedding.simplex()
            face_index_in_tet = face_embedding.vertices()[3] # which vertex of the tet is not on the face
            coorientation = coorientations[tet.index()][face_index_in_tet]
            does_face_agree_with_tet = face_embedding.vertices().sign() # +1 if they agree on orientation (?)
            does_face_agree_with_coorientation.append(coorientation * does_face_agree_with_tet)

        assert does_face_agree_with_coorientation[0] == does_face_agree_with_coorientation[1]
        out.append(does_face_agree_with_coorientation[0])
    return out

def get_tet_top_vert_nums(tet_vert_coorientations, tet_num):
    vert_coorientations = tet_vert_coorientations[tet_num]
    top_vert_nums = []
    for i in range(4):
        if vert_coorientations[i] == -1:
            top_vert_nums.append(i)
    assert len(top_vert_nums) == 2
    return top_vert_nums

def get_tet_top_and_bottom_edges(tet_vert_coorientations, tet):
    tet_num = tet.index()
    top_vert_nums = get_tet_top_vert_nums(tet_vert_coorientations, tet_num)
    bottom_vert_nums = list(set([0,1,2,3]) - set(top_vert_nums))
    top_edge_num = vert_pair_to_edge_num[tuple(top_vert_nums)]
    bottom_edge_num = vert_pair_to_edge_num[tuple(bottom_vert_nums)]
    return [tet.edge(top_edge_num), tet.edge(bottom_edge_num)]

def get_top_and_bottom_nums(tet_vert_coors, tet_num):
    t0, t1 = get_tet_top_vert_nums(tet_vert_coors, tet_num)
    bottom_vert_nums = [0,1,2,3]
    bottom_vert_nums.remove(t0)
    bottom_vert_nums.remove(t1)
    b0, b1 = bottom_vert_nums
    return [(t0,t1), (b0,b1)]


### This should be "taut_symmetry_group"

@liberal
def symmetry_group_size(tri, angle, return_isoms = False):
    isoms = tri.findAllIsomorphisms(tri)
    count = 0
    taut_isoms = []
    for isom in isoms:
        if isom.facePerm(0).sign() == 1:  ## fixes orientation of the triangulation
            if angle == apply_isom_to_angle_struct_list(angle, isom):  ## fixes taut angle structure
                coorientations = is_transverse_taut(tri, angle, return_type = 'tet_vert_coorientations')
                if coorientations[0][0] == coorientations[isom.tetImage(0)][isom.facetPerm(0)[0]]:  ## fixes transverse structure
                    count += 1
                    if return_isoms:
                        taut_isoms.append(isom)
    if return_isoms:
        return taut_isoms
    return count

def get_tet_above_edge(tri, angle, edge, tet_vert_coorientations = None, get_tet_below_edge = False):
    """
    Find the tetrahedron above (or below) an edge
    """
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    embeddings = edge.embeddings()
    for embed in embeddings:
        tet = embed.tetrahedron()
        top_edge, bottom_edge = get_tet_top_and_bottom_edges(tet_vert_coorientations, tet)
        if get_tet_below_edge: 
            if top_edge == edge:
                return tet
        else:
            if bottom_edge == edge:
                return tet

def top_bottom_embeddings_of_faces(tri, angle, tet_vert_coorientations = None):
    """
    returns two lists: one whose ith entry is the top embedding of face i, another whose ith entry is the bottom embedding of face i
    (top embedding = embedding as a top face)
    """
    
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    
    n = tri.countTetrahedra()
    
    top_embeddings = [None]*2*n
    bottom_embeddings = [None]*2*n
    
    for face in tri.triangles():
        embed0 = face.embedding(0)
        embed1 = face.embedding(1)
        tet0_index = embed0.simplex().index()
        tet1_index = embed1.simplex().index()
        index_of_face_in_tet0 = embed0.vertices()[3]
        index_of_face_in_tet1 = embed1.vertices()[3]
        if tet_vert_coorientations[tet0_index][index_of_face_in_tet0] == 1: # coorientation points out of tet0, i.e. tet0 is below f
            top_embedding = embed0
            bottom_embedding = embed1
        else:
            top_embedding = embed1
            bottom_embedding = embed0
        top_embeddings[face.index()] = top_embedding
        bottom_embeddings[face.index()] = bottom_embedding
        
    for i in range(2*n):
        assert top_embeddings[i] != None and bottom_embeddings[i] != None
    
    return top_embeddings, bottom_embeddings

def edge_side_face_collections(triangulation, angle_struct, tet_vert_coorientations = None, return_tets = False, order_bottom_to_top = True):
    """
    For each edge e, find the two collections of (face numbers, edge_index in that face corresponding to e)
    or also (tet numbers, edge index in that tet corresponding to e)
    on either side of the pis, ordered from bottom to top if we have coorientations
    """

    out_triangles = []
    out_tets = []
    for edge_num in range(triangulation.countEdges()):
        edge = triangulation.edge(edge_num)
        embeddings = edge.embeddings()

        triangle_num_sets = [[],[],[]]
        tet_num_sets = [[],[],[]]
        which_set_to_add_to = 0
        for embed in embeddings:
            tet = embed.tetrahedron()
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            # as we walk around the edge, leading is in front of us, trailing is behind us
            # FIX - the following link is broken. 
            # see http://regina.sourceforge.net/engine-docs/classregina_1_1NTetrahedron.html#a54d99721b2ab2a0a0a72b6216b436440
            f = tet.triangle(leading_vert_num)  ### this is the face behind the tetrahedron as we walk around
            fmapping = tet.faceMapping(2, leading_vert_num)
            index_of_opposite_vert_in_f = fmapping.inverse()[trailing_vert_num]
            assert f.edge(index_of_opposite_vert_in_f) == edge
            triangle_num_sets[which_set_to_add_to].append((f.index(), index_of_opposite_vert_in_f))
            if there_is_a_pi_here(angle_struct, embed):
                which_set_to_add_to += 1
            else:
                tet_num_sets[which_set_to_add_to].append((tet.index(), embed.face())) ### embed.face is the edge number of the tet corresponding to the edge
        triangle_num_sets = [triangle_num_sets[2] + triangle_num_sets[0], triangle_num_sets[1]]
        tet_num_sets = [tet_num_sets[2] + tet_num_sets[0], tet_num_sets[1]]
        ## we wrap around, have to combine the two lists on the side we started and ended on

        triangle_num_sets[1].reverse() ## flip one so they are going the same way.
        tet_num_sets[1].reverse()
        ## Now if we have coorientations, make them both go up
        if tet_vert_coorientations == None: ### try to install them
            tet_vert_coorientations = is_transverse_taut(triangulation, angle_struct, return_type = "tet_vert_coorientations")
        if tet_vert_coorientations != False:
            embed = embeddings[0]
            tet = embed.tetrahedron()
            vert_perm = embed.vertices()
            leading_vert_num = vert_perm[3]
            if (tet_vert_coorientations[tet.index()][leading_vert_num] == +1) != (not order_bottom_to_top): 
            ### then coorientation points out of this tetrahedron through this face,
            ### which is behind the tetrahedron as we walk around
            ### so we are the wrong way round (if we are ordering from bottom to top)
                triangle_num_sets[0].reverse()
                triangle_num_sets[1].reverse()
                tet_num_sets[0].reverse()
                tet_num_sets[1].reverse()
        out_triangles.append(triangle_num_sets)
        out_tets.append(tet_num_sets)
    if not return_tets:
        return out_triangles
    else:
        return out_triangles, out_tets


