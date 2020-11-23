#
# transverse_taut.py
#

# functions for working with transverse  taut ideal triangulations.

import regina # needed inside of imported files
from taut import liberal, apply_isom_to_angle_struct_list, vert_pair_to_edge_num

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

@liberal
def symmetry_group_size(tri, angle):
    isoms = tri.findAllIsomorphisms(tri)
    count = 0
    for isom in isoms:
        if isom.facePerm(0).sign() == 1:  ## fixes orientation of the triangulation
            if angle == apply_isom_to_angle_struct_list(angle, isom):  ## fixes taut angle structure
                coorientations = is_transverse_taut(tri, angle, return_type = 'tet_vert_coorientations')
                if coorientations[0][0] == coorientations[isom.tetImage(0)][isom.facetPerm(0)[0]]:  ## fixes transverse structure
                    count += 1
    return str(count)
