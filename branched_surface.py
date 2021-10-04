#
# branched_surface.py
#

# Code for working with branched surfaces.

from taut import isosig_to_tri_angle, isosig_from_tri_angle, edge_number_to_edge_pair, edge_num_to_vert_pair
from taut import vert_pair_to_edge_num, unsorted_vert_pair_to_edge_pair, edge_pair_to_edge_numbers
from itertools import product

def isosig_to_tri_angle_branch(isosig):
    """
    Given a taut branched isosig, returns an oriented regina triangulation,
    the list of angles for the taut angle structure, and the branched surface
    for the new labelling.
    """
    sig_parts = isosig.split("_")
    tri, angle = isosig_to_tri_angle(sig_parts[0] + '_' + sig_parts[1])
    branch = list(sig_parts[2])
    branch = [ord(letter) - 97 for letter in branch]
    assert all([0 <= l and l <= 11 for l in branch])
    return tri, angle, branch


def isosig_from_tri_angle_branch(tri, angle, branch):
    """
    Given a triangulation and taut angle structure, generate the taut
    isosig.
    """
    pass

# dealing with relabelling triangulations

### Given vertex labels of a tet, specify the branched surface in that tet as follows
### Record the large edge of the branched surface as an edge number for this tet (0 through 5)
### This gives an edge pair for the (large, tiny) ((L,T) on the diagram below)
### We also need to know the pair of edges that are mixed (M). The remaining pair are small (S).
### Mixed pair is either +1 or -1 mod 3 in the edge pair numbering convention.
###   

### Train tracks:
###                      M
###             0 +------------+ 3
###               |     |    ,'|
###               |   ,'|  ,'  |
###   top faces  S|--"  |,'    | 
###    of tet     |   T,'|  _--|S
###               |  ,'  |,'   |
###               |,'    |     |
###             2 +------------+ 1
###                      M                       
###
###             0 +------------+ 3
###               |`.    \     |
###  bottom faces |  `.   |,-. |
###     of tet    |    `.,'   "| 
###               |_  _,'`.L   |
###               | "" |   `.  |
###               |     \    `.|
###             2 +------------+ 1

### This is branched surface with edge 0 and mixed "offset" -1. This gets branched surface "b" (in the ordering "a" to "l")

def branch_num_to_large_edge_and_mixed_edge_pair_num(branch_num):
    large_edge_num = branch_num // 2 ### divide by two, round down
    handedness_bit = branch_num % 2
    large_tiny_edge_pair_num = edge_number_to_edge_pair(large_edge_num)
    if handedness_bit == 0:
        mixed_edge_pair_num = (large_tiny_edge_pair_num + 1) % 3
    else:
        mixed_edge_pair_num = (large_tiny_edge_pair_num - 1) % 3
    return large_edge_num, mixed_edge_pair_num

def branch_num_from_large_edge_and_mixed_edge_pair_num(large_edge_num, mixed_edge_pair_num):
    large_tiny_edge_pair_num = edge_number_to_edge_pair(large_edge_num)
    if mixed_edge_pair_num == (large_tiny_edge_pair_num + 1) % 3:
        handedness_bit = 0
    else:
        assert mixed_edge_pair_num == (large_tiny_edge_pair_num - 1) % 3
        handedness_bit = 1
    return 2 * large_edge_num + handedness_bit

def apply_isom_to_branched_surface(branch, isom):
    """
    Given a branched surface and an isomorphism of a triangulation,
    return the branched surface relative to the new triangulation.
    """
    new_branch = [None] * len(branch)
    for i in range(len(branch)):
        mapped_tet_index = isom.tetImage(i)
        tetPerm = isom.facetPerm(i)
        large_edge_num, mixed_edge_pair_num = branch_num_to_large_edge_and_mixed_edge_pair_num(branch[i])
        large_edge_verts = edge_num_to_vert_pair[large_edge_num]
        new_large_edge_verts = tuple([tetPerm[v] for v in large_edge_verts])
        new_large_edge_num = vert_pair_to_edge_num[new_large_edge_verts]

        mixed_edge_verts = (0, mixed_edge_pair_num + 1) ### one of the mixed edges
        new_mixed_edge_verts = [tetPerm[v] for v in mixed_edge_verts]
        new_mixed_edge_pair_num = unsorted_vert_pair_to_edge_pair(new_mixed_edge_verts)
        # new_mixed_edge_verts.sort()
        # new_mixed_edge_pair_num = vert_pair_to_edge_pair[tuple(new_mixed_edge_verts)]

        new_branch[mapped_tet_index] = branch_num_from_large_edge_and_mixed_edge_pair_num(new_large_edge_num, new_mixed_edge_pair_num)
    return new_branch

def lex_smallest_branched_surface(tri, branch):   
    """
    Finds the lexicographically smallest branched surface among
    symmetries of the one we have.
    """

    ### Maybe we dont want to use this because we always have an angle structure:
    ### The lex smallest angle structure determines the numbering of tetrahedra and vertices.

    all_isoms = tri.findAllIsomorphisms(tri)
    all_branches = []
    for isom in all_isoms:
        all_branches.append( apply_isom_to_branched_surface(branch, isom) )
    all_branches.sort()
    print(all_branches)
    return all_branches[0]

def large_edge_of_face(branch_num, face_num):
    """
    Which edge of the face in this tet is large for the traintrack in that face?
    Returns the vertex number opposite that edge, a number between 0 and 3 for this tet
    """
    large_edge_num, mixed_edge_pair_num = branch_num_to_large_edge_and_mixed_edge_pair_num(branch_num)
    large_edge_verts = edge_num_to_vert_pair[large_edge_num]
    tiny_edge_verts = [n for n in [0,1,2,3] if n not in large_edge_verts]

    if face_num in tiny_edge_verts:
        return tiny_edge_verts[1 - tiny_edge_verts.index(face_num)]  ## the other one
    else:
        assert face_num in large_edge_verts
        mixed_edge_nums = list(edge_pair_to_edge_numbers(mixed_edge_pair_num))
        mixed_edges_vert_pairs = [edge_num_to_vert_pair[edge_num] for edge_num in mixed_edge_nums]
        for pair in mixed_edges_vert_pairs:
            if face_num in pair:
                return pair[1 - pair.index(face_num)]

def is_branched(tri, branch):
    """
    Determine whether or not the branched surfaces in each tetrahedron
    glue up consistently to give a branched surface in the manifold.
    """

    for face_num in range(tri.countTriangles()):
        face = tri.triangle(face_num)

        embed0 = face.embedding(0)
        tet_num0 = embed0.simplex().index()
        tet_0_face_num = embed0.face()
        vertices0 = embed0.vertices() # Maps vertices (0,1,2) of face to the corresponding vertex numbers of tet0

        large_edge_in_face0 = vertices0.inverse()[large_edge_of_face(branch[tet_num0], tet_0_face_num)]

        embed1 = face.embedding(1)
        tet_num1 = embed1.simplex().index()
        tet_1_face_num = embed1.face()
        vertices1 = embed1.vertices() # Maps vertices (0,1,2) of face to the corresponding vertex numbers of tet1

        large_edge_in_face1 = vertices1.inverse()[large_edge_of_face(branch[tet_num1], tet_1_face_num)]

        if large_edge_in_face0 != large_edge_in_face1:
            return False
    return True

def all_branched_surfaces(tri):
    n = tri.countTetrahedra()
    candidates = product([0,1,2,3,4,5,6,7,8,9,10,11], repeat = n)
    out = []
    for cand in candidates:
        if is_branched(tri, cand):
            out.append(cand)
    return out

def determine_possible_branch_given_two_faces(faces_info):
    """
    Given train tracks on two faces of a tet, what are the possible branch surfaces inside it
    """
    f0, f1 = faces_info
    face0, large_vert0 = f0   # large vert is vertex number in the tet that in the face is opposite large edge of that face
    face1, large_vert1 = f1

    possible_branches = []

    ## 4 cases:

    ### large edges for the two faces are at the shared edge: get two branched surfaces - choice of which is mixed edge pair
    ### one large edge is on shared edge, the other is not: shared edge is mixed. other points at large edge for tet
    ### neither large edge is on shared edge but they are not opposite in the tet: get two branched surfaces - either could be large
    ### neither large edge is on shared edge and they are opposite in the tet: shared edge is tiny

    if face0 == large_vert1 and face1 == large_vert0: ### case 1
        ### large edges for the two faces are at the shared edge: get two branched surfaces - choice of which is mixed edge pair
        large_edge_num = 5 - vert_pair_to_edge_num[(face0, face1)]
        for handedness_bit in range(2):
            possible_branches.append(2 * large_edge_num + handedness_bit)
    elif (face0 == large_vert1) ^ (face1 == large_vert0): # xor ### case 2
        if face1 == large_vert0:
            face_pointing_at_shared = face0
            face_pointing_away = face1
            large_vert_pointing_at_shared = large_vert0
            large_vert_pointing_away = large_vert1
        else:
            face_pointing_at_shared = face1
            face_pointing_away = face0
            large_vert_pointing_at_shared = large_vert1
            large_vert_pointing_away = large_vert0
        large_edge_num = 5 - vert_pair_to_edge_num[(face_pointing_away, large_vert_pointing_away)]
        mixed_edge_pair_num = unsorted_vert_pair_to_edge_pair[(face_pointing_at_shared, face_pointing_away)]
        possible_branches.append( branch_num_from_large_edge_and_mixed_edge_pair_num(large_edge_num, mixed_edge_pair_num) )
    else:  ### neither face0 == large_vert1 nor face1 == large_vert0
        if large_vert0 == large_vert1: ### case 3
            large_edge_for_face0 = 5 - vert_pair_to_edge_num[(face0, large_vert0)]
            large_edge_for_face1 = 5 - vert_pair_to_edge_num[(face1, large_vert1)]
            possible_branches.append( branch_num_from_large_edge_and_mixed_edge_pair_num(large_edge_for_face0, edge_number_to_edge_pair(large_edge_for_face1) ))
            possible_branches.append( branch_num_from_large_edge_and_mixed_edge_pair_num(large_edge_for_face1, edge_number_to_edge_pair(large_edge_for_face0) ))
        else: ### case 4
            large_edge_num = vert_pair_to_edge_num[(face0, face1)]
            mixed_edge_pair_num = unsorted_vert_pair_to_edge_pair[(face0, large_vert0)]
            possible_branches.append( branch_num_from_large_edge_and_mixed_edge_pair_num(large_edge_num, mixed_edge_pair_num) )
    return possible_branches

def check_consistency():
    for branch_num in range(12):
        for i in range(4):
            large_verti = large_edge_of_face(branch_num, i)
            for j in range(i):
                large_vertj = large_edge_of_face(branch_num, j)
                assert branch_num in determine_possible_branch_given_two_faces([(i, large_verti),(j, large_vertj)])








