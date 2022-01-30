#
# branched_surface.py
#

# Code for working with branched surfaces.

from taut import isosig_to_tri_angle, isosig_from_tri_angle, edge_number_to_edge_pair, edge_num_to_vert_pair
from taut import vert_pair_to_edge_num, unsorted_vert_pair_to_edge_pair, edge_pair_to_edge_numbers
from itertools import product
import string

from veering import is_veering
from transverse_taut import is_transverse_taut

def isosig_to_tri_angle_branch(isosig):
    """
    Given a taut branched isosig, returns an oriented regina triangulation,
    the list of angles for the taut angle structure, and the branched surface
    for the new labelling.
    """
    sig_parts = isosig.split("_")
    tri, angle, isom = isosig_to_tri_angle(sig_parts[0] + '_' + sig_parts[1], return_isom = True)
    branch = list(sig_parts[2])
    branch = [ord(letter) - 97 for letter in branch]
    assert all([0 <= l and l <= 11 for l in branch])
    branch = apply_isom_to_branched_surface(branch, isom)
    assert is_branched(tri, branch) 
    return tri, angle, branch


def isosig_from_tri_angle_branch(tri, angle, branch):
    """
    Given a triangulation and taut angle structure and a branched surface, generate the taut branched
    isosig.
    """
    taut_isoSig, isom, regina_tri = isosig_from_tri_angle(tri, angle, return_isom = True, return_Regina_tri = True)
    branch = apply_isom_to_branched_surface(branch, isom)
    assert is_branched(regina_tri, branch)
    branch_sig = "".join([string.ascii_lowercase[b] for b in branch])
    return taut_isoSig + '_' + branch_sig
    

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

def branch_num_to_large_edge_and_mixed_edge_pair_num(branch_num, return_small = False):
    large_edge_num = branch_num // 2 ### divide by two, round down
    handedness_bit = branch_num % 2
    large_tiny_edge_pair_num = edge_number_to_edge_pair(large_edge_num)
    if handedness_bit == 0:
        mixed_edge_pair_num = (large_tiny_edge_pair_num + 1) % 3
        small_edge_pair_num = (large_tiny_edge_pair_num - 1) % 3
    else:
        mixed_edge_pair_num = (large_tiny_edge_pair_num - 1) % 3
        small_edge_pair_num = (large_tiny_edge_pair_num + 1) % 3
    if return_small:
        return large_edge_num, mixed_edge_pair_num, small_edge_pair_num
    else:
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
        new_mixed_edge_pair_num = unsorted_vert_pair_to_edge_pair[tuple(new_mixed_edge_verts)]
        # new_mixed_edge_verts.sort()
        # new_mixed_edge_pair_num = vert_pair_to_edge_pair[tuple(new_mixed_edge_verts)]

        new_branch[mapped_tet_index] = branch_num_from_large_edge_and_mixed_edge_pair_num(new_large_edge_num, new_mixed_edge_pair_num)
    return new_branch

def apply_swaps_to_branched_surface(branch, swaps): ### workaround since regina 6 didnt support constructing our own isomorphisms (could fix in regina 7)
    """
    Given a branched surface and a list of Perm4s, one for each tetrahedron,
    apply the permutations to the branched surface
    """
    for i in range(len(branch)):
        tetPerm = swaps[i]
        large_edge_num, mixed_edge_pair_num = branch_num_to_large_edge_and_mixed_edge_pair_num(branch[i])
        large_edge_verts = edge_num_to_vert_pair[large_edge_num]
        new_large_edge_verts = tuple([tetPerm[v] for v in large_edge_verts])
        new_large_edge_num = vert_pair_to_edge_num[new_large_edge_verts]

        mixed_edge_verts = (0, mixed_edge_pair_num + 1) ### one of the mixed edges
        new_mixed_edge_verts = [tetPerm[v] for v in mixed_edge_verts]
        new_mixed_edge_pair_num = unsorted_vert_pair_to_edge_pair[tuple(new_mixed_edge_verts)]
        # new_mixed_edge_verts.sort()
        # new_mixed_edge_pair_num = vert_pair_to_edge_pair[tuple(new_mixed_edge_verts)]

        branch[i] = branch_num_from_large_edge_and_mixed_edge_pair_num(new_large_edge_num, new_mixed_edge_pair_num)
    # return new_branch

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
    # print(all_branches)
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

def has_non_sing_semiflow(tri, branch):
    if not is_branched(tri, branch):
        return False
    for edge in tri.edges():
        mixed_count = 0
        for embed in edge.embeddings():
            b = branch[embed.simplex().index()]
            large_edge_num, mixed_edge_pair_num = branch_num_to_large_edge_and_mixed_edge_pair_num(b)
            verts = embed.vertices()
            edge_pair = unsorted_vert_pair_to_edge_pair[ (verts[0], verts[1]) ]
            if edge_pair == mixed_edge_pair_num:
                mixed_count += 1
                if mixed_count > 2:
                    return False
        if mixed_count != 2:
            return False
    return True

def all_branched_surfaces(tri):
    ### warning, this will be exponentially slow for not tiny triangulations
    n = tri.countTetrahedra()
    candidates = product([0,1,2,3,4,5,6,7,8,9,10,11], repeat = n)
    out = []
    for cand in candidates:
        if is_branched(tri, cand):
            out.append(cand)
    return out

def determine_possible_branch_given_two_faces(faces, large_verts):
    """
    Given train tracks on two faces of a tet, what are the possible branch surfaces inside it
    """
    face0, face1 = faces   # large vert is vertex number in the tet that in the face is opposite large edge of that face
    large_vert0, large_vert1 = large_verts
    assert face0 != face1
    assert face0 != large_vert0
    assert face1 != large_vert1

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

def determine_possible_branch_given_three_faces(faces, large_verts):
    """
    Given train tracks on three faces of a tet, what are the possible branch surfaces inside it
    """
    assert len(set(faces)) == 3
    assert faces[0] != large_verts[0]
    assert faces[1] != large_verts[1]
    assert faces[2] != large_verts[2]
    fourth_face = [i for i in [0,1,2,3] if i not in faces][0]
    ### find mixed edges
    mixed_indices = [] ### face indices that point at mixed edges
    for i in range(3):
        if large_verts[i] in faces:
            j = faces.index(large_verts[i])
            if large_verts[j] != faces[i]:
                mixed_indices.append(i)
    if len(mixed_indices) != 1:
        return None
    mixed_index = mixed_indices[0]
    face_pointing_at_mixed_edge = faces[mixed_index]
    face_pointing_away_from_mixed_edge = large_verts[mixed_index]
    mixed_edge_pair_num = unsorted_vert_pair_to_edge_pair[(face_pointing_at_mixed_edge, face_pointing_away_from_mixed_edge)]

    j = faces.index(large_verts[mixed_index])
    face_pointed_at_by_face_pointing_away_from_mixed_edge = large_verts[j]
    if face_pointed_at_by_face_pointing_away_from_mixed_edge == fourth_face:  ### then (fourth_face, face_pointing_away_from_mixed_edge) is tiny
        large_edge_num = 5 - vert_pair_to_edge_num[(fourth_face, face_pointing_away_from_mixed_edge)]
        return branch_num_from_large_edge_and_mixed_edge_pair_num(large_edge_num, mixed_edge_pair_num)
    else:
        large_edge_num = vert_pair_to_edge_num[(fourth_face, face_pointing_at_mixed_edge)]
        return branch_num_from_large_edge_and_mixed_edge_pair_num(large_edge_num, mixed_edge_pair_num)

def check_consistency():
    for branch_num in range(12):
        for i in range(4):
            large_verti = large_edge_of_face(branch_num, i)
            for j in range(i):
                large_vertj = large_edge_of_face(branch_num, j)
                faces = (i, j)
                large_verts = (large_verti, large_vertj)
                assert branch_num in determine_possible_branch_given_two_faces(faces, large_verts)

    for branch_num in range(12):
        for i in range(4):
            large_verti = large_edge_of_face(branch_num, i)
            for j in range(i):
                large_vertj = large_edge_of_face(branch_num, j)
                for k in range(j):
                    large_vertk = large_edge_of_face(branch_num, k)
                    faces = [i,j,k]
                    large_verts = [large_verti, large_vertj, large_vertk]
                    # print(faces, large_verts, branch_num)
                    # print(determine_possible_branch_given_three_faces(faces, large_verts))
                    assert branch_num == determine_possible_branch_given_three_faces(faces, large_verts)

def upper_branched_surface(tri, angle, return_lower = False):
    """Returns the upper branched surface for a veering triangulation"""
    veering_colours = is_veering(tri, angle, return_type = "veering_colours")
    tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    assert veering_colours != False and tet_vert_coorientations != False
    branch = []
    for i, a in enumerate(angle):
        assert tet_vert_coorientations[i][0] == tet_vert_coorientations[i][a+1] ### the other end of the pi edge at 0
        large_edge_num = vert_pair_to_edge_num[(0, a+1)]
        tiny_edge_num = 5 - large_edge_num
        if (tet_vert_coorientations[i][0] == -1) != return_lower: ### into the tetrahedron through face 0, xor get lower branched surface
            large_edge_num, tiny_edge_num = tiny_edge_num, large_edge_num # swap
        tiny_edge_colour = veering_colours[tri.tetrahedron(i).edge(tiny_edge_num).index()]
        mixed_edge_pair_num = (vert_pair_to_edge_num[(0, a+1)] + 1) % 3  ### can view as an edge num too...
        if veering_colours[tri.tetrahedron(i).edge(mixed_edge_pair_num).index()] != tiny_edge_colour:
            mixed_edge_pair_num = (mixed_edge_pair_num + 1) % 3 ### then its the other one
        branch.append( branch_num_from_large_edge_and_mixed_edge_pair_num(large_edge_num, mixed_edge_pair_num) )
    # assert is_branched(tri, branch)
    assert has_non_sing_semiflow(tri, branch)
    return branch

### The following function doesnt care about the transverse taut structure, and probably it should.
### We should also have the pachner moves carry a transverse taut structure with them.
def branch_veeringness(tri, angle, branch):
    """Count the number of tetrahedra for which the branch looks like it should in a veering triangulation"""
    count = 0
    for i in range(tri.countTetrahedra()):
        a = angle[i]
        b = branch[i]
        large_edge_num, _ = branch_num_to_large_edge_and_mixed_edge_pair_num(b)
        if large_edge_num in [a, 5 - a]:
            count += 1
    return count










