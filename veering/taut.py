#
# taut.py
#

# functions for building, recording, and working with (transverse)
# taut ideal triangulations.  We assume that the given manifolds are
# orientable.


import regina
from functools import wraps


# convention - A "regina isosig" is a string specifying a labelled,
# triangulated three-manifold.  An "angle string" specifies the
# location of the pi-angles of a (presumably) taut structure on such a
# triangulation.  A "taut isosig" is a string of the form
#
# regina_isosig + "_" + angle_string


# "Be liberal in what you accept" - Postel's law
def liberal(func):
    @wraps(func)
    def liberal_wrapper(*args, **kwargs):
        if type(args[0]) is str:
            sig = args[0]
            if "_" in sig:
                tri, angle = isosig_to_tri_angle(sig)
                args = (tri, angle) + args[1:]
            else:
                tri = regina.Triangulation3.fromIsoSig(sig)
                args = (tri,) + args[1:]
        return func(*args, **kwargs)

    return liberal_wrapper


# converting between vertices, edges, and edge pairs.


edge_num_to_vert_pair  = {0: (0, 1), 1: (0, 2), 2: (0, 3), 3: (1, 2), 4: (1, 3), 5: (2, 3)}    
vert_pair_to_edge_num = {(0, 1):0, (1, 0):0, (0, 2):1, (2, 0):1, (0, 3):2, (3, 0):2, (1, 2):3, (2, 1):3, (1, 3):4, (3, 1):4, (2, 3):5, (3, 2):5}
        
vert_pair_to_edge_pair = {(0, 1): 0, (2, 3): 0, (0, 2): 1, (1, 3): 1, (0, 3): 2, (1, 2): 2}
unsorted_vert_pair_to_edge_pair = {(0, 1): 0, (1, 0): 0, (2, 3): 0, (3, 2): 0, (0, 2): 1, (2, 0): 1, (1, 3): 1, (3, 1): 1, (0, 3): 2, (3, 0): 2, (1, 2): 2, (2, 1): 2}

def edge_pair_to_edge_numbers(e):
    return (e, 5 - e)

def edge_number_to_edge_pair(n):
    if n < 3:
        return n
    else:
        return 5 - n

# converting a taut isosig to (tri,angle) pair and back again


def isosig_to_tri_angle(isosig, return_isom = False):
    """
    Given a taut isosig, returns an oriented regina triangulation and
    the list of angles for the taut angle structure, for the new
    labelling.
    If return_isom, returns the isomorphism from the regina isosig triang to the oriented triangulation
    """
    data = isosig.split("_")
    isosig, angle = data[0], data[1]  ## we don't care if there is extra data in the sig, such as a branched surface
    tri = regina.Triangulation3.fromIsoSig(isosig)
    angle = [int(i) for i in list(angle)]
    isom = fix_orientations(tri, angle, return_isom = return_isom)  # this does not alter what the angles_string should be
    assert tri.isOriented()
    assert is_taut(tri, angle)
    if return_isom:
        return tri, angle, isom
    else:
        return tri, angle

def isoms_move_tetrahedra_to_same_tetrahedra(isom1, isom2):
    assert isom1.size() == isom2.size()
    for i in range(isom1.size()):
        if isom1.simpImage(i) != isom2.simpImage(i):
            return False
    return True

def isosig_from_tri_angle(tri, angle, return_isom = False, return_Regina_tri = False):
    """
    Given a triangulation and taut angle structure, generate the taut
    isosig. If return_isom, give the isom from the original triang to 
    the symmetry of the isosig triang with lex smallest angle struct.
    """
    isosig, isosig_isom = tri.isoSigDetail()  # isom is the mapping from the original triangulation to the isosig triangulation
    isosig_tri = regina.Triangulation3.fromIsoSig(isosig)
    angle = apply_isom_to_angle_struct_list(angle, isosig_isom)

    smallest_angle, isom2 = lex_smallest_angle_structure(isosig_tri, angle, return_isom = True)
    smallest_angle_string = "".join([str(num) for num in smallest_angle])

    out = isosig + "_" + smallest_angle_string
    if return_isom or return_Regina_tri:
        out = [out]
    if return_isom:
        out.append(isom2 * isosig_isom)   ### should work in regina 7
    if return_Regina_tri:
        out.append(isosig_tri)

    return out

    # if return_isom:
    #     return isosig + "_" + smallest_angle_string, isom2 * isosig_isom  ### should work in regina 7
    #     # ### regina's implementation of isomorphisms doesn't give group operations, so we cannot write isom2*isom... instead we have to do this:
    #     # lex_smallest_tri = isom2.apply(isosig_tri)
    #     # all_isoms = tri.findAllIsomorphisms(lex_smallest_tri)
    #     # for isom in all_isoms:
    #     #     if isoms_move_tetrahedra_to_same_tetrahedra(isom, isosig_isom): ### uses fact that isom2 doesn't reorder the tetrahedra, only their vertices inside
    #     #         return isosig + "_" + smallest_angle_string, isom
    #     # assert False ### should never get here
    # else:
    #     return isosig + "_" + smallest_angle_string


# checking tautness

@liberal
def is_taut(tri, angle, return_totals = False):
    totals = [0] * tri.countEdges()
    for i, tet in enumerate(tri.tetrahedra()):
        edge_nums = edge_pair_to_edge_numbers(angle[i])
        for e in edge_nums:
            totals[tet.edge(e).index()] += 1
    if return_totals == True:
        return totals
    
    return all(total == 2 for total in totals)


# dealing with relabelling triangulations


def apply_isom_to_angle_struct_list(original_angle_struct_list, isom, return_edge_pair = True):
    """
    Given a taut angle structure and an isomorphism of a triangulation,
    return the taut angle structure relative to the new triangulation.
    """
    new_angle_struct_list = [None] * len(original_angle_struct_list)
    for i in range(len(original_angle_struct_list)):
        mapped_tet_index = isom.tetImage(i)
        original_triang_pi_edge = edge_num_to_vert_pair[original_angle_struct_list[i]]
        pi_edge = [isom.facetPerm(i)[original_triang_pi_edge[0]], isom.facetPerm(i)[original_triang_pi_edge[1]]]
        pi_edge.sort()
        pi_edge = tuple(pi_edge)
        if return_edge_pair:
            pi_number = vert_pair_to_edge_pair[pi_edge]
        else:
            pi_number = vert_pair_to_edge_num[pi_edge]
        new_angle_struct_list[mapped_tet_index] = pi_number
    return new_angle_struct_list

def lex_smallest_angle_structure(tri, angle, return_isom = False):
    """
    Finds the lexicographically smallest angle structure among
    symmetries of the one we have.
    """
    all_isoms = tri.findAllIsomorphisms(tri)
    # print('num isoms', len(all_isoms))
    all_angles = []
    for isom in all_isoms:
        all_angles.append( (apply_isom_to_angle_struct_list(angle, isom), isom) )
    all_angles.sort(key = lambda pair: pair[0])
    # print('all angles', all_angles)
    if return_isom:
        return all_angles[0]
    else:
        return all_angles[0][0]

def num_taut_automorphisms(tri, angle):
    """
    Finds the number of taut automorphisms
    """
    all_isoms = tri.findAllIsomorphisms(tri)
    all_angles = []
    for isom in all_isoms:
        all_angles.append( apply_isom_to_angle_struct_list(angle, isom) )
    return all_angles.count(all_angles[0])


# functions to deal with orientations

# Note that regina's native .orient doesn't tell you what labelling
# changes occurred, which we need


def find_orientations(triangulation):
    """
    Assuming the triangulation is orientable, finds the current
    orientations of the tetrahedra.
    """
    orientations = [0] * triangulation.countTetrahedra()
    orientations[0] = +1  # fix the orientation for tet 0
    frontier = [(0, 0), (0, 1), (0, 2), (0, 3)]
    while len(frontier) > 0:
        tet_num, face_num = frontier.pop()
        tet = triangulation.tetrahedron(tet_num)
        adjtet, adjgluing = tet.adjacentTetrahedron(face_num), tet.adjacentGluing(face_num)
        if orientations[adjtet.index()] != 0:  # we have already been to that tetrahedron
            assert orientations[adjtet.index()] == adjgluing.sign() * -1 * orientations[tet.index()], "triangulation not orientable!"
            frontier.remove((adjtet.index(), adjgluing[face_num]))  # remove the reverse gluing from the frontier
        else:
            orientations[adjtet.index()] = adjgluing.sign() * -1 * orientations[tet.index()]
            for i in range(4):
                if i != adjgluing[face_num]:  # don't go backwards
                    frontier.append((adjtet.index(), i))
    for orient in orientations:
        assert orient != 0, "triangulation not connected!"
    return orientations


def reverse_tet_orientation(triangulation, tet, pi_location):
    # only need triangulation here for diagnostics
    """
    Reglues tet into triangulation with reversed orientation.  Returns
    which permutation of the tetrahedron was used.
    """
    swaps = {0: regina.Perm4(1, 0, 2, 3), 1: regina.Perm4(2, 1, 0, 3), 2: regina.Perm4(3, 1, 2, 0)}
    swap = swaps[pi_location]
    adjtets = []
    oldadjgluings = []
    adjgluings = []
    self_faces = []
    other_faces = []
    for face in range(4):
        adjtet = tet.adjacentTetrahedron(face)
        adjtets.append(adjtet)
        oldadjgluing = tet.adjacentGluing(face)
        oldadjgluings.append(tet.adjacentGluing(face))
        if adjtet != tet:
            adjgluings.append(oldadjgluing * swap)  # swaps 0 with 1 in this tet, not the other tet
            other_faces.append(face)
        else:
            adjgluings.append(swap * oldadjgluing * swap)
            self_faces.append(face)
    tet.isolate()  # undo all gluings
    for face in other_faces:
        tet.join(swap[face], adjtets[face], adjgluings[face])
    if len(self_faces) > 0:  # assume it must be two
        assert len(self_faces) == 2, "4 self faces, must be an error..."
        face = self_faces[0]  # only need to make one of the identifications
        tet.join(swap[face], adjtets[face], adjgluings[face])
    return swap

def moves_tetrahedra(isom):
    for i in range(isom.size()):
        if isom.simpImage(i) != i:
            return True
    return False

def fix_orientations(tri, angle, return_isom = False):
    """
    Fix the orientations of the tetrahedra in triangulation so that
    they are consistently oriented.  We choose how to flip each
    tetrahedron so that the angle structure list does not change.
    If return_isom, return the isomorphism from the original tri to the fixed orientations tri
    """
    orig_tri = regina.Triangulation3(tri)

    old_orientations = find_orientations(tri)
    swaps = []
    for i, orientation in enumerate(old_orientations):
        if orientation == -1:
            swaps.append( reverse_tet_orientation(tri, tri.tetrahedron(i), angle[i]) )
        else:
            swaps.append( regina.Perm4() ) ## identity

    if return_isom:  
        out_isom = regina.Isomorphism3.identity(len(swaps))
        for i, p in enumerate(swaps):
            out_isom.setFacetPerm(i, p)
        return out_isom

        ### Regina 6 didn't let us build the isom directly from perms in python. Regina 7 does, using setFacetPerm
        ### old:
        # all_isoms = orig_tri.findAllIsomorphisms(tri) ### we will be order two, so we dont care which way this goes
        # for isom in all_isoms:
        #     if not moves_tetrahedra(isom):
        #         for i in range(tri.countTetrahedra()):
        #             assert swaps[i] == isom.facetPerm(i)
        #             assert isom == out_isom
        #         return isom
        assert False ## should never get here




# functions for converting regina's angle structure format to ours


def pi_edgepair(regina_angle_struct, tet_num):
    """
    Given a regina angle structure that is taut, tells us which pair
    of edges have the pi angle.
    """
    # 0 is edge pair 01|23, 1 is 02|13, 2 is 03|12
    for edgepair in range(3):
        if regina_angle_struct.angle(tet_num, edgepair) > 0:
            return edgepair
    assert False # we shouldn't be able to get here.


def taut_regina_angle_struct_to_taut_struct(regina_angle_struct):
    """
    Convert a taut regina angle structure to list of which edge pair
    is the pi pair.
    """
    out = []
    for tet_num in range(regina_angle_struct.triangulation().countTetrahedra()):
        out.append(pi_edgepair(regina_angle_struct, tet_num))
    return out


pi_edgepair_dict = { (1,0,0) : 0, (0,1,0) : 1, (0,0,1) : 2 }


def charge_to_angle(charge):
    """
    Given a list of 3*n integers with each triple of the form (1,0,0),
    (0,1,0), or (0,0,1), convert to our angle structure format.
    """
    assert len(charge) % 3 == 0
    n = int(round(len(charge)/3))
    out = []
    for i in range(n):
        tet = charge[ 3*i : 3*i+3 ]
        out.append( pi_edgepair_dict[tuple(tet)] )
    return out

def angle_to_charge(angle, flipper_format = False):
    """
    Given a list of n integers in [0,2], convert to charge format.
    """
    if flipper_format:
        # The veering code uses "vertex with 0 (minus one)". On the
        # other hand, flipper and t3m use "vertex with 3".
        # See line 25 of
        # https://github.com/MarkCBell/flipper/blob/master/flipper/kernel/taut.py
        angle = [2 - a for a in angle]
    out = [0] * (3*len(angle))
    for i, a in enumerate(angle):
        out[3*i + a] = 1
    if flipper_format:
        # flipper adds a variable to homogenise, so we do the same.
        out.append(-1)
    return out

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
