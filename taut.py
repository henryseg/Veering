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


edge_num_to_vert_pair  = {0: (0, 1), 1: (0, 2), 2: (0, 3), 3: (2, 3), 4: (1, 3), 5: (1, 2)}    
vert_pair_to_edge_num = {(0, 1):0, (1, 0):0, (0, 2):1, (2, 0):1, (0, 3):2, (3, 0):2, (1, 2):3, (2, 1):3, (1, 3):4, (3, 1):4, (2, 3):5, (3, 2):5}
        
vert_pair_to_edge_pair = {(0, 1): 0, (2, 3): 0, (0, 2): 1, (1, 3): 1, (0, 3): 2, (1, 2): 2}
unsorted_vert_pair_to_edge_pair = {(0, 1): 0, (1, 0): 0, (2, 3): 0, (3, 2): 0, (0, 2): 1, (2, 0): 1, (1, 3): 1, (3, 1): 1, (0, 3): 2, (3, 0): 2, (1, 2): 2, (2, 1): 2}

def edge_pair_to_edge_numbers(e):
    return (e, 5 - e)


# converting a taut isosig to (tri,angle) pair and back again


def isosig_to_tri_angle(isosig):
    """
    Given a taut isosig, returns an oriented regina triangulation and
    the list of angles for the taut angle structure, for the new
    labelling.
    """
    isosig, angle = isosig.split("_")
    tri = regina.Triangulation3.fromIsoSig(isosig)
    angle = [int(i) for i in list(angle)]
    fix_orientations(tri, angle)  # this does not alter what the angles_string should be
    assert tri.isOriented()
    return tri, angle


def isosig_from_tri_angle(tri, angle):
    """
    Given a triangulation and taut angle structure, generate the taut
    isosig.
    """
    isosig, isom = tri.isoSigDetail()  # isom is the mapping between the original triangulation and the isosig triangulation
    isosig_tri = regina.Triangulation3.fromIsoSig(isosig)
    # isosig_tri_angle_struct_list = apply_isom_to_angle_struct_list(angle, isom)
    angle = apply_isom_to_angle_struct_list(angle, isom)

    # # Now find the lexicographically smallest angle structure string of symmetries of the one we have
    # all_isoms = isosig_tri.findAllIsomorphisms(isosig_tri)
    # all_angles_lists = []
    # for isom in all_isoms:
    #     all_angles_lists.append( apply_isom_to_angle_struct_list(isosig_tri_angle_struct_list, isom) )
    # all_angles_strings = ["".join([str(num) for num in angles_list]) for angles_list in all_angles_lists]
    # all_angles_strings.sort()

    smallest_angle = lex_smallest_angle_structure(isosig_tri, angle)
    smallest_angle_string = "".join([str(num) for num in smallest_angle])

    return isosig + "_" + smallest_angle_string


# checking tautness


@liberal
def is_taut(tri, angle):
    totals = [0] * tri.countEdges()
    for i, tet in enumerate(tri.tetrahedra()):
        edge_nums = edge_pair_to_edge_numbers(angle[i])
        for e in edge_nums:
            totals[tet.edge(e).index()] += 1
    # print totals
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

def lex_smallest_angle_structure(tri, angle):
    """find the lexicographically smallest angle structure among symmetries of the one we have"""
    all_isoms = tri.findAllIsomorphisms(tri)
    all_angles = []
    for isom in all_isoms:
        all_angles.append( apply_isom_to_angle_struct_list(angle, isom) )
    all_angles.sort()
    return all_angles[0]

def num_taut_automorphisms(tri, angle):
    """find the number of taut automorphisms"""
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


def fix_orientations(tri, angle):
    """
    Fix the orientations of the tetrahedra in triangulation so that
    they are consistently oriented.  We choose how to flip each
    tetrahedron so that the angle structure list does not change.
    """
    old_orientations = find_orientations(tri)
    for i, orientation in enumerate(old_orientations):
        if orientation == -1:
            reverse_tet_orientation(tri, tri.tetrahedron(i), angle[i])


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




