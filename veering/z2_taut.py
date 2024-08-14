#
# taut.py
#

# functions for working with ZZ2 taut ideal tris.  

import regina
from copy import deepcopy

from .taut import liberal, edge_pair_to_edge_numbers, vert_pair_to_edge_pair
from .fundamental_domain import non_tree_face_loops


# check ZZ2 tautness


@liberal
def is_z2_taut(tri, angle):
    totals = [0] * tri.countEdges()
    for i, tet in enumerate(tri.tetrahedra()):
        edge_nums = edge_pair_to_edge_numbers(angle[i])
        for e in edge_nums:
            totals[tet.edge(e).index()] += 1
    # print totals
    return all(total%2 == 0 for total in totals)


def find_z2_taut_structures(tri):
    """try all possibilities of flat tetrahedra, looking for a z2 taut structure on tri"""

    num_edges = tri.countEdges()
    remaining_edge_degrees = [e.degree() for e in tri.edges()]
    pass_back_structures = []
    try_build_z2_taut_struct(tri, [], [0]*num_edges, remaining_edge_degrees, pass_back_structures)
    return pass_back_structures


def try_build_z2_taut_struct(tri, partial_taut_structure, num_pis_at_edges, remaining_edge_degrees, pass_back_structures):
    """recursive function explores tree of partial taut structures looking for z2 taut. builds taut structure
    by adding tetrahedra in sequentially by tetrahedron number"""

    if len(partial_taut_structure) == tri.countTetrahedra():
        # print partial_taut_structure, num_pis_at_edges, veering_directions
        assert is_z2_taut(tri, partial_taut_structure)
        pass_back_structures.append(partial_taut_structure)
        return True
    #try to add next tetrahedron in
    tet_index = len(partial_taut_structure)
    tet = tri.tetrahedron(tet_index)
    for edgepair in range(3):
        try_add_z2_taut_tet(tri, partial_taut_structure, num_pis_at_edges, remaining_edge_degrees, tet, edgepair, pass_back_structures)


def try_add_z2_taut_tet(tri, partial_taut_structure, num_pis_at_edges, remaining_edge_degrees, tet, edgepair, pass_back_structures):
    """try to add the given tet, tet_index to the structure with given edgepair being pi"""
    #first see if num_pis_at_edges is ok
    new_num_pis_at_edges = deepcopy(num_pis_at_edges)
    new_remaining_edge_degrees = deepcopy(remaining_edge_degrees)
    
    pi_edge_1, pi_edge_2 = edgepair, 5 - edgepair
    edge_1_index = tet.edge(pi_edge_1).index()
    edge_2_index = tet.edge(pi_edge_2).index()
    new_num_pis_at_edges[edge_1_index] += 1
    new_num_pis_at_edges[edge_2_index] += 1
    for e in range(6):
        e_index = tet.edge(e).index()
        new_remaining_edge_degrees[e_index] -= 1
        assert new_remaining_edge_degrees[e_index] >= 0
        if new_remaining_edge_degrees[e_index] == 0:
            if new_num_pis_at_edges[e_index] % 2 != 0:
                return False

    #survived this far, so go deeper into the tree
    new_partial_taut_structure = deepcopy(partial_taut_structure)
    new_partial_taut_structure.append(edgepair)
    #print new_partial_taut_structure, new_num_pis_at_edges, new_veering_directions
    try_build_z2_taut_struct(tri, new_partial_taut_structure, new_num_pis_at_edges, new_remaining_edge_degrees, pass_back_structures)


def is_trivial_in_cohomology(tri, angle):
    """
    Test all loops in this triangulation, do any of them pass an odd number
    of pi's. If so, return False
    """
    # print('angle', angle)
    loops = non_tree_face_loops(tri, include_tetrahedra = True)
    for (face_inds, tets) in loops:
        # print('face_inds', face_inds)
        n = len(tets)
        pi_count = 0
        for i in range(n):
            tet = tets[i]
            face_ind0 = face_inds[i]
            face_ind1 = face_inds[(i+1)%n]
            tet_face_indices = [tet.triangle(j).index() for j in range(4)]
            k0 = tet_face_indices.index(face_ind0)
            k1 = tet_face_indices.index(face_ind1)
            if k0 == k1:
                assert n == 1
                k1 = tet_face_indices.index(face_ind1, k0 + 1)  ## start looking from index k0 + 1
            edge_pair = vert_pair_to_edge_pair[tuple(sorted((k0, k1)))]
            # print('k0, tet, k1, edge_pair', k0, tet.index(), k1, edge_pair, angle[tet.index()] == edge_pair)
            if angle[tet.index()] == edge_pair:
                pi_count += 1
        if pi_count%2 != 0:
            return False
    return True


def cohomology_loops(tri):
    loops = non_tree_face_loops(tri, include_tetrahedra = True)
    num_tet = tri.countTetrahedra()
    out = []
    for (face_inds, tets) in loops:
        equ = [0] * (3 * num_tet)
        n = len(tets)
        for i in range(n):
            tet = tets[i]
            face_ind0 = face_inds[i]
            face_ind1 = face_inds[(i+1)%n]
            tet_face_indices = [tet.triangle(j).index() for j in range(4)]
            k0 = tet_face_indices.index(face_ind0)
            k1 = tet_face_indices.index(face_ind1)
            if k0 == k1:
                assert n == 1
                k1 = tet_face_indices.index(face_ind1, k0 + 1)  ## start looking from index k0 + 1
            edge_pair = vert_pair_to_edge_pair[tuple(sorted((k0, k1)))]
            equ[3*tet.index() + edge_pair] += 1  ### everything is mod 2, so can't tell the difference between + or - anyway
        equ = [x % 2 for x in equ]
        out.append(equ)
    return out


def find_cohomology_trivial_z2_taut_structures(tri):
    structs = find_z2_taut_structures(tri)
    return [s for s in structs if is_trivial_in_cohomology(tri, s)]


def test():
    # t = regina.Triangulation3.fromIsoSig('cPcbbbiht')
    # t = regina.Triangulation3.fromIsoSig('dLQbcccdero')
    # t = regina.SnapPeaTriangulation('/Users/segerman/Dropbox/Schleimer-Segerman/Veering/Finiteness/NotesAndPictures/m019.tri')
    ### t = regina.SnapPeaCensusManifold(regina.SnapPeaCensusManifold.SEC_5, 19).construct() ### seems to be broken
    # t = regina.SnapPeaTriangulation('/Users/segerman/Dropbox/Schleimer-Segerman/Veering/Finiteness/NotesAndPictures/m015.tri')
    # t = regina.SnapPeaTriangulation('/Users/segerman/Dropbox/Schleimer-Segerman/Veering/Finiteness/NotesAndPictures/m011.tri')
    #t = regina.SnapPeaTriangulation('/Users/segerman/Dropbox/Schleimer-Segerman/Veering/Finiteness/NotesAndPictures/m129.tri')
    t = regina.SnapPeaTriangulation('/Users/segerman/Dropbox/Schleimer-Segerman/Veering/Finiteness/NotesAndPictures/m129_v2.tri')
    # t.orient()
    # structs = find_z2_taut_structures(t)
    # for s in structs:
    #     if is_trivial_in_cohomology(t,s):
    #         print(s, is_trivial_in_cohomology(t,s))

    structs = find_cohomology_trivial_z2_taut_structures(t)
    print(structs)
    #for s in structs:
    #    print(s)


if __name__ == '__main__':
    test()
