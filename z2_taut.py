#
# taut.py
#

# functions for working with Z2 taut ideal tris.  

import regina
from copy import deepcopy
from taut import liberal, edge_pair_to_edge_numbers

# check Z2 tautness


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

def test():
    t = regina.Triangulation3.fromIsoSig('cPcbbbiht')
    # t = regina.Triangulation3.fromIsoSig('dLQbcccdero')
    t.orient()
    structs = find_z2_taut_structures(t)
    for s in structs:
        print(s)
        print(is_z2_taut(t, s))

if __name__ == '__main__':
    test()
