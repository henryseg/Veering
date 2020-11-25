import regina
from copy import deepcopy
from taut import isosig_from_tri_angle
from file_io import parse_data_file, write_data_file
from transverse_taut import is_transverse_taut

def try_add_veering_tet(triangulation, partial_taut_structure, num_pis_at_edges, veering_directions, tet, edgepair, pass_back_structures):
    """try to add the given tet, tet_index to the structure with given edgepair being pi"""
    #first see if num_pis_at_edges is ok
    new_num_pis_at_edges = deepcopy(num_pis_at_edges)
    
    pi_edge_1, pi_edge_2 = edgepair, 5 - edgepair
    edge_1_index = tet.edge(pi_edge_1).index()
    edge_2_index = tet.edge(pi_edge_2).index()
    new_num_pis_at_edges[edge_1_index] += 1
    new_num_pis_at_edges[edge_2_index] += 1
    if new_num_pis_at_edges[edge_1_index] > 2 or new_num_pis_at_edges[edge_2_index] > 2:
        return False
    
    #now label edges L or R
    new_veering_directions = deepcopy(veering_directions)
    
    L_edgepair = (edgepair + 1) % 3 #L or R depends on orientation, some convention here should work
    L_edge_1, L_edge_2 = L_edgepair, 5 - L_edgepair
    L_edge_1_index = tet.edge(L_edge_1).index()
    L_edge_2_index = tet.edge(L_edge_2).index()
    #print 'L edge indices',L_edge_1_index, L_edge_2_index
    if new_veering_directions[L_edge_1_index] == 'R':
        return False
    else:
        new_veering_directions[L_edge_1_index] = 'L'
    if new_veering_directions[L_edge_2_index] == 'R':
        return False
    else:
        new_veering_directions[L_edge_2_index] = 'L'
    R_edgepair = (edgepair + 2) % 3
    R_edge_1, R_edge_2 = R_edgepair, 5 - R_edgepair 
    R_edge_1_index = tet.edge(R_edge_1).index()
    R_edge_2_index = tet.edge(R_edge_2).index()
    #print 'R edge indices',R_edge_1_index, R_edge_2_index
    if new_veering_directions[R_edge_1_index] == 'L':
        return False
    else:
        new_veering_directions[R_edge_1_index] = 'R'
    if new_veering_directions[R_edge_2_index] == 'L':
        return False
    else:
        new_veering_directions[R_edge_2_index] = 'R'
    #survived this far, so go deeper into the tree
    new_partial_taut_structure = deepcopy(partial_taut_structure)
    new_partial_taut_structure.append(edgepair)
    #print new_partial_taut_structure, new_num_pis_at_edges, new_veering_directions
    try_build_veering_struct(triangulation, new_partial_taut_structure, new_num_pis_at_edges, new_veering_directions, pass_back_structures)


def try_build_veering_struct(triangulation, partial_taut_structure, num_pis_at_edges, veering_directions, pass_back_structures):
    """recursive function explores tree of partial taut structures looking for veering. builds taut structure
    by adding tetrahedra in sequentially by tetrahedron number"""
    # print partial_taut_structure, num_pis_at_edges, veering_directions
    if len(partial_taut_structure) == triangulation.countTetrahedra():
        #veering directions must all check out, all num_pis_at_edges must be 2 (by counting argument)
        # print partial_taut_structure, num_pis_at_edges, veering_directions
        pass_back_structures.append(partial_taut_structure)
        return True
    #try to add next tetrahedron in
    tet_index = len(partial_taut_structure)
    tet = triangulation.tetrahedron(tet_index)
    for edgepair in range(3):
        try_add_veering_tet(triangulation, partial_taut_structure, num_pis_at_edges, veering_directions, tet, edgepair, pass_back_structures)

def find_veering_structures(triangulation):
    """try all possibilities of flat tetrahedra, looking for a veering structure on triangulation"""
    # first do quick test on edge degrees
    for e in triangulation.edges():
        if e.degree() < 4:
            return []

    num_tet = triangulation.countTetrahedra()
    num_edges = triangulation.countEdges()
    pass_back_structures = []
    try_build_veering_struct(triangulation, [], [0]*num_tet, ['-']*num_edges, pass_back_structures)
    return pass_back_structures

def find_veering_isosigs_from_regina_isosig(regina_isosig, transverse_taut_only = True):
    t = regina.Triangulation3.fromIsoSig(regina_isosig)
    t.orient()
    veer_structs = find_veering_structures(t)
    if transverse_taut_only:
        veer_structs = [struct for struct in veer_structs if is_transverse_taut(t, struct)]
    veering_isosigs = [isosig_from_tri_angle(t, alpha) for alpha in veer_structs]
    veering_isosigs = list(set(veering_isosigs))
    veering_isosigs.sort()
    return veering_isosigs

def find_veering_isosigs_from_regina_isosig_file(filename):
    data = parse_data_file(filename)
    out = []
    for i, regina_isosig in enumerate(data):
        if i%500 == 0:
            print float(i)/float(len(data))
        out.extend( find_veering_isosigs_from_regina_isosig(regina_isosig) )
    out = [[foo] for foo in out]
    write_data_file(out, filename.split('.')[0] + '_veering_isosigs.txt')

if __name__ == '__main__':
    # t = regina.readSnapPea('m016_canonical_022.tri')
    # t = regina.SnapPeaTriangulation('m016_canonical_022.tri')
    # t = regina.Triangulation3(t)
    # print t.edges()

    regina_isosig = 'iLLLQPcbeefeghhhhhhqhqxxa'
    print find_veering_isosigs_from_regina_isosig(regina_isosig)
    # find_veering_isosigs_from_regina_isosig_file('Veering_census/veering_covers.txt')