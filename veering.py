#
# veering.py
#

# Check if a taut triangulation is veering

import regina # needed inside of imported files
from taut import liberal

def tet_type(triangulation, tet_num, veering_colours):
    num_L = [ veering_colours[triangulation.tetrahedron(tet_num).edge(i).index()] for i in range(6) ].count("L")
    if num_L == 2:
        return "right"
    elif num_L == 3:
        return "toggle"
    elif num_L == 4:
        return "left"
    assert False

@liberal
def is_veering(tri, angle, return_type = "boolean"):
    """
    checks to see if this triangulation with taut angle structure is
    veering
    """
    # return type can be "boolean", "veering_colours", "veering_partition" or "num_toggles"
    # all return False if it is not veering. If it is veering:
    # "veering_colours" returns the data of which edges are "L" (bLue) or "R" (Right) veering,
    # "veering_partition" rewrites this as a list of the "L" edges followed by a list of the "R" edges
    # "num_toggles" returns the number of toggle tetrahedra
    # "tet_types" returns list of tets, either "toggle", "left" or "right" veering

    # first do quick test on edge degrees
    for e in tri.edges():
        if e.degree() < 4:
            return False

    # now do proper check
    veering_colours = ["-"] * tri.countEdges()
    for tet_num in range(tri.countTetrahedra()):
        tet = tri.tetrahedron(tet_num)
        edgepair = angle[tet_num]
        # now label edges L or R
        L_edgepair = (edgepair + 1) % 3 # L or R depends on orientation, this takes us from pi edge pair to L edge pair
        L_edge_1, L_edge_2 = L_edgepair, 5 - L_edgepair
        L_edge_1_index = tet.edge(L_edge_1).index()
        L_edge_2_index = tet.edge(L_edge_2).index()
        if veering_colours[L_edge_1_index] == "R" or veering_colours[L_edge_2_index] == "R":
            return False
        else:
            veering_colours[L_edge_1_index] = "L"
            veering_colours[L_edge_2_index] = "L"
        R_edgepair = (edgepair + 2) % 3
        R_edge_1, R_edge_2 = R_edgepair, 5 - R_edgepair
        R_edge_1_index = tet.edge(R_edge_1).index()
        R_edge_2_index = tet.edge(R_edge_2).index()
        if veering_colours[R_edge_1_index] == "L" or veering_colours[R_edge_2_index] == "L":
            return False
        else:
            veering_colours[R_edge_1_index] = "R"
            veering_colours[R_edge_2_index] = "R"
    if return_type == "veering_colours":
        return veering_colours
    elif return_type == "tet_types" or return_type == "num_toggles":
        tet_types = []
        for tet_num in range(tri.countTetrahedra()):
            tet_types.append(tet_type(tri, tet_num, veering_colours))
        if return_type == "tet_types":
            return tet_types
        else:
            return tet_types.count("toggle")
    elif return_type == "veering_partition":
        Llist = []
        Rlist = []
        for i, direction in enumerate(veering_colours):
            if direction == "L":
                Llist.append(i)
            else:
                Rlist.append(i)
        return (Llist, Rlist)
    else:
        return True
