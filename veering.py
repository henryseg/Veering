#
# veering.py
#

# Check if a taut triangulation is veering

import regina #needed inside of imported files
from transverse_taut import is_transverse_taut, get_tet_top_and_bottom_edges
from taut import liberal, is_taut, vert_pair_to_edge_num, edge_num_to_vert_pair

def tet_type(triangulation, tet_num, veering_colours):
    num_L = [ veering_colours[triangulation.tetrahedron(tet_num).edge(i).index()] for i in range(6) ].count("blue")
    if num_L == 2:
        return "red"
    elif num_L == 3:
        return "toggle"
    elif num_L == 4:
        return "blue"
    assert False

@liberal
def is_veering(tri, angle, return_type = "boolean"):
    """
    checks to see if this triangulation with taut angle structure is
    veering
    """
    # return type can be "boolean", "veering_colours", "veering_partition" or "num_toggles"
    # all return False if it is not veering. If it is veering:
    # "veering_colours" returns the data of which edges are "blue" (bLue - Left veering) or "red" (Right veering) veering,
    # "veering_partition" rewrites this as a list of the "blue" edges followed by a list of the "red" edges
    # "num_toggles" returns the number of toggle tetrahedra
    # "tet_types" returns list of tets, either "toggle", "red" or "blue" veering

    assert is_taut(tri, angle)

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
        if veering_colours[L_edge_1_index] == "red" or veering_colours[L_edge_2_index] == "red":
            return False
        else:
            veering_colours[L_edge_1_index] = "blue"
            veering_colours[L_edge_2_index] = "blue"
        R_edgepair = (edgepair + 2) % 3
        R_edge_1, R_edge_2 = R_edgepair, 5 - R_edgepair
        R_edge_1_index = tet.edge(R_edge_1).index()
        R_edge_2_index = tet.edge(R_edge_2).index()
        if veering_colours[R_edge_1_index] == "blue" or veering_colours[R_edge_2_index] == "blue":
            return False
        else:
            veering_colours[R_edge_1_index] = "red"
            veering_colours[R_edge_2_index] = "red"
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
            if direction == "blue":
                Llist.append(i)
            else:
                Rlist.append(i)
        return (Llist, Rlist)
    else:
        assert return_type == 'boolean'
        return True

class veering_triangulation():
    """Container class for a triangulation with transverse veering structure, possibly with hyperbolic shapes"""
    def __init__(self, tri, angle, tet_shapes = None):
        self.tri = tri
        self.angle = angle
        assert is_taut(tri, angle)
        self.coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
        assert self.coorientations != False
        self.veering_colours = is_veering(tri, angle, return_type = "veering_colours")
        assert self.veering_colours != False
        self.tet_types = is_veering(tri, angle, return_type = "tet_types")
        self.tet_shapes = tet_shapes

    def get_edge_between_verts_index(self, tet_num, verts):
        ### the following dict turns a vert pair into index of edge within a tetrahedron
        edge_num = vert_pair_to_edge_num[tuple(verts)]
        edge = self.tri.tetrahedron(tet_num).edge(edge_num)
        return edge.index()

    def get_edge_between_verts_colour(self, tet_num, verts):
        """returns the veering direction (colour) for the given edge of tetrahedron"""
        return self.veering_colours[self.get_edge_between_verts_index(tet_num, verts)]

    def exotic_angles(self):
        """Given a transverse taut veering triangulation, find the exotic taut structures (as defined in Futer-Gueritaud)"""

        ### From Futer-Gueritaud's leading-trailing deformations, the exotic angles should be derivable even without
        ### a transverse-taut structure, but it seems easier with.

        exotic_upper = []
        exotic_lower = []
        for i, tet in enumerate(self.tri.tetrahedra()):
            top_edge, bottom_edge = get_tet_top_and_bottom_edges(self.coorientations, tet)
            top_colour, bottom_colour = self.veering_colours[top_edge.index()], self.veering_colours[bottom_edge.index()]
            tet_angle = self.angle[i]

            equator = (tet_angle + 1) % 3
            equator_colour = self.veering_colours[ tet.edge(equator).index() ]
            assert self.veering_colours[ tet.edge(5 - equator).index() ] == equator_colour
            if equator_colour == top_colour:
                exotic_upper.append(equator)
            else:
                exotic_upper.append((equator + 1)%3)

            if equator_colour == bottom_colour:
                exotic_lower.append(equator)
            else:
                exotic_lower.append((equator + 1)%3)
        assert is_taut(self.tri, exotic_upper) and is_taut(self.tri, exotic_lower)
        return [exotic_upper, exotic_lower]

    ### imported methods

    def is_edge_orientable(self, return_type = "boolean"):
        from edge_orientability import is_edge_orientable as is_eo
        return is_eo(self.tri, self.angle, return_type = return_type)





