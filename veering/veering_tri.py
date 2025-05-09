#
# veering_tri.py
#

# Provides a checker for veering isosigs and a veering_triangulations
# class.

import regina #needed inside of imported files

from .fundamental_domain import non_tree_face_loops_oriented
from .taut import liberal, is_taut, vert_pair_to_edge_num, edge_num_to_vert_pair
from .transverse_taut import (is_transverse_taut, get_tet_top_and_bottom_edges,
                              top_bottom_embeddings_of_faces, edge_side_face_collections,
                              get_tet_top_vert_nums)


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
    """
    Container class for a triangulation with transverse veering
    structure, possibly with hyperbolic shapes
    """
    def __init__(self, tri, angle, tet_shapes = None, field = None):
        self.tri = tri
        self.angle = angle
        assert is_taut(tri, angle)
        self.coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
        assert self.coorientations != False
        self.veering_colours = is_veering(tri, angle, return_type = "veering_colours")
        assert self.veering_colours != False
        self.tet_types = is_veering(tri, angle, return_type = "tet_types")
        self.tet_shapes = tet_shapes
        self.field = field

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
        from .edge_orientability import is_edge_orientable as is_eo
        return is_eo(self.tri, self.angle, return_type = return_type)


def is_AB_turn(vt, face0, face1, face0_dir, face1_dir):
    """
    For the "turn" (face0, face1, face0_dir, face1_dir) in the veering
    triangulation vt, we decide if it is an "anti-branching" (AB) turn
    as defined at the top of page 16 of arxiv:2008.04836.  In more
    detail: we travel along face0 in direction face0_dir (+1 if with
    the coorientation) into a tet t. We then leave through face1 in
    direction face1_dir.  We return True if this turn is an AB turn:
    the triangles are adjacent along an equatorial edge of t of the
    same colour as the top diagonal of the edge.
    """
    top, bottom = top_bottom_embeddings_of_faces(vt.tri, vt.angle)
    if face0_dir == 1:
        embed0 = bottom[face0]
    else:
        embed0 = top[face0]
    if face1_dir == -1:
        embed1 = bottom[face1]
    else:
        embed1 = top[face1]
    t0 = embed0.tetrahedron()
    t1 = embed1.tetrahedron()
    # print(t0.index(), t1.index())
    assert(t0 == t1)
    t = t0

    if face0_dir != face1_dir:
        return False
    f0 = embed0.face()
    f1 = embed1.face()
    equatorial_nums = [0,1,2,3]
    equatorial_nums.remove(f0)
    equatorial_nums.remove(f1)
    equatorial_colour = vt.get_edge_between_verts_colour(t.index(), equatorial_nums)
    top_vert_nums = get_tet_top_vert_nums(vt.coorientations, t.index())
    top_colour = vt.get_edge_between_verts_colour(t.index(), top_vert_nums)
    return top_colour == equatorial_colour


def loop_twistednesses(tri, angle):
    """
    Returns the (list of the) images of the face generators under the
    edge-orientation homomorphism.  For each face of the triangulation
    we take its canonical loop (using the canonical spanning tree).
    We then compute the image (either +1 or -1) by counting the parity
    of number of AB turns along the loop.  We return these as a list.
    (Note that the homomorphism is trivial on tree faces.)  See
    Proposition 5.7 of arxiv:2101.12162.
    """
    vt = veering_triangulation(tri, angle)
    twistednesses_dict = {}
    oriented_loops, all_signs = non_tree_face_loops_oriented(tri, angle)
    for i in range(len(oriented_loops)):
        loop = oriented_loops[i]
        signs = all_signs[i]
        count = 0
        for j in range(len(loop)):
            f0, f1 = loop[j], loop[(j+1)%len(loop)]
            f0d, f1d = signs[j], signs[(j+1)%len(loop)]
            if is_AB_turn(vt, f0, f1, f0d, f1d):
                count += 1
        twistednesses_dict[loop[0]] = (-1)**(count % 2)  # first in loop is the non-tree face
    for i in range(tri.countTriangles()):
        if i not in twistednesses_dict:
            twistednesses_dict[i] = 1
    
    return [twistednesses_dict[i] for i in range(tri.countTriangles())]
