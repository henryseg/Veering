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
    checks to see if this triangulation with this taut angle
    structure is veering
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

### Code for getting nice positions of tetrahedra, primarily for drawing midsurfaces

def get_edge_num_between_verts(v0, v1):
    edge_vert_index_map = {(0,1):0, (0,2):1, (0,3):2, (1,2): 3, (1,3):4, (2,3):5 }  ## should import
    edge = [v0,v1]
    edge.sort()
    return edge_vert_index_map[tuple(edge)]

def get_edge_between_verts(tetrahedron, v0, v1):
    edge_num_in_tet = get_edge_num_between_verts(v0, v1)
    edge = tetrahedron.edge(edge_num_in_tet)
    return edge

def get_edge_between_verts_index(tetrahedron, v0, v1):
    edge = get_edge_between_verts(tetrahedron, v0, v1)
    return edge.index()

def get_edge_between_verts_oriented(tetrahedron, v0, v1, edge_orientations_relative_to_regina = None):
    edge_num_in_tet = get_edge_num_between_verts(v0, v1)
    perm = tetrahedron.faceMapping(1, edge_num_in_tet)
    edge_num = get_edge_between_verts_index(tetrahedron, v0, v1)
    if (edge_orientations_relative_to_regina == None) or (edge_orientations_relative_to_regina[edge_num] == 1):
        return (perm[0], perm[1])  ### these are v0, v1 in order given by underlying edge
    else:
        return (perm[1], perm[0]) 

def is_same_orientation_as_regina(tetrahedron, v0, v1, edge_orientations_relative_to_regina = None):
    return get_edge_between_verts_oriented(tetrahedron, v0, v1, edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)[0] == v0

def get_edge_between_verts_colour(veering_colours, tetrahedron, v0, v1):
    """returns the veering direction (colour) for the given edge, between v0, v1"""
    edge_index = get_edge_between_verts_index( tetrahedron, v0, v1)
    return veering_colours[edge_index]

def get_vert_locations(triangulation, tet_num, veering_colours, coorientations):
    coors = coorientations[tet_num]
    top_vertices = []
    bottom_vertices = []
    #coor for a tet looks like [1, -1, 1, -1], is 1 for pointing out of tet, -1 for in
    for i in range(4):
        if coors[i] == 1:
            bottom_vertices.append(i) 
            # Assume coorientation is upwards. If == 1, face coorientation is outwards, 
            # means this face is on top, means corresponding vertex is on bottom
        else:
            top_vertices.append(i) 
    col = get_edge_between_verts_colour(veering_colours, triangulation.tetrahedron(tet_num), top_vertices[0], bottom_vertices[0])
    if col != "red":
        bottom_vertices = [bottom_vertices[1], bottom_vertices[0]]
    return (top_vertices, bottom_vertices)

def rotate_vertices(vert_posns):
    top_verts, bottom_verts = vert_posns
    return ([top_verts[1], top_verts[0]], [bottom_verts[1], bottom_verts[0]])

###           top[0]
###          /   |   \
### bottom[0]--- | ---bottom[1]
###          \   |   /
###           top[1]


def get_edge_to_left_of_tet(triangulation, veering_colours, tet_vert_posns_below, tet_vert_posns_above, tet_num_right, side_of_zigzag = 'upper'):
    # print tet_vert_posns_below 
    # print tet_vert_posns_above
    # print tet_num_right
    tet_right = triangulation.tetrahedron(tet_num_right)  
    if side_of_zigzag == 'upper':          ### which is the edge to the left of tet_right depends on the colour of the annulus we are in
        tet_right_vert_posns = tet_vert_posns_above[tet_num_right]
        edge_of_tet_right = tet_right_vert_posns[1]  ### bottom edge of the tet, is in our annulus
    else:
        tet_right_vert_posns = tet_vert_posns_below[tet_num_right]
        edge_of_tet_right = tet_right_vert_posns[0]  ### top edge of the tet, is in our annulus
    annulus_colour = get_edge_between_verts_colour(veering_colours, tet_right, edge_of_tet_right[0], edge_of_tet_right[1]) 

    if annulus_colour == "red":
        edge_num = get_edge_between_verts_index(tet_right, tet_right_vert_posns[0][0], tet_right_vert_posns[1][0])
        # print 'to left of ', tet_num_right, 'is edge', edge_num, 'red', side_of_zigzag
    else:
        edge_num = get_edge_between_verts_index(tet_right, tet_right_vert_posns[0][1], tet_right_vert_posns[1][0])
        # print 'to left of ', tet_num_right, 'is edge', edge_num, 'blue', side_of_zigzag
    return edge_num

def orientations_agree(triangulation, veering_colours, tet_num_below_edge_num, tet_vert_posns_below, tet_vert_posns_above, tet_num_right):
    #### compare orientation of a tet in lower half of annulus to the tet below the edge to left of original tet
    # print 'orientations agree function'
    edge_num = get_edge_to_left_of_tet(triangulation, veering_colours, tet_vert_posns_below, tet_vert_posns_above, tet_num_right, side_of_zigzag = 'lower') 
    tet_num_below = tet_num_below_edge_num[edge_num]
    ### check for consistency with the tet below an edge against the tet to the right of the edge
    if veering_colours[edge_num] == "red":  
        front_vert_above = tet_vert_posns_below[tet_num_right][1][0]  ## bottom[0]
        back_vert_above = tet_vert_posns_below[tet_num_right][0][0]   ## top[0]
    else:
        front_vert_above = tet_vert_posns_below[tet_num_right][0][1]  ## top[1]
        back_vert_above = tet_vert_posns_below[tet_num_right][1][0]   ## bottom[0]
    front_vert_below = tet_vert_posns_below[tet_num_below][0][1]  # top[1]
    back_vert_below = tet_vert_posns_below[tet_num_below][0][0]   # top[0]

    v0, v1 = get_edge_between_verts_oriented(triangulation.tetrahedron(tet_num_right), front_vert_above, back_vert_above)
    w0, w1 = get_edge_between_verts_oriented(triangulation.tetrahedron(tet_num_below), front_vert_below, back_vert_below)

    return ((v0 == front_vert_above) and (w0 == front_vert_below)) or ((v1 == front_vert_above) and (w1 == front_vert_below))

def find_cut_edges(triangulation, veering_colours, tet_num_below_edge_num, tet_vert_posns_below, tet_vert_posns_above, zigzags):
    cut_edges = []
    for zigzag in zigzags:
        tet_class_below, tet_class_above = zigzag
        for tet_num in tet_class_below:
            if not orientations_agree(triangulation, veering_colours, tet_num_below_edge_num, tet_vert_posns_below, tet_vert_posns_above, tet_num):
                # print 'find cut edges'
                cut_edges.append(get_edge_to_left_of_tet(triangulation, veering_colours, tet_vert_posns_below, tet_vert_posns_above, tet_num, side_of_zigzag = 'lower')) 
    return cut_edges

def get_nice_edge_orientations_relative_to_regina(triangulation, veering_colours, tet_vert_posns_below, tet_vert_posns_above, zigzags):
    ## all edges on bottom side of an annulus should point towards us (but we draw them pointing away??)
    edge_orientations_relative_to_regina = [1] * triangulation.countEdges()
    for zigzag in zigzags:
        tet_class_below, tet_class_above = zigzag
        for tet_num in tet_class_below:
            edge_num = get_edge_to_left_of_tet(triangulation, veering_colours, tet_vert_posns_below, tet_vert_posns_above, tet_num, side_of_zigzag = 'lower')

            top_vertices, bottom_vertices = tet_vert_posns_below[tet_num]
            if veering_colours[edge_num] == "red":
                verts = [top_vertices[0], bottom_vertices[0]] #left
            else:
                verts = [bottom_vertices[0], top_vertices[1]] #left
            if is_same_orientation_as_regina(triangulation.tetrahedron(tet_num), verts[0], verts[1]):
                edge_orientations_relative_to_regina[edge_num] = -1
    return edge_orientations_relative_to_regina

def get_consistent_tet_vert_posns(triangulation, angle, tet_types, coorientations):
    veering_colours = is_veering(triangulation, angle, return_type = 'veering_colours')
    # coorientations = is_transverse_taut(triangulation, angle_structure, return_type = 'tet_vert_coorientations')
    tet_vert_posns_below = []  ## positions of the tetrahedron vertices in the top down view. These determine which way round to draw diamonds
    tet_vert_posns_above = []  ## there are two for tracking the upper and lower pairs of triangles on the tetrahedra
    ### but above is talking about triangles in the top half of a zigzag, so the triangles at the bottom of tetrahedra, and vice versa. Confusing, I know...

    for tet_num in range(triangulation.countTetrahedra()):
        tet_vert_posns_below.append( get_vert_locations(triangulation, tet_num, veering_colours, coorientations) )
        tet_vert_posns_above.append( get_vert_locations(triangulation, tet_num, veering_colours, coorientations) )

    ### get map from edges to the tetrahedra below them
    edge_num_below_tet_num = []
    tet_num_below_edge_num = [''] * triangulation.countTetrahedra()
    for tet_num in range(triangulation.countTetrahedra()):
        top_vertices, bottom_vertices = tet_vert_posns_below[tet_num]
        edge_num = get_edge_between_verts_index(triangulation.tetrahedron(tet_num), top_vertices[0], top_vertices[1])
        edge_num_below_tet_num.append(edge_num)
        tet_num_below_edge_num[edge_num] = tet_num

    ### Next make tet_vert_posns consistent ...
    ### find zigzags, then find cuts. These are edges where we have to twist as we glue

    blue_zigzags = []
    red_zigzags = []
    tet_nums_to_visit = list(range(triangulation.countTetrahedra())) ## these will all be on the bottom of zigzags
    while len(tet_nums_to_visit) > 0: # find all components of all colour complex regions
        ### first time around, assume that the last tet is good, move from there.
        first_tet_num = tet_nums_to_visit.pop() ## start new component
        frontier_tet_nums = [first_tet_num]  ### these are possible tets to start exploring from
        while len(frontier_tet_nums) > 0:
            tet_num = frontier_tet_nums.pop()  # build lower half of zigzag starting from this tet
            edge_class_below = []  
            tet_class_below = []
            edge_class_above = []
            tet_class_above = []
            while True:    # go along the zigzag
                tet = triangulation.tetrahedron(tet_num)
                top_vertices, bottom_vertices = tet_vert_posns_below[tet_num]
                tet_above_right = tet.adjacentSimplex( bottom_vertices[0] )  #adjacent tet to the tet opposite bottom left vert (so in top half of zigzag)
                tet_above_right_from_tet_perm = tet.adjacentGluing( bottom_vertices[0] )

                if get_edge_between_verts_colour(veering_colours, tet, top_vertices[0], top_vertices[1]) == "red": ## top edge of orig tet
                    ### then if consistent orientation, top[0] of orig tet should map to top[0] of tet_above_right, and 
                    ### tet to right of orig tet is through tet_above_right opposite this vertex
                    colour_constant = 0
                else:
                    colour_constant = 1
                
###           top[0]
###          /   |   \
### bottom[0]--- | ---bottom[1]
###          \   |   /
###           top[1]

                edge_class_below.append( get_edge_between_verts_index(tet, bottom_vertices[0], top_vertices[colour_constant]) ) 
                # edge to the left of this tet
                edge_class_above.append( get_edge_between_verts_index(tet, top_vertices[0], top_vertices[1]) )
                # edge to the left of the above right tet
                tet_class_below.append( tet_num )
                tet_class_above.append( tet_above_right.index() )

                if tet_above_right_from_tet_perm[ bottom_vertices[1] ] != tet_vert_posns_above[tet_above_right.index()][1][1]: # bottom_vertices[1] should match whether red or blue
                    ### then orientation is wrong
                    # print 'rotate tet above right:', tet_above_right.index()
                    tet_vert_posns_above[tet_above_right.index()] = rotate_vertices( tet_vert_posns_above[tet_above_right.index()] )

                tet_right_from_tet_above_right_perm = tet_above_right.adjacentGluing( tet_above_right_from_tet_perm[top_vertices[colour_constant]] )
                tet_right = tet_above_right.adjacentSimplex( tet_above_right_from_tet_perm[top_vertices[colour_constant]] )

                if tet_right_from_tet_above_right_perm[tet_above_right_from_tet_perm[ bottom_vertices[1] ]] != tet_vert_posns_below[tet_right.index()][0][colour_constant]:
                    ### then orientation is wrong
                    # print 'rotate tet right:', tet_right.index()
                    tet_vert_posns_below[tet_right.index()] = rotate_vertices( tet_vert_posns_below[tet_right.index()] )

                if tet_right.index() in tet_nums_to_visit:  #we have not looped yet
                    tet_nums_to_visit.remove(tet_right.index())
                    if tet_right.index() in frontier_tet_nums:  # if we hit it while exploring horizontally, we are already building the loop its in, 
                        frontier_tet_nums.remove(tet_right.index()) # so we don't need to explore from there

                    tet_num = tet_right.index() # do this next
                else:
                    break # we have looped horizontally: stop
            zigzag = [tet_class_below, tet_class_above]
            if colour_constant == 0:
                red_zigzags.append(zigzag)
            else:
                blue_zigzags.append(zigzag)

            # now go downwards from the edges in the edge class below
            for tet_num in tet_class_below:
                # print 'go downwards'
                edge_num = get_edge_to_left_of_tet(triangulation, veering_colours, tet_vert_posns_below, tet_vert_posns_above, tet_num, side_of_zigzag = 'lower')
                tet_num_below = tet_num_below_edge_num[edge_num]
                if tet_num_below in tet_nums_to_visit: 
                    if tet_num_below not in frontier_tet_nums:  # just in case, don't want duplicates in the list
                        frontier_tet_nums.append(tet_num)

                    if not orientations_agree(triangulation, veering_colours, tet_num_below_edge_num, tet_vert_posns_below, tet_vert_posns_above, tet_num):
                        # print 'rotate tet below:', tet_num_below
                        tet_vert_posns_below[tet_num_below] = rotate_vertices(tet_vert_posns_below[tet_num_below])
                # else: # we have looped vertically. 
            break

    red_zigzags.sort(key = lambda x: len(x[0]))   # sorted short to long
    blue_zigzags.sort(key = lambda x: len(x[0]))
    zigzags = red_zigzags + blue_zigzags

    while True:
        made_change = False
        for i, zigzag in enumerate(zigzags):
            tet_class_below, tet_class_above = zigzag
            if len([tet_num for tet_num in tet_class_below if tet_types[tet_num] == 'toggle']) == 0:  #all fan tetrahedra
                # print 'all fan tetrahedra'
                if not orientations_agree(triangulation, veering_colours, tet_num_below_edge_num, tet_vert_posns_below, tet_vert_posns_above, tet_class_below[0]):
                    for tet_num in tet_class_below:
                        assert not orientations_agree(triangulation, veering_colours, tet_num_below_edge_num, tet_vert_posns_below, tet_vert_posns_above, tet_num)
                    # print 'and orientations dont match below: rotate this entire annulus'  ### also rotates some other stuff for toggles...
                    for tet_num_above in tet_class_above:
                        tet_vert_posns_above[tet_num_above] = rotate_vertices(tet_vert_posns_above[tet_num_above])
                    for tet_num_below in tet_class_below:
                        tet_vert_posns_below[tet_num_below] = rotate_vertices(tet_vert_posns_below[tet_num_below])
                    tet_class_below.append(tet_class_below.pop(0)) # fence post stuff
                    tet_class_below.reverse() 
                    tet_class_above.reverse()
                    zigzags[i] = [tet_class_below, tet_class_above]
                    made_change = True
                    break # out of for loop
        if made_change == False:  #otherwise, have to keep shifting down until we make a pass through and make no changes
            break

    cut_edges = find_cut_edges(triangulation, veering_colours, tet_num_below_edge_num, tet_vert_posns_below, tet_vert_posns_above, zigzags)
    # print tet_vert_posns_below
    # print zigzags

    return tet_vert_posns_below, tet_vert_posns_above, zigzags, cut_edges
