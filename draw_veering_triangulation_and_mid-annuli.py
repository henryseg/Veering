import sys
import os

import veering
import transverse_taut
import veering_detect
import taut
import pyx
import math

# pyx.text.set(mode="latex")
# pyx.text.preamble(r"\usepackage{times}")

red = pyx.color.rgb.red
blue = pyx.color.rgb.blue
dark_red = pyx.color.rgb(0.5,0.0,0.0)
dark_blue = pyx.color.rgb(0.0,0.0,0.5)
light_red = pyx.color.rgb(1.0,0.8,0.8)
light_blue = pyx.color.rgb(0.8,0.8,1.0)
direction_to_col = {'blue':blue, 'red':red}

green = pyx.color.rgb(0.0,0.5,0.0)
purple = pyx.color.rgb(0.5,0.0,0.5)
grey = pyx.color.rgb(0.5,0.5,0.5)

tet_width = 5.0
tet_gap = 1.0
lw1 = 0.1  # line width
lw2 = 0.05  # thinner line width for edges in front of and behind diamonds
rw = 0.25  # rectangle width
hlw = 0.2 # hinge cut line width
as1 = 0.7
as2 = 0.35
dot_rad = 0.35  #0.1
# dot = '$\odot$'
# cross = '$\otimes$'
dot = '.'
cross = 'x'

green_purple_offset = 0.06
# green_purple_offset = 0.0
perp_offset = 0.12

### some of these things exist elsewhere in the code - should be factored

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

def get_edge_between_verts_colour(veering_colours, tetrahedron, v0, v1):
    """returns the veering direction (colour) for the given edge, between v0, v1"""
    edge_index = get_edge_between_verts_index( tetrahedron, v0, v1)
    return direction_to_col[ veering_colours[edge_index] ]

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

def draw_edge(canvas, triangulation, veering_colours, tet_num, vert_posns, v0, v1, context = 'tetrahedron', dotted = False, starred = False, edge_orientations_relative_to_regina = None):
    col = get_edge_between_verts_colour(veering_colours, triangulation.tetrahedron(tet_num), v0, v1)
    v0, v1 = get_edge_between_verts_oriented(triangulation.tetrahedron(tet_num), v0, v1, edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
    w0 = vert_posns[v0]
    w1 = vert_posns[v1]

    if context == 'tetrahedron':
        style = [pyx.deco.stroked([col]), pyx.style.linewidth(lw1), pyx.deco.earrow([pyx.deco.filled([col])], size = as1)]
    else:
        style = [pyx.deco.stroked([col]), pyx.style.linewidth(lw2), pyx.deco.earrow([pyx.deco.filled([col])], size = as2)]
        if dotted:
            style.append(pyx.style.linestyle.dotted)

    canvas.stroke(pyx.path.line( w0.real, w0.imag, w1.real, w1.imag), style )
    
    edge_center = 0.5*(w0+w1)
    diff = w1-w0
    perp = diff*1j
    # perp = complex(-diff[1], diff[0])

    if context == 'tetrahedron':
        center = (1.0/4.0)*(vert_posns[0]+vert_posns[1]+vert_posns[2]+vert_posns[3])
        if abs(diff.real) < 0.01:  #edge is vertical
            pos = edge_center + 0.26*diff + 0.06*perp
        elif abs(diff.imag) < 0.01:  #edge is horizontal
            pos = edge_center + 0.3*diff + 0.05*perp
        else:
            pos = edge_center + 0.17*diff + 0.08*perp
    else:
        pos = edge_center + 0.35j  
    edge_string = str(get_edge_between_verts_index(triangulation.tetrahedron(tet_num), v0, v1))
    # if starred:  # the stars died, 2018-09-10 RIP
    #     edge_string = edge_string + '*'
    canvas.text(pos.real, pos.imag, edge_string, textattrs=[pyx.text.size(sizename="normalsize"), pyx.text.halign.center, pyx.text.vshift.middlezero, col])

def draw_vertex_num(canvas, vert_posns, vert_num): 
    center = (1.0/4.0)*(vert_posns[0]+vert_posns[1]+vert_posns[2]+vert_posns[3])
    pos = -0.1 * center + 1.1 * vert_posns[vert_num]
    canvas.text(pos.real, pos.imag, str(vert_num), textattrs=[pyx.text.size(sizename="small"), pyx.text.halign.center, pyx.text.vshift.middlezero])

def draw_triangle_num(canvas, veering_colours, vert_posns, tetrahedron, tri_num): #, offset_down = False):
    tri_verts = range(4)
    tri_verts.pop(tri_num)
    cols = [get_edge_between_verts_colour(veering_colours, tetrahedron, tri_verts[0], tri_verts[1]),
            get_edge_between_verts_colour(veering_colours, tetrahedron, tri_verts[1], tri_verts[2]),
            get_edge_between_verts_colour(veering_colours, tetrahedron, tri_verts[2], tri_verts[0])]
    reds = [c for c in cols if c == red]
    if len(reds) == 2:
        tri_col = dark_red
    else:
        tri_col = dark_blue
    center = (1.0/4.0)*(vert_posns[0]+vert_posns[1]+vert_posns[2]+vert_posns[3])
    pos = (1.0/3.0) * (vert_posns[0]+vert_posns[1]+vert_posns[2]+vert_posns[3] - vert_posns[tri_num]) 
    pos = -1.3 * center + 2.3 * pos
    tri_index = tetrahedron.triangle(tri_num).index()
    canvas.text(pos.real, pos.imag, str(tri_index), textattrs=[pyx.text.size(sizename="normalsize"), pyx.text.halign.center, pyx.text.vshift.middlezero, tri_col])

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
    if col != red:
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

    if annulus_colour == red:
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
    if veering_colours[edge_num] == 'red':  
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

def get_consistent_tet_vert_posns(triangulation, angle, tet_types, coorientations):
    veering_colours = veering.is_veering(triangulation, angle, return_type = 'veering_colours')
    # coorientations = veering.is_transverse_taut(triangulation, angle_structure, return_type = 'tet_vert_coorientations')
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
    tet_nums_to_visit = range(triangulation.countTetrahedra()) ## these will all be on the bottom of zigzags
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

                if get_edge_between_verts_colour(veering_colours, tet, top_vertices[0], top_vertices[1]) == red: ## top edge of orig tet
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

def get_nice_edge_orientations_relative_to_regina(triangulation, veering_colours, tet_vert_posns_below, tet_vert_posns_above, zigzags):
    ## all edges on bottom side of an annulus should point towards us (but we draw them pointing away??)
    edge_orientations_relative_to_regina = [1] * triangulation.countEdges()
    for zigzag in zigzags:
        tet_class_below, tet_class_above = zigzag
        for tet_num in tet_class_below:
            edge_num = get_edge_to_left_of_tet(triangulation, veering_colours, tet_vert_posns_below, tet_vert_posns_above, tet_num, side_of_zigzag = 'lower')

            top_vertices, bottom_vertices = tet_vert_posns_below[tet_num]
            if veering_colours[edge_num] == 'red':
                verts = [top_vertices[0], bottom_vertices[0]] #left
            else:
                verts = [bottom_vertices[0], top_vertices[1]] #left
            if is_same_orientation_as_regina(triangulation.tetrahedron(tet_num), verts[0], verts[1]):
                edge_orientations_relative_to_regina[edge_num] = -1
    return edge_orientations_relative_to_regina

class midannulus:
    def __init__(self, tet_class_below, tet_class_above, top_edge_indices, bottom_edge_indices):
        self.tet_class_below = tet_class_below
        self.tet_class_above = tet_class_above
        self.top_edge_indices = top_edge_indices
        self.bottom_edge_indices = bottom_edge_indices
    def rotate(self, k):
        self.tet_class_below = self.tet_class_below[k:] + self.tet_class_below[:k]
        self.tet_class_above = self.tet_class_above[k:] + self.tet_class_above[:k]
        self.top_edge_indices = self.top_edge_indices[k:] + self.top_edge_indices[:k]
        self.bottom_edge_indices = self.bottom_edge_indices[k:] + self.bottom_edge_indices[:k]
    def flip(self, tet_vert_posns_above, tet_vert_posns_below):
        # print 'flipping:', self.tet_class_above, self.tet_class_below, self.top_edge_indices, self.bottom_edge_indices
        for tet_num_above in self.tet_class_above:
            tet_vert_posns_above[tet_num_above] = rotate_vertices(tet_vert_posns_above[tet_num_above])
        for tet_num_below in self.tet_class_below:
            tet_vert_posns_below[tet_num_below] = rotate_vertices(tet_vert_posns_below[tet_num_below])
        # self.tet_class_below = self.tet_class_below[1:] + self.tet_class_below[:1]  # fence post stuff
        self.tet_class_below.reverse() 
        self.tet_class_above.reverse()
        self.tet_class_above = self.tet_class_above[1:] + self.tet_class_above[:1] # fence post stuff
        self.top_edge_indices.reverse()
        self.bottom_edge_indices.reverse()
        self.bottom_edge_indices = self.bottom_edge_indices[-1:] + self.bottom_edge_indices[:-1] # fence post stuff
        # print 'flip result:', self.tet_class_above, self.tet_class_below, self.top_edge_indices, self.bottom_edge_indices
        

class edge_between_midannuli:
    def __init__(self, index, is_cut):
        self.index = index
        self.midannulus_below = None
        self.midannulus_above = None
        self.is_cut = is_cut

def has_orientable_mid_surfaces(edges_between_midannuli):
    return all([not e.is_cut for e in edges_between_midannuli])

def make_midannuli_and_edges(triangulation, tet_vert_posns_below, tet_vert_posns_above, zigzags, cut_edges):
    ### go through the midannuli looking for remaining twisted gluings. Flip the annuli pictures over when we can improve. Also reorder zigzags so that they print in a good order.
    zigzags = zigzags[:] ## make shallow copy
    zigzags_top_edge_indices = []
    zigzags_bottom_edge_indices = []
    for i, zigzag in enumerate(zigzags):
        top_edge_indices = []
        bottom_edge_indices = []
        tet_class_below, tet_class_above = zigzag
        # print 'zigzag', zigzag
        for j in range(len(tet_class_below)):
            tet_num_below = tet_class_below[j]
            tet_num_above = tet_class_above[j]
            # print tet_num_below, tet_vert_posns_below[j]
            top_edge_indices.append( get_edge_between_verts_index(triangulation.tetrahedron(tet_num_below), tet_vert_posns_below[tet_num_below][0][0], tet_vert_posns_below[tet_num_below][0][1]) )
            bottom_edge_indices.append( get_edge_between_verts_index(triangulation.tetrahedron(tet_num_above), tet_vert_posns_above[tet_num_above][1][0], tet_vert_posns_above[tet_num_above][1][1]) )
        zigzags_top_edge_indices.append( top_edge_indices )
        bottom_edge_indices = bottom_edge_indices[-1:] + bottom_edge_indices[:-1]  # off by one issue
        zigzags_bottom_edge_indices.append( bottom_edge_indices )

    edges_between_midannuli = []
    for i in range(triangulation.countEdges()):
        edges_between_midannuli.append( edge_between_midannuli(i, i in cut_edges) )
    midannuli = []
    for i, zigzag in enumerate(zigzags):
        this_midannulus = midannulus( zigzag[0], zigzag[1], zigzags_top_edge_indices[i], zigzags_bottom_edge_indices[i] )
        midannuli.append( this_midannulus )
        for edge_index in zigzags_top_edge_indices[i]:
            edges_between_midannuli[edge_index].midannulus_below = this_midannulus
        for edge_index in zigzags_bottom_edge_indices[i]:
            edges_between_midannuli[edge_index].midannulus_above = this_midannulus
    return midannuli, edges_between_midannuli

def flip_and_order_midannuli(midannuli, edges_between_midannuli, tet_types, tet_vert_posns_above, tet_vert_posns_below):
    #### find connectivity between zigzags and edges.
    #### go along first zigzag Z, add zigzags above edges above Z, flip when edge is a twist and zigzag on other side is new...

    midannuli_to_do = midannuli[:] # shallow copy
    new_midannuli = [] # reordered list
    midannulus_hopper = [] # next ones to do go in the hopper
    while len(midannuli_to_do) + len(midannulus_hopper) > 0:
        if len(midannulus_hopper) > 0: 
            this_midannulus = midannulus_hopper.pop(0)
        else: ### don't start with a mid-annulus that has only fans on its bottom
            for i in range(len(midannuli_to_do)):
                found_toggle = False
                # print midannuli_to_do[i].tet_class_below
                for tet_index in midannuli_to_do[i].tet_class_below:
                    if tet_types[tet_index] == 'toggle':
                        found_toggle = True
                        break
                if found_toggle:
                    this_midannulus = midannuli_to_do.pop(i)  # new component, this midannulus will not be flipped
                    break
        new_midannuli.append(this_midannulus)
        for ei in this_midannulus.top_edge_indices:
            above_midannulus = edges_between_midannuli[ei].midannulus_above
            if above_midannulus in midannuli_to_do:
                midannuli_to_do.remove(above_midannulus)
                midannulus_hopper.append(above_midannulus)
                # now decide whether or not to flip it: is this edge in cut_edges?
                if edges_between_midannuli[ei].is_cut:
                    above_midannulus.flip(tet_vert_posns_above, tet_vert_posns_below)
                    for ej in above_midannulus.bottom_edge_indices:
                        edges_between_midannuli[ej].is_cut = not edges_between_midannuli[ej].is_cut
                    for ej in above_midannulus.top_edge_indices:
                        edges_between_midannuli[ej].is_cut = not edges_between_midannuli[ej].is_cut
                # rotate it around so that the edge we are glued at is first
                # print ei, above_midannulus.bottom_edge_indices.index(ei)
                offset = above_midannulus.bottom_edge_indices.index(ei)
                # print above_midannulus.bottom_edge_indices, offset
                above_midannulus.rotate(offset)
        ## do nothing if above_midannulus already in the hopper
    # for mid in new_midannuli:
    #     print mid.tet_class_above, mid.tet_class_below, mid.top_edge_indices, mid.bottom_edge_indices

    return new_midannuli
                
def draw_tetrahedra(c, triangulation, veering_colours, tet_vert_posns_below, edge_orientations_relative_to_regina = None, draw_green = True, draw_purple = True):
    for tet_num in range(triangulation.countTetrahedra()):
        tetrahedron = triangulation.tetrahedron(tet_num)

        tet_center = complex(tet_num * (tet_width + tet_gap), 0.0)
        tet_n = tet_center + complex(0.0, tet_width*0.5)
        tet_s = tet_center + complex(0.0, -tet_width*0.5)
        tet_w = tet_center + complex(-tet_width*0.5, 0.0)
        tet_e = tet_center + complex(tet_width*0.5, 0.0)

        ## draw tet number
        c.text(tet_center.real, 0.7*tet_width, str(tet_num), textattrs=[pyx.text.size(sizename="Huge"), pyx.text.halign.center, pyx.text.vshift.middlezero])

        top_vertices, bottom_vertices = tet_vert_posns_below[tet_num]  # just use the below orientations to draw the tetrahedra

        ### assume that tet_n is the first top vertex. All other positions determined by this:
        ### tet_s is second top vertex
        ### tet_w is connected to tet_n by a red edge - make tet_w be bottom_vertices[0] if necessary
        vert_posns = ['','','','']
        vert_posns[top_vertices[0]] = tet_n
        vert_posns[top_vertices[1]] = tet_s
        vert_posns[bottom_vertices[0]] = tet_w
        vert_posns[bottom_vertices[1]] = tet_e

        vert_posns2 = [v + complex(0.0, -(tet_width + tet_gap)) for v in vert_posns] 
        
        ### top triangles
        ## draw top edge
        v0, v1 = top_vertices
        top_colour = get_edge_between_verts_colour(veering_colours, tetrahedron, v0, v1)

        draw_edge(c, triangulation, veering_colours, tet_num, vert_posns, v0, v1, edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
        
        ## draw 0 angle edges, upper figure
        #red
        draw_edge(c, triangulation, veering_colours, tet_num, vert_posns, top_vertices[0], bottom_vertices[0], edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
        draw_edge(c, triangulation, veering_colours, tet_num, vert_posns, top_vertices[1], bottom_vertices[1], edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
        #blue
        draw_edge(c, triangulation, veering_colours, tet_num, vert_posns, top_vertices[0], bottom_vertices[1], edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
        draw_edge(c, triangulation, veering_colours, tet_num, vert_posns, top_vertices[1], bottom_vertices[0], edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)            
        
###           top[0]
###          /   |   \
### bottom[0]--- | ---bottom[1]
###          \   |   /
###           top[1]

        # west triangle then east_triangle
        top_triangle_vertices = [[top_vertices[0], bottom_vertices[0], top_vertices[1]], [top_vertices[1], bottom_vertices[1], top_vertices[0]]]
        #assume top_colour is red, rotate if not... 
        
        for tv in top_triangle_vertices:
            triangle_corners = [vert_posns[vert_num] for vert_num in tv]

            edge_orderings = [1,1,1]   #edges are sw, top, nw
            for i in range(3):
                v0, v1 = tv[(i+1)%3], tv[(i+2)%3]
                w0, w1 = get_edge_between_verts_oriented(tetrahedron, v0, v1, edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
                if w0 != v0:
                    edge_orderings[i] = -1
            if top_colour == blue: # have to rotate [a,b,c] to [c,a,b]
                triangle_corners = triangle_corners[2:] + triangle_corners[:2]
                edge_orderings = edge_orderings[2:] + edge_orderings[:2]
            draw_triangle_green_and_purple(c, triangle_corners, top_colour, edge_orderings, draw_green = draw_green, draw_purple = draw_purple)

        draw_triangle_num(c, veering_colours, vert_posns, tetrahedron, bottom_vertices[0])
        draw_triangle_num(c, veering_colours, vert_posns, tetrahedron, bottom_vertices[1])

        ### bottom triangles
        ## draw bottom edge
        v0, v1 = bottom_vertices
        draw_edge(c, triangulation, veering_colours, tet_num, vert_posns2, v0, v1, edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
        bottom_colour = get_edge_between_verts_colour(veering_colours, tetrahedron, v0, v1)

        ## draw 0 angle edges, lower figure
        #nw red

        draw_edge(c, triangulation, veering_colours, tet_num, vert_posns2, top_vertices[0], bottom_vertices[0], edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
        #se red
        draw_edge(c, triangulation, veering_colours, tet_num, vert_posns2, top_vertices[1], bottom_vertices[1], edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
        #ne blue
        draw_edge(c, triangulation, veering_colours, tet_num, vert_posns2, top_vertices[0], bottom_vertices[1], edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
        #sw blue     
        draw_edge(c, triangulation, veering_colours, tet_num, vert_posns2, top_vertices[1], bottom_vertices[0], edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
        
        # north triangle then south_triangle
        top_triangle_vertices = [[bottom_vertices[0], bottom_vertices[1], top_vertices[0]], [bottom_vertices[1], bottom_vertices[0], top_vertices[1]]]
        #assume bottom_colour is red, rotate if not... 
        
        for j, tv in enumerate(top_triangle_vertices):
            triangle_corners = [vert_posns2[vert_num] for vert_num in tv]

            edge_orderings = [1,1,1]   #edges are sw, top, nw
            for i in range(3):
                v0, v1 = tv[(i+1)%3], tv[(i+2)%3]
                w0, w1 = get_edge_between_verts_oriented(tetrahedron, v0, v1, edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
                if w0 != v0:
                    edge_orderings[i] = -1
            if bottom_colour == blue: # have to rotate [a,b,c] to [b,c,a]
                triangle_corners = triangle_corners[1:] + triangle_corners[:1]
                edge_orderings = edge_orderings[1:] + edge_orderings[:1]
            draw_triangle_green_and_purple(c, triangle_corners, bottom_colour, edge_orderings, draw_green = draw_green, draw_purple = draw_purple)
        
        draw_triangle_num(c, veering_colours, vert_posns2, tetrahedron, top_vertices[0])
        draw_triangle_num(c, veering_colours, vert_posns2, tetrahedron, top_vertices[1])


        for i in range(4):
            draw_vertex_num(c, vert_posns, i)
            draw_vertex_num(c, vert_posns2, i)


def bezier_path_from_4_vecs(v0,v1,v2,v3):
    return pyx.path.curve( v0.real, v0.imag, v1.real, v1.imag, v2.real, v2.imag, v3.real, v3.imag)

def polygonal_path(points, closed = False):
    out = []
    p = points[0]
    out.append(pyx.path.moveto(p[0],p[1]))
    for p in points[1:]:
        out.append(pyx.path.lineto(p[0],p[1]))
    if closed:
        out.append(pyx.path.closepath())
    return pyx.path.path(*out)

def draw_triangle_green_and_purple(c, triangle_corners, majority_colour, edge_orderings, draw_green = True, draw_purple = True):
    ### triangle corners is list of three vector positions, going anticlockwise around the triangle, starting opposite the singleton colour edge
    ### majority_colour says which colour has two out of the three
    ### edge_orderings is a list of three booleans, for whether the (green, purple) order is the same as anticlockwise around the triangle or not

    singleton_edge_direction = triangle_corners[2] - triangle_corners[1]
    # perp = perp_offset*complex(-singleton_edge_direction[1], singleton_edge_direction[0])
    perp = perp_offset*singleton_edge_direction*1j

    mps_green = []
    mps_purple = []
    for i in range(3):
        # back_pt = 0.5-edge_orderings[i]*green_purple_offset)*triangle_corners[(i+1)%3] + (0.5+edge_orderings[i]*green_purple_offset)*triangle_corners[(i+2)%3])
        # mps_green.append( (0.5+edge_orderings[i]*green_purple_offset)*triangle_corners[(i+1)%3] + (0.5-edge_orderings[i]*green_purple_offset)*triangle_corners[(i+2)%3] )
        # mps_purple.append( (0.5-edge_orderings[i]*green_purple_offset)*triangle_corners[(i+1)%3] + (0.5+edge_orderings[i]*green_purple_offset)*triangle_corners[(i+2)%3] )
        mps_green.append( 0.5*triangle_corners[(i+1)%3] + 0.5*triangle_corners[(i+2)%3] )
        mps_purple.append( 0.5*triangle_corners[(i+1)%3] + 0.5*triangle_corners[(i+2)%3] )

    green_mid = 0.5*(mps_green[1]+mps_green[2])
    purple_mid = 0.5*(mps_purple[1]+mps_purple[2])

    straight_path_green = bezier_path_from_4_vecs(mps_green[1], green_mid, green_mid, mps_green[2])
    straight_path_purple = bezier_path_from_4_vecs(mps_purple[1], purple_mid, purple_mid, mps_purple[2])
    if majority_colour == red:
        # step_1_green = (0.5-green_purple_offset)*triangle_corners[2] + (0.5+green_purple_offset)*triangle_corners[1] + perp
        # step_1_purple = (0.5+green_purple_offset)*triangle_corners[2] + (0.5-green_purple_offset)*triangle_corners[1] + perp
        step_1_green = 0.5*triangle_corners[2] + 0.5*triangle_corners[1] + perp
        step_1_purple = 0.5*triangle_corners[2] + 0.5*triangle_corners[1] + perp
        step_2_green = (0.5-green_purple_offset)*triangle_corners[2] + (0.5+green_purple_offset)*triangle_corners[1] + 2*perp
        step_2_purple = (0.5+green_purple_offset)*triangle_corners[2] + (0.5-green_purple_offset)*triangle_corners[1] + 2*perp
        green_end = mps_green[2]
        purple_end = mps_purple[1]
    else:
        step_1_green = 0.5*triangle_corners[2] + 0.5*triangle_corners[1] + perp
        step_1_purple = 0.5*triangle_corners[2] + 0.5*triangle_corners[1] + perp
        step_2_green = (0.5+green_purple_offset)*triangle_corners[2] + (0.5-green_purple_offset)*triangle_corners[1] + 2*perp
        step_2_purple = (0.5-green_purple_offset)*triangle_corners[2] + (0.5+green_purple_offset)*triangle_corners[1] + 2*perp
        green_end = mps_green[1]
        purple_end = mps_purple[2]

    # connector_1_green = bezier_path_from_4_vecs(mps_green[0], mps_green[0] + 0.5*perp, step_1_green - 0.5*perp, step_1_green)
    # connector_2_green = bezier_path_from_4_vecs(step_1_green, step_1_green + 0.5*perp, step_2_green - 0.5*perp, step_2_green)  
    # connector_1_purple = bezier_path_from_4_vecs(mps_purple[0], mps_purple[0] + 0.5*perp, step_1_purple - 0.5*perp, step_1_purple)
    # connector_2_purple = bezier_path_from_4_vecs(step_1_purple, step_1_purple + 0.5*perp, step_2_purple - 0.5*perp, step_2_purple)  
    
    # arc_green = bezier_path_from_4_vecs(step_2_green, step_2_green + 0.5*perp, 1.2*green_mid - 0.2*green_end, 0.2*green_mid + 0.8*green_end )
    # arc_purple = bezier_path_from_4_vecs(step_2_purple, step_2_purple + 0.5*perp, 1.2*purple_mid - 0.2*purple_end, 0.2*purple_mid + 0.8*purple_end)

    arc_green = bezier_path_from_4_vecs(mps_green[0], mps_green[0] + 0.25*perp, 0.8*green_mid + 0.2*green_end, 0.0*green_mid + 1.0*green_end )
    arc_purple = bezier_path_from_4_vecs(mps_purple[0], mps_purple[0] + 0.25*perp, 0.8*purple_mid + 0.2*purple_end, 0.0*purple_mid + 1.0*purple_end)

    draw_grey = draw_green and draw_purple
    if draw_purple:
        # if not draw_grey:
            # c.stroke(straight_path_purple,  [pyx.style.linewidth(lw2), pyx.deco.stroked([purple])])
            # c.stroke(connector_1_purple, [pyx.style.linewidth(lw2), pyx.deco.stroked([purple])])
        # c.stroke(connector_2_purple, [pyx.style.linewidth(lw2), pyx.deco.stroked([purple])])
        c.stroke(arc_purple, [pyx.style.linewidth(lw2), pyx.deco.stroked([purple])])
    if draw_green:
        # if not draw_grey:
        #     c.stroke(straight_path_green,  [pyx.style.linewidth(lw2), pyx.deco.stroked([green])])
        #     c.stroke(connector_1_green, [pyx.style.linewidth(lw2), pyx.deco.stroked([green])])
        # c.stroke(connector_2_green, [pyx.style.linewidth(lw2), pyx.deco.stroked([green])])
        c.stroke(arc_green, [pyx.style.linewidth(lw2), pyx.deco.stroked([green])])
    if draw_grey:
        c.stroke(straight_path_green,  [pyx.style.linewidth(1.5*lw2), pyx.deco.stroked([grey])])
        # c.stroke(connector_1_green, [pyx.style.linewidth(1.5*lw2), pyx.deco.stroked([grey])])

def draw_half_diamond(c, diamond_pts, triangulation, edges_between_midannuli, veering_colours, tet_types, tet_vert_posns, tet_num, is_upper_half = True, edge_orientations_relative_to_regina = None):
    top_vertices, bottom_vertices = tet_vert_posns[tet_num]

    red_vertex_pairs = [[top_vertices[0], top_vertices[1]], #top
                          [top_vertices[0], bottom_vertices[0]], #left
                          [bottom_vertices[1], bottom_vertices[0]], #bottom
                          [bottom_vertices[1], top_vertices[1]], # right
                          [bottom_vertices[0], top_vertices[1]],  # front, read from left to right, is on top left and bottom right of diamond
                          [top_vertices[0], bottom_vertices[1]]   # back, read from left to right, is on top right and bottom left of diamond
                         ]
    red_tri_nums = [ bottom_vertices[1], top_vertices[1], top_vertices[0], bottom_vertices[0] ]
    blue_vertex_pairs = [[top_vertices[0], top_vertices[1]], #top
                          [bottom_vertices[0], top_vertices[1]], #left
                          [bottom_vertices[0], bottom_vertices[1]], #bottom
                          [top_vertices[0], bottom_vertices[1]], # right
                          [bottom_vertices[0], top_vertices[0]], # back, read from left to right, is on top left and bottom right of diamond (slope = 1)
                          [top_vertices[1], bottom_vertices[1]] # front, read from left to right, is on top right and bottom left of diamond (slope = -1)
                         ]
    blue_tri_nums = [ bottom_vertices[1], top_vertices[0], top_vertices[1], bottom_vertices[0] ]

###           top[0]
###          /   |   \
### bottom[0]--- | ---bottom[1]
###          \   |   /
###           top[1]

    diamond_center = 0.25*(diamond_pts[0] + diamond_pts[1] + diamond_pts[2] + diamond_pts[3])
    diamond_midpts = [0.5*(diamond_pts[i] + diamond_pts[(i+1)%4]) for i in range(4)]

    scale_factor = 1.0  ### use 0.8 to move them off
    ### scale diamond_pts and diamond_midpts for labels
    label_diamond_pts = [scale_factor*(p - diamond_center)+diamond_center for p in diamond_pts]
    label_diamond_midpts = [scale_factor*(p - diamond_center)+diamond_center for p in diamond_midpts]

    tetrahedron = triangulation.tetrahedron(tet_num)
    if is_upper_half:
        col = get_edge_between_verts_colour(veering_colours, tetrahedron, top_vertices[0], top_vertices[1])
        corner_indices = [3,0,1]
        tet_label_shift = 0.12*tet_width
    else:
        col = get_edge_between_verts_colour(veering_colours, tetrahedron, bottom_vertices[0], bottom_vertices[1])
        corner_indices = [1,2,3]
        tet_label_shift = -0.12*tet_width
    if col == red:
        light_col = light_red
        dark_col = dark_red
        vertex_pairs = red_vertex_pairs
        tri_nums = red_tri_nums
    else:
        light_col = light_blue
        dark_col = dark_blue
        vertex_pairs = blue_vertex_pairs
        tri_nums = blue_tri_nums

    
    polygon = pyx.box.polygon( [[diamond_pts[i].real, diamond_pts[i].imag] for i in corner_indices] )
    c.fill(polygon.path(), [light_col])

    if tet_types[tet_num] == 'toggle':
        p0, p1 = [0.7*(p - diamond_center) + diamond_center for p in [diamond_pts[1], diamond_pts[3]]]
        c.stroke(pyx.path.line( p0.real, p0.imag, p1.real, p1.imag), [pyx.style.linewidth(hlw), pyx.style.linecap.round] )

    ## front and back edges
    for j in range(2):
        edge_verts = vertex_pairs[j+4]
        if j == 0:
            if is_upper_half:
                diamond_pt_indices = [1,0]
            else: 
                diamond_pt_indices = [2,3]
        else:
            if is_upper_half:
                diamond_pt_indices = [0,3]
            else: 
                diamond_pt_indices = [1,2]
        edge_endpt_posns_pair = [diamond_pts[diamond_pt_indices[k]] for k in range(2)]
        edge_center = 0.5*(edge_endpt_posns_pair[0] + edge_endpt_posns_pair[1])
        edge_endpt_posns_pair = [0.5*p + 0.5*edge_center for p in edge_endpt_posns_pair]
        edge_endpt_posns = ['','','','']
        for k in range(2):    ### hack so we can reuse tetrahedron drawing code
            edge_endpt_posns[edge_verts[k]] = edge_endpt_posns_pair[k]  
        if ((j == 0 and col == blue) or (j==1 and col == red)):  
            dotted = True
        else:
            dotted = False
        draw_edge(c, triangulation, veering_colours, tet_num, edge_endpt_posns, edge_verts[0], edge_verts[1], context = 'diamond', dotted = dotted, edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)

    for j,i in enumerate(corner_indices):
        verts = vertex_pairs[i]
        edge_index = get_edge_between_verts_index(tetrahedron, verts[0], verts[1])
        if edges_between_midannuli[edge_index].is_cut:
            c.stroke(pyx.path.circle(diamond_pts[i].real, diamond_pts[i].imag, dot_rad))
        edge_string = str(edge_index)
        if is_same_orientation_as_regina(tetrahedron, verts[0], verts[1], edge_orientations_relative_to_regina = edge_orientations_relative_to_regina):  
            edge_string = edge_string + dot
        else:
            edge_string = edge_string + cross
        ### edge labels
        c.text(label_diamond_pts[i].real, label_diamond_pts[i].imag, edge_string, textattrs=[pyx.text.size(sizename="large"), pyx.text.halign.center, pyx.text.vshift.middlezero, col])
        
        ## triangle labels
        if j != 2:  ### last one doesn't get a triangle label
            ### triangle labels
            tri_string = tetrahedron.triangle(tri_nums[i]).index()
            c.text(label_diamond_midpts[i].real, label_diamond_midpts[i].imag, tri_string, textattrs=[pyx.text.size(sizename="large"), pyx.text.halign.center, pyx.text.vshift.middlezero, dark_col])

    # draw tet_num label
    c.text(diamond_center.real, diamond_center.imag + tet_label_shift, str(tet_num), textattrs=[pyx.text.size(sizename="Huge"), pyx.text.halign.center, pyx.text.vshift.middlezero])

def draw_midannuli(c, triangulation, midannuli, edges_between_midannuli, veering_colours, tet_types, tet_vert_posns_below, tet_vert_posns_above, edge_orientations_relative_to_regina):
    for j,this_midannulus in enumerate(midannuli):
        tet_class_below, tet_class_above = this_midannulus.tet_class_below, this_midannulus.tet_class_above

        for k in range(len(tet_class_below)): 
            diamond_pts = [complex(k * (tet_width) -tet_width*0.5*math.sin(i*2*math.pi/4), 
                    tet_width*0.5*math.cos(i*2*math.pi/4) + j*(tet_gap + 0.5*tet_width) + tet_width) for i in range(4)]
            diamond_pts_up = [v + complex(0.5*tet_width, 0.5*tet_width) for v in diamond_pts]

            draw_half_diamond(c, diamond_pts, triangulation, edges_between_midannuli, veering_colours, tet_types, tet_vert_posns_below, tet_class_below[k], is_upper_half = True, edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)
            draw_half_diamond(c, diamond_pts_up, triangulation, edges_between_midannuli, veering_colours, tet_types, tet_vert_posns_above, tet_class_above[k], is_upper_half = False, edge_orientations_relative_to_regina = edge_orientations_relative_to_regina)


def draw_triangulation(triangulation, midannuli_filename, tetrahedra_filename, angle = None, draw_stuff = True, draw_green = True, draw_purple = True):
    """make a picture of the triangulation, save to output_filename. Assumes that filename is of form '*_xxxx.tri' where xxxx is the angle structure for veering, unless is input in angle_structure_str"""
    # print 'drawing:', midannuli_filename
    if angle == None:
        veering_structures = veering_detect.find_veering_structures(triangulation)
        assert len(veering_structures) > 0
        if len(veering_structures) > 1:
            print 'warning, not drawing all veering structures!'
        angle_structure = veering_structures[0]
        print 'angle structure:', str(angle)

    veering_colours = veering.is_veering(triangulation, angle, return_type = 'veering_colours')
    assert veering_colours != False

    coorientations = transverse_taut.is_transverse_taut(triangulation, angle, return_type = 'tet_vert_coorientations')
    if coorientations == False:
        # print 'not transverse taut! Not drawing triangulation' 
        return False
    # coorientations = [[-c for c in t] for t in coorientations] # flip all

    tet_types = veering.is_veering(triangulation, angle, return_type = 'tet_types')
    tet_vert_posns_below, tet_vert_posns_above, zigzags, cut_edges = get_consistent_tet_vert_posns(triangulation, angle, tet_types, coorientations)
    edge_orientations_relative_to_regina = get_nice_edge_orientations_relative_to_regina(triangulation, veering_colours, tet_vert_posns_below, tet_vert_posns_above, zigzags)

    midannuli, edges_between_midannuli = make_midannuli_and_edges(triangulation, tet_vert_posns_below, tet_vert_posns_above, zigzags, cut_edges)
    midannuli = flip_and_order_midannuli(midannuli, edges_between_midannuli, tet_types, tet_vert_posns_above, tet_vert_posns_below) ### adds flips to midannuli, also reorders midannuli list
   
    if draw_stuff:
        ct = pyx.canvas.canvas() 
        draw_tetrahedra(ct, triangulation, veering_colours, tet_vert_posns_below, edge_orientations_relative_to_regina = edge_orientations_relative_to_regina, draw_green = draw_green, draw_purple = draw_purple)
        ct.writePDFfile(tetrahedra_filename)

        cm = pyx.canvas.canvas() 
        draw_midannuli(cm, triangulation, midannuli, edges_between_midannuli, veering_colours, tet_types, tet_vert_posns_below, tet_vert_posns_above, edge_orientations_relative_to_regina)
        cm.writePDFfile(midannuli_filename)

    return has_orientable_mid_surfaces(edges_between_midannuli)

def draw_triangulation_from_veering_isosig(veering_isosig, midannuli_filename = None, tetrahedra_filename = None, draw_stuff = True):
    tri, angle = taut.isosig_to_tri_angle(veering_isosig)
    if midannuli_filename == None:
        midannuli_filename = 'Images/Mid-annuli/' + veering_isosig + '_mid-annuli.pdf'
    if tetrahedra_filename == None:
        tetrahedra_filename = 'Images/Triangulations/' + veering_isosig + '_tetrahedra.pdf'
    return draw_triangulation(tri, midannuli_filename, tetrahedra_filename, angle = angle, draw_stuff = draw_stuff, draw_green = True, draw_purple = True)

def draw_triangulations_from_veering_isosigs_file(veering_isosigs_filename, midannuli_dirname, tetrahedra_dirname):
    veering_isosigs_list = veering_isosigs.parse_data_file(veering_isosigs_filename)
    for veering_isosig in veering_isosigs_list:
        draw_triangulation_from_veering_isosig(veering_isosig, 
            midannuli_filename = midannuli_dirname + '/' + veering_isosig + '_mid-annuli.pdf',
            tetrahedra_filename = tetrahedra_dirname + '/' + veering_isosig + '_tetrahedra.pdf')

def calculate_orientable_mid_surfaces_from_veering_isosigs_file(veering_isosigs_filename, file_out):
    veering_isosigs_list = veering_isosigs.parse_data_file(veering_isosigs_filename)
    out_data = []
    for i, veering_isosig in enumerate(veering_isosigs_list):
        if i%1000 == 0:
            print float(i)/float(len(veering_isosigs_list))
        orientable_mid_surfaces = draw_triangulation_from_veering_isosig(veering_isosig, draw_stuff = False)
        out_data.append(veering_isosig + ' ' + str(orientable_mid_surfaces))
    file_out = open('Veering_census/' + file_out, 'w')
    for line in out_data:
        file_out.write(line + '\n')
    file_out.close()

if __name__ == '__main__':
    # draw_triangulations_from_veering_isosigs_file('Veering_census/veering_census_up_to_12.txt', '../Veering_mid-annuli', '../Veering_tetrahedra')

    # calculate_well_framed_from_veering_isosigs_file('Veering_census/veering_census.txt', 'veering_census_with_well_framed.txt')
    # calculate_well_framed_from_veering_isosigs_file('Veering_census/veering_euler.txt', 'veering_euler_with_well_framed.txt')
    # calculate_well_framed_from_veering_isosigs_file('Veering_census/veering_covers_veering_isosigs.txt', 'veering_covers_with_well_framed.txt')


    # draw_triangulation_from_veering_isosig('dLQacccjsnk_200') 
    # draw_triangulation_from_veering_isosig('iLMzMPcbcdefghhhhhqqqqhxq_12200221')
    # draw_triangulation_from_veering_isosig('lLLAvPAQcbcdejhihkjjktsfxhxhqsjfw_20102120121')
    # draw_triangulation_from_veering_isosig('lLLvLQAQcbefgigijkikkxxxgbrgrqraq_01110220202')
    # draw_triangulation_from_veering_isosig('lLLvLQPQcbeffhgjikjkkxxxfjssfxfhq_10221010102')
    # draw_triangulation_from_veering_isosig('mvvLPQwQQfjfihgfkllklkhahaaaaqqxxaa_210012222211')
   
    # draw_triangulation_from_veering_isosig('gLLAQcdecfffhsermws_122201') 
   
    # draw_triangulation_from_veering_isosig('iLLLQPcbeefeghhhhhhqhqxxa_21112001')
    # draw_triangulation_from_veering_isosig('iLLLQPcbeefeghhhhhhqhqxxa_01110221')

    # draw_triangulations_from_veering_isosigs_file('Gadgets/double_gadget_veering_isosigs.txt', 'Gadgets/Veering_mid-annuli', 'Gadgets/Veering_tetrahedra' )

    # draw_triangulations_from_veering_isosigs_file('iLLLQPcbeefeghhhhhhqhqxxa_veer_isosigs.txt', 'iLLLQPcbeefeghhhhhhqhqxxa_mid-annuli', 'iLLLQPcbeefeghhhhhhqhqxxa_veering_tetrahedra' )

    # draw_triangulations_from_veering_isosigs_file('m016_surgeries.txt', 'm016_surgeries_mid-annuli', 'm016_surgeries_veering_tetrahedra' )



    # sigs = ['ivvPQQcfhghgfghfaaaaaaaaa_01122000', 
    #         'mvvLPQwQQfghffjikllklkaaaaaaaaaaaaa_102021111100', 
    #         'mvvLPQwQQhjhihgflkkllkaaaaaaaaaaaaa_022220001122',
    #         'oLvLAvwQQQccfghkmnkjklnmmnhialaatiqqqffff_20102220210100',
    #         'oLvLAwzPQQccfghhilkjnmmnmnhialgoiqqffffff_20102220200011']

    flat_toggles = ['qLLvAvAMQLQkbeehklmnjnnppopooxxxahahxxxaxqxxxq_2111200221111100',
                    'qLLvAvPMQLQkbeehlkmnjnnpoopopxxxaaxxxxxaxaxxax_2111200221111100',
                    'qvvLPQMvQLQkfgfhhgfknlmoppopohahhaaahaqqaqqaaa_1222211100222200']
    for sig in flat_toggles:
        draw_triangulation_from_veering_isosig(sig)

    # draw_triangulation_from_veering_isosig('mvvLPQwQQfghffjikllklkaaaaaaaaaaaaa_102021111100')

    # draw_triangulation_from_veering_isosig('cPcbbbiht_12')

    # from file_io import parse_data_file
    # # census = parse_data_file('Data/veering_census.txt')

    # # for sig in census[:20]:
    # #     draw_triangulation_from_veering_isosig(sig)

    # # vol_inc = parse_data_file('Data/volume_rise.txt')
    # # for line in vol_inc:
    # #     sig, tet_num = line.split(" ")
    # #     draw_triangulation_from_veering_isosig(sig, midannuli_filename = 'Images/Mid-annuli/Volume_rise/' + sig + '_tet_' + tet_num +'_mid-annuli.pdf')

    # data = parse_data_file('Data/low_volume_stacks.txt')
    # for line in data:
    #     split_data = line.split(" ")
    #     sig, vol = split_data[:2]
    #     stack = "".join(split_data[2:])
    #     draw_triangulation_from_veering_isosig(sig, midannuli_filename = 'Images/Mid-annuli/Low_volume_stacks/' + sig + stack +'_mid-annuli.pdf')




