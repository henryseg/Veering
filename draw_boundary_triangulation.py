
from file_io import parse_data_file
from taut import isosig_to_tri_angle
from veering_triangulation import veering_triangulation

import pyx ### vector graphics 

class vector(list):
  def __add__(self, other):
    return self.__class__(map(lambda x,y: x+y, self, other)) #self.__class__ is vector, unless i am a polynomial or something!

  def __neg__(self):
    return self.__class__(map(lambda x: -x, self))

  def __sub__(self, other):
    return self.__class__(map(lambda x,y: x-y, self, other))

  def __mul__(self, other): #mult by scalar
    return self.__class__(map(lambda x: x*other, self))

  def __rmul__(self, other):
    return (self*other)

anglechoice_face2vert = {(0,0):1, (0,1):0, (0,2):3, (0,3):2,
                         (1,0):2, (1,1):3, (1,2):0, (1,3):1,
                         (2,0):3, (2,1):2, (2,2):1, (2,3):0}

### example usage: next_pi_vertex = anglechoice_face2vert[ (self.angle_structure[next_tet.index()], inf_vert) ]

class ladder_unit:
    """a triangle in a ladder, together with associated data"""

    def __init__(self, vt, tet_face):
        self.vt = vt
        self.tet_face = tet_face
        self.left_vertices = []  # of the ladder, not in terms of veering colours
        self.right_vertices = []
        self.calculate_left_right_vertices()
        self.vertex_posns = {}

    def __str__(self):
        return str(self.tet_face)
    def __repr__(self):
        return str(self.tet_face)

    def calculate_left_right_vertices(self):
        tet_num, inf_vert = self.tet_face
        verts = get_triangle_vertex_order(inf_vert)
        red_vertices = []
        blue_vertices = []
        for v in verts:
            if get_edge_between_verts_colour(self.vt, tet_num, inf_vert, v) == 'R':
                red_vertices.append(v)
            else:
                blue_vertices.append(v)
        ## make pi_vertex first on the side that it is on
        pi_vertex = anglechoice_face2vert[(self.vt.angle[tet_num], inf_vert)]
        if pi_vertex in red_vertices:
            i = red_vertices.index(pi_vertex)
            red_vertices = red_vertices[i:] + red_vertices[:i]
        else:
            i = blue_vertices.index(pi_vertex)
            blue_vertices = blue_vertices[i:] + blue_vertices[:i]
        
        if self.vt.coorientations[tet_num][inf_vert] == -1:  # flip ladder in this case
            self.left_vertices = blue_vertices
            self.right_vertices = red_vertices
            self.left_vertices.reverse()
            self.right_vertices.reverse()  # get the pi vertex above the other one on the same side
        else:
            self.left_vertices = red_vertices
            self.right_vertices = blue_vertices

    def add_vertex_posns(self, posns_dict):
        self.vertex_posns = posns_dict

    def is_on_left(self):
        return len(self.left_vertices) == 2

    def is_on_right(self):
        return len(self.right_vertices) == 2

    def draw_triangle_label(self, my_canvas, ladder_width = None, delta = 0.2):
        if self.is_on_left():
            rungs = [self.draw_triangle_curvy_edge(my_canvas, self.right_vertices[0], vert, ladder_width, delta = delta, draw = False) for vert in self.left_vertices]
        else:
            rungs = [self.draw_triangle_curvy_edge(my_canvas, self.left_vertices[0], vert, ladder_width, delta = delta, draw = False) for vert in self.right_vertices]
        points = [rung.at(0.7*rung.arclen()) for rung in rungs]
        pos = 0.5*(vector(points[0]) + vector(points[1]))

        my_canvas.text(pos[0], pos[1], "$"+str(self.tet_face[0])+"_"+str(self.tet_face[1])+"$", textattrs=[pyx.text.halign.center, pyx.text.vshift.middlezero])

    def draw_corner_symbol(self, my_canvas, vertex, symbol):
        vertex_names = self.vertex_posns.keys()
        posns = [vector(self.vertex_posns[vertex_name]) for vertex_name in vertex_names]
        center = (1.0/3.0)*(posns[0]+posns[1]+posns[2])
        pos = (1.0/3.0)*(center + 2.0*vector(self.vertex_posns[vertex]))
        my_canvas.text(pos[0], pos[1], symbol)

    def draw_vertex_dot(self, my_canvas, vertex):
        x,y = self.vertex_posns[vertex][0], self.vertex_posns[vertex][1]
        colour = get_edge_between_verts_colour(self.vt, self.tet_face[0], self.tet_face[1], vertex)
        draw_vertex_colour(my_canvas, (x,y), colour)

    def draw_vertex_dots(self, my_canvas):
        for v in self.left_vertices + self.right_vertices:
            self.draw_vertex_dot(my_canvas, v)

    def draw_triangle_curvy_edge(self, my_canvas, v0, v1, ladder_width, delta = 0.2, draw=True):
        colours = {'L':pyx.color.rgb.blue, 'R':pyx.color.rgb.red}
        vp0 = self.vertex_posns[v0]
        vp1 = self.vertex_posns[v1]
        veering_dir = get_edge_between_verts_colour(self.vt, self.tet_face[0], v0, v1)
        colour = colours[ veering_dir ]
        if abs(vp0[0] - vp1[0]) < 0.01: ###hack, this is a vertical edge, directions depend on colour
            out_path = pyx.path.line( vp0[0], vp0[1], vp1[0], vp1[1])
            if draw: my_canvas.stroke(out_path, [pyx.deco.stroked([colour])])
            return out_path
        elif (vp0[0] < vp1[0]) != (get_edge_between_verts_colour(self.vt, self.tet_face[0], self.tet_face[1], v0) == 'R'):
            sign = +1
        else: ### if we are at a red dot and going right, we should be concave up
            sign = -1
        shift = sign * delta * ladder_width
        out_path = pyx.path.curve( vp0[0], vp0[1], vp0[0], vp0[1] + shift, vp1[0], vp1[1] + shift, vp1[0], vp1[1])
        if draw: my_canvas.stroke(out_path,  [pyx.deco.stroked([colour])])
        return out_path

    def draw_face_label(self, my_canvas, face_num, ladder_width = None, delta = 0.2):  
        tet = self.vt.tri.tetrahedron( self.tet_face[0] )
        triangle = tet.triangle(face_num)
        triangle_num = triangle.index()

        vertex_names = self.vertex_posns.keys()
        edge_verts = vertex_names[:]
        edge_verts.remove(face_num)

        edge_path = self.draw_triangle_curvy_edge(my_canvas, edge_verts[0], edge_verts[1], ladder_width, delta = 0.2, draw=False)
        pos = edge_path.at(0.5*edge_path.arclen())

        my_canvas.text(pos[0], pos[1], "$"+str(triangle_num)+"$", textattrs=[pyx.text.halign.center, pyx.text.vshift.middlezero, pyx.color.rgb(0,0.5,0)])

    def draw_labels_curvy(self, my_canvas, ladder_width, delta = 0.2):
        vertex_names = self.vertex_posns.keys()
        tet_index, face = self.tet_face
        angle_choice = self.vt.angle[tet_index]
        pi_vert = anglechoice_face2vert[(angle_choice, face)]  
        if pi_vert in self.left_vertices:
            singleton = self.right_vertices[0]
            third = self.left_vertices[ (self.left_vertices.index(pi_vert) + 1) % 2 ]
        else:
            singleton = self.left_vertices[0]
            third = self.right_vertices[ (self.right_vertices.index(pi_vert) + 1) % 2 ]

        posns = [vector(self.vertex_posns[vertex_name]) for vertex_name in vertex_names]
        center = (1.0/3.0)*(posns[0]+posns[1]+posns[2])

        magic_number = 0.6
        for i in range(3):
            vname = vertex_names[i]

            self.draw_face_label(my_canvas, vname, ladder_width = ladder_width, delta=delta)
            
            if vname != pi_vert:
                other_verts = vertex_names[:]
                other_verts.remove(vname)
                paths = [self.draw_triangle_curvy_edge(my_canvas, vname, other_verts[j], ladder_width, delta = delta, draw=False) for j in range(2)]
                # points = [path.at(0.2*path.arclen()) for path in paths]
                if vname == singleton:
                    amount = 0.2
                else:
                    amount = 0.1
                points = [path.at(amount*ladder_width) for path in paths]
                pos = 0.5*(vector(points[0]) + vector(points[1]))
            else:
                if self.vertex_posns[vname][0] > center[0]:
                    sign = -1
                else:
                    sign = +1
                pos = vector( self.vertex_posns[vname] ) + vector([sign*ladder_width * 0.03, 0])

            my_canvas.text(pos[0], pos[1], "$"+str(vertex_names[i])+"$", textattrs=[pyx.text.size(sizename="scriptsize"), pyx.text.halign.center, pyx.text.vshift.middlezero])

    def draw_triangle_edges(self, my_canvas, ladder_width = None, delta = 0.2):
        vertex_names = self.vertex_posns.keys()
        for lv in self.left_vertices:
            for rv in self.right_vertices:
                self.draw_triangle_curvy_edge(my_canvas, lv, rv, ladder_width, delta = delta)
        if self.is_on_left():
            s0,s1 = self.left_vertices
        else:
            s0,s1 = self.right_vertices
        self.draw_triangle_curvy_edge(my_canvas, s0,s1, ladder_width, delta = delta)  

class ladder:
    """ladder of triangles in cusp triangulation of a veering triangulation"""
    def __init__(self, vt, start_tet_face):
        self.ladder_unit_list = []
        self.vt = vt
        self.make_ladder(start_tet_face)
    
    def __str__(self):
        return '[' + ','.join([str(lu) for lu in self.ladder_unit_list]) + ']'

    def __repr__(self):
        return str(self)

    def add_ladder_unit(self, tet_face):
        """add a tet_face to the ladder"""
        self.ladder_unit_list.append(ladder_unit(self.vt, tet_face))

    def count_vertices(self):
        left_length = sum([len(tri.left_vertices) for tri in self.ladder_unit_list]) - len(self.ladder_unit_list)
        right_length = sum([len(tri.right_vertices) for tri in self.ladder_unit_list]) - len(self.ladder_unit_list)
        return left_length, right_length

    def word(self):
        out = []
        for unit in self.ladder_unit_list:
            if unit.is_on_right():
                out.append('l') #left turn up the ladder
            else:
                out.append('r')
        return out
        
    def draw(self, my_canvas, width = 5.0, height = 10.0, originx = 0.0, originy = 0.0, delta = 0.2):
        left_len, right_len = self.count_vertices()
        first_ladder_unit = self.ladder_unit_list[0]
 
        left_pos = 0
        right_pos = 0
        for ladder_unit in self.ladder_unit_list:
            if ladder_unit.is_on_right():
                back_coords = (originx + width, originy + height*float(right_pos)/float(right_len))
                back_label = ladder_unit.right_vertices[0]
                right_pos += 1
            else:
                back_coords = (originx, originy + height*float(left_pos)/float(left_len))
                back_label = ladder_unit.left_vertices[0]
                left_pos += 1
                    
            left_front_coords = (originx, originy + height*float(left_pos)/float(left_len))
            right_front_coords = (originx + width, originy + height*float(right_pos)/float(right_len))

            ladder_unit.add_vertex_posns({back_label:back_coords, ladder_unit.left_vertices[-1]:left_front_coords, ladder_unit.right_vertices[-1]:right_front_coords})

            ladder_unit.draw_triangle_label(my_canvas, ladder_width = width, delta = delta)
            
            ladder_unit.draw_triangle_edges(my_canvas, ladder_width = width, delta = delta)
            ladder_unit.draw_vertex_dots(my_canvas)
            ladder_unit.draw_labels_curvy(my_canvas, width, delta = delta)

    def make_ladder(self, start_tet_face):
        """build a ladder starting from start_tet_face, going along the ladder with the pi vertex as the trailing vertex"""
        
        (tet_index, inf_vert) = start_tet_face

        current_tet = self.vt.tri.tetrahedron(tet_index)
        current_inf_vert = inf_vert
        current_tet_face = start_tet_face  
        while True: 
            self.add_ladder_unit(current_tet_face)        
            current_pi_vertex = anglechoice_face2vert[ (self.vt.angle[current_tet.index()], current_inf_vert) ]
            gluing = current_tet.adjacentGluing(current_pi_vertex)

            current_tet = current_tet.adjacentTetrahedron(current_pi_vertex)
            current_inf_vert = gluing[current_inf_vert]
            current_tet_face = (current_tet.index(),current_inf_vert)
            if current_tet_face == start_tet_face:
                break
        ### if needed, flip ladder
        tet_num, inf_vert = self.ladder_unit_list[0].tet_face
        if self.vt.coorientations[tet_num][inf_vert] == -1: 
            self.ladder_unit_list.reverse()

    def get_rung_holonomy(self):
        """returns vector of form [a,b,c,d,e,f,...] where a,b,c are the powers of z0, 1/(1-z0) and (z0-1)/z0 appearing in the holonomy going up the ladder"""
        out = [0]*3*self.vt.tri.countTetrahedra()
        for unit in self.ladder_unit_list:
            tet, face = unit.tet_face
            if unit.is_on_left(): #we are turning clockwise about right vertex
                sign = -1
                edge = tuple([face, unit.right_vertices[0]].sort())
            else:
                sign = 1
                edge = tuple([face, unit.left_vertices[0]].sort())
            edge_to_angle = {(0,1):0, (2,3):0, (0,2):1, (1,3):1, (0,3):2, (1,2):2}
            out[3*unit.tet + edge_to_angle[edge]] += sign
        return out

def draw_edge_label(my_canvas, edge_endpoint_coords):
    endpts = map(vector, edge_endpoint_coords)
    center = 0.5*(endpts[0] + endpts[1])
    square = pyx.path.rect(center[0] - 0.1, center[1] - 0.1, 0.2, 0.2)
    my_canvas.fill(square, [pyx.deco.filled([pyx.color.rgb.black])])

def draw_vertex_colour(my_canvas, coords, veering_direction):
    colours = {'L':pyx.color.rgb.blue, 'R':pyx.color.rgb.red}
    circ = pyx.path.circle(coords[0], coords[1] ,0.1)
    my_canvas.fill(circ, [pyx.deco.filled([colours[veering_direction]])])

def draw_symmetry_symbol(my_canvas, coords):
    circ = pyx.path.circle(coords[0], coords[1] ,0.18)
    my_canvas.stroke(circ, [pyx.color.rgb.green])

def get_edge_between_verts_colour(vt, tet_num, v1, v2):
    """returns the veering direction (colour) for the given edge of tetrahedron"""
    ### the following dict turns a vert pair into index of edge within a tetrahedron
    vert_pair_to_edge_index = {(0,1):0, (1,0):0, (0,2):1, (2,0):1, (0,3):2, (3,0):2, (1,2):3, (2,1):3, (1,3):4, (3,1):4, (2,3):5, (3,2):5}
    edge_num = vert_pair_to_edge_index[(v1,v2)]
    edge = vt.tri.tetrahedron(tet_num).edge(edge_num)
    return vt.veering_colours[edge.index()]

def get_triangle_vertex_order(face):
    """returns an anticlockwise ordering on the vertices of a triangle 
    consistent with our orientation convention"""
    if face % 2 == 0:
        triangle_verts = [0,1,2,3]
    else:
        triangle_verts = [0,3,2,1]
    triangle_verts.remove(face)
    
    return triangle_verts

def find_unit_index(new_tet_face, my_ladder):
    """find the index of this tet_face in the ladder"""
    unit_index = 0
    ladder_length = len(my_ladder.ladder_unit_list)
    while my_ladder.ladder_unit_list[unit_index].tet_face != new_tet_face:
        unit_index += 1
        if unit_index >= ladder_length:
            print 'error, couldnt find matching unit'
            return None
    return unit_index

class torus_triangulation:
    """list of ladders stacked next to each other"""

    def __init__(self, vt, start_tet_face):
        self.ladder_list = []
        self.vt = vt
        self.make_torus_triangulation(start_tet_face)
        self.tet_faces = []
        for l in self.ladder_list:
            for lu in l.ladder_unit_list:
                self.tet_faces.append(lu.tet_face)
        # print self.tet_faces

    def add_ladder(self, tet_face):
        self.ladder_list.append(ladder(self.vt, tet_face))

    def is_tet_face_in_ladders(self, tet_face):
        """Checks to see if a tet_face is already been used"""
        for L in self.ladder_list:
            for lu in L.ladder_unit_list:
                if lu.tet_face == tet_face:
                    return True
        return False

    def words(self):
        out = []
        for ladder in self.ladder_list:
            out.append(ladder.word())
        return out

    def draw(self, my_canvas, ladder_width = 5.0, height = 10.0, originx = 0.0, originy = 0.0):
        num_ladders = len(self.ladder_list)
        #ladder_width = width / float(num_ladders)
        for i,L in enumerate(self.ladder_list):
            L.draw(my_canvas, width = ladder_width, height = height, originx = originx + ladder_width * i, originy = originy)
        self.draw_symmetries(my_canvas)

    def make_torus_triangulation(self, start_tet_face):
        """build a torus triangulation by building multiple ladders"""
        tet_num, face = start_tet_face

        ### find a loop with canonical horizontal slope
        rho = []
        current_tet_face = start_tet_face
        while current_tet_face not in rho:
            rho.append(current_tet_face)
            tet_num, inf_vert = current_tet_face
            pi_vertex = anglechoice_face2vert[ (self.vt.angle[tet_num], inf_vert) ]
            verts = get_triangle_vertex_order(inf_vert)
            ## find back vertex, is the other end of edge of this triangle of the minority colour from the pi vertex
            verts.remove(pi_vertex)
            u,v = verts

            pi_colour = get_edge_between_verts_colour(self.vt, tet_num, inf_vert, pi_vertex)
            for w in verts:
                if get_edge_between_verts_colour(self.vt, tet_num, inf_vert, w) == pi_colour:
                    trailing_vertex = w
                    verts.remove(w)
                    pivot_vertex = verts[0]
                    break
            leading_vertex = pi_vertex
            ### now forget the pis, we pivot around pivot_vertex until the leading edge (opposite trailing vertex) colour changes:
            start_leading_edge_colour = get_edge_between_verts_colour(self.vt, tet_num, leading_vertex, pivot_vertex)
            assert start_leading_edge_colour == pi_colour # for some reason...
            tet = self.vt.tri.tetrahedron(tet_num)
            while get_edge_between_verts_colour(self.vt, tet.index(), leading_vertex, pivot_vertex) == start_leading_edge_colour:
                gluing = tet.adjacentGluing(trailing_vertex)
                new_tet = tet.adjacentTetrahedron(trailing_vertex)
                inf_vert = gluing[inf_vert]
                pivot_vertex = gluing[pivot_vertex]
                leading_vertex, trailing_vertex = gluing[trailing_vertex], gluing[leading_vertex]
                tet = new_tet
            current_tet_face = (tet.index(), inf_vert)
        sideways = rho[rho.index(current_tet_face):]   ### remove initial tail

        tet_num, inf_vert = sideways[0]
        if self.vt.coorientations[tet_num][inf_vert] == 1: 
            sideways = sideways[1:] + sideways[:1] ## convention: first ladder is convex down, aka green

        ## make sideways go to the right rather than to the left
        tet_num, inf_vert = sideways[0]
        verts = get_triangle_vertex_order(inf_vert)
        cols = [get_edge_between_verts_colour(self.vt, tet_num, inf_vert, v) for v in verts]
        if cols.count('R') == 2:  # pi is on the right, so sideways must be going left
            sideways.reverse()
            sideways = sideways[1:] + sideways[:1] # maintain first ladder convex down

        # print 'len sideways', len(sideways)
        # print 'sideways', sideways
        for tet_face in sideways:
            if self.is_tet_face_in_ladders(tet_face): ### sideways may wrap multiple times around the torus
                break
            self.add_ladder(tet_face)
        # for l in self.ladder_list:
        #     print l

    def ct_start_edges(self):
        out = []
        for i,l in enumerate(self.ladder_list):
            if i%2==0:   ### get both ladder poles from the even ladders to be initially moving the correct way around the ladder pole edges
            ### need input for get_ct_edge_above(current_tet, vert_posns, edge_vertex, face_vertex, old_ct_edge, depth_increment = 1, verbose = 0.0)
                pass
        return out

    def draw_symmetries(self, my_canvas, draw=True):
        count = 0
        dict_of_tet_pairings = {}
        for ladder_index, ladder in enumerate(self.ladder_list):
            for unit_index, unit in enumerate(ladder.ladder_unit_list):
                if unit.is_on_left():
                    edge_axis_unit_pair = self.find_edge_axis_unit_pair(ladder_index, unit_index) #find the pair of triangles to start developing from, checking symmetry
                    vertex_axis_unit_pair = self.find_vertex_axis_unit_pair(ladder_index, unit_index)

                    if self.is_symmetric_torus(edge_axis_unit_pair, dict_of_tet_pairings = dict_of_tet_pairings): #axis is midpoint of the edge to the left of this triangle
                        v0, v1 = vector(unit.vertex_posns[unit.left_vertices[0]]), vector(unit.vertex_posns[unit.left_vertices[1]])
                        draw_symmetry_symbol(my_canvas, 0.5 * (v0+v1))
                        count += 1
                    if self.is_symmetric_torus(vertex_axis_unit_pair, dict_of_tet_pairings = dict_of_tet_pairings): #axis is the back left vertex from this triangle
                        draw_symmetry_symbol(my_canvas, unit.vertex_posns[unit.left_vertices[0]])
                        count += 1
        #check if tet_pairings is 'order 2'
        for tet_pairings in dict_of_tet_pairings.keys():
            used_numbers = []
            for tet_pairing in tet_pairings:
                if tet_pairing[0] == tet_pairing[1]:
                    if tet_pairing[0] in used_numbers:
                        #print 'not order 2!'
                        del dict_of_tet_pairings[tet_pairings]
                        break #break out of for loop
                    used_numbers.append(tet_pairing[0])
                else:
                    if tet_pairing[0] in used_numbers or tet_pairing[1] in used_numbers:
                        #print 'not order 2!'
                        del dict_of_tet_pairings[tet_pairings]
                        break #break out of for loop
                    used_numbers.append(tet_pairing[0])
                    used_numbers.append(tet_pairing[1])
        #print 'total number of symmetries: ', count, ' dict: ', dict_of_tet_pairings
        print '| symmetries, num of: ', dict_of_tet_pairings 

    def find_edge_axis_unit_pair(self, ladder_index, unit_index): #explore the triangulation, see if the midpoint of the edge to the left of this triangle is a 180 deg rotation axis
        unit_posn_A = (ladder_index, unit_index)
        unit_A = self.ladder_list[ladder_index].ladder_unit_list[unit_index]
        tet_B = self.vt.tri.tetrahedron(unit_A.tet_face[0]).adjacentTetrahedron(unit_A.right_vertices[0])
        tet_B_index = tet_B.index()
        gluing = self.vt.tri.tetrahedron(unit_A.tet_face[0]).adjacentGluing(unit_A.right_vertices[0])
        face_B = gluing[unit_A.tet_face[1]]
        ladder_B_index = (ladder_index - 1) % len(self.ladder_list)
        unit_posn_B = (ladder_B_index, find_unit_index((tet_B_index, face_B), self.ladder_list[ladder_B_index]) )
        #print unit_posn_A, unit_A.tet_face, tet_B_index, face_B
        return (unit_posn_A, unit_posn_B)

    def find_vertex_axis_unit_pair(self, ladder_index, unit_index):
        unit_posn_A = (ladder_index, unit_index)
        ladder_A = self.ladder_list[ladder_index]
        unit_A = ladder_A.ladder_unit_list[unit_index]
        #now move down the ladder to find the symmetric triangle on the ladder to the left
        len_A = len(self.ladder_list[ladder_index].ladder_unit_list)
        temp_index = (unit_index - 1) % len_A
        while ladder_A.ladder_unit_list[temp_index].is_on_right():
            temp_index = (temp_index - 1) % len_A
        #have found the unit opposite the unit we want...
        unit_temp = ladder_A.ladder_unit_list[temp_index]

        ####edit next section to change unit_A to unit_temp
        
        tet_B = self.vt.tri.tetrahedron(unit_temp.tet_face[0]).adjacentTetrahedron(unit_temp.right_vertices[0])
        tet_B_index = tet_B.index()
        gluing = self.vt.tri.tetrahedron(unit_temp.tet_face[0]).adjacentGluing(unit_temp.right_vertices[0])
        face_B = gluing[unit_temp.tet_face[1]]
        ladder_B_index = (ladder_index - 1) % len(self.ladder_list)
        unit_posn_B = (ladder_B_index, find_unit_index((tet_B_index, face_B), self.ladder_list[ladder_B_index]) )
        #print unit_posn_A, unit_A.tet_face, tet_B_index, face_B
        return (unit_posn_A, unit_posn_B)
        
    def is_symmetric_torus(self, axis_unit_pair, dict_of_tet_pairings = {}):
        #go right from ladder_A, left from ladder_B until we cover whole torus

        #also record pairings of tetrahedron numbers
        tet_pairings = []
        
        num_ladders = len(self.ladder_list)
        for i in range(num_ladders / 2):
            if not self.is_symmetric_ladder_pair(axis_unit_pair, tet_pairings = tet_pairings):
                return False
            axis_unit_pair = self.find_next_axis_unit_pair(axis_unit_pair)
        #so it is a symmetry point
        #print tet_pairings
        tet_pairings = tuple(tet_pairings) #so hashable, so it can be put in the dict
        if tet_pairings in dict_of_tet_pairings:
            dict_of_tet_pairings[tet_pairings] = dict_of_tet_pairings[tet_pairings] + 1
        else:
            dict_of_tet_pairings[tet_pairings] = 1
        return True

    def is_symmetric_ladder_pair(self, axis_unit_pair, tet_pairings = []):
        # print axis_unit_pair
        unit_posn_A, unit_posn_B = axis_unit_pair
        ladder_A_index, unit_A_index = unit_posn_A
        ladder_B_index, unit_B_index = unit_posn_B
        ladder_A, ladder_B = self.ladder_list[ladder_A_index], self.ladder_list[ladder_B_index]
        len_A, len_B = len(ladder_A.ladder_unit_list), len(ladder_B.ladder_unit_list)
        if len_A != len_B:
            return False
        for j in range(len_A):
            if not ( ladder_A.ladder_unit_list[ (unit_A_index + j) % len_A ].is_on_left() != ladder_B.ladder_unit_list[ (unit_B_index - j) % len_B ].is_on_left() ):
                return False  #this is if both are on left or both are on right, meaning they are not symmetrical
            else: #add this pair to the list of tet_pairings
                tet_pairing = [ladder_A.ladder_unit_list[ (unit_A_index + j) % len_A ].tet_face[0], ladder_B.ladder_unit_list[ (unit_B_index - j) % len_B ].tet_face[0] ]
                tet_pairing.sort()
                tet_pairing = tuple(tet_pairing)
                if tet_pairing not in tet_pairings:
                    tet_pairings.append(tet_pairing)
                    tet_pairings.sort()
        return True
    
    def find_next_axis_unit_pair(self, axis_unit_pair):
        #go left from ladder_A, right from ladder_B
        unit_posn_A, unit_posn_B = axis_unit_pair
        ladder_A_index, unit_A_index = unit_posn_A
        ladder_B_index, unit_B_index = unit_posn_B
        ladder_A, ladder_B = self.ladder_list[ladder_A_index], self.ladder_list[ladder_B_index]
        len_A, len_B = len(ladder_A.ladder_unit_list), len(ladder_B.ladder_unit_list)
        while ladder_A.ladder_unit_list[unit_A_index].is_on_left():
            unit_A_index = (unit_A_index + 1) % len_A
            unit_B_index = (unit_B_index - 1) % len_B
        #now have the tris that we can move through to the next ladders
        unit_A = ladder_A.ladder_unit_list[unit_A_index]
        new_tet_A = self.vt.tri.tetrahedron(unit_A.tet_face[0]).adjacentTetrahedron(unit_A.left_vertices[0])
        new_tet_A_index = new_tet_A.index()
        gluing = self.vt.tri.tetrahedron(unit_A.tet_face[0]).adjacentGluing(unit_A.left_vertices[0])
        new_face_A = gluing[unit_A.tet_face[1]]
        new_ladder_A_index = (ladder_A_index + 1) % len(self.ladder_list)
        new_unit_posn_A = (new_ladder_A_index, find_unit_index((new_tet_A_index, new_face_A), self.ladder_list[new_ladder_A_index]) )

        unit_B = ladder_B.ladder_unit_list[unit_B_index]
        new_tet_B = self.vt.tri.tetrahedron(unit_B.tet_face[0]).adjacentTetrahedron(unit_B.right_vertices[0])
        new_tet_B_index = new_tet_B.index()
        gluing = self.vt.tri.tetrahedron(unit_B.tet_face[0]).adjacentGluing(unit_B.right_vertices[0])
        new_face_B = gluing[unit_B.tet_face[1]]
        new_ladder_B_index = (ladder_B_index - 1) % len(self.ladder_list)
        new_unit_posn_B = (new_ladder_B_index, find_unit_index((new_tet_B_index, new_face_B), self.ladder_list[new_ladder_B_index]) )

        return (new_unit_posn_A, new_unit_posn_B)     

class boundary_triangulation:
    """list of torus_triangulations for all boundary components of the manifold"""

    def __init__(self,vt):
        self.torus_triangulation_list = []
        self.vt = vt
        self.tet_face_list = self.generate_tet_face_list()
        while self.tet_face_list != []:
            start_tet_face = self.tet_face_list[0]
            self.add_torus_triangulation(start_tet_face)
        
    def add_torus_triangulation(self, start_tet_face):
        new_torus_triangulation = torus_triangulation(self.vt, start_tet_face)
        for tet_face in new_torus_triangulation.tet_faces:
            self.tet_face_list.remove(tet_face)
        self.torus_triangulation_list.append(new_torus_triangulation)

    def generate_tet_face_list(self):
        out = []
        for i in range(self.vt.tri.countTetrahedra()):
            for j in range(4):
                out.append((i,j))
        return out

    def words(self):
        out = []
        for torus_triang in self.torus_triangulation_list:
            out.append(torus_triang.words())
        return out

    def draw(self, my_canvas, ladder_width = 5.0, torus_triangulation_height = 10.0, originx = 0.0, originy = 0.0):
        num_torus_triangulations = len(self.torus_triangulation_list)
        for i,T in enumerate(self.torus_triangulation_list):
            print 'cusp:', i, '| num ladders:', len(T.ladder_list),
            T.draw(my_canvas, ladder_width = ladder_width, height = torus_triangulation_height, originx = originx, originy = originy + i * (torus_triangulation_height * 1.1))
            
def generate_boundary_triangulation(tri, angle, output_filename, draw = True):
    """make a picture of the boundary triangulation, save to output_filename. Assumes that filename is of form '*_xxxx.tri' where xxxx is the angle structure for veering, unless is input in angle_structure_str"""
    vt = veering_triangulation(tri, angle)
    B = boundary_triangulation(vt)

    if draw:
        c = pyx.canvas.canvas()
        B.draw(c, ladder_width = 10.0, torus_triangulation_height = 20.0)   #large!
        c.writePDFfile(output_filename)

def draw_triangulation_boundary_from_veering_isosig(veering_isosig, output_filename = None):
    tri, angle = isosig_to_tri_angle(veering_isosig)
    if output_filename == None:
        output_filename = veering_isosig + '.pdf'
    generate_boundary_triangulation(tri, angle, output_filename)

def draw_triangulations_from_veering_isosigs_file(veering_isosigs_filename, output_dirname, num_to_draw = None):
    veering_isosigs_list = parse_data_file(veering_isosigs_filename)
    if num_to_draw != None:
        to_draw = veering_isosigs_list[:num_to_draw]
    else:
        to_draw = veering_isosigs_list
    for veering_isosig in to_draw:
        print veering_isosig
        draw_triangulation_boundary_from_veering_isosig(veering_isosig, output_filename = output_dirname + '/' + veering_isosig + '.pdf')

draw_triangulations_from_veering_isosigs_file('Data/veering_census.txt', 'Images/Boundary_triangulations', num_to_draw = 5)

# draw_triangulation_boundary_from_veering_isosig('cPcbbbiht_12')
# draw_triangulation_boundary_from_veering_isosig('gLLAQbecdfffhhnkqnc_120012')


