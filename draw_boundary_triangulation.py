
from file_io import parse_data_file, read_from_pickle
from taut import isosig_to_tri_angle
from veering_triangulation import veering_triangulation
from develop_ideal_hyperbolic_tetrahedra import developed_position, develop_verts_CP1, unknown_vert_to_known_verts_ordering, convert_to_complex
from veering_cannon_thurston import ct_edge, get_ct_edge_above, develop_cannon_thurston

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

class tet_face:
    def __init__(self, tet_num, face, verts_CP1 = None):
        self.tet_num = tet_num
        self.face = face
        self.verts_CP1 = verts_CP1
        self.origin_in_C = complex(0,0)  ### used to move things around for drawing

    def __eq__(self, other): ## don't care about verts_CP1
        return self.tet_num == other.tet_num and self.face == other.face 

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return '(' + str(self.tet_num) + ',' + str(self.face) + ')'

class ladder_unit(tet_face):
    """a triangle in a ladder, together with associated data"""

    def __init__(self, vt, tf):
        self.vt = vt
        # self.tet_face = tf
        tet_face.__init__(self, tf.tet_num, tf.face, verts_CP1 = tf.verts_CP1)
        self.left_vertices = []  # of the ladder, not in terms of veering colours
        self.right_vertices = []
        self.calculate_left_right_vertices()
        self.verts_C = {}

    # def __str__(self):
    #     return str(self.tet_face)
    # def __repr__(self):
    #     return str(self.tet_face)

    def calculate_left_right_vertices(self):
        tet_num, inf_vert = self.tet_num, self.face
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

    def add_verts_C(self, posns_dict):
        self.verts_C = posns_dict

    def is_on_left(self):
        return len(self.left_vertices) == 2

    def is_on_right(self):
        return len(self.right_vertices) == 2

    def draw_triangle_label(self, my_canvas, curvy = False, args = {}):
        #print self.verts_C
        delta = args['delta']
        if not curvy:
            vertex_names = self.verts_C.keys()
            posns = [self.verts_C[vertex_name] for vertex_name in vertex_names]
            pos = (1.0/3.0)*(posns[0]+posns[1]+posns[2])
            pos = [pos.real, pos.imag]
        else:
            ladder_width = args['ladder_width']
            if self.is_on_left():
                rungs = [self.draw_triangle_curvy_edge(my_canvas, self.right_vertices[0], vert, ladder_width, delta = delta, draw = False) for vert in self.left_vertices]
            else:
                rungs = [self.draw_triangle_curvy_edge(my_canvas, self.left_vertices[0], vert, ladder_width, delta = delta, draw = False) for vert in self.right_vertices]
            points = [rung.at(0.7*rung.arclen()) for rung in rungs]
            pos = [0.5*(points[0][0] + points[1][0]), 0.5*(points[0][1] + points[1][1])]

        my_canvas.text(pos[0], pos[1], "$"+str(self.tet_num)+"_"+str(self.face)+"$", textattrs=[pyx.text.halign.center, pyx.text.vshift.middlezero])

    def draw_corner_and_face_labels(self, my_canvas):
        vertex_names = self.verts_C.keys()
        posns = [self.verts_C[vertex_name] for vertex_name in vertex_names]
        center = (1.0/3.0)*(posns[0]+posns[1]+posns[2])
        for i in range(3):
            self.draw_face_label(my_canvas, vertex_names[i], curvy = False)
            pos = (1.0/3.0)*(center + 2.0*posns[i])
            my_canvas.text(pos.real, pos.imag, "$"+str(vertex_names[i])+"$", textattrs=[pyx.text.size(sizename="scriptsize"), pyx.text.halign.center, pyx.text.vshift.middlezero])

    def draw_vertex_dot(self, my_canvas, vertex):
        x,y = self.verts_C[vertex].real, self.verts_C[vertex].imag
        colour = get_edge_between_verts_colour(self.vt, self.tet_num, self.face, vertex)
        draw_vertex_colour(my_canvas, (x,y), colour)

    def draw_vertex_dots(self, my_canvas):
        for v in self.left_vertices + self.right_vertices:
            self.draw_vertex_dot(my_canvas, v)

    def draw_triangle_curvy_edge(self, my_canvas, v0, v1, ladder_width, delta = 0.2, draw=True):
        colours = {'L':pyx.color.rgb.blue, 'R':pyx.color.rgb.red}
        vp0 = self.verts_C[v0]
        vp1 = self.verts_C[v1]
        colour = colours[ get_edge_between_verts_colour(self.vt, self.tet_num, v0, v1) ]
        if abs(vp0.real - vp1.real) < 0.01: ###hack, this is a vertical edge, directions depend on colour
            out_path = pyx.path.line( vp0.real, vp0.imag, vp1.real, vp1.imag)
            if draw: my_canvas.stroke(out_path, [pyx.deco.stroked([colour])])
            return out_path
        elif (vp0.real < vp1.real) != (get_edge_between_verts_colour(self.vt, self.tet_num, self.face, v0) == 'R'):
            sign = +1
        else: ### if we are at a red dot and going right, we should be concave up
            sign = -1
        shift = sign * delta * ladder_width
        out_path = pyx.path.curve( vp0.real, vp0.imag, vp0.real, vp0.imag + shift, vp1.real, vp1.imag + shift, vp1.real, vp1.imag)
        if draw: my_canvas.stroke(out_path,  [pyx.deco.stroked([colour])])
        return out_path

    def draw_triangle_edge(self, my_canvas, v0, v1):
        colours = {'L':pyx.color.rgb.blue, 'R':pyx.color.rgb.red}
        vp0 = self.verts_C[v0]
        vp1 = self.verts_C[v1]
        colour = colours[ get_edge_between_verts_colour(self.vt, self.tet_num, v0, v1) ]
        my_canvas.stroke(pyx.path.line( vp0.real, vp0.imag, vp1.real, vp1.imag),  [pyx.deco.stroked([colour])]  )

    def draw_face_label(self, my_canvas, face_num, ladder_width = None, curvy = True, delta = 0.2):  
        tet = self.vt.tri.tetrahedron( self.tet_num )
        triangle = tet.triangle(face_num)
        triangle_num = triangle.index()

        vertex_names = self.verts_C.keys()
        edge_verts = vertex_names[:]
        edge_verts.remove(face_num)

        if curvy:
            edge_path = self.draw_triangle_curvy_edge(my_canvas, edge_verts[0], edge_verts[1], ladder_width, delta = 0.2, draw=False)
            pos = edge_path.at(0.5*edge_path.arclen())
        else:
            pos_C = 0.5*(self.verts_C[edge_verts[0]] + self.verts_C[edge_verts[1]])
            pos = [pos_C.real, pos_C.imag]
        my_canvas.text(pos[0], pos[1], "$"+str(triangle_num)+"$", textattrs=[pyx.text.halign.center, pyx.text.vshift.middlezero, pyx.color.rgb(0,0.5,0)])

    def draw_labels_curvy(self, my_canvas, ladder_width, delta = 0.2):
        vertex_names = self.verts_C.keys()
        tet_index, face = self.tet_num, self.face
        angle_choice = self.vt.angle[tet_index]
        pi_vert = anglechoice_face2vert[(angle_choice, face)]  
        if pi_vert in self.left_vertices:
            singleton = self.right_vertices[0]
            third = self.left_vertices[ (self.left_vertices.index(pi_vert) + 1) % 2 ]
        else:
            singleton = self.left_vertices[0]
            third = self.right_vertices[ (self.right_vertices.index(pi_vert) + 1) % 2 ]

        posns = [self.verts_C[vertex_name] for vertex_name in vertex_names]
        center = (1.0/3.0)*(posns[0]+posns[1]+posns[2])

        magic_number = 0.6
        for i in range(3):
            vname = vertex_names[i]

            self.draw_face_label(my_canvas, vname, ladder_width = ladder_width, delta = delta)
            
            if vname != pi_vert:
                other_verts = vertex_names[:]
                other_verts.remove(vname)
                paths = [self.draw_triangle_curvy_edge(my_canvas, vname, other_verts[j], ladder_width, delta = delta, draw=False) for j in range(2)]
                if vname == singleton:
                    amount = 0.2
                else:
                    amount = 0.1
                points = [path.at(amount*ladder_width) for path in paths]
                pos = [0.5*(points[0][0] + points[1][0]), 0.5*(points[0][1] + points[1][1])]
            else:
                if self.verts_C[vname].real > center.real:
                    sign = -1
                else:
                    sign = +1
                pos = self.verts_C[vname] + complex(sign*ladder_width * 0.03, 0)
                pos = [pos.real, pos.imag]

            my_canvas.text(pos[0], pos[1], "$"+str(vertex_names[i])+"$", textattrs=[pyx.text.size(sizename="scriptsize"), pyx.text.halign.center, pyx.text.vshift.middlezero])

    def draw_triangle_edges(self, my_canvas, curvy = True, args = {}):
        #ladder_width = None, delta = 0.2):
        vertex_names = self.verts_C.keys()
        for lv in self.left_vertices:
            for rv in self.right_vertices:
                if curvy:
                    self.draw_triangle_curvy_edge(my_canvas, lv, rv, args['ladder_width'], delta = args['delta'])
                else:
                    self.draw_triangle_edge(my_canvas, lv, rv)
        if self.is_on_left():
            s0,s1 = self.left_vertices
        else:
            s0,s1 = self.right_vertices
        if curvy:
            self.draw_triangle_curvy_edge(my_canvas, s0, s1, args['ladder_width'], delta = args['delta'])  
        else:
            self.draw_triangle_edge(my_canvas, s0, s1)

    def get_ct_edge(self, ladder_is_even):
        # assert self.is_on_left()
        # start_ct_edge = ct_edge(self.tet_num, None, self.verts_CP1, 0, self.vt)
        ### when on an even ladder (concave down in our pictures) the edge_vertex is the inf vertex, the face_vertex is the singleton
        ### vice versa for odd ladders  
        if self.is_on_left():
            singleton = self.right_vertices[0]
        else:
            singleton = self.left_vertices[0]
        if ladder_is_even:
            edge_vertex = self.face # inf vertex
            face_vertex = singleton
        else:
            edge_vertex = singleton
            face_vertex = self.face # inf vertex
        self.ct_edge = get_ct_edge_above(self.vt.tri.tetrahedron(self.tet_num), self.verts_CP1, edge_vertex, face_vertex, self.vt, 0, depth_increment = 0, verbose = 0.0)

    def generate_ct(self, ladder_is_even = True, args = {}):
        #max_depth = 1, epsilon = 0.02, verbose = 0.0):
        self.get_ct_edge(ladder_is_even)
        self.ct_developed_edges = develop_cannon_thurston([self.ct_edge], max_depth = args['ct_depth'], epsilon = args['ct_epsilon'], verbose = 0.1)

    # def draw_ct(self, canv, origin, geom_complex_scale, veering_colour = None, lw = 0.005):
    """draw edges one by one"""
    #     for edge in self.ct_developed_edges:
    #         s, e = edge.start_complex, edge.end_complex
    #         s = geom_complex_scale * ( s + origin ) 
    #         e = geom_complex_scale * ( e + origin ) 
    #         canv.stroke(pyx.path.line(s.real, s.imag, e.real, e.imag), [pyx.style.linewidth(lw), colour])

    def draw_ct_path(self, canv, origin, geom_complex_scale, veering_colour = None, lw = 0.005):
        """draw path of edges"""
        # if veering_colour == 'R':
        #     colour = pyx.color.rgb(0.5,0.0,0.0)  # dark red
        #     grad = pyx.color.gradient.RedBlack
        # else:
        #     colour = pyx.color.rgb(0.0,0.0,0.5)  # dark blue
        #     grad = pyx.color.gradient.BlueBlack
        grad = pyx.color.gradient.Hue
        edges = self.ct_developed_edges
        comp_coords = [edges[0].start_complex] + [edge.end_complex for edge in edges]
        comp_coords = [geom_complex_scale * (c + origin) for c in comp_coords]
        p = pyx.path.path( pyx.path.moveto(comp_coords[0].real, comp_coords[0].imag) )
        for coord in comp_coords[1:]: 
          p.append( pyx.path.lineto(coord.real, coord.imag) )
        # canv.stroke(p, [pyx.style.linewidth(lw), colour])
        canv.stroke(p, [pyx.style.linewidth(lw), pyx.deco.colorgradient(grad)])


class ladder:
    """ladder of triangles in cusp triangulation of a veering triangulation"""
    def __init__(self, torus_triang, start_tf):
        self.torus_triang = torus_triang
        self.ladder_unit_list = []
        self.vt = torus_triang.vt
        self.holonomy = None
        self.is_even = None
        self.ladder_origin = start_tf.origin_in_C
        self.make_ladder(start_tf)

    def __str__(self):
        return '[' + ','.join([str(lu) for lu in self.ladder_unit_list]) + ']'

    def __repr__(self):
        return str(self)

    def count_vertices(self):
        left_length = sum([len(tri.left_vertices) for tri in self.ladder_unit_list]) - len(self.ladder_unit_list)
        right_length = sum([len(tri.right_vertices) for tri in self.ladder_unit_list]) - len(self.ladder_unit_list)
        return left_length, right_length
        
    def calc_verts_C(self, args = {}):
        origin = self.ladder_origin
        left_len, right_len = self.count_vertices()
        left_pos = 0
        right_pos = 0
        for ladder_unit in self.ladder_unit_list:
            if args['style'] == 'ladders':
                width = args['ladder_width']
                height = args['ladder_height']
                if ladder_unit.is_on_right():
                    back_coords = origin + complex(width, height*float(right_pos)/float(right_len))
                    back_label = ladder_unit.right_vertices[0]
                    right_pos += 1
                else:
                    back_coords = origin + complex(0, height*float(left_pos)/float(left_len))
                    back_label = ladder_unit.left_vertices[0]
                    left_pos += 1
                        
                left_front_coords = origin + complex(0, height*float(left_pos)/float(left_len))
                right_front_coords = origin + complex(width, height*float(right_pos)/float(right_len))

                ladder_unit.add_verts_C({back_label:back_coords, ladder_unit.left_vertices[-1]:left_front_coords, ladder_unit.right_vertices[-1]:right_front_coords})
            else:
                assert args['style'] == 'geometric'
                posns_dict = {}
                temp_origin = origin
                if args['ct_depth'] >= 0 and self == self.torus_triang.ladder_list[-1] and ladder_unit.is_on_right():
                    ### when drawing cannon-thurston, put these triangles on the left, so the cannon-thurston paths have triangles on both sides
                    tet = self.vt.tri.tetrahedron(ladder_unit.tet_num)
                    left_vert = ladder_unit.left_vertices[0]
                    neighbour_tet = tet.adjacentTetrahedron(left_vert)
                    gluing = tet.adjacentGluing(left_vert)
                    neighbour_inf_vert = gluing[ladder_unit.face]
                    neighbour_index_in_ladder_list = self.torus_triang.ladder_list[0].ladder_unit_list.index( tet_face(neighbour_tet.index(), neighbour_inf_vert) )
                    neighbour = self.torus_triang.ladder_list[0].ladder_unit_list[neighbour_index_in_ladder_list]
                    back_vert = ladder_unit.right_vertices[0]
                    neighbour_back_vert = gluing[back_vert]
                    offset = convert_to_complex(neighbour.verts_CP1[neighbour_back_vert]) - convert_to_complex(ladder_unit.verts_CP1[back_vert]) 
                    temp_origin = self.torus_triang.ladder_list[0].ladder_origin + offset
                for i in range(4):
                    if ladder_unit.face != i:  # don't include infinity vertex
                        c = convert_to_complex(ladder_unit.verts_CP1[i]) 

                        c = args['geom_complex_scale'] * ( c + temp_origin ) 
                        posns_dict[i] = c
                ladder_unit.add_verts_C(posns_dict)

    def draw(self, my_canvas, args = {}, origin = complex(0,0)):
        delta = args['delta']
        origin = self.ladder_origin
        for ladder_unit in self.ladder_unit_list:
            if args['style'] == 'ladders':
                width = args['ladder_width']
                ladder_unit.draw_triangle_label(my_canvas, curvy = True, args = args)
                ladder_unit.draw_triangle_edges(my_canvas, curvy = True, args = args)
                ladder_unit.draw_labels_curvy(my_canvas, width, delta = delta)
            else:                
                ladder_unit.draw_triangle_label(my_canvas, curvy = False, args = args)
                ladder_unit.draw_triangle_edges(my_canvas, curvy = False, args = args)
                ladder_unit.draw_corner_and_face_labels(my_canvas)
                if args['ct_depth'] >= 0:
                    if ladder_unit.is_on_left():
                        veering_colour = get_edge_between_verts_colour(self.vt, ladder_unit.tet_num, ladder_unit.left_vertices[0], ladder_unit.left_vertices[1])
                        ladder_unit.generate_ct(ladder_is_even = self.is_even, args = args)
                        # print len(ladder_unit.ct_developed_edges)
                        # ladder_unit.draw_ct(my_canvas, origin, args['geom_complex_scale'], colour = colour)
                        ladder_unit.draw_ct_path(my_canvas, origin, args['geom_complex_scale'], veering_colour = veering_colour)
            ladder_unit.draw_vertex_dots(my_canvas)

    def make_ladder(self, start_tf):
        """build a ladder starting from start_tet_face, going along the ladder with the pi vertex as the trailing vertex"""
        
        (tet_index, inf_vert) = start_tf.tet_num, start_tf.face

        current_tet = self.vt.tri.tetrahedron(tet_index)
        current_inf_vert = inf_vert
        current_tf = start_tf
        while True: 
            self.ladder_unit_list.append(ladder_unit(self.vt, current_tf))   
            current_pi_vertex = anglechoice_face2vert[ (self.vt.angle[current_tet.index()], current_inf_vert) ]
            gluing = current_tet.adjacentGluing(current_pi_vertex)
            current_tet = current_tet.adjacentTetrahedron(current_pi_vertex)
            verts_CP1 = None
            if current_tf.verts_CP1 != None:
                verts_CP1 = develop_verts_CP1(current_tf.verts_CP1, gluing, current_pi_vertex, self.vt.tet_shapes[current_tet.index()])
            current_inf_vert = gluing[current_inf_vert]
            current_tf = tet_face( current_tet.index(), current_inf_vert, verts_CP1 = verts_CP1 )
            if current_tf == start_tf:
                not_inf_vert = (start_tf.face + 1) % 4
                if current_tf.verts_CP1 != None:
                    self.holonomy = convert_to_complex(current_tf.verts_CP1[not_inf_vert]) - convert_to_complex(start_tf.verts_CP1[not_inf_vert]) 
                break
        ### if needed, flip ladder
        tet_num, inf_vert = self.ladder_unit_list[0].tet_num, self.ladder_unit_list[0].face
        if self.vt.coorientations[tet_num][inf_vert] == -1: 
            self.ladder_unit_list.reverse()

def draw_vertex_colour(my_canvas, coords, veering_direction):
    colours = {'L':pyx.color.rgb.blue, 'R':pyx.color.rgb.red}
    circ = pyx.path.circle(coords[0], coords[1] ,0.1)
    my_canvas.fill(circ, [pyx.deco.filled([colours[veering_direction]])])

def draw_symmetry_symbol(my_canvas, coords):
    circ = pyx.path.circle(coords.real, coords.imag ,0.18)
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
    while my_ladder.ladder_unit_list[unit_index] != new_tet_face:  ### because of how equality is defined for tet_faces
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
        self.sideways_holonomy = None
        self.make_torus_triangulation(start_tet_face)
        self.tet_faces = []
        for i,l in enumerate(self.ladder_list):
            # print 'ladder', i
            for lu in l.ladder_unit_list:
                # print lu
                self.tet_faces.append(lu)

    def draw(self, my_canvas, args = {}): #style = 'ladders', ladder_width = 5.0, height = 10.0, , ct_depth = 5):
        ### geometric_scale_factor is just so that the text looks good
        for i,L in enumerate(self.ladder_list):
            if args['style'] == 'ladders':
                L.ladder_origin = complex(args['ladder_width'] * i, 0)  ## ignore any stuff already in ladder_origin
            elif args['style'] == 'geometric':
                holonomy = self.ladder_list[0].holonomy
                L.ladder_origin = L.ladder_origin + (i%2) * holonomy 
                geom_complex_scale = args['geometric_scale_factor']*len(self.ladder_list[0].ladder_unit_list) * complex(0,-1) / self.ladder_list[0].holonomy ## rotate and scale
                args['geom_complex_scale'] = geom_complex_scale
            L.calc_verts_C(args = args)
        
        self.draw_symmetries(my_canvas)

        for L in self.ladder_list:
            L.draw(my_canvas, args = args)

    def find_sideways(self, start_tet_face):
        tet_num, face = start_tet_face.tet_num, start_tet_face.face
        """find a loop with canonical horizontal slope"""
        rho = []
        current_tet_face = start_tet_face

        verts_CP1 = None
        if self.vt.tet_shapes != None:  ### set up vertex positions on CP1 for first tetrahedron, so face is at infinity
            verts_CP1 = [None, None, None, None]
            verts_CP1[face] = [1,0]
            verts_CP1[3-face] = [0,1]
            verts_CP1[(face+2)%4] = [1,1]
            last_vert = 3 - ((face+2)%4)
            ordering = unknown_vert_to_known_verts_ordering[last_vert] 
            verts_CP1[last_vert] = developed_position(verts_CP1[ordering[0]], verts_CP1[ordering[1]], verts_CP1[ordering[2]], self.vt.tet_shapes[tet_num])
            # print 'tet, inf, shapes', tet_num, face, verts_CP1
        current_tf = tet_face(current_tet_face.tet_num, current_tet_face.face, verts_CP1 = verts_CP1) 
        
        while current_tf not in rho:  ### ignoring the verts_CP1 information because of our equality definition for tet_faces
            rho.append(current_tf)
            tet_num, inf_vert = current_tf.tet_num, current_tf.face
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
            if self.vt.tet_shapes != None:
                verts_CP1 = current_tf.verts_CP1
            while get_edge_between_verts_colour(self.vt, tet.index(), leading_vertex, pivot_vertex) == start_leading_edge_colour:
                # print tet.index(), '|iplt|', inf_vert, pivot_vertex, leading_vertex, trailing_vertex
                gluing = tet.adjacentGluing(trailing_vertex)
                new_tet = tet.adjacentTetrahedron(trailing_vertex)
                new_inf_vert = gluing[inf_vert]
                new_pivot_vertex = gluing[pivot_vertex]
                new_leading_vertex, new_trailing_vertex = gluing[trailing_vertex], gluing[leading_vertex]
                if self.vt.tet_shapes != None:
                    verts_CP1 = develop_verts_CP1(verts_CP1, gluing, trailing_vertex, self.vt.tet_shapes[new_tet.index()])
                    # print 'tet, inf, shapes', new_tet.index(), new_inf_vert, verts_CP1
                    assert verts_CP1[new_inf_vert] == [1,0]
                tet = new_tet
                inf_vert, pivot_vertex, leading_vertex, trailing_vertex = new_inf_vert, new_pivot_vertex, new_leading_vertex, new_trailing_vertex
            current_tf = tet_face( tet.index(), inf_vert, verts_CP1 = verts_CP1 )
        sideways = rho[rho.index(current_tf):]   ### remove initial tail. Here .index ignores the verts_CP1 data of a tet_face

        not_inf_vert = (current_tf.face + 1) % 4
        self.sideways_holonomy = convert_to_complex(current_tf.verts_CP1[not_inf_vert]) - convert_to_complex(sideways[0].verts_CP1[not_inf_vert]) 
        ### in cases when i(sideways, ladderpole slope) > 1, this isn't what we want for moving things around

        tet_num, inf_vert = sideways[0].tet_num, sideways[0].face
        if self.vt.coorientations[tet_num][inf_vert] == 1: 
            # print 'choice 1'
            sideways[0].origin_in_C += self.sideways_holonomy
            sideways = sideways[1:] + sideways[:1] ## convention: first ladder is convex down, aka green

        ## make sideways go to the right rather than to the left
        tet_num, inf_vert = sideways[0].tet_num, sideways[0].face
        verts = get_triangle_vertex_order(inf_vert)
        cols = [get_edge_between_verts_colour(self.vt, tet_num, inf_vert, v) for v in verts]
        if cols.count('R') == 2:  # pi is on the right, so sideways must be going left
            # print 'choice 2'
            sideways.reverse()
            sideways[0].origin_in_C -= self.sideways_holonomy # put it back where we found it
            sideways = sideways[1:] + sideways[:1] # maintain first ladder convex down
        return sideways

    def make_torus_triangulation(self, start_tet_face):
        """build a torus triangulation by building multiple ladders"""
        
        sideways = self.find_sideways(start_tet_face)

        for tf in sideways:
            if len(self.ladder_list) > 0 and tf in self.ladder_list[0].ladder_unit_list:
            # if self.is_tet_face_in_ladders(tf): ### sideways may wrap multiple times around the torus
                break
            self.ladder_list.append(ladder(self, tf))
        if self.vt.tet_shapes != None:
            for i, L in enumerate(self.ladder_list):
                L.is_even = (i%2 == 0)
                assert abs( (-1)**(i%2) * L.holonomy - self.ladder_list[0].holonomy ) < 0.001 ## all ladder holonomies the same

    def draw_symmetries(self, my_canvas, draw=True):
        count = 0
        dict_of_tet_pairings = {}
        for ladder_index, ladder in enumerate(self.ladder_list):
            for unit_index, unit in enumerate(ladder.ladder_unit_list):
                if unit.is_on_left():
                    edge_axis_unit_pair = self.find_edge_axis_unit_pair(ladder_index, unit_index) #find the pair of triangles to start developing from, checking symmetry
                    vertex_axis_unit_pair = self.find_vertex_axis_unit_pair(ladder_index, unit_index)

                    if self.is_symmetric_torus(edge_axis_unit_pair, dict_of_tet_pairings = dict_of_tet_pairings): #axis is midpoint of the edge to the left of this triangle
                        v0, v1 = unit.verts_C[unit.left_vertices[0]], unit.verts_C[unit.left_vertices[1]]
                        draw_symmetry_symbol(my_canvas, 0.5 * (v0+v1))
                        count += 1
                    if self.is_symmetric_torus(vertex_axis_unit_pair, dict_of_tet_pairings = dict_of_tet_pairings): #axis is the back left vertex from this triangle
                        draw_symmetry_symbol(my_canvas, unit.verts_C[unit.left_vertices[0]])
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
        # print '| symmetries, num of: ', dict_of_tet_pairings 

    def find_edge_axis_unit_pair(self, ladder_index, unit_index): #explore the triangulation, see if the midpoint of the edge to the left of this triangle is a 180 deg rotation axis
        unit_posn_A = (ladder_index, unit_index)
        unit_A = self.ladder_list[ladder_index].ladder_unit_list[unit_index]
        tet_B = self.vt.tri.tetrahedron(unit_A.tet_num).adjacentTetrahedron(unit_A.right_vertices[0])
        tet_B_index = tet_B.index()
        gluing = self.vt.tri.tetrahedron(unit_A.tet_num).adjacentGluing(unit_A.right_vertices[0])
        face_B = gluing[unit_A.face]
        ladder_B_index = (ladder_index - 1) % len(self.ladder_list)
        unit_posn_B = (ladder_B_index, find_unit_index(tet_face(tet_B_index, face_B), self.ladder_list[ladder_B_index]) )
        #print unit_posn_A, unit_A, tet_B_index, face_B
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
        
        tet_B = self.vt.tri.tetrahedron(unit_temp.tet_num).adjacentTetrahedron(unit_temp.right_vertices[0])
        tet_B_index = tet_B.index()
        gluing = self.vt.tri.tetrahedron(unit_temp.tet_num).adjacentGluing(unit_temp.right_vertices[0])
        face_B = gluing[unit_temp.face]
        ladder_B_index = (ladder_index - 1) % len(self.ladder_list)
        unit_posn_B = (ladder_B_index, find_unit_index(tet_face(tet_B_index, face_B), self.ladder_list[ladder_B_index]) )
        #print unit_posn_A, unit_A, tet_B_index, face_B
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
                tet_pairing = [ladder_A.ladder_unit_list[ (unit_A_index + j) % len_A ].tet_num, ladder_B.ladder_unit_list[ (unit_B_index - j) % len_B ].tet_num ]
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
        new_tet_A = self.vt.tri.tetrahedron(unit_A.tet_num).adjacentTetrahedron(unit_A.left_vertices[0])
        new_tet_A_index = new_tet_A.index()
        gluing = self.vt.tri.tetrahedron(unit_A.tet_num).adjacentGluing(unit_A.left_vertices[0])
        new_face_A = gluing[unit_A.face]
        new_ladder_A_index = (ladder_A_index + 1) % len(self.ladder_list)
        new_unit_posn_A = (new_ladder_A_index, find_unit_index(tet_face(new_tet_A_index, new_face_A), self.ladder_list[new_ladder_A_index]) )

        unit_B = ladder_B.ladder_unit_list[unit_B_index]
        new_tet_B = self.vt.tri.tetrahedron(unit_B.tet_num).adjacentTetrahedron(unit_B.right_vertices[0])
        new_tet_B_index = new_tet_B.index()
        gluing = self.vt.tri.tetrahedron(unit_B.tet_num).adjacentGluing(unit_B.right_vertices[0])
        new_face_B = gluing[unit_B.face]
        new_ladder_B_index = (ladder_B_index - 1) % len(self.ladder_list)
        new_unit_posn_B = (new_ladder_B_index, find_unit_index(tet_face(new_tet_B_index, new_face_B), self.ladder_list[new_ladder_B_index]) )

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
                out.append( tet_face(i,j) )
        return out

    def draw(self, output_filename, args = {}):
        canvases = []
        
        for i,T in enumerate(self.torus_triangulation_list):
            print 'cusp:', i, '| num ladders:', len(T.ladder_list)
            c = pyx.canvas.canvas()
            T.draw(c, args = args)
            canvases.append(c)

        out_canvas = pyx.canvas.canvas()
        height_offset = 0.0
        for i,c in enumerate(canvases):
            out_canvas.insert(c, attrs=[pyx.trafo.translate(-c.bbox().left(), height_offset - c.bbox().bottom())])
            height_offset += c.bbox().height() + 0.05 ### add a tiny bit to stop crashes due to line width
        out_canvas.writePDFfile(output_filename)

def generate_boundary_triangulation(tri, angle, args = {}, output_filename = None, draw = True, ct_depth = 20):
    """make a picture of the boundary triangulation, save to output_filename. Assumes that filename is of form '*_xxxx.tri' where xxxx is the angle structure for veering, unless is input in angle_structure_str"""
    vt = veering_triangulation(tri, angle, tet_shapes = args['tet_shapes'])
    B = boundary_triangulation(vt)
    if args['draw']:
        B.draw(output_filename, args = args)  
        
def draw_triangulation_boundary_from_veering_isosig(veering_isosig, args = {}, output_filename = None, verbose = 0.0):
    if verbose > 0.0: print args
    tri, angle = isosig_to_tri_angle(veering_isosig)
    if output_filename == None:
        output_filename = veering_isosig + '.pdf'
    generate_boundary_triangulation(tri, angle, args = args, output_filename = output_filename)

def draw_triangulations_from_veering_isosigs_file(veering_isosigs_filename, output_dirname, args = {}, num_to_draw = None):
    veering_isosigs_list = parse_data_file(veering_isosigs_filename)
    if num_to_draw != None:
        to_draw = veering_isosigs_list[:num_to_draw]
    else:
        to_draw = veering_isosigs_list

    shapes_data = read_from_pickle('Data/veering_shapes_up_to_ten_tetrahedra.pkl')
    names = shapes_data.keys()
    for veering_isosig in to_draw:
        print veering_isosig
        args['tet_shapes'] = shapes_data[veering_isosig]
        draw_triangulation_boundary_from_veering_isosig(veering_isosig, args = args, output_filename = output_dirname + '/' + veering_isosig + '.pdf')

if __name__ == "__main__":

    # Set 'ct_depth': <some non-negative integer> to do cannon-thurston
    args = {'draw':True, 'ct_depth':500, 'ct_epsilon':0.01, 'geometric_scale_factor': 1.5, 'delta': 0.2, 'ladder_width': 10.0, 'ladder_height': 20.0}

    # args['style'] = 'ladders'
    # draw_triangulations_from_veering_isosigs_file('Data/veering_census.txt', 'Images/Boundary_triangulations/Ladders', args = args, num_to_draw = 20)
    args['style'] = 'geometric'
    # draw_triangulations_from_veering_isosigs_file('Data/veering_census.txt', 'Images/Boundary_triangulations/Geometric', args = args, num_to_draw = 50)

    

    # # name = 'cPcbbbdxm_10'
    # name = 'cPcbbbiht_12'
    # # # # name = 'eLMkbcddddedde_2100'
    # # name = 'eLAkbccddhhsqs_1220'
    # name = 'fLAMcaccdeejsnaxk_20010'
    # # # # name = 'fLLQcbecdeepuwsua_20102'
    # # # # name = 'fLLQcbeddeehhbghh_01110'
    # # # # name = 'jLAwwAQbcbdfghihihhwhnaaxrn_211211021'
    # name = 'eLAkaccddjsnak_2001'
    # name = 'eLAkbbcdddhwqj_2102'
    name = 'dLQacccjsnk_200'

    shapes_data = read_from_pickle('Data/veering_shapes_up_to_ten_tetrahedra.pkl')
    args['tet_shapes'] = shapes_data[name]
    draw_triangulation_boundary_from_veering_isosig(name, args = args, output_filename = name + '_' + str(args['ct_depth']) + '_' + str(args['ct_epsilon']) + '.pdf', verbose = 1.0)


    # names = ['kLALPPzkbcbefghgijjxxnsaaqkqqs_0110021020',
    # 'kLALPPzkcbbegfhgijjhhrwaaxnxxn_1221100101',
    # 'kLAMLLAkcbbdeghihjjhhrhhkaarxn_1221211201']
    # for name in names:
    #     print name
    #     draw_triangulation_boundary_from_veering_isosig(name, style = 'ladders', tet_shapes = None) 
        # draw_triangulation_boundary_from_veering_isosig(name, style = 'geometric', tet_shapes = shapes_data[name])


