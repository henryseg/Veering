from file_io import parse_data_file, read_from_pickle
from taut import isosig_to_tri_angle
from veering import veering_triangulation
from develop_ideal_hyperbolic_tetrahedra import developed_position, unknown_vert_to_known_verts_ordering
from basic_math import sign

class vertex:
    def __init__(self, pos):
        self.pos = pos
        self.anticlockwise_edge_is_interesting = None ## is the anticlockwise adjacent coastal edge one we want to draw?
        self.anticlockwise_edge_is_long = None
    
    def __repr__(self):
        if not self.pos.is_infinity():
            return str(self.pos.complex())
        else:
            return "inf"
    
    def __str__(self):
        return self.__repr__()

class landscape_triangle:
    def __init__(self, continent, face_index, is_upper, is_red, vertices, neighbours):
        self.continent = continent
        self.index = face_index ## in the quotient manifold
        self.is_upper = is_upper ## true if we are on the upper landscape of the continent (or were, before we got buried) 
        self.is_red = is_red ## if self has two red edges
        self.vertices = vertices  ## list of three vertices, oriented anticlockwise as viewed from above, with the special vertex first
        ### the special vertex is incident to two edges of the same colour
        self.neighbours = neighbours ## list of three triangles incident to this one, opposite corresponding vertices
        self.is_buried = False  ## either it is inside the continent, or we aren't interested in it for drawing purposes

    def __str__(self):
        return 'continent_ind,triang_ind,upper,red,vertices,buried ' + str([self.continent.triangles.index(self), self.index, self.is_upper, self.is_red, self.vertices, self.is_buried])

    def downriver_index(self):
        """index of triangle in self.neighbours that is downriver of self"""
        if self.is_upper == self.is_red: return 2
        else: return 1

    def downriver(self):
        """neighbour that is downriver of self"""
        return self.neighbours[self.downriver_index()]

    def not_downriver(self):
        """the other two neighbours, in anticlockwise order"""
        i = (self.downriver_index() + 1) % 3
        j = (self.downriver_index() + 2) % 3
        return [ self.neighbours[i], self.neighbours[j] ]

    def update_contacts(self, old_neighbour, new_neighbour):
        ind = self.neighbours.index(old_neighbour)
        self.neighbours[ind] = new_neighbour

    def outwards_tet(self):
        face = self.continent.vt.tri.face(2,self.index)
        embed0 = face.embedding(0)
        embed1 = face.embedding(1)
        tet0 = embed0.simplex()
        tet1 = embed1.simplex()
        if (self.continent.vt.coorientations[tet0.index()][embed0.face()] == -1) == (self.is_upper): ## -1 means into the tetrahedron, so this is above triangle
            return (tet0, embed0)
        else:
            return (tet1, embed1)

    def edge_indices(self):
        tet, embed = self.outwards_tet()  
        face_a = embed.face()  ### our number in the regina tet

        ###   c---R----b
        ###   |`d    ,'|     faces a, b on bottom, c, d on top
        ###   L  ` ,' a| 
        ###   |b ,' .  L 
        ###   |,'    c.| 
        ###   a----R---d 

        ###   b---R----d
        ###   |`a    ,'|     faces c, d on bottom, a, b on top
        ###   L  ` ,' c| 
        ###   |d ,' .  L 
        ###   |,'    b.| 
        ###   c----R---a 

        for i in range(4):
            if i != face_a:
                if self.continent.vt.coorientations[tet.index()][face_a] == self.continent.vt.coorientations[tet.index()][i]:
                    face_b = i
        verts = range(4)
        verts.remove(face_a)
        verts.remove(face_b)
        if (sign([verts[0], face_b, verts[1], face_a]) == 1) == (not self.is_upper): ## get orientation correct  
            face_d, face_c = verts  
        else:
            face_c, face_d = verts
        edge_b_index = self.continent.vt.get_edge_between_verts_index(tet.index(), (face_c, face_d))
        edge_c_index = self.continent.vt.get_edge_between_verts_index(tet.index(), (face_d, face_b))
        edge_d_index = self.continent.vt.get_edge_between_verts_index(tet.index(), (face_b, face_c))
        if self.is_red == self.is_upper:
            return [edge_c_index, edge_d_index, edge_b_index]
        else:
            return [edge_d_index, edge_b_index, edge_c_index]

    def check_against_neighbours(self):
        """sanity checks"""
        if not self.is_buried:

            for i in range(3):
                neighbour = self.neighbours[i]
                assert neighbour.neighbours.count(self) == 1
                assert self.neighbours.count(neighbour) == 1
                j = neighbour.neighbours.index(self)
                self_edge_indices = self.edge_indices()
                neighbour_edge_indices = neighbour.edge_indices()
                ## check vertices
                try:
                    if self.is_upper == neighbour.is_upper:
                        assert self.vertices[(i+1)%3] == neighbour.vertices[(j-1)%3] 
                        assert self.vertices[(i-1)%3] == neighbour.vertices[(j+1)%3]
                    else:
                        assert self.vertices[(i+1)%3] == neighbour.vertices[(j+1)%3] 
                        assert self.vertices[(i-1)%3] == neighbour.vertices[(j-1)%3] 
                except AssertionError:
                    print 'failed vertex check'
                    print i, j
                    print self
                    print neighbour
                    raise
                ## check our colouring for consistency
                try:
                    self_col = ((i == 0) == self.is_red)  ### if True then shared edge is blue
                    neighbour_col = ((j == 0) == neighbour.is_red)
                    assert self_col == neighbour_col
                except AssertionError:
                    print 'failed colour check'
                    print i, j
                    print self
                    print neighbour
                    raise
                ## check edge indices
                try:
                    assert self_edge_indices[i] == neighbour_edge_indices[j]
                except AssertionError:
                    print 'failed edge index check'
                    print i, j
                    print self
                    print self_edge_indices
                    print neighbour
                    print neighbour_edge_indices
                    raise

    def coastal_indices(self):
        return [ i for i in range(3) if self.is_upper != self.neighbours[i].is_upper ] 

        ##  i+1
        ##   *-------* i  
        ##    \     /
        ##     \   / 
        ##      \ /
        ##       * i-1 

    def interesting_coastal_indices(self):
        return [ i for i in self.coastal_indices() if self.vertices[(i+1)%3].anticlockwise_edge_is_interesting ]

    def interesting_long_coastal_indices(self):
        out = []
        ci = self.interesting_coastal_indices()
        for i in ci:
            u = self.vertices[(i+1)%3]
            v = self.vertices[(i+2)%3]
            if abs(u.pos.complex() - v.pos.complex()) > self.continent.max_interesting_edge_length:
                out.append(i)
        return out

class continent:
    def __init__(self, vt, initial_tet_face, desired_vertices = []):
        # print 'initializing continent'
        self.vt = vt
        self.triangles = []
        self.first_non_buried_index = None
        self.num_tetrahedra = 1

        self.tet_face = initial_tet_face
        self.desired_vertices = desired_vertices
        self.boundary_triangulation_vertices = set()

        self.vertices = [vertex(v) for v in self.tet_face.verts_pos]
        self.infinity = self.vertices[self.tet_face.face]
        assert self.infinity.pos.is_infinity()
        self.coast = None
        self.num_long_edges = None
        self.max_interesting_edge_length = None

        # self.infinity = vertex([1,0])
        # self.vertices = [ self.infinity, vertex([0,1]), vertex([1,1]), vertex([self.vt.tet_shapes[initial_tet_num],1]) ]
        
        # self.vertices = [None, None, None, None]  ### came from draw_boundary_triangulation
        # self.vertices[0] = self.infinity
        # self.vertices[3-0] = vertex([0,1])
        # self.vertices[(0+2)%4] = vertex([1,1])
        # last_vert = 3 - ((0+2)%4)
        # ordering = unknown_vert_to_known_verts_ordering[last_vert] 
        # last_vert_CP1 = developed_position(self.vertices[ordering[0]].CP1, self.vertices[ordering[1]].CP1, self.vertices[ordering[2]].CP1, self.vt.tet_shapes[initial_tet_num])
        # self.vertices[last_vert] = vertex(last_vert_CP1)


        # self.vertices_adjacent_to_infinity = []   ## list of (vertex, colour of edge from vertex to infinity)

        # for i in range(4):
        #     if i != self.tet_face.face:
        #         col = self.vt.get_edge_between_verts_colour(self.tet_face.tet_num, (self.tet_face.face, i))
        #         self.vertices_adjacent_to_infinity.append( (self.vertices[i], col) )

        ###   c---R----b
        ###   |`d    ,'|     faces a, b on bottom, c, d on top
        ###   L  ` ,' a| 
        ###   |b ,' .  L 
        ###   |,'    c.| 
        ###   a----R---d 

        upper_face_nums = []
        lower_face_nums = []
        for i in range(4):
            if self.vt.coorientations[self.tet_face.tet_num][i] == +1:
                upper_face_nums.append(i)
            else:
                lower_face_nums.append(i)

        face_a, face_b = lower_face_nums ## decree that face_a is the smaller index
        if sign([face_a, upper_face_nums[1], face_b, upper_face_nums[0]]) == 1: ## get orientation correct
            face_c, face_d = upper_face_nums
        else:
            face_d, face_c = upper_face_nums

        tet = vt.tri.tetrahedron(self.tet_face.tet_num)
        face_a_index = tet.face(2,face_a).index()
        face_b_index = tet.face(2,face_b).index()
        face_c_index = tet.face(2,face_c).index()
        face_d_index = tet.face(2,face_d).index()

        upper_edge_colour = self.vt.get_edge_between_verts_colour(self.tet_face.tet_num, lower_face_nums)
        lower_edge_colour = self.vt.get_edge_between_verts_colour(self.tet_face.tet_num, upper_face_nums)

        ab_is_red = ( lower_edge_colour == 'R' )
        ab_is_upper = False
        cd_is_red = ( upper_edge_colour == 'R' )
        cd_is_upper = True

        triangle_a = landscape_triangle(self, face_a_index, ab_is_upper, ab_is_red, None, None)
        triangle_b = landscape_triangle(self, face_b_index, ab_is_upper, ab_is_red, None, None)
        triangle_c = landscape_triangle(self, face_c_index, cd_is_upper, cd_is_red, None, None)
        triangle_d = landscape_triangle(self, face_d_index, cd_is_upper, cd_is_red, None, None)

        self.triangles.extend([triangle_a, triangle_b, triangle_c, triangle_d])

        ## now for the vertices

        vert_a = self.vertices[face_a]
        vert_b = self.vertices[face_b]
        vert_c = self.vertices[face_c]
        vert_d = self.vertices[face_d]

        if ab_is_red:
            triangle_a.vertices = [vert_c, vert_d, vert_b]
            triangle_b.vertices = [vert_d, vert_c, vert_a]
        else:
            triangle_a.vertices = [vert_d, vert_b, vert_c]
            triangle_b.vertices = [vert_c, vert_a, vert_d]

        if cd_is_red:
            triangle_c.vertices = [vert_a, vert_d, vert_b]
            triangle_d.vertices = [vert_b, vert_c, vert_a]
        else:
            triangle_c.vertices = [vert_b, vert_a, vert_d]
            triangle_d.vertices = [vert_a, vert_b, vert_c]

        ## now for the neighbours

        if ab_is_red:
            triangle_a.neighbours = [triangle_c, triangle_d, triangle_b]
            triangle_b.neighbours = [triangle_d, triangle_c, triangle_a]
        else:
            triangle_a.neighbours = [triangle_d, triangle_b, triangle_c]
            triangle_b.neighbours = [triangle_c, triangle_a, triangle_d]

        if cd_is_red:
            triangle_c.neighbours = [triangle_a, triangle_d, triangle_b]
            triangle_d.neighbours = [triangle_b, triangle_c, triangle_a]
        else:
            triangle_c.neighbours = [triangle_b, triangle_a, triangle_d]
            triangle_d.neighbours = [triangle_a, triangle_b, triangle_c]
        # self.sanity_check()

        if self.desired_vertices != []:
            for v in self.vertices:
                if v != self.infinity:
                    self.check_vertex_desired(v) 

    def check_vertex_desired(self, v, epsilon = 0.001):
        v_in_C = v.pos.complex()
        for w in self.desired_vertices:
            if abs(v_in_C - w) < epsilon:
                self.desired_vertices.remove(w)
                self.boundary_triangulation_vertices.add(v)
                break

    def sanity_check(self):
        for tri in self.triangles:
            tri.check_against_neighbours()

    def bury(self, triangle):
        while not triangle.is_buried:
            self.silt(triangle)
            
    def flow(self, triangle):
        """returns the triangle all the way downriver, and whether it is coastal"""
        while True:
            neighbour = triangle.downriver()
            if neighbour.is_upper != triangle.is_upper:
                return (triangle, True)
            else:
                neighbour2 = neighbour.downriver()  
                if triangle == neighbour2:
                    return (triangle, False)
                else:
                    triangle = neighbour

    def silt(self, triangle):
        """flow downriver from this triangle until we hit either a sink, or the coast, add one tetrahedron there"""
        downriver_triangle, is_coastal = self.flow(triangle)
        if is_coastal:
            self.coastal_fill(downriver_triangle)
        else:
            self.in_fill(downriver_triangle)
        self.num_tetrahedra += 1
        # if self.num_tetrahedra % 1000 == 0:
        #     print self.num_long_edges
        # self.sanity_check()

    def build_fundamental_domain(self, max_num_tetrahedra = 50000):
        self.first_non_buried_index = 0
        while len(self.desired_vertices) > 0 and self.num_tetrahedra < max_num_tetrahedra:  # will go a little over because we check after each bury, which adds many tetrahedra
            tri = self.triangles[self.first_non_buried_index]  
            self.bury(tri)
            self.first_non_buried_index += 1
            while self.triangles[self.first_non_buried_index].is_buried:
            # while self.triangles[first_non_buried_index].is_buried or self.triangles[first_non_buried_index].is_upper:
                self.first_non_buried_index += 1
        self.update_coast()  ## we don't update this as we build

    def build(self, max_interesting_edge_length = 0.1, max_num_tetrahedra = 50000):  # build until all edges we want to draw are short
        self.max_interesting_edge_length = max_interesting_edge_length
        print 'max_interesting_edge_length', max_interesting_edge_length
        self.update_coast()
        
        ## count number of long edges, mark vertices as long
        self.num_long_edges = 0
        for i,v in enumerate(self.coast):
            if v.anticlockwise_edge_is_interesting:
                w = self.coast[i+1]  ### if we have looped then you are Icarus: close to infinity so surely not interesting 
                if abs(v.pos.complex() - w.pos.complex()) > self.max_interesting_edge_length:
                    v.anticlockwise_edge_is_long = True
                    self.num_long_edges += 1

        ## now build

        while self.num_long_edges > 0 and self.num_tetrahedra < max_num_tetrahedra: 
            tri = self.triangles[self.first_non_buried_index]  

            if tri.interesting_long_coastal_indices() != []:
                self.bury(tri)
            self.first_non_buried_index += 1
            while self.triangles[self.first_non_buried_index].is_buried:
                self.first_non_buried_index += 1
        print 'num_long_edges', self.num_long_edges, 'num_tetrahedra', self.num_tetrahedra
        self.update_coast()
        print 'num_long_edges_direct_count', self.count_long_edges()
        print 'max_length', self.calculate_max_interesting_coast_edge_length()



### this plan doesnt seem to work...
    # def bury_uninteresting_triangles(self, interesting_segments):
    #     """We only care about certain parts of the coast for drawing purposes, bury triangles that flow to other parts of the coast."""
    #     for triangle in self.triangles:
    #         coastal_triangle, is_coastal = self.flow(triangle)
    #         assert is_coastal
    #         ind = coastal_triangle.downriver_index()
    #         u = coastal_triangle.vertices[(ind - 1)%3]
    #         v = coastal_triangle.vertices[(ind + 1)%3]  
    #         u_ind, v_ind = self.coast.index(u), self.coast.index(v)
    #         assert u_ind < v_ind  

    #         found_segment = False
    #         for segment in interesting_segments:
    #             if segment[0] <= u_ind and v_ind <= segment[1] or segment[1] <= u_ind and v_ind <= segment[0]:  ### then we are interesting
    #                 found_segment = True
    #                 break
    #
    #         triangle.is_buried = not found_segment ## if we never found it, bury it because we don't care about it or its potential descendents 

    def in_fill(self, triangle):
        # print 'in fill'
        
        neighbour = triangle.downriver()
        assert not triangle.is_buried
        assert not neighbour.is_buried
        assert triangle == neighbour.downriver() and triangle.is_upper == neighbour.is_upper
     
        ###   b---R----t
        ###   |`a    ,'|     is_upper, so faces t and n are below, a and b are new triangles above
        ###   L  ` ,' n|
        ###   |t ,' .  L     
        ###   |,'    b.|
        ###   n----R---a 
 
        ###   n---R----b
        ###   |`t    ,'|     not is_upper, so t and n are above, a and b are new triangles below
        ###   L  ` ,' a| 
        ###   |b ,' .  L     
        ###   |,'    n.|
        ###   a----R---t       

        ## We find the easy information so we can build the new triangles
        tet, embed = triangle.outwards_tet()  
        face_t = embed.face()

        for i in range(4):
            if i != face_t:
                if self.vt.coorientations[tet.index()][face_t] == self.vt.coorientations[tet.index()][i]:
                    face_n = i
        verts = range(4)
        verts.remove(face_t)
        verts.remove(face_n)

        if (sign([verts[0], face_n, verts[1], face_t]) == 1) == (not triangle.is_upper): ## get orientation correct  
            face_b, face_a = verts  
        else:
            face_a, face_b = verts

        face_a_index = tet.face(2,face_a).index()
        face_b_index = tet.face(2,face_b).index()

        far_edge_colour = self.vt.get_edge_between_verts_colour(tet.index(), (face_t, face_n))
        ab_is_red = ( far_edge_colour == 'R' )
        ab_is_upper = triangle.is_upper

        triangle_a = landscape_triangle(self, face_a_index, ab_is_upper, ab_is_red, None, None)
        triangle_b = landscape_triangle(self, face_b_index, ab_is_upper, ab_is_red, None, None)
        self.triangles.extend([triangle_a, triangle_b])

        ## now for the vertices

        if triangle.is_red == triangle.is_upper:
            vert_a, vert_b, vert_n = triangle.vertices
            vert_bn, vert_an, vert_t = neighbour.vertices
        else:
            vert_b, vert_n, vert_a = triangle.vertices
            vert_an, vert_t, vert_bn = neighbour.vertices

        assert vert_b == vert_bn and vert_a == vert_an

        if ab_is_red == ab_is_upper:
            triangle_a.vertices = [vert_t, vert_b, vert_n]
            triangle_b.vertices = [vert_n, vert_a, vert_t]
        else:
            triangle_a.vertices = [vert_n, vert_t, vert_b]
            triangle_b.vertices = [vert_t, vert_n, vert_a]
        
        ## now for the neighbours

        if ab_is_red == ab_is_upper:
            triangle_a.neighbours = [triangle.not_downriver()[0], triangle_b, neighbour.not_downriver()[1]]
            triangle_b.neighbours = [neighbour.not_downriver()[0], triangle_a, triangle.not_downriver()[1]]
        else:
            triangle_a.neighbours = [neighbour.not_downriver()[1], triangle.not_downriver()[0], triangle_b]
            triangle_b.neighbours = [triangle.not_downriver()[1], neighbour.not_downriver()[0], triangle_a]

        triangle.not_downriver()[0].update_contacts(triangle, triangle_a)
        triangle.not_downriver()[1].update_contacts(triangle, triangle_b)
        neighbour.not_downriver()[0].update_contacts(neighbour, triangle_b)
        neighbour.not_downriver()[1].update_contacts(neighbour, triangle_a)

        if self.desired_vertices != []:
            new_edge = [vert_t, vert_n]
            if self.infinity in new_edge:
                new_edge.remove(self.infinity)
                self.check_vertex_desired(new_edge.pop()) 

        triangle.is_buried = True
        neighbour.is_buried = True

    def coastal_fill(self, triangle):
        # print 'coastal fill ' + str(self.triangles.index(triangle)) + ' triang ind ' + str(triangle.index) 
        neighbour_australian = triangle.downriver() ## it's upside down from us
        assert triangle.is_upper != neighbour_australian.is_upper
        assert triangle in neighbour_australian.neighbours ## don't know which one
        assert not triangle.is_buried
        assert not neighbour_australian.is_buried

        ###   b---R----t
        ###   |`a    ,'|     is_upper, so faces t and c are below, a and b are new triangles above
        ###   L  ` ,' c|
        ###   |t ,' .  L     
        ###   |,'    b.|
        ###   c----R---a 
 
        ###   c---R----b
        ###   |`t    ,'|     not is_upper, so t and c are above, a and b are new triangles below
        ###   L  ` ,' a| 
        ###   |b ,' .  L     
        ###   |,'    c.|
        ###   a----R---t   

        ## We find the easy information so we can build the new triangles
        tet, embed = triangle.outwards_tet()  
        face_t = embed.face()
        # print 'tet.index(), face_t', tet.index(), face_t

        for i in range(4):
            if i != face_t:
                if self.vt.coorientations[tet.index()][face_t] == self.vt.coorientations[tet.index()][i]:
                    face_c = i

        verts = range(4)
        verts.remove(face_t)
        verts.remove(face_c)

        if (sign([verts[0], face_c, verts[1], face_t]) == 1) == (not triangle.is_upper): ## get orientation correct  
            face_b, face_a = verts  
        else:
            face_a, face_b = verts
        # print 'face_a, face_b, face_c', face_a, face_b, face_c
        face_a_index = tet.face(2,face_a).index()
        face_b_index = tet.face(2,face_b).index()
        face_c_index = tet.face(2,face_c).index()
        # print 'face triangulation indices: a b c', face_a_index, face_b_index, face_c_index

        far_edge_colour = self.vt.get_edge_between_verts_colour(tet.index(), (face_c, face_t))
        # print 'far edge colour', far_edge_colour
        ab_is_red = ( far_edge_colour == 'R' )
        c_is_red = triangle.is_red
        ab_is_upper = triangle.is_upper
        c_is_upper = not triangle.is_upper

        triangle_a = landscape_triangle(self, face_a_index, ab_is_upper, ab_is_red, None, None)
        triangle_b = landscape_triangle(self, face_b_index, ab_is_upper, ab_is_red, None, None)
        triangle_c = landscape_triangle(self, face_c_index,  c_is_upper,  c_is_red, None, None)

        self.triangles.extend([triangle_a, triangle_b, triangle_c])

        ## now for the vertices

        ### let's find out vert_a, vert_b, vert_c, which are the vertices (with CP1 data) opposite faces.
        if triangle.is_red == triangle.is_upper:
            vert_a, vert_b, vert_c = triangle.vertices
        else:
            vert_b, vert_c, vert_a = triangle.vertices

        tet_vert_posns = [None, None, None, None]
        tet_vert_posns[face_a] = vert_a.pos
        tet_vert_posns[face_b] = vert_b.pos
        tet_vert_posns[face_c] = vert_c.pos

        ### next: permute the triangle verts in CP1 into tet order. Then plug through tet_ordering so we can develop
        tet_shape = self.vt.tet_shapes[tet.index()]
        # print tet_shape
        tet_ordering = unknown_vert_to_known_verts_ordering[face_t]
        vert_t = vertex( developed_position(tet_vert_posns[tet_ordering[0]], tet_vert_posns[tet_ordering[1]], tet_vert_posns[tet_ordering[2]], tet_shape) )
        # print [ vertex(p) for p in [tet_vert_posns[tet_ordering[0]], tet_vert_posns[tet_ordering[1]], tet_vert_posns[tet_ordering[2]]] ]
        # print vert_t
        self.vertices.append(vert_t)

        if self.desired_vertices != []:
            if self.infinity in triangle.vertices:
                self.check_vertex_desired(vert_t) 

        vert_t.anticlockwise_edge_is_interesting = vert_a.anticlockwise_edge_is_interesting

        if vert_a.anticlockwise_edge_is_interesting:
            if vert_a.anticlockwise_edge_is_long:
                self.num_long_edges -= 1

            vert_a.anticlockwise_edge_is_long = (abs(vert_a.pos.complex() - vert_t.pos.complex()) > self.max_interesting_edge_length)
            if vert_a.anticlockwise_edge_is_long:
                self.num_long_edges += 1

            vert_t.anticlockwise_edge_is_long = (abs(vert_t.pos.complex() - vert_b.pos.complex()) > self.max_interesting_edge_length)
            if vert_t.anticlockwise_edge_is_long:
                self.num_long_edges += 1   

        if ab_is_red == ab_is_upper:
            triangle_a.vertices = [vert_t, vert_b, vert_c]
            triangle_b.vertices = [vert_c, vert_a, vert_t]
        else:
            triangle_a.vertices = [vert_c, vert_t, vert_b]
            triangle_b.vertices = [vert_t, vert_c, vert_a]

        if c_is_red != c_is_upper:
            triangle_c.vertices = [vert_b, vert_a, vert_t]
        else:
            triangle_c.vertices = [vert_a, vert_t, vert_b]

        ### now for the neighbours

        if ab_is_red == ab_is_upper:
            triangle_a.neighbours = [triangle.not_downriver()[0], triangle_b, triangle_c]
            triangle_b.neighbours = [triangle_c, triangle_a, triangle.not_downriver()[1]]
        else:
            triangle_a.neighbours = [triangle_c, triangle.not_downriver()[0], triangle_b]
            triangle_b.neighbours = [triangle.not_downriver()[1], triangle_c, triangle_a]

        if c_is_red != c_is_upper:
            triangle_c.neighbours = [triangle_b, triangle_a, neighbour_australian]
        else:
            triangle_c.neighbours = [triangle_a, neighbour_australian, triangle_b]

        triangle.not_downriver()[0].update_contacts(triangle, triangle_a)
        triangle.not_downriver()[1].update_contacts(triangle, triangle_b)
        neighbour_australian.update_contacts(triangle, triangle_c)

        triangle.is_buried = True
  
    def update_coast(self):
        """Returns vertices in anticlockwise order as viewed from above"""  
        out = []
        for tri in self.triangles:
            if not tri.is_buried:
                break  ## found an initial unburied tri 
        vert_index = 0  
        initial_vert = tri.vertices[vert_index]
        vert = initial_vert   
        while out == [] or vert != initial_vert:
            if tri.neighbours[(vert_index - 1) % 3].is_upper != tri.is_upper: ## we are coastal
                vert_index = (vert_index + 1) % 3   
                vert = tri.vertices[vert_index]
                out.append(vert)
            else: # walk to the next triangle
                tri = tri.neighbours[(vert_index - 1) % 3]  
                vert_index = tri.vertices.index(vert) 

        ##  i+1
        ##   *-------* i  
        ##    \     /
        ##     \   / 
        ##      \ /
        ##       * i-1 

        ## now rotate to put infinity first
        inf_vert_index = out.index( self.infinity )
        out = out[inf_vert_index:] + out[:inf_vert_index]
        self.coast = out

    def mark_interesting_segments(self, interesting_segments):
        for i,v in enumerate(self.coast):
            for segment in interesting_segments:
                if segment[0] <= i < segment[1]:
                    v.anticlockwise_edge_is_interesting = True

    def segment_between(self, u, v):
        """return the segment of the coast between u and v inclusive"""
        assert u != v
        u_ind = self.coast.index(u)
        v_ind = self.coast.index(v)
        if u_ind < v_ind:
            return self.coast[u_ind : v_ind + 1]
        else:
            return self.coast[v_ind : u_ind + 1]
        
    def vertices_and_edges_adjacent_to_infinity(self):
        ## vertices is set of (vertex, colour of edge from vertex to infinity)
        ## edges is list of ( (u,v), colour ) where u, v are vertices at either end of an edge opposite infinity in a triangle
        vertices = set()
        edges = []
        for triangle in self.triangles:
            if self.infinity in triangle.vertices:
                ind = triangle.vertices.index(self.infinity)
                endpoints = ( triangle.vertices[(ind+1)%3], triangle.vertices[(ind+2)%3] )
                edge_veering_colour = self.vt.veering_colours[ triangle.edge_indices()[ind] ]
                edges.append( (endpoints, edge_veering_colour) ) 

                for i in range(1,3):
                    vert = triangle.vertices[(ind+i)%3]
                    vertex_veering_colour = self.vt.veering_colours[ triangle.edge_indices()[(ind-i)%3] ]
                    vertices.add( (vert, vertex_veering_colour) )
        return vertices, edges

    def calculate_max_interesting_coast_edge_length(self):
        max_length = 0.0
        for i,v in enumerate(self.coast):
            if v.anticlockwise_edge_is_interesting:
                edge_length = abs( v.pos.complex() - self.coast[i+1].pos.complex() )
                if edge_length > max_length:
                    max_length = edge_length
        return max_length

    def count_long_edges(self):
        out = 0
        for i,v in enumerate(self.coast):
            if v.anticlockwise_edge_is_interesting:
                edge_length = abs( v.pos.complex() - self.coast[i+1].pos.complex() )
                if edge_length > self.max_interesting_edge_length:
                    out += 1
        return out


if __name__ == '__main__':

    veering_isosig = 'dLQacccjsnk_200'

    shapes_data = read_from_pickle('Data/veering_shapes_up_to_ten_tetrahedra.pkl')

    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle, tet_shapes = shapes_data[veering_isosig])
    # con = continent( vt )
    # con.build(5)
    # print con.coast()





