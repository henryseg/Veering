from file_io import parse_data_file, read_from_pickle
from taut import isosig_to_tri_angle
from veering import veering_triangulation
from develop_ideal_hyperbolic_tetrahedra import developed_position, unknown_vert_to_known_verts_ordering
from basic_math import sign

class vertex:
    def __init__(self, continent, pos):
        self.continent = continent
        self.edges = []
        self.pos = pos
        self.ladderpole_ancestors = set() ## which ladderpole edges did you come from
        self.continent.vertices.append(self)
        self.coastal_index = None
        self.circle_pos = None
    
    def is_ladderpole_descendant(self):
        return len(self.ladderpole_ancestors) != 0

    def edge_between(self, other):
        for e in self.edges:
            if other in e.vertices:
                return e  ### note, no parallel edges in univ cover of a veering triangulation
        assert False

    def __repr__(self):
        if self.pos == None:
            return str(coastal_index)
        elif not self.pos.is_infinity():
            return str(self.pos.complex())
        else:
            return "inf"
    
    def __str__(self):
        return self.__repr__()

class landscape_edge:
    def __init__(self, continent, vertices, is_red): ## could add: is_red, edge_index
        self.continent = continent
        self.vertices = vertices
        for v in self.vertices:
            v.edges.append(self)
        self.continent.edges.append(self)
        self.is_red = is_red
        self.boundary_triangles = [] ### updated in update_boundary. If this edge is on boundary of the continent then there will be two triangles in here.
        # try:
        #     assert self.length() > 0.0001
        # except:
        #     print 'edge too short', self.vertices
        #     raise

    def __repr__(self):
        u, v = self.vertices
        return ' '.join( [str(self.continent.edges.index(self)), 'edge', str(u), str(v), str(self.length())] )

    def length(self):
        u, v = self.vertices
        return abs(u.pos.complex() - v.pos.complex())

    def is_long(self):
        return self.length() > self.continent.max_length

    def is_under_ladderpole(self):
        return any( v.is_ladderpole_descendant() for v in self.vertices )

    def midpoint(self):
        u, v = self.vertices
        return 0.5*(u.pos.complex() + v.pos.complex())

    def shared_vertex(self, other):
        intersection = set(self.vertices) & set(other.vertices) 
        assert len(intersection) == 1
        return intersection.pop()

    def is_coastal(self):
        return self.boundary_triangles[0].is_upper != self.boundary_triangles[1].is_upper

    def is_watershed(self):
        if self.is_coastal():
            return False
        else:
            return (self.boundary_triangles[0].downriver() != self.boundary_triangles[1]) and (self.boundary_triangles[1].downriver() != self.boundary_triangles[0]) 
                   ### fall if one of these is equal, sink if both are equal
    
    def is_coastal_sink(self, upper = True):  ### have to say if we are interested in the upper or lower landscape
        if not self.is_coastal():
            return False
        else:
            for tri in self.boundary_triangles:
                if tri.is_upper == upper: ### we are looking at the correct triangle
                    return tri.edges[tri.downriver_index()] == self

class landscape_triangle:
    def __init__(self, continent, face_index, is_upper, is_red, vertices, edges, neighbours):
        self.continent = continent
        self.index = face_index ## in the quotient manifold
        self.is_upper = is_upper ## true if we are on the upper landscape of the continent (or were, before we got buried) 
        self.is_red = is_red ## if self has two red edges
        self.vertices = vertices  ## list of three vertices, oriented anticlockwise as viewed from above, with the special vertex first
        ### the special vertex is incident to two edges of the same colour

        ###    *                *
        ###   / \              / \
        ###  R   R     or     B   B  
        ### /     \          /     \
        ### ---B---          ---R---

        self.edges = edges
        self.upper_tet = None
        self.lower_tet = None
        self.neighbours = neighbours ## list of three triangles incident to this one, opposite corresponding vertices
        self.is_buried = False  ## either it is inside the continent, or we aren't interested in it for drawing purposes
        self.continent.triangles.append(self)

    def __str__(self):
        return 'continent_ind,triang_ind,upper,red,vertices,buried ' + str([self.continent.triangles.index(self), self.index, self.is_upper, self.is_red, self.vertices, self.is_buried])

    def set_upper_tet(self, con_tet):
        self.upper_tet = con_tet
        con_tet.lower_triangles.append(self)
        assert len(con_tet.lower_triangles) <= 2

    def set_lower_tet(self, con_tet):
        self.lower_tet = con_tet
        con_tet.upper_triangles.append(self)
        assert len(con_tet.upper_triangles) <= 2

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

    def outwards_tet_all_vert_indices(self):
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
        verts = list(range(4))
        verts.remove(face_a)
        verts.remove(face_b)
        if (sign([verts[0], face_b, verts[1], face_a]) == 1) == (not self.is_upper): ## get orientation correct  
            face_d, face_c = verts  
        else:
            face_c, face_d = verts
        return (face_a, face_b, face_c, face_d)      

    def outwards_tet_our_vert_indices(self):
        face_a, face_b, face_c, face_d = self.outwards_tet_all_vert_indices()
        if self.is_red == self.is_upper:
            return [face_c, face_d, face_b]
        else:
            return [face_d, face_b, face_c]

    def edge_indices(self):
        face_a, face_b, face_c, face_d = self.outwards_tet_all_vert_indices()
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
                    print('failed vertex check')
                    print((i, j))
                    print(self)
                    print(neighbour)
                    raise
                ## check our colouring for consistency
                try:
                    self_col = ((i == 0) == self.is_red)  ### if True then shared edge is blue
                    neighbour_col = ((j == 0) == neighbour.is_red)
                    assert self_col == neighbour_col
                except AssertionError:
                    print('failed colour check')
                    print((i, j))
                    print(self)
                    print(neighbour)
                    raise
                ## check edge indices
                try:
                    assert self_edge_indices[i] == neighbour_edge_indices[j]
                except AssertionError:
                    print('failed edge index check')
                    print((i, j))
                    print(self)
                    print(self_edge_indices)
                    print(neighbour)
                    print(neighbour_edge_indices)
                    raise

    def coastal_indices(self):
        assert not self.is_buried
        return [ i for i in range(3) if self.is_upper != self.neighbours[i].is_upper ] 

        ##  i+1
        ##   *-------* i  
        ##    \     /
        ##     \   / 
        ##      \ /
        ##       * i-1 

    def ladderpole_descendant_coastal_indices(self):
        return [ i for i in self.coastal_indices() if self.edges[i].is_under_ladderpole() ]

    def ladderpole_descendant_long_coastal_indices(self):
        return [ i for i in self.ladderpole_descendant_coastal_indices() if self.edges[i].is_long() ]

    # def distance_to_prong(self):
    def end_of_river_edge(self):
        assert not self.is_buried  
        tri, is_coastal = self.continent.flow(self)
        return tri.edges[tri.downriver_index()]

    def all_river_crossing_edges(self):
        assert not self.is_buried  
        tris, is_coastal = self.continent.flow(self, return_all = True)
        return [tri.edges[tri.downriver_index()] for tri in tris]

    def shared_edge(self, other):
        intersection = set(self.edges) & set(other.edges) 
        assert len(intersection) == 1
        return intersection.pop()

    # def dist_from_lox(self, edge):
    #     vt = self.continent.vt
    #     edge_index = self.edges.index(edge)
    #     tet, embed = self.outwards_tet()

    #     ### get the tetrahedron outwards of the edge 
    #     our_regina_vert_indices = self.outwards_tet_our_vert_indices()

    #     init_verts_pos = [None, None, None, None]
    #     for i, ind in enumerate(our_regina_vert_indices):
    #         init_verts_pos[ind] = self.vertices[i].pos

    #     tet_shape = vt.tet_shapes[tet.index()]
    #     tet_ordering = unknown_vert_to_known_verts_ordering[embed.face]
    #     init_verts_pos[embed.face] = developed_position(tet_vert_posns[tet_ordering[0]], tet_vert_posns[tet_ordering[1]], tet_vert_posns[tet_ordering[2]], tet_shape) 

    #     lox = loxodromic_from_flag(vt, tet.index(), embed.face, our_regina_vert_indices[edge_index], init_verts_pos = init_verts_pos)

    #     ### has two fixed points... get dist of midpoint of the edge to the correct fixed point...

class continent_tetrahedron:
    def __init__(self, continent, tet_index):
        self.continent = continent
        self.index = tet_index ## in the quotient manifold
        self.upper_triangles = [] 
        self.lower_triangles = []
        self.continent.tetrahedra.append(self)

    def upper_edge(self):
        return self.upper_triangles[0].shared_edge(self.upper_triangles[1])

    def lower_edge(self):
        return self.lower_triangles[0].shared_edge(self.lower_triangles[1])

    def edges(self):
        out = self.upper_triangles[0].edges[:] ### copy of list
        out.remove(self.upper_edge())
        return out + self.upper_triangles[1].edges + [self.lower_edge()]

    def vertices(self):
        out = set([])
        for tri in self.upper_triangles:
            out = out.union(set(tri.vertices))
        return out

class continent:
    def __init__(self, vt, initial_tet_face, desired_vertices = []):
        # print 'initializing continent'
        self.vt = vt
        self.triangles = [] 
        self.first_non_buried_index = None
        self.num_tetrahedra = 1

        self.tet_face = initial_tet_face
        self.desired_vertices = desired_vertices[:]  ## copy
        self.boundary_triangulation_vertices = set()

        self.edges = []
        self.vertices = []
        self.tetrahedra = []
        for v in self.tet_face.verts_pos:
            vertex(self, v)  ## creates and adds to the list of vertices
        self.infinity = self.vertices[self.tet_face.face]
        if self.infinity.pos != None:
            assert self.infinity.pos.is_infinity() 
        self.coast = None
        self.max_length = None

        self.upper_landscape_triangles = set([]) ### updates in update_boundary
        self.lower_landscape_triangles = set([]) 
        self.upper_landscape_edges = set([])
        self.lower_landscape_edges = set([])
        self.coastal_edges = []

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

        ab_is_red = ( lower_edge_colour == "red" )  ## the triangles, not the edge
        ab_is_upper = False
        cd_is_red = ( upper_edge_colour == "red" )
        cd_is_upper = True

        triangle_a = landscape_triangle(self, face_a_index, ab_is_upper, ab_is_red, None, None, None)
        triangle_b = landscape_triangle(self, face_b_index, ab_is_upper, ab_is_red, None, None, None)
        triangle_c = landscape_triangle(self, face_c_index, cd_is_upper, cd_is_red, None, None, None)
        triangle_d = landscape_triangle(self, face_d_index, cd_is_upper, cd_is_red, None, None, None)

        ## add the tetrahedron

        con_tet = continent_tetrahedron(self, self.tet_face.tet_num)
        triangle_a.set_upper_tet(con_tet)
        triangle_b.set_upper_tet(con_tet)
        triangle_c.set_lower_tet(con_tet)
        triangle_d.set_lower_tet(con_tet)

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

        ## now for the edges

        ###   c---R----b
        ###   |`d    ,'|     faces a, b on bottom, c, d on top
        ###   L  ` ,' a| 
        ###   |b ,' .  L 
        ###   |,'    c.| 
        ###   a----R---d 

        edge_ab = landscape_edge(self, [vert_a, vert_b], upper_edge_colour == "red")
        edge_ac = landscape_edge(self, [vert_a, vert_c], False)
        edge_ad = landscape_edge(self, [vert_a, vert_d], True)
        edge_bc = landscape_edge(self, [vert_b, vert_c], True)
        edge_bd = landscape_edge(self, [vert_b, vert_d], False)
        edge_cd = landscape_edge(self, [vert_c, vert_d], lower_edge_colour == "red")

        if ab_is_red: ## the triangles a and b, not the edge
            triangle_a.edges = [edge_bd, edge_bc, edge_cd]
            triangle_b.edges = [edge_ac, edge_ad, edge_cd]
        else:
            triangle_a.edges = [edge_bc, edge_cd, edge_bd]
            triangle_b.edges = [edge_ad, edge_cd, edge_ac]

        if cd_is_red: 
            triangle_c.edges = [edge_bd, edge_ab, edge_ad]
            triangle_d.edges = [edge_ac, edge_ab, edge_bc]
        else:
            triangle_c.edges = [edge_ad, edge_bd, edge_ab]
            triangle_d.edges = [edge_bc, edge_ac, edge_ab]

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

        # self.update_boundary()

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
            
    def flow(self, triangle, return_all = False):
        """returns the triangle all the way downriver, and whether it is coastal"""
        all_out = []
        while True:
            all_out.append(triangle)
            neighbour = triangle.downriver()
            if neighbour.is_upper != triangle.is_upper:
                if return_all:
                    return (all_out, True)
                else:
                    return (triangle, True)
            else:
                neighbour2 = neighbour.downriver()  
                if triangle == neighbour2:
                    if return_all:
                        return (all_out, False)
                    else:
                        return (triangle, False)
                else:
                    triangle = neighbour

    def silt(self, triangle):
        """flow downriver from this triangle until we hit either a sink, or the coast, add one tetrahedron there"""
        downriver_triangle, is_coastal = self.flow(triangle)
        if is_coastal:
            new_tet = self.coastal_fill(downriver_triangle)
        else:
            new_tet = self.in_fill(downriver_triangle)
        # self.sanity_check()
        return new_tet

    def bury(self, triangle):
        while not triangle.is_buried:
            new_tet = self.silt(triangle)
        if triangle.is_upper:
            assert triangle in new_tet.lower_triangles
        else:
            assert triangle in new_tet.upper_triangles
        return new_tet ### return the last added tetrahedron
  
    def update_boundary(self):
        """Installs coastal vertices in anticlockwise order as viewed from above
           Also for all boundary edges tells them neighbouring faces"""  

        self.upper_landscape_triangles = set([])
        self.lower_landscape_triangles = set([])
        for tri in self.triangles:
            if not tri.is_buried:
                if tri.is_upper:
                    self.upper_landscape_triangles.add(tri)
                else:
                    self.lower_landscape_triangles.add(tri)

        self.upper_landscape_edges = set([])
        self.lower_landscape_edges = set([])

        for e in self.edges:
            e.boundary_triangles = []

        for tri in self.upper_landscape_triangles:
            for e in tri.edges:
                e.boundary_triangles.append(tri)
                self.upper_landscape_edges.add(e)
        for tri in self.lower_landscape_triangles:
            for e in tri.edges:
                e.boundary_triangles.append(tri)
                self.lower_landscape_edges.add(e)

        self.coast = []
        for tri in self.triangles:
            if not tri.is_buried:
                break  ## found an initial unburied tri 
        vert_index = 0  
        initial_vert = tri.vertices[vert_index]
        vert = initial_vert   
        while self.coast == [] or vert != initial_vert:
            if tri.neighbours[(vert_index - 1) % 3].is_upper != tri.is_upper: ## we are coastal
                vert_index = (vert_index + 1) % 3   
                vert = tri.vertices[vert_index]
                self.coast.append(vert)
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
        inf_vert_index = self.coast.index( self.infinity )
        self.coast = self.coast[inf_vert_index:] + self.coast[:inf_vert_index]

        ## update coastal_edges
        self.coastal_edges = []
        for i in range(len(self.coast)):
            u = self.coast[i]
            v = self.coast[(i+1)%len(self.coast)]
            for e in u.edges:
                if v in e.vertices:
                    e.vertices = [u, v]
                    self.coastal_edges.append(e)
                    break


    def make_convex(self):
        ### new triangles are added to the end of the list so this is safe.
        ### new sinks created by in-fills have new triangles that point at them, so we get everything.
        for tri in self.triangles:
            if not tri.is_buried:
                downriver_triangle, is_coastal = self.flow(tri)
                if not is_coastal:
                    self.in_fill(downriver_triangle)
        self.update_boundary()
                        
    def build_fundamental_domain_old(self, max_num_tetrahedra = 50000):
        self.first_non_buried_index = 0
        while len(self.desired_vertices) > 0 and self.num_tetrahedra < max_num_tetrahedra:  # will go a little over because we check after each bury, which adds many tetrahedra
            tri = self.triangles[self.first_non_buried_index]  
            self.bury(tri)
            self.first_non_buried_index += 1
            while self.triangles[self.first_non_buried_index].is_buried:
            # while self.triangles[first_non_buried_index].is_buried or self.triangles[first_non_buried_index].is_upper:
                self.first_non_buried_index += 1
        self.update_boundary()  ## we don't update this as we build

    ### old version builds lots of things we dont care about, this is much faster.
    def build_fundamental_domain(self, max_num_tetrahedra = 50000):
        self.first_non_buried_index = 0
        while len(self.desired_vertices) > 0 and self.num_tetrahedra < max_num_tetrahedra:  # will go a little over because we check after each bury, which adds many tetrahedra
            tri = self.triangles[self.first_non_buried_index] 
             ### if this tri is incident to infinity, bury it
            if self.infinity in tri.vertices:
                self.bury(tri)
            self.first_non_buried_index += 1
            while self.triangles[self.first_non_buried_index].is_buried:
            # while self.triangles[first_non_buried_index].is_buried or self.triangles[first_non_buried_index].is_upper:
                self.first_non_buried_index += 1
        self.update_boundary()  ## we don't update this as we build

    def build_naive(self, max_num_tetrahedra = 50000):  ### just keep building until we hit max tetrahedra
        self.first_non_buried_index = 0
        while self.num_tetrahedra < max_num_tetrahedra:  
            tri = self.triangles[self.first_non_buried_index]  
            self.bury(tri)
            self.first_non_buried_index += 1
            while self.triangles[self.first_non_buried_index].is_buried:
                self.first_non_buried_index += 1
        self.update_boundary()  ## we don't update this as we build

    def build_on_coast(self, max_length = 0.1, max_num_tetrahedra = 50000):  # build until all edges we want to draw are short
        self.max_length = max_length
        # print(('max_length', max_length))
        self.update_boundary()

        ## now build

        # while self.num_long_edges > 0 and self.num_tetrahedra < max_num_tetrahedra: 
        while self.num_tetrahedra < max_num_tetrahedra and self.first_non_buried_index < len(self.triangles): 
            tri = self.triangles[self.first_non_buried_index]  

            if tri.ladderpole_descendant_long_coastal_indices() != []:
                self.bury(tri)
            self.first_non_buried_index += 1
            while self.first_non_buried_index < len(self.triangles) and self.triangles[self.first_non_buried_index].is_buried:
                self.first_non_buried_index += 1
        # print(('num_tetrahedra', self.num_tetrahedra))
        hit_max_tetrahedra = self.num_tetrahedra >= max_num_tetrahedra
        # print(('hit max tetrahedra', hit_max_tetrahedra))
        self.update_boundary()
        # print(('num_long_edges_direct_count', self.count_long_edges()))
        # print(('max_coastal_edge_length', self.calculate_max_ladderpole_descendant_coast_edge_length()))
        return hit_max_tetrahedra

    def build_make_long_descendant_edges_internal(self, max_length = 0.1, max_num_tetrahedra = 50000):  # build until all edges we want to draw are short
        self.max_length = max_length
        # print(('max_length', max_length))

        ## now build

        while self.first_non_buried_index < len(self.triangles) and self.num_tetrahedra < max_num_tetrahedra: 
            tri = self.triangles[self.first_non_buried_index]  

            if any( [(edge.is_under_ladderpole() and edge.is_long()) for edge in tri.edges] ):
                self.bury(tri)
            self.first_non_buried_index += 1
            while self.first_non_buried_index < len(self.triangles) and self.triangles[self.first_non_buried_index].is_buried:
                self.first_non_buried_index += 1
        # print(('num_tetrahedra', self.num_tetrahedra))
        hit_max_tetrahedra = self.num_tetrahedra >= max_num_tetrahedra
        # print(('hit max tetrahedra', hit_max_tetrahedra))
        self.update_boundary()
        # print(('num_long_edges_direct_count', self.count_long_edges()))
        # print(('max_coastal_edge_length', self.calculate_max_ladderpole_descendant_coast_edge_length()))
        return hit_max_tetrahedra

    def build_explore_prongs(self, max_length = 0.1, max_num_tetrahedra = 50000):  # build until all edges we want to draw are short
        self.max_length = max_length
        # print(('max_length', max_length))

        ## now build

        while self.first_non_buried_index < len(self.triangles) and self.num_tetrahedra < max_num_tetrahedra: 
            tri = self.triangles[self.first_non_buried_index]  

            crossing_edges = tri.all_river_crossing_edges()
            for e in crossing_edges:   ### previously only used the last edge
            # for e in crossing_edges[-1:]:   ### only the last edge
                if e.is_under_ladderpole():
                    u = tri.vertices[tri.downriver_index()]
                    m = e.midpoint()
                    distance_to_prong = abs(u.pos.complex() - m)
                    if distance_to_prong > max_length:
                        self.bury(tri)
                        break
                # dist_to_v = abs(u.pos.complex() - v.pos.complex())
                # dist_to_w = abs(u.pos.complex() - w.pos.complex())
                # if dist_to_v > self.max_length or dist_to_w > self.max_length:
                #     if 0.01 < dist_to_v/dist_to_w < 100:  ## ratio is not too weird
                #         self.bury(tri)
                #     else:
                #         print "big ratio..."
            self.first_non_buried_index += 1
            while self.first_non_buried_index < len(self.triangles) and self.triangles[self.first_non_buried_index].is_buried:
                self.first_non_buried_index += 1
        # print(('num_tetrahedra', self.num_tetrahedra))
        hit_max_tetrahedra = self.num_tetrahedra >= max_num_tetrahedra
        # print(('hit max tetrahedra', hit_max_tetrahedra))
        self.update_boundary()
        # print(('num_long_edges_direct_count', self.count_long_edges()))
        # print(('max_coastal_edge_length', self.calculate_max_ladderpole_descendant_coast_edge_length()))
        return hit_max_tetrahedra

    # def build_loxodromics(self, max_length = 0.1, max_num_tetrahedra = 50000):
    #     self.max_length = max_length
    #     print 'max_length', max_length
    #     self.update_boundary()

    #     ## now build

    #     while self.first_non_buried_index < len(self.triangles) and self.num_tetrahedra < max_num_tetrahedra: 
    #         tri = self.triangles[self.first_non_buried_index]  

    #         if any( [(edge.is_ladderpole_descendant and tri.dist_from_lox(edge) > self.max_length) for edge in tri.edges] ):
    #             self.bury(tri)
    #         self.first_non_buried_index += 1
    #         while self.first_non_buried_index < len(self.triangles) and self.triangles[self.first_non_buried_index].is_buried:
    #             self.first_non_buried_index += 1
    #     print 'num_tetrahedra', self.num_tetrahedra
    #     self.update_boundary()
    #     print 'num_long_edges_direct_count', self.count_long_edges()
    #     print 'max_coastal_edge_length', self.calculate_max_ladderpole_descendant_coast_edge_length()

    def build_long_and_mid(self, max_length = 0.1, max_num_tetrahedra = 50000):  
        self.max_length = max_length
        # print(('max_length', max_length))
        self.update_boundary()

        ## now build

        while self.first_non_buried_index < len(self.triangles) and self.num_tetrahedra < max_num_tetrahedra: 
            tri = self.triangles[self.first_non_buried_index]  

            is_long = any( (edge.is_under_ladderpole() and edge.is_long()) for edge in tri.edges )
            
            mid_is_far = False
            crossing_edges = tri.all_river_crossing_edges()
            for e in crossing_edges:   ### previously only used the last edge
                if e.is_under_ladderpole():
                    u = tri.vertices[tri.downriver_index()]
                    m = e.midpoint()
                    distance_to_prong = abs(u.pos.complex() - m)
                    if distance_to_prong > max_length: 
                        mid_is_far = True
                        break

            if is_long or mid_is_far:
                self.bury(tri)

            self.first_non_buried_index += 1
            while self.first_non_buried_index < len(self.triangles) and self.triangles[self.first_non_buried_index].is_buried:
                self.first_non_buried_index += 1
        # print(('num_tetrahedra', self.num_tetrahedra))
        hit_max_tetrahedra = self.num_tetrahedra >= max_num_tetrahedra
        # print(('hit max tetrahedra', hit_max_tetrahedra))
        self.update_boundary()
        # print(('num_long_edges_direct_count', self.count_long_edges()))
        # print(('max_coastal_edge_length', self.calculate_max_ladderpole_descendant_coast_edge_length()))
        return hit_max_tetrahedra


    def in_fill(self, triangle):
        # print 'in fill'
        
        neighbour = triangle.downriver()
        assert not triangle.is_buried
        assert not neighbour.is_buried
        assert triangle == neighbour.downriver() and triangle.is_upper == neighbour.is_upper and triangle.is_red == neighbour.is_red
     
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
        verts = list(range(4))
        verts.remove(face_t)
        verts.remove(face_n)

        if (sign([verts[0], face_n, verts[1], face_t]) == 1) == (not triangle.is_upper): ## get orientation correct  
            face_b, face_a = verts  
        else:
            face_a, face_b = verts

        face_a_index = tet.face(2,face_a).index()
        face_b_index = tet.face(2,face_b).index()

        far_edge_colour = self.vt.get_edge_between_verts_colour(tet.index(), (face_t, face_n))
        ab_is_red = ( far_edge_colour == "red" )
        ab_is_upper = triangle.is_upper

        triangle_a = landscape_triangle(self, face_a_index, ab_is_upper, ab_is_red, None, None, None)
        triangle_b = landscape_triangle(self, face_b_index, ab_is_upper, ab_is_red, None, None, None)

        ## add the tetrahedron 

        con_tet = continent_tetrahedron(self, tet.index())
        if triangle.is_upper:
            triangle.set_upper_tet(con_tet)
            neighbour.set_upper_tet(con_tet)
            triangle_a.set_lower_tet(con_tet)
            triangle_b.set_lower_tet(con_tet)
        else:
            triangle.set_lower_tet(con_tet)
            neighbour.set_lower_tet(con_tet)
            triangle_a.set_upper_tet(con_tet)
            triangle_b.set_upper_tet(con_tet)            

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
        
        ## now for the edges

        edge_tn = landscape_edge(self, [vert_t, vert_n], far_edge_colour == "red") ## never coastal

        if triangle.is_red == triangle.is_upper:
            edge_bn, edge_na, edge_ab = triangle.edges 
            edge_at, edge_tb, edge_ba = neighbour.edges
        else:
            edge_na, edge_ab, edge_bn = triangle.edges 
            edge_tb, edge_ba, edge_at = neighbour.edges
        assert edge_ab == edge_ba

        if ab_is_red == ab_is_upper:
            triangle_a.edges = [edge_bn, edge_tn, edge_tb]
            triangle_b.edges = [edge_at, edge_tn, edge_na]
        else:
            triangle_a.edges = [edge_tb, edge_bn, edge_tn]
            triangle_b.edges = [edge_na, edge_at, edge_tn]

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
        self.num_tetrahedra += 1
        return con_tet ### pass back the tetrahedron we added

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

        ### at is red iff not triangle.is_upper
        ### bt is red iff triangle.is_upper

        ## We find the easy information so we can build the new triangles
        tet, embed = triangle.outwards_tet()  
        face_t = embed.face()
        # print 'tet.index(), face_t', tet.index(), face_t

        for i in range(4):
            if i != face_t:
                if self.vt.coorientations[tet.index()][face_t] == self.vt.coorientations[tet.index()][i]:
                    face_c = i

        verts = list(range(4))
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
        ab_is_red = ( far_edge_colour == "red" )
        c_is_red = triangle.is_red
        ab_is_upper = triangle.is_upper
        c_is_upper = not triangle.is_upper

        triangle_a = landscape_triangle(self, face_a_index, ab_is_upper, ab_is_red, None, None, None)
        triangle_b = landscape_triangle(self, face_b_index, ab_is_upper, ab_is_red, None, None, None)
        triangle_c = landscape_triangle(self, face_c_index,  c_is_upper,  c_is_red, None, None, None)

        ## add the tetrahedron 

        con_tet = continent_tetrahedron(self, tet.index())
        if triangle.is_upper:
            triangle.set_upper_tet(con_tet)
            triangle_c.set_upper_tet(con_tet)
            triangle_a.set_lower_tet(con_tet)
            triangle_b.set_lower_tet(con_tet)
        else:
            triangle.set_lower_tet(con_tet)
            triangle_c.set_lower_tet(con_tet)
            triangle_a.set_upper_tet(con_tet)
            triangle_b.set_upper_tet(con_tet)  

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

        if self.vt.tet_shapes != None:
            ### next: permute the triangle verts in CP1 into tet order. Then plug through tet_ordering so we can develop
            tet_shape = self.vt.tet_shapes[tet.index()]
            # print tet_shape
            tet_ordering = unknown_vert_to_known_verts_ordering[face_t]
            pos = developed_position(tet_vert_posns[tet_ordering[0]], tet_vert_posns[tet_ordering[1]], tet_vert_posns[tet_ordering[2]], tet_shape)
        
            # ancestors = [ a for a in vert_a.ladderpole_ancestors if a in vert_b.ladderpole_ancestors ]
            ancestors = vert_a.ladderpole_ancestors.intersection(vert_b.ladderpole_ancestors)

            vert_t = vertex( self, pos)

            vert_t.ladderpole_ancestors = ancestors
        else:
            vert_t = vertex( self, None)

        if self.desired_vertices != []:
            if self.infinity in triangle.vertices:
                self.check_vertex_desired(vert_t) 

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

        ### now for the edges

        edge_at = landscape_edge(self, [vert_a, vert_t], not triangle.is_upper) ## coastal
        edge_bt = landscape_edge(self, [vert_b, vert_t], triangle.is_upper) ## coastal
        edge_ct = landscape_edge(self, [vert_c, vert_t], far_edge_colour == "red") ## never coastal

        if triangle.is_red == triangle.is_upper:
            edge_bc, edge_ca, edge_ab = triangle.edges 
        else:
            edge_ca, edge_ab, edge_bc = triangle.edges

        if ab_is_red == ab_is_upper:
            triangle_a.edges = [edge_bc, edge_ct, edge_bt]
            triangle_b.edges = [edge_at, edge_ct, edge_ca]
        else:
            triangle_a.edges = [edge_bt, edge_bc, edge_ct]
            triangle_b.edges = [edge_ca, edge_at, edge_ct]

        if c_is_red != c_is_upper:
            triangle_c.edges = [edge_at, edge_bt, edge_ab]
        else:
            triangle_c.edges = [edge_bt, edge_ab, edge_at]

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
        self.num_tetrahedra += 1
        return con_tet ### pass back the tetrahedron we added

    def mark_ladderpole_descendants(self, ladderpole_descendant_segments):
        for i, v in enumerate(self.coast):
            for j, segment in enumerate(ladderpole_descendant_segments):
                if segment[0] <= i <= segment[1]:
                    v.ladderpole_ancestors.add(j)

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

    def calculate_max_ladderpole_descendant_coast_edge_length(self):
        m = 0.0
        for i,v in enumerate(self.coast):
            w = self.coast[(i+1)%len(self.coast)]
            if len(v.ladderpole_ancestors.intersection(w.ladderpole_ancestors)) > 0:
                edge_length = abs( v.pos.complex() - w.pos.complex() )
                if edge_length > m:
                    m = edge_length
        return m

    def count_long_edges(self):
        out = 0
        for i,v in enumerate(self.coast):
            w = self.coast[(i+1)%len(self.coast)]
            if len(v.ladderpole_ancestors.intersection(w.ladderpole_ancestors)) > 0:
                edge_length = abs( v.pos.complex() - w.pos.complex() )
                if edge_length > self.max_length:
                    out += 1
        return out

    def lightning_dividers(self, special_vertices):

        ## for every special vertex v, do the following for upper and lower landscapes. For every river R with source v,
        ## follow it to the coast, let e be the mouth. On the lower landscape, we define a set of edges, the "dividers", for the pair
        ## (v,e). These are the edges that link (v,e) in the circular ordering. We return all lists of dividers, for upper/lower, and for each river.

        ## remark: there are other lightning curves for, say, edges, or more generally, for an arbitrary pair of points on S^1(alpha)
        ## A more general algorithm: lets do it for an edge on the upper landscape, say e. Take its endpoints p, q.
        ## Let e_i be the edges on the lower landscape that link e. We calculate this from the coast.
        ## We then put these in order: we say that e_i < e_j if e_i separates p from e_j.
        ## Then clean up as in the following function.

        ## first, find rivers starting from infinity on the upper landscape
        dividers = [[],[]]  ### first list is upper, second list is lower
        for tri in self.triangles:
            if not tri.is_buried:
                river_source = tri.vertices[tri.downriver_index()]
                if river_source == self.infinity or river_source in special_vertices: 
                    river, got_to_coast = self.flow(tri, return_all = True)
                    assert got_to_coast
                    river_edges = [tri.edges[tri.downriver_index()] for tri in river]

                    if tri.is_upper:                        
                        ### now push river down
                        ### deal with tri next to river_source first. It must have its track pointing to a 0 angle 
                        ### edge of the tet below because upper track cusps dont point towards the top
                        ### edge of a tetrahedron. So push down, repeat until cannot push down any more
                        
                        while river[0].lower_tet != None:
                            new_river_tris = river[0].lower_tet.lower_triangles
                            if river_source in new_river_tris[1].vertices:
                                new_river_tris.reverse()
                            assert river_source not in new_river_tris[1].vertices
                            river = new_river_tris + river[1:]
                            river_edges = [river[0].shared_edge(river[1])] + river_edges

                        ### now we only have to worry about pushing down through tetrahedra not incident to river_source
                        while True:  ## this is quadratic time, could be linear... but who cares
                            restart_while_loop = False
                            for i in range(1, len(river)):
                                the_lower_tet = river[i].lower_tet
                                if the_lower_tet != None:
                                    upper_tris = the_lower_tet.upper_triangles
                                    other_tri = upper_tris[(upper_tris.index(river[i]) + 1) %2]
                                    a, b = river_edges[i-1], river_edges[i]
                                    assert a in river[i].edges and b in river[i].edges 
                                    if a not in other_tri.edges and b not in other_tri.edges: 
                                        ## so the path doesnt touch the other upper triangle,
                                        ## so we can push down through a normal triangle
                                        lower_tris = the_lower_tet.lower_triangles
                                        if a in lower_tris[1].edges:
                                            lower_tris.reverse()
                                        assert a in lower_tris[0].edges and b in lower_tris[1].edges
                                        river = river[:i] + lower_tris + river[i+1:]
                                        river_edges = river_edges[:i] + [the_lower_tet.lower_edge()] + river_edges[i:]
                                        restart_while_loop = True
                                        break # out of for loop
                            if restart_while_loop:
                                continue
                            for i in range(1, len(river) - 1):
                                if river[i].lower_tet != None and river[i].lower_tet == river[i+1].lower_tet: ## can push down through this tet
                                    the_lower_tet = river[i].lower_tet
                                    a, b = river_edges[i-1], river_edges[i+1]
                                    assert a in river[i].edges and b in river[i+1].edges 
                                    lower_tris = the_lower_tet.lower_triangles
                                    if a in lower_tris[1].edges:
                                        lower_tris.reverse()
                                    assert a in lower_tris[0].edges and b in lower_tris[1].edges
                                    river = river[:i] + lower_tris + river[i+2:]
                                    river_edges = river_edges[:i] + [the_lower_tet.lower_edge()] + river_edges[i+1:]
                                    restart_while_loop = True
                                    break # out of for loop
                            if restart_while_loop:
                                continue
                        # now we exit the while loop
                            break
                        ### finally, assert that there is nothing below any tri in the river at the end 
                        for tri in river:
                            assert tri.lower_tet == None
                        # curve = [e.midpoint() for e in river_edges]
                        dividers[0].append(river_edges)   ## midpoints of these are our lightning curve

                    else:  ### its a triangle on the bottom... swap all the uppers with lowers... 
                        ### there is some serious code factoring that could be done here!
                    
                        ### now push river up
                        
                        while river[0].upper_tet != None:
                            new_river_tris = river[0].upper_tet.upper_triangles
                            if river_source in new_river_tris[1].vertices:
                                new_river_tris.reverse()
                            assert river_source not in new_river_tris[1].vertices
                            river = new_river_tris + river[1:]
                            river_edges = [river[0].shared_edge(river[1])] + river_edges

                        ### now we only have to worry about pushing up through tetrahedra not incident to river_source
                        while True:
                            restart_while_loop = False
                            for i in range(1, len(river)):
                                the_upper_tet = river[i].upper_tet
                                if the_upper_tet != None:
                                    lower_tris = the_upper_tet.lower_triangles
                                    other_tri = lower_tris[(lower_tris.index(river[i]) + 1) %2]
                                    a, b = river_edges[i-1], river_edges[i]
                                    assert a in river[i].edges and b in river[i].edges 
                                    if a not in other_tri.edges and b not in other_tri.edges: 
                                        ## so the path doesnt touch the other lower triangle,
                                        ## so we can push up through a normal triangle
                                        upper_tris = the_upper_tet.upper_triangles
                                        if a in upper_tris[1].edges:
                                            upper_tris.reverse()
                                        assert a in upper_tris[0].edges and b in upper_tris[1].edges
                                        river = river[:i] + upper_tris + river[i+1:]
                                        river_edges = river_edges[:i] + [the_upper_tet.upper_edge()] + river_edges[i:]
                                        restart_while_loop = True
                                        break # out of for loop
                            if restart_while_loop:
                                continue
                            for i in range(1, len(river) - 1):
                                if river[i].upper_tet != None and river[i].upper_tet == river[i+1].upper_tet: ## can push down through this tet
                                    the_upper_tet = river[i].upper_tet
                                    a, b = river_edges[i-1], river_edges[i+1]
                                    assert a in river[i].edges and b in river[i+1].edges 
                                    upper_tris = the_upper_tet.upper_triangles
                                    if a in upper_tris[1].edges:
                                        upper_tris.reverse()
                                    assert a in upper_tris[0].edges and b in upper_tris[1].edges
                                    river = river[:i] + upper_tris + river[i+2:]
                                    river_edges = river_edges[:i] + [the_upper_tet.upper_edge()] + river_edges[i+1:]
                                    restart_while_loop = True
                                    break # out of for loop
                            if restart_while_loop:
                                continue
                        # now we exit the while loop
                            break
                        ### finally, assert that there is nothing below any tri in the river at the end 
                        for tri in river:
                            assert tri.upper_tet == None
                        # curve = [e.midpoint() for e in river_edges]
                        dividers[1].append(river_edges)   ## midpoints of these are our lightning curve
        
        ### now clean up each divider by removing "fan edges"
        # out = []
        # for divider in dividers:
        #     if len(divider) <= 2:
        #         out.append(divider)
        #         continue
        #     clean_divider = []
        #     a, b = divider[:2]
        #     x = a.shared_vertex(b)
        #     clean_divider.append(a)
        #     for c in divider[2:]:
        #         y = b.shared_vertex(c)
        #         if x == y:
        #             a, b = a, c
        #         else:
        #             clean_divider.append(b)
        #             a, b = b, c
        #             x = y
        #     clean_divider.append(divider[-1])
        #     out.append(clean_divider)

        return dividers

if __name__ == '__main__':

    veering_isosig = 'dLQacccjsnk_200'

    shapes_data = read_from_pickle('Data/veering_shapes_up_to_ten_tetrahedra.pkl')

    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle, tet_shapes = shapes_data[veering_isosig])
    # con = continent( vt )
    # con.build(5)
    # print con.coast()





