from file_io import parse_data_file, read_from_pickle
from taut import isosig_to_tri_angle
from veering import veering_triangulation
from develop_ideal_hyperbolic_tetrahedra import developed_position, convert_to_complex, unknown_vert_to_known_verts_ordering
from basic_math import sign

class vertex:
    def __init__(self, CP1):
        self.CP1 = CP1
        ### maybe an age
    def __repr__(self):
        if self.CP1 != [1,0]:
            return str(convert_to_complex(self.CP1))
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
        self.is_buried = False

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

class continent:
    def __init__(self, vt, initial_tet_num = 0):
        # print 'initializing continent'
        self.vt = vt
        self.triangles = []
        self.infinity = vertex([1,0])
        self.vertices = [ self.infinity, vertex([0,1]), vertex([1,1]), vertex([self.vt.tet_shapes[initial_tet_num],1]) ]
        self.vertices_adjacent_to_infinity = []   ## list of (vertex, colour of edge from vertex to infinity)

        self.num_tetrahedra = 1

        for i in range(1,4):
            col = self.vt.get_edge_between_verts_colour(initial_tet_num, (0, i))
            self.vertices_adjacent_to_infinity.append( (self.vertices[i], col) )

        ###   c---R----b
        ###   |`d    ,'|     faces a, b on bottom, c, d on top
        ###   L  ` ,' a| 
        ###   |b ,' .  L 
        ###   |,'    c.| 
        ###   a----R---d 

        upper_face_nums = []
        lower_face_nums = []
        for i in range(4):
            if self.vt.coorientations[initial_tet_num][i] == +1:
                upper_face_nums.append(i)
            else:
                lower_face_nums.append(i)

        face_a, face_b = lower_face_nums ## decree that face_a is the smaller index
        if sign([face_a, upper_face_nums[1], face_b, upper_face_nums[0]]) == 1: ## get orientation correct
            face_c, face_d = upper_face_nums
        else:
            face_d, face_c = upper_face_nums

        tet = vt.tri.tetrahedron(initial_tet_num)
        face_a_index = tet.face(2,face_a).index()
        face_b_index = tet.face(2,face_b).index()
        face_c_index = tet.face(2,face_c).index()
        face_d_index = tet.face(2,face_d).index()

        upper_edge_colour = self.vt.get_edge_between_verts_colour(initial_tet_num, lower_face_nums)
        lower_edge_colour = self.vt.get_edge_between_verts_colour(initial_tet_num, upper_face_nums)

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

    def sanity_check(self):
        for tri in self.triangles:
            tri.check_against_neighbours()

    def bury(self, triangle):
        while not triangle.is_buried:
            self.silt(triangle)
            
    def silt(self, triangle):
        """flow downriver from this triangle until we hit either a sink, or the coast, add one tetrahedron there"""
        while True:
            neighbour = triangle.downriver()
            if neighbour.is_upper != triangle.is_upper:
                self.coastal_fill(triangle)
                break
            else:
                neighbour2 = neighbour.downriver()
                if triangle == neighbour2:
                    self.in_fill(triangle)
                    break
                else:
                    triangle = neighbour
        self.num_tetrahedra += 1
        # self.sanity_check()

    def build(self, max_num_tetrahedra):
        first_non_buried_index = 0
        while self.num_tetrahedra < max_num_tetrahedra:  # will go a little over because we check after each bury, which adds many tetrahedra
            tri = self.triangles[first_non_buried_index]
            self.bury(tri)
            first_non_buried_index += 1
            while self.triangles[first_non_buried_index].is_buried:
            # while self.triangles[first_non_buried_index].is_buried or self.triangles[first_non_buried_index].is_upper:
                first_non_buried_index += 1

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
        tet_vert_posns[face_a] = vert_a.CP1
        tet_vert_posns[face_b] = vert_b.CP1
        tet_vert_posns[face_c] = vert_c.CP1

        ### next: permute the triangle verts in CP1 into tet order. Then plug through tet_ordering so we can develop
        tet_shape = self.vt.tet_shapes[tet.index()]
        # print tet_shape
        tet_ordering = unknown_vert_to_known_verts_ordering[face_t]
        vert_t = vertex( developed_position(tet_vert_posns[tet_ordering[0]], tet_vert_posns[tet_ordering[1]], tet_vert_posns[tet_ordering[2]], tet_shape) )
        # print [ vertex(p) for p in [tet_vert_posns[tet_ordering[0]], tet_vert_posns[tet_ordering[1]], tet_vert_posns[tet_ordering[2]]] ]
        # print vert_t
        self.vertices.append(vert_t)
        if self.infinity in triangle.vertices:
            which_one = [vert_a, vert_b, vert_c].index(self.infinity)
            face_inf = [face_a, face_b, face_c][which_one]
            col = self.vt.get_edge_between_verts_colour(tet.index(), (face_t, face_inf))
            self.vertices_adjacent_to_infinity.append( (vert_t, col) )

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

    def equator(self):
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

        ## now rotate to put infinity first
        inf_vert_index = out.index( self.vertices[0] )
        out = out[inf_vert_index:] + out[:inf_vert_index]
        return out

    def edges_adjacent_to_infinity(self):
        ## list of ( (u,v), colour ) where u, v are vertices at either end of an edge opposite infinity in a triangle
        out = []
        for triangle in self.triangles:
            if self.infinity in triangle.vertices:
                ind = triangle.vertices.index(self.infinity)
                endpoints = ( triangle.vertices[(ind+1)%3], triangle.vertices[(ind+2)%3] )
                veering_colour = self.vt.veering_colours[ triangle.edge_indices()[ind] ]
                out.append( (endpoints, veering_colour) ) 
        return out

if __name__ == '__main__':

    veering_isosig = 'dLQacccjsnk_200'

    shapes_data = read_from_pickle('Data/veering_shapes_up_to_ten_tetrahedra.pkl')

    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle, tet_shapes = shapes_data[veering_isosig])
    con = continent( vt )
    con.build(5)
    print con.equator()





