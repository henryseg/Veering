from file_io import parse_data_file, read_from_pickle
from taut import isosig_to_tri_angle
from veering import veering_triangulation
from develop_ideal_hyperbolic_tetrahedra import developed_position, develop_verts_CP1, convert_to_complex, unknown_vert_to_known_verts_ordering
from basic_math import sign

class vertex:
    def __init__(self, CP1_pos):
        self.CP1 = CP1_pos
        ### maybe an age

class landscape_triangle:
    def __init__(self, face_index, is_upper, is_red, vertices, neighbours):
        self.index = face_index ## in the quotient manifold
        self.is_upper = is_upper ## true if we are on the upper landscape of the continent (or were, before we got buried) 
        self.is_red = is_red ## if self has two red edges
        self.vertices = vertices  ## list of three vertices, oriented anticlockwise as viewed from above, with the special vertex first
        ### the special vertex is incident to two edges of the same colour
        self.neighbours = neighbours ## list of three triangles incident to this one, opposite corresponding vertices
        self.is_buried = False

    def downriver_index(self):
        """index of triangle in self.neighbours that is downriver of self"""
        if is_upper == is_red: return 2
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

class continent:
    def __init__(self, vt, initial_tet_num = 0):
        self.vt = vt
        self.triangles = []
        self.vertices = [ [1,0],[0,1],[1,1],[self.vt.tet_shapes[initial_tet_num],1] ]
     

        ###   c---R----b
        ###   |`d    ,'|     faces a, b on bottom, c, d on top
        ###   L  ` ,' a| 
        ###   |b ,' .  L     
        ###   |,'    c.|
        ###   a----R---d  

        #### HERE ... go look at triangles creation in coastal_fill

        upper_face_nums = []
        lower_face_nums = []
        for i in range(4):
            if self.vt.coorientations[initial_tet_num][i] == +1:
                upper_face_nums.append(i)
            else:
                lower_face_nums.append(i)
        upper_edge_colour = self.vt.get_edge_between_verts_colour(initial_tet_num, lower_face_nums)
        lower_edge_colour = self.vt.get_edge_between_verts_colour(initial_tet_num, upper_face_nums)
        face_a, face_c = lower_face_nums ## decree that face_a is the smaller index
        if sign([face_a, upper_face_nums[0], face_c, upper_face_nums[1]]) == 1: ## get orientation correct
            face_b, face_d = upper_face_nums
        else:
            face_d, face_b = upper_face_nums

    def incident_tet(self, triangle):
        assert not triangle.is_buried

        face = self.vt.tri.face(2,triangle.index)
        embed0 = face.embedding(0)
        embed1 = face.embedding(1)
        tet0 = embed0.simplex()
        tet1 = embed1.simplex()
        if (self.vt.coorientations[tet0.index()][embed0.face()] == -1) == (triangle.is_upper): ## -1 means into the tetrahedron, so this is above triangle
            return (tet0, embed0)
        else:
            return (tet1, embed1)

    def bury(self, triangle):
        while not triangle.is_buried:
            self.silt(triangle)
            
    def silt(self, triangle):
        """flow downriver from this triangle until we hit either a sink, or the coast, add one tetrahedron there"""
        while True:
            neighbour = triangle.downriver()
            if neighbour.is_upper != triangle.is_upper:
                self.coastal_fill(triangle)
            else:
                neighbour2 = neighbour.downriver()
                if triangle == neighbour2:
                    self.in_fill(triangle)
                else:
                    triangle = neighbour

    def in_fill(self, triangle):
        neighbour = triangle.downriver()
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
        tet, embed = self.incident_tet(triangle)  
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

        far_edge_colour = self.vt.get_edge_between_verts_colour(tet.index(), (face_a, face_b))
        ab_is_red = ( far_edge_colour == 'R' )
        ab_is_upper = triangle.is_upper

        triangle_a = landscape_triangle(face_a_index, ab_is_upper, ab_is_red, None, None)
        triangle_b = landscape_triangle(face_b_index, ab_is_upper, ab_is_red, None, None)
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
        neighbour_australian = triangle.downriver() ## it's upside down from us
        assert triangle.is_upper != neighbour_australian.is_upper
        assert triangle in neighbour_australian.neighbours ## don't know which one
        assert triangle != neighbour_australian.downriver() ## although we know its not this one

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
        tet, embed = self.incident_tet(triangle)  
        face_t = embed.face()

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

        face_a_index = tet.face(2,face_a).index()
        face_b_index = tet.face(2,face_b).index()
        face_c_index = tet.face(2,face_c).index()

        far_edge_colour = self.vt.get_edge_between_verts_colour(tet.index(), (face_a, face_b))
        ab_is_red = ( far_edge_colour == 'R' )
        c_is_red = triangle.is_red
        ab_is_upper = triangle.is_upper
        c_is_upper = not triangle.is_upper

        triangle_a = landscape_triangle(face_a_index, ab_is_upper, ab_is_red, None, None)
        triangle_b = landscape_triangle(face_b_index, ab_is_upper, ab_is_red, None, None)
        triangle_c = landscape_triangle(face_c_index,  c_is_upper,  c_is_red, None, None)

        self.triangles.extend([triangle_a, triangle_b, triangle_c])

        ## now for the vertices

        ### let's find out vert_a, vert_b, vert_c, which are the vertices (with CP1 data) opposite faces.
        if triangle.is_red == triangle.is_upper:
            vert_a, vert_b, vert_c = triangle.vertices
        else:
            vert_b, vert_c, vert_a = triangle.vertices

        tet_vert_posns = [None, None, None, None]
        tet_vert_posns[face_a] = vert_a.CP1_pos
        tet_vert_posns[face_b] = vert_b.CP1_pos
        tet_vert_posns[face_c] = vert_c.CP1_pos

        ### next: permute the triangle verts in CP1 into tet order. Then plug through tet_ordering so we can develop
        tet_shape = self.vt.tet_shapes[tet.index()]
        tet_ordering = unknown_vert_to_known_verts_ordering(face_t)
        vert_t = vertex( developed_position(tet_vertex_posns[tet_ordering[0]], tet_vertex_posns[tet_ordering[1]], tet_vertex_posns[tet_ordering[2]], tet_shape) )

        self.vertices.append(vert_t)

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


if __name__ == '__main__':
    pass

