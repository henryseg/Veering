from file_io import parse_data_file, read_from_pickle
from taut import isosig_to_tri_angle
from veering_triangulation import veering_triangulation
from develop_ideal_hyperbolic_tetrahedra import developed_position, develop_verts_CP1, unknown_vert_to_known_verts_ordering, convert_to_complex

edge_vert_index_map = {(0,1):0, (0,2):1, (0,3):2, (1,2):3, (1,3):4, (2,3):5}

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
    def __init__(self, vt):
        self.vt = vt
        self.triangles = []
        self.vertices = [ [1,0],[0,1],[1,1],[vt.tet_shapes[0],1] ]
        ### make first tetrahedron triangles...

    def incident_tet(self, triangle):
        assert not triangle.is_buried

        face = self.vt.tri.face(2,triangle.index)
        embed0 = face.embedding(0)
        embed1 = face.embedding(1)
        tet0 = embed0.simplex()
        tet1 = embed1.simplex()
        if (vt.coorientations[tet0.index()][embed0.face()] == -1) == (triangle.is_upper): ## -1 means into the tetrahedron, so this is above triangle
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

        ###       3
        ###   0--------3
        ###   |`.a   ,'|     is_upper, so t and n are below, a and b are new triangles above
        ### 0 |  ` ,' n|
        ###   |t ,' .  | 2   the tet labels here orient this tetrahedron in a way INconsistent with our veering convention
        ###   |,'   b`.|
        ###   1--------2 
        ###        1

        ###       0
        ###   1--------0
        ###   |`.t   ,'|     not is upper, so t and n are above, a and b are new triangles below
        ### 1 |  ` ,' a| 
        ###   |b ,' .  | 3   the tet labels here orient this tetrahedron in a way consistent with our veering convention
        ###   |,'   n`.|
        ###   2--------3 
        ###        2

        ind_t = triangle.downriver_index()
        ind_n = neighbour.downriver_index()
        assert triangle.vertices[(ind_t - 1) % 3] == neighbour.vertices[(ind_n + 1) % 3]
        assert triangle.vertices[(ind_t + 1) % 3] == neighbour.vertices[(ind_n - 1) % 3]

        vertices = [ triangle.vertices[(ind_t - 1) % 3], 
                     triangle.vertices[ind_t],
                     triangle.vertices[(ind_n - 1) % 3], 
                     triangle.vertices[ind_n] ]

        neighbours = triangle.not_downriver() + neighbour.not_downriver()

        tet, embed = self.incident_tet(triangle)
        face_t = embed.face()

        ordering = unknown_vert_to_known_verts_ordering[face_t]  ## we must use orientation of tetrahedron to distinguish a from b
        if triangle.is_upper:
            ordering.reverse()

        for i in range(4):
            if i != face_t:
                if vt.coorientations[tet.index()][face_t] == vt.coorientations[tet.index()][i]:
                    face_n = i
                else:
                    ind = ordering.index(face_n)
                    face_a = ordering[(ind + 1) % 3]
                    face_b = ordering[(ind - 1) % 3]

        face_a_index = tet.face(2,face_a).index()
        face_b_index = tet.face(2,face_b).index()

        far_edge_vert_nums = [face_a, face_b]
        far_edge_vert_nums.sort()
        far_edge_num = edge_vert_index_map[tuple(far_edge_vert_nums)]
        far_edge = tet.edge(far_edge_num)
        is_red = ( vt.veering_colours[far_edge.index()] == 'R' )

        is_upper = triangle.is_upper

        if is_red:
            order_a = [3,0,1]
            order_b = [1,2,3]
        else:
            order_a = [1,3,0]
            order_b = [3,1,2]

        vertices_a = [vertices[i] for i in order_a]
        neighbours_a = [neighbours[i] for i in order_a]
        vertices_b = [vertices[i] for i in order_b]
        neighbours_b = [neighbours[i] for i in order_b]       

        triangle_a = landscape_triangle(face_a_index, is_upper, is_red, vertices_a, neighbours_a)
        triangle_b = landscape_triangle(face_b_index, is_upper, is_red, vertices_b, neighbours_b)

        neighbours[0].update_contacts(triangle, triangle_a)
        neighbours[1].update_contacts(triangle, triangle_b)
        neighbours[2].update_contacts(neighbour, triangle_b)
        neighbours[3].update_contacts(neighbour, triangle_a)

        self.triangles.extend([triangle_a, triangle_b])
        triangle.is_buried = True
        neighbour.is_buried = True

    def coastal_fill(self, triangle):
        pass

