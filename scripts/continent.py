#
# continent.py
#

from veering.file_io import parse_data_file, read_from_pickle
from veering.basic_math import sign
from veering.taut import isosig_to_tri_angle, vert_pair_to_edge_num, edge_num_to_vert_pair
from veering.veering_tri import veering_triangulation
from veering.fundamental_domain import spanning_dual_tree

from develop_ideal_hyperbolic_tetrahedra import developed_position, unknown_vert_to_known_verts_ordering
from circular_order import are_anticlockwise, are_linking
from sage.rings.rational_field import QQ 


class vertex:
    def __init__(self, continent, pos):
        self.continent = continent
        self.edges = []
        self.triangles_with_special_vertex_here = [] ### the edges of the triangle meeting here have the same colour
        self.triangles_without_special_vertex_here = [] ### the other two kinds of triangle corners
        self.tetrahedra = []
        self.pos = pos  ### complex position, not used for circular pictures
        self.ladderpole_ancestors = set() ## which ladderpole edges did you come from
        self.continent.vertices.append(self)
        self.circle_pos = None
        self.cusp_leaves = []
    
    def coastal_index(self):
        return self.continent.coast.index(self)

    def chronological_index(self):
        return self.continent.vertices.index(self)

    def is_ladderpole_descendant(self):
        return len(self.ladderpole_ancestors) != 0

    def edge_between(self, other):
        for e in self.edges:
            if other in e.vertices:
                return e  ### note, no parallel edges in univ cover of a veering triangulation
        assert False

    def __repr__(self):
        if self.pos == None:
            return 'chron' + str(self.chronological_index())
        elif not self.pos.is_infinity():
            return str(self.pos.complex())
        else:
            return "inf"
    
    def __str__(self):
        return self.__repr__()

class landscape_edge:
    def __init__(self, continent, vertices, edge_index, is_red): 
        self.continent = continent
        self.vertices = vertices
        for v in self.vertices:
            v.edges.append(self)
        self.triangles = []
        self.upper_tet = None
        self.lower_tet = None
        self.side_tetrahedra = []
        self.other_end = {self.vertices[0]:self.vertices[1], self.vertices[1]:self.vertices[0]}
        self.continent.edges.append(self)
        self.index = edge_index
        self.is_red = is_red
        self.cusp_leaves = [] ### should be empty if not a coastal edge
        self.rectangle_sides_data = [None, None, None, None]
        # try:
        #     assert self.length() > 0.0001
        # except:
        #     print 'edge too short', self.vertices
        #     raise

    def __repr__(self):
        u, v = self.vertices
        # return ' '.join( [str(self.continent.edges.index(self)), 'edge', str(u), str(v), str(self.length())] )
        return '_'.join( [str(self.continent.edges.index(self)), 'edge', str(u), str(v)] )

    def coastal_index(self):
        return self.continent.coastal_edges.index(self)

    def end_positions(self):
        return [v.coastal_index() for v in self.vertices]

    def links(self, other):
        """Returns True if this edge links other (could be leaf or edge)"""
        a1, a2 = self.end_positions()
        b1, b2 = other.end_positions()
        return are_linking(a1, a2, b1, b2)

    def is_coastal(self):
        return self in self.continent.coastal_edges

    def boundary_triangles(self):
        out = []
        for tri in self.triangles:
            if tri.is_boundary():
                out.append(tri)
        assert len(out) == 2 or len(out) == 0
        return out

    def rectangle_sides(self):
        if not None in self.rectangle_sides_data:
            return self.rectangle_sides_data
        ### otherwise update things
        #     -3-v
        #    |  /|
        #    0 e 2
        #    |/  |  ### e is red <=> 1, 3 are purple and 0, 2 are green
        #    u-1-

        u, v = self.vertices
        leaves = u.cusp_leaves
        first_leaf = leaves[0]
        last_leaf = leaves[-1]
        if are_anticlockwise(u.coastal_index(), first_leaf.end_positions()[1], v.coastal_index()):
            self.rectangle_sides_data[1] = first_leaf
        elif are_anticlockwise(u.coastal_index(), v.coastal_index(), last_leaf.end_positions()[1]):
            self.rectangle_sides_data[0] = last_leaf
        else:      
            for i, leaf in enumerate(leaves):
                if are_anticlockwise(u.coastal_index(), leaf.end_positions()[1], v.coastal_index()):
                    self.rectangle_sides_data[0] = leaves[i-1]
                    self.rectangle_sides_data[1] = leaf
                    break
        leaves = v.cusp_leaves
        first_leaf = leaves[0]
        last_leaf = leaves[-1]
        if are_anticlockwise(v.coastal_index(), first_leaf.end_positions()[1], u.coastal_index()):
            self.rectangle_sides_data[3] = first_leaf
        elif are_anticlockwise(v.coastal_index(), u.coastal_index(), last_leaf.end_positions()[1]):
            self.rectangle_sides_data[2] = last_leaf
        else:      
            for i, leaf in enumerate(leaves):
                if are_anticlockwise(v.coastal_index(), leaf.end_positions()[1], u.coastal_index()):
                    self.rectangle_sides_data[2] = leaves[i-1]
                    self.rectangle_sides_data[3] = leaf
                    break
        if not None in self.rectangle_sides_data:
            a, b, c, d = self.rectangle_sides_data
            # if not (a.links(d) and b.links(c)):
            #     print('rectangle sides linking failed')
            #     print(self.vertices)
            assert a.links(d) and b.links(c)
            assert a.is_upper == c.is_upper
            assert b.is_upper == d.is_upper
            assert a.is_upper != b.is_upper
            assert self.is_red == a.is_upper
        return self.rectangle_sides_data

    def green_purple_rectangle_sides(self):
        if None in self.rectangle_sides_data: ### either the continent doesnt have the leaves or the list needs updating
            self.rectangle_sides() ### update 
        assert not None in self.rectangle_sides_data ### if still not working then should run ensure_continent_contains_rectangle first
        a, b, c, d = self.rectangle_sides_data
        if self.is_red:
            return [[a, c], [b, d]]
        else:
            return [[b, d], [a, c]]

    def rectangle_sides_by_vertex(self):
        if None in self.rectangle_sides_data: ### either the continent doesnt have the leaves or the list needs updating
            self.rectangle_sides() ### update 
        assert not None in self.rectangle_sides_data ### if still not working then should run ensure_continent_contains_rectangle first
        a, b, c, d = self.rectangle_sides_data
        u, v = self.vertices
        return {u : [a, b], v : [c, d]}

    def ensure_continent_contains_rectangle(self):
        """Add to the continent to ensure that this edge has all its cusp leaves in the continent"""
        while None in self.rectangle_sides():
            # print('edge has this many triangles', len(self.triangles))
            found_a_tri_to_bury = False
            for tri in self.triangles:
                # print('tri.index', tri.index)
                if not tri.is_buried():
                    found_a_tri_to_bury = True
                    # print('bury triangle with index', tri.index)
                    self.continent.bury(tri)
                    assert tri.is_buried()
                    break
            assert found_a_tri_to_bury ### We checked all triangles and all are buried so we should be done
        assert not None in self.rectangle_sides()

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

class cusp_leaf:
    def __init__(self, continent, cusp, coastal_edge, is_upper):
        self.continent = continent
        self.cusp = cusp
        self.coastal_edge = coastal_edge
        self.is_upper = is_upper 

    def end_positions(self):
        start = self.cusp.coastal_index()
        index_on_edge = self.coastal_edge.cusp_leaves.index(self)
        end = self.coastal_edge.coastal_index() + QQ( ((index_on_edge) + 1)/(len(self.coastal_edge.cusp_leaves) + 1) )
        return (start, end)

    def links(self, other):
        """Returns True if this cusp leaf links other (could be leaf or edge)"""
        a1, a2 = self.end_positions()
        b1, b2 = other.end_positions()
        return are_linking(a1, a2, b1, b2)

        ### have faces know about the cusp leaves that start there, and vice versa?
        ### leave until we need it

    def sees_to_its_left(self, cusp): 
        ### if the cusp is to the left of the oriented leaf, return true
        ### In particular, a leaf does _not_ see its own cusp to its left 

        start, end = self.end_positions()
        return are_anticlockwise(start, end, cusp.coastal_index())

    # def is_entirely_to_one_side_of(self, tet): ### are all of the cusps of tet on one or the other side of self?
    #     if self.cusp in tet.vertices:
    #         return False
    #     are_to_left = [self.sees_to_its_left(v) for v in tet.vertices]
    #     return are_to_left.count(True) % 4 == 0 ### True if all are on one side or the other

    def separates(self, cusp_list_1, cusp_list_2): ### returns True if and only if all of cusp_list_1 is on one side of self, and all of cusp_list_2 is on the other side
        if self.cusp in cusp_list_1 or self.cusp in cusp_list_2:
            return False
        are_to_left_1 = [self.sees_to_its_left(v) for v in cusp_list_1]
        are_to_left_2 = [self.sees_to_its_left(v) for v in cusp_list_2]
        if all(are_to_left_1) and not any(are_to_left_2):
            return True
        if all(are_to_left_2) and not any(are_to_left_1):
            return True   
        return False

    def weakly_separates(self, cusp_list_1, cusp_list_2): 
        ### like separates, but ignores if self.cusp is in either list
        cusp_list_1b = [c for c in cusp_list_1 if c != self.cusp]
        cusp_list_2b = [c for c in cusp_list_2 if c != self.cusp]
        return self.separates(cusp_list_1b, cusp_list_2b)
        
    def is_entirely_to_one_side_of(self, tet):
        return self.separates(tet.vertices, [])

class landscape_triangle:
    def __init__(self, continent, face_index, is_upper, is_red):
        self.continent = continent
        self.index = face_index ## in the quotient manifold
        self.is_upper = is_upper ## true if we are on the upper landscape of the continent (or were, before we got buried) 
        self.is_red = is_red ## if self has two red edges
        self.vertices = []  ## list of three vertices, oriented anticlockwise as viewed from above, with the special vertex first
        ### the special vertex is incident to two edges of the same colour

        ###    *                *
        ###   / \              / \
        ###  R   R     or     B   B  
        ### /     \          /     \
        ### ---B---          ---R---

        self.edges = [] ## opposite corresponding vertices
        self.upper_tet = None
        self.lower_tet = None
        self.neighbours = [] ## list of three triangles incident to this one, opposite corresponding vertices
        self.continent.triangles.append(self)

    def __str__(self):
        return 'continent_ind,triang_ind,upper,red,vertices,buried ' + str([self.continent.triangles.index(self), self.index, self.is_upper, self.is_red, self.vertices, self.is_buried()])

    def update_vertices(self, vertices):
        self.vertices = vertices
        vertices[0].triangles_with_special_vertex_here.append(self)
        vertices[1].triangles_without_special_vertex_here.append(self)
        vertices[2].triangles_without_special_vertex_here.append(self)

    def update_edges(self, edges):
        self.edges = edges
        for e in self.edges:
            e.triangles.append(self)
        for i in range(3):  ### edge i should be opposite vertex i
            assert not self.vertices[i] in self.edges[i].vertices
            assert self.vertices[(i+1)%3] in self.edges[i].vertices
            assert self.vertices[(i+2)%3] in self.edges[i].vertices

    def is_buried(self):
        return self.upper_tet != None and self.lower_tet != None

    def is_boundary(self): ### for convenience
        return not self.is_buried()

    def set_upper_tet(self, con_tet):
        self.upper_tet = con_tet
        con_tet.lower_triangles.append(self)
        assert len(con_tet.lower_triangles) <= 2

    def set_lower_tet(self, con_tet):
        self.lower_tet = con_tet
        con_tet.upper_triangles.append(self)
        assert len(con_tet.upper_triangles) <= 2

    def inwards_con_tet(self):
        assert self.is_boundary()
        if self.is_upper:
            return self.lower_tet
        else:
            return self.upper_tet

    def outwards_con_tet(self):
        assert self.is_buried()
        if self.is_upper:
            return self.upper_tet
        else:
            return self.lower_tet

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

    def outwards_tet(self, return_other_tet = False, return_other_tet_and_embed = False): ### for the universal cover tet outside the continent, next to this triangle, get the corresponding quotient tet
        face = self.continent.vt.tri.face(2,self.index)
        embed0 = face.embedding(0)
        embed1 = face.embedding(1)
        tet0 = embed0.simplex()
        tet1 = embed1.simplex()
        if (self.continent.vt.coorientations[tet0.index()][embed0.face()] == -1) == (self.is_upper): ## -1 means into the tetrahedron, so this is above triangle
            if return_other_tet_and_embed:
                return (tet0, embed0, tet1, embed1)
            elif return_other_tet:
                return (tet0, embed0, tet1)
            else:
                return (tet0, embed0)
        else:
            if return_other_tet_and_embed:
                return (tet1, embed1, tet0, embed0)
            elif return_other_tet:
                return (tet1, embed1, tet0)
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
        if self.is_boundary():

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
        assert self.is_boundary()
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
        assert self.is_boundary()  
        tri, is_coastal = self.continent.flow(self)
        return tri.edges[tri.downriver_index()]

    def all_river_crossing_edges(self):
        assert self.is_boundary()  
        tris, is_coastal = self.continent.flow(self, return_all = True)
        return [tri.edges[tri.downriver_index()] for tri in tris]

    def shared_edge(self, other):
        intersection = set(self.edges) & set(other.edges) 
        assert len(intersection) == 1
        return intersection.pop()

    def tet_on_other_side(self, tet):
        if tet == self.upper_tet:
            return self.lower_tet
        if tet == self.lower_tet:
            return self.upper_tet
        assert False

    def get_non_special_vertex_cusp_leaves(self): ### returns dict from non-special vertices to the [internal, boundary] leaves for self at the vertex
        out = {}
        for v in self.vertices[1:]:
            v_edges = [e for e in self.edges if v in e.vertices]
            assert len(v_edges) == 2
            for e in v_edges:
                e.ensure_continent_contains_rectangle()
            rect_sides_at_v = [e.rectangle_sides_by_vertex()[v] for e in v_edges]
            internal_leaf_at_v = list(set(rect_sides_at_v[0]).intersection(set(rect_sides_at_v[1]))) 
            assert len(internal_leaf_at_v) == 1
            boundary_leaves_at_v = list(set(rect_sides_at_v[0]).symmetric_difference(set(rect_sides_at_v[1])))
            assert len(boundary_leaves_at_v) == 2
            out[v] = {'internal':internal_leaf_at_v[0], 'boundary':boundary_leaves_at_v}
        return out


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
    def __init__(self, continent, tet_index, came_from):
        self.continent = continent
        self.index = tet_index ## in the quotient manifold
        self.coorientations = continent.vt.coorientations[tet_index]
        self.upper_triangles = []  ### this gets set when we set this tet as lower_tet for some triangle. No particular meaning to the order of this list.
        self.lower_triangles = []
        self.upper_edge = None
        self.lower_edge = None
        self.equatorial_edges = []
        self.continent.tetrahedra.append(self)
        self.vertices = [None, None, None, None] ## ordered according to the indices of the vertices downstairs in the manifold
        self.gluings = [None, None, None, None] ## a gluing specifies another tetrahedron and the Perm4 from the downstairs manifold
        self.came_from = came_from ### the landscape_triangle that we were glued onto. Initial continent_tetrahedron has None

    def set_upper_edge(self, e):
        self.upper_edge = e
        e.lower_tet = self

    def set_lower_edge(self, e):
        self.lower_edge = e
        e.upper_tet = self

    def set_equatorial_edges(self, edges):
        self.equatorial_edges = edges
        for e in edges:
            e.side_tetrahedra.append(self)

    def set_vertex(self, vert_num, vert):
        self.vertices[vert_num] = vert
        vert.tetrahedra.append(self)

    # def upper_edge(self):
    #     return self.upper_triangles[0].shared_edge(self.upper_triangles[1])

    # def lower_edge(self):
    #     return self.lower_triangles[0].shared_edge(self.lower_triangles[1])

    def edges(self):
        return self.equatorial_edges + [self.upper_edge, self.lower_edge]
        # out = self.upper_triangles[0].edges[:] ### copy of list
        # out.remove(self.upper_edge())
        # return out + self.upper_triangles[1].edges + [self.lower_edge()]

    # def vertices(self):
    #     out = set([])
    #     for tri in self.upper_triangles:
    #         out = out.union(set(tri.vertices))
    #     return out

    def ordered_faces(self):
        out = []
        for i in range(4):
            face_index = self.continent.vt.tri.tetrahedron(self.index).triangle(i).index()
            if self.continent.vt.coorientations[self.index][i] == +1: ## is upper triangle
                for triangle in self.upper_triangles:
                    if triangle.index == face_index:
                        out.append(triangle)
                        break
            else:
                for triangle in self.lower_triangles:
                    if triangle.index == face_index:
                        out.append(triangle)
                        break
        # print(out)
        assert len(out) == 4
        return out

    def edge(self, ind):
        """given index in this tetrahedron of the edges downstairs in the manifold, return the corresponding continent edge"""
        faces = self.ordered_faces()
        vert_pair = edge_num_to_vert_pair[ind]
        face_pair = [0,1,2,3]
        face_pair.remove(vert_pair[0])
        face_pair.remove(vert_pair[1])
        return faces[face_pair[0]].shared_edge(faces[face_pair[1]])

    def ordered_edges(self):
        out = []
        for i in range(6):
            out.append(self.edge(i))
        return out

    # def ordered_vertices(self):
    #     out = []
    #     tet_faces = self.ordered_faces()
    #     tet_vertices = self.vertices
    #     for i in range(4):
    #         face_vertices = tet_faces[i].vertices
    #         diff = set(tet_vertices) - set(face_vertices)
    #         assert len(diff) == 1
    #         vert = diff.pop()
    #         out.append(vert)
    #     return out
        
    def path_to_init_tet(self):
        """Returns list of landscape_triangles which connect self back to the initial tet"""
        path = []
        current_tet = self
        while current_tet.came_from != None:
            path.append(current_tet.came_from)
            current_tet = current_tet.came_from.tet_on_other_side(current_tet)
        return path

    def path_to_other_tet(self, other_tet):
        """Returns list of landscape_triangles which connect self to other"""
        path1 = self.path_to_init_tet()
        path2 = other_tet.path_to_init_tet()
        ### Now trim off the parts of the path they agree on
        path1.reverse()
        path2.reverse()
        if len(path1) > 0 and len(path2) > 0:
            i = 0
            while path1[i] == path2[i]:
                i += 1
                if len(path1) <= i or len(path2) <= i:
                    break
            path1 = path1[i:]
            path2 = path2[i:]

        path1.reverse()        
        return path1 + path2

    def face_num_path_to_other_tet(self, other_tet):
        path = self.path_to_other_tet(other_tet)
        face_num_path = []
        current_tet = self
        for triangle in path:
            face_num_path.append( current_tet.ordered_faces().index(triangle) )
            current_tet = triangle.tet_on_other_side(current_tet)
        assert current_tet == other_tet
        return face_num_path

    def follow_face_num_path(self, path):
        """From self, follow the path of face nums, return the continent tet we get to. Will build more continent if necessary"""
        current_tet = self
        for face_num in path:
            triangle = current_tet.ordered_faces()[face_num]
            if triangle.is_boundary(): ### then build more continent so we can go to the other side of face
                current_tet = self.continent.bury(triangle) ## bury the face and return the tet on the other side
            else:
                current_tet = triangle.tet_on_other_side(current_tet)
        return current_tet



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

        self.upper_landscape_triangles = set([]) 
        self.lower_landscape_triangles = set([]) 
        self.upper_landscape_edges = set([])
        self.lower_landscape_edges = set([])
        self.coastal_edges = None

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
        face_a_index = tet.face(2,face_a).index() ### index for the triangle, face_a is the index for the tet
        face_b_index = tet.face(2,face_b).index()
        face_c_index = tet.face(2,face_c).index()
        face_d_index = tet.face(2,face_d).index()

        upper_edge_colour = self.vt.get_edge_between_verts_colour(self.tet_face.tet_num, lower_face_nums)
        lower_edge_colour = self.vt.get_edge_between_verts_colour(self.tet_face.tet_num, upper_face_nums)

        tris_ab_are_red = ( lower_edge_colour == "red" )  ## the triangles, not the edge
        tris_ab_are_upper = False
        tris_cd_are_red = ( upper_edge_colour == "red" )
        tris_cd_are_upper = True

        triangle_a = landscape_triangle(self, face_a_index, tris_ab_are_upper, tris_ab_are_red)
        triangle_b = landscape_triangle(self, face_b_index, tris_ab_are_upper, tris_ab_are_red)
        triangle_c = landscape_triangle(self, face_c_index, tris_cd_are_upper, tris_cd_are_red)
        triangle_d = landscape_triangle(self, face_d_index, tris_cd_are_upper, tris_cd_are_red)

        ## now for the vertices

        vert_a = self.vertices[face_a]
        vert_b = self.vertices[face_b]
        vert_c = self.vertices[face_c]
        vert_d = self.vertices[face_d]

        if tris_ab_are_red:
            triangle_a.update_vertices( [vert_c, vert_d, vert_b] )
            triangle_b.update_vertices( [vert_d, vert_c, vert_a] )
        else:
            triangle_a.update_vertices( [vert_d, vert_b, vert_c] )
            triangle_b.update_vertices( [vert_c, vert_a, vert_d] )

        if tris_cd_are_red:
            triangle_c.update_vertices( [vert_a, vert_d, vert_b] )
            triangle_d.update_vertices( [vert_b, vert_c, vert_a] )
        else:
            triangle_c.update_vertices( [vert_b, vert_a, vert_d] )
            triangle_d.update_vertices( [vert_a, vert_b, vert_c] )

        ## add the edges

        ###   c---R----b
        ###   |`d    ,'|     faces a, b on bottom, c, d on top
        ###   L  ` ,' a| 
        ###   |b ,' .  L 
        ###   |,'    c.| 
        ###   a----R---d 

        edge_ab_index = tet.face(1,vert_pair_to_edge_num[(face_a, face_b)]).index()  ### face_a is the same as the vert number in the tet
        edge_ac_index = tet.face(1,vert_pair_to_edge_num[(face_a, face_c)]).index()
        edge_ad_index = tet.face(1,vert_pair_to_edge_num[(face_a, face_d)]).index()
        edge_bc_index = tet.face(1,vert_pair_to_edge_num[(face_b, face_c)]).index()
        edge_bd_index = tet.face(1,vert_pair_to_edge_num[(face_b, face_d)]).index()
        edge_cd_index = tet.face(1,vert_pair_to_edge_num[(face_c, face_d)]).index()

        edge_ab = landscape_edge(self, [vert_a, vert_b], edge_ab_index, upper_edge_colour == "red")
        edge_ac = landscape_edge(self, [vert_a, vert_c], edge_ac_index, False)
        edge_ad = landscape_edge(self, [vert_a, vert_d], edge_ad_index, True)
        edge_bc = landscape_edge(self, [vert_b, vert_c], edge_bc_index, True)
        edge_bd = landscape_edge(self, [vert_b, vert_d], edge_bd_index, False)
        edge_cd = landscape_edge(self, [vert_c, vert_d], edge_cd_index, lower_edge_colour == "red")

        ### make the coast (vertices) and coastal edges
        self.coast = [vert_a, vert_d, vert_b, vert_c]
        self.coastal_edges = [edge_ad, edge_bd, edge_bc, edge_ac]
        ## now rotate to put infinity first
        inf_vert_index = self.coast.index( self.infinity )
        self.coast = self.coast[inf_vert_index:] + self.coast[:inf_vert_index]
        self.coastal_edges = self.coastal_edges[inf_vert_index:] + self.coastal_edges[:inf_vert_index]

        if tris_ab_are_red: ## the triangles a and b, not the edge
            triangle_a.update_edges( [edge_bd, edge_bc, edge_cd] )
            triangle_b.update_edges( [edge_ac, edge_ad, edge_cd] )
        else:
            triangle_a.update_edges( [edge_bc, edge_cd, edge_bd] )
            triangle_b.update_edges( [edge_ad, edge_cd, edge_ac] )

        if tris_cd_are_red: 
            triangle_c.update_edges( [edge_bd, edge_ab, edge_ad] )
            triangle_d.update_edges( [edge_ac, edge_ab, edge_bc] )
        else:
            triangle_c.update_edges( [edge_ad, edge_bd, edge_ab] )
            triangle_d.update_edges( [edge_bc, edge_ac, edge_ab] )

        ## now for the neighbours

        if tris_ab_are_red:
            triangle_a.neighbours = [triangle_c, triangle_d, triangle_b]
            triangle_b.neighbours = [triangle_d, triangle_c, triangle_a]
        else:
            triangle_a.neighbours = [triangle_d, triangle_b, triangle_c]
            triangle_b.neighbours = [triangle_c, triangle_a, triangle_d]

        if tris_cd_are_red:
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

        ## add the tetrahedron

        con_tet = continent_tetrahedron(self, self.tet_face.tet_num, None)
        for i in range(4):
            con_tet.set_vertex(i, self.vertices[i]) ## copy of the vertices of the continent
        con_tet.set_upper_edge(edge_ab)
        con_tet.set_lower_edge(edge_cd) 
        con_tet.set_equatorial_edges([edge_ad, edge_bd, edge_bc, edge_ac])
        triangle_a.set_upper_tet(con_tet)
        triangle_b.set_upper_tet(con_tet)
        triangle_c.set_lower_tet(con_tet)
        triangle_d.set_lower_tet(con_tet)

        ## now cusp leaves
        ###   c---R----b
        ###   |`d    ,'|     faces a, b on bottom, c, d on top
        ###   L  ` ,' a| 
        ###   |b ,' .  L 
        ###   |,'    c.| 
        ###   a----R---d 

        if tris_cd_are_red: ### upper edge
            a_leaf = cusp_leaf(self, vert_a, edge_bc, True)
            edge_bc.cusp_leaves.append(a_leaf)
            b_leaf = cusp_leaf(self, vert_b, edge_ad, True)
            edge_ad.cusp_leaves.append(b_leaf)
        else:
            a_leaf = cusp_leaf(self, vert_a, edge_bd, True)
            edge_bd.cusp_leaves.append(a_leaf)
            b_leaf = cusp_leaf(self, vert_b, edge_ac, True)
            edge_ac.cusp_leaves.append(b_leaf)
        vert_a.cusp_leaves.append(a_leaf)
        vert_b.cusp_leaves.append(b_leaf)

        if tris_ab_are_red: ### lower edge
            c_leaf = cusp_leaf(self, vert_c, edge_ad, False)
            edge_ad.cusp_leaves.append(c_leaf)
            d_leaf = cusp_leaf(self, vert_d, edge_bc, False)
            edge_bc.cusp_leaves.append(d_leaf)
        else:
            c_leaf = cusp_leaf(self, vert_c, edge_bd, False)
            edge_bd.cusp_leaves.insert(0, c_leaf)  # prepend
            d_leaf = cusp_leaf(self, vert_d, edge_ac, False)
            edge_ac.cusp_leaves.insert(0, d_leaf)
        vert_c.cusp_leaves.append(c_leaf)
        vert_d.cusp_leaves.append(d_leaf)

        # self.build_boundary_data()
        assert self.is_convex()

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
        """flow downriver from this triangle until we hit the coast, do a coastal fill there"""
        downriver_triangle, is_coastal = self.flow(triangle)
        self.coastal_fill(downriver_triangle) ## also in_fills everything possible from that coastal_fill
        assert self.is_convex()

    def bury(self, triangle):
        ### assume continent is convex before we start
        assert self.is_convex()
        assert not triangle.is_buried()
        while not triangle.is_buried():
            self.silt(triangle)
        return triangle.outwards_con_tet()

    def in_fill(self, triangle):
        # print('in fill')
        
        neighbour = triangle.downriver()
        assert not triangle.is_buried()
        assert not neighbour.is_buried()
        assert triangle == neighbour.downriver() 
        assert triangle.is_upper == neighbour.is_upper 
        assert triangle.is_red == neighbour.is_red
     
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
        tet, embed, other_tet = triangle.outwards_tet(return_other_tet = True)  
        face_t = embed.face()
        
        n_tet, n_embed, n_other_tet = neighbour.outwards_tet(return_other_tet = True)  
        face_n = n_embed.face()

        assert tet == n_tet

        other_con_tet = triangle.inwards_con_tet()
        n_other_con_tet = neighbour.inwards_con_tet()

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
        tris_ab_are_red = ( far_edge_colour == "red" )
        tris_ab_are_upper = triangle.is_upper

        triangle_a = landscape_triangle(self, face_a_index, tris_ab_are_upper, tris_ab_are_red)
        triangle_b = landscape_triangle(self, face_b_index, tris_ab_are_upper, tris_ab_are_red)

        ## add the tetrahedron 

        con_tet = continent_tetrahedron(self, tet.index(), triangle)

        ### do the other triangle as well...

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

        if tris_ab_are_red == tris_ab_are_upper:
            triangle_a.update_vertices( [vert_t, vert_b, vert_n] )
            triangle_b.update_vertices( [vert_n, vert_a, vert_t] )
        else:
            triangle_a.update_vertices( [vert_n, vert_t, vert_b] )
            triangle_b.update_vertices( [vert_t, vert_n, vert_a] )
        
        ## add gluings and vertices to con_tet

        gluing = tet.adjacentGluing(face_t)
        con_tet.gluings[face_t] = (other_tet, gluing)
        other_face_t = gluing[face_t]
        assert other_con_tet.gluings[other_face_t] == None ### make sure it is not already glued to something
        other_con_tet.gluings[other_face_t] = (tet, other_tet.adjacentGluing(other_face_t))

        n_gluing = n_tet.adjacentGluing(face_n)
        con_tet.gluings[face_n] = (n_other_tet, n_gluing)
        other_face_n = n_gluing[face_n]
        assert n_other_con_tet.gluings[other_face_n] == None ### make sure it is not already glued to something
        n_other_con_tet.gluings[other_face_n] = (n_tet, n_other_tet.adjacentGluing(other_face_n))

        for vert_num in range(4):
            if vert_num != face_t and vert_num != face_n:
                vt = other_con_tet.vertices[gluing[vert_num]]
                vn = n_other_con_tet.vertices[n_gluing[vert_num]]
                assert vt == vn
                # con_tet.vertices[vert_num] = vt
                con_tet.set_vertex(vert_num, vt)
            elif vert_num == face_t:
                # con_tet.vertices[vert_num] = n_other_con_tet.vertices[n_gluing[vert_num]] 
                con_tet.set_vertex(vert_num, n_other_con_tet.vertices[n_gluing[vert_num]])
            else:
                assert vert_num == face_n
                # con_tet.vertices[vert_num] = other_con_tet.vertices[gluing[vert_num]]
                con_tet.set_vertex(vert_num, other_con_tet.vertices[gluing[vert_num]])
                
        ## now for the edges

        edge_tn_index = tet.face(1, vert_pair_to_edge_num[(face_t, face_n)]).index()
        edge_tn = landscape_edge(self, [vert_t, vert_n], edge_tn_index, far_edge_colour == "red") ## never coastal

        if triangle.is_red == triangle.is_upper:
            edge_bn, edge_na, edge_ab = triangle.edges 
            edge_at, edge_tb, edge_ba = neighbour.edges
        else:
            edge_na, edge_ab, edge_bn = triangle.edges 
            edge_tb, edge_ba, edge_at = neighbour.edges
        assert edge_ab == edge_ba

        if tris_ab_are_red == tris_ab_are_upper:
            triangle_a.update_edges( [edge_bn, edge_tn, edge_tb] )
            triangle_b.update_edges( [edge_at, edge_tn, edge_na] )
        else:
            triangle_a.update_edges( [edge_tb, edge_bn, edge_tn] )
            triangle_b.update_edges( [edge_na, edge_at, edge_tn] )

        assert edge_tn.links(edge_ab)
        assert not edge_tn.links(edge_bn)
        assert not edge_tn.links(edge_na)
        assert not edge_tn.links(edge_at)
        assert not edge_tn.links(edge_tb)

        if triangle.is_upper:
            con_tet.set_upper_edge(edge_tn)
            con_tet.set_lower_edge(edge_ab)
        else: 
            con_tet.set_upper_edge(edge_ab)
            con_tet.set_lower_edge(edge_tn)
        con_tet.set_equatorial_edges([edge_at, edge_tb, edge_bn, edge_na])

        ## now for the neighbours

        if tris_ab_are_red == tris_ab_are_upper:
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

        self.num_tetrahedra += 1
        return (triangle_a, triangle_b) ### pass back the new landscape triangles

    def coastal_fill(self, triangle):
        # print('coastal fill ' + str(self.triangles.index(triangle)) + ' triang ind ' + str(triangle.index)) 
        neighbour_australian = triangle.downriver() ## it's upside down from us
        assert triangle.is_upper != neighbour_australian.is_upper
        assert triangle in neighbour_australian.neighbours ## don't know which one
        assert not triangle.is_buried()
        assert not neighbour_australian.is_buried()

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
        tet, embed, other_tet = triangle.outwards_tet(return_other_tet = True)  
        face_t = embed.face()
        # print 'tet.index(), face_t', tet.index(), face_t

        for i in range(4):
            if i != face_t:
                if self.vt.coorientations[tet.index()][face_t] == self.vt.coorientations[tet.index()][i]:
                    face_c = i

        other_con_tet = triangle.inwards_con_tet()

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
        tris_ab_are_red = ( far_edge_colour == "red" ) ## faces, not talking about an edge
        c_is_red = triangle.is_red
        tris_ab_are_upper = triangle.is_upper
        c_is_upper = not triangle.is_upper

        triangle_a = landscape_triangle(self, face_a_index, tris_ab_are_upper, tris_ab_are_red)
        triangle_b = landscape_triangle(self, face_b_index, tris_ab_are_upper, tris_ab_are_red)
        triangle_c = landscape_triangle(self, face_c_index,  c_is_upper,  c_is_red)

        ## add the tetrahedron 

        con_tet = continent_tetrahedron(self, tet.index(), triangle)

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

        ##### insert vert_t into self.coast
        if vert_b == self.infinity:
            self.coast.append( vert_t )
        else:
            self.coast.insert( self.coast.index(vert_b), vert_t )

        if self.desired_vertices != []:
            if self.infinity in triangle.vertices:
                self.check_vertex_desired(vert_t) 

        if tris_ab_are_red == tris_ab_are_upper:
            triangle_a.update_vertices( [vert_t, vert_b, vert_c] )
            triangle_b.update_vertices( [vert_c, vert_a, vert_t] )
        else:
            triangle_a.update_vertices( [vert_c, vert_t, vert_b] )
            triangle_b.update_vertices( [vert_t, vert_c, vert_a] )

        if c_is_red != c_is_upper:
            triangle_c.update_vertices( [vert_b, vert_a, vert_t] )
        else:
            triangle_c.update_vertices( [vert_a, vert_t, vert_b] )

        ## add gluings and vertices to con_tet

        gluing = tet.adjacentGluing(face_t)
        con_tet.gluings[face_t] = (other_tet, gluing)
        other_face_t = gluing[face_t]
        assert other_con_tet.gluings[other_face_t] == None ### make sure it is not already glued to something
        other_con_tet.gluings[other_face_t] = (tet, other_tet.adjacentGluing(other_face_t))

        for vert_num in range(4):
            if vert_num != face_t:
                # con_tet.vertices[vert_num] = other_con_tet.vertices[gluing[vert_num]]
                con_tet.set_vertex(vert_num, other_con_tet.vertices[gluing[vert_num]])
            else:
                assert vert_num == face_t
                # con_tet.vertices[vert_num] = vert_t
                con_tet.set_vertex(vert_num, vert_t)

        ### now for the edges

        edge_at_index = tet.face(1, vert_pair_to_edge_num[(face_a, face_t)]).index()
        edge_bt_index = tet.face(1, vert_pair_to_edge_num[(face_b, face_t)]).index()
        edge_ct_index = tet.face(1, vert_pair_to_edge_num[(face_c, face_t)]).index()

        edge_at = landscape_edge(self, [vert_a, vert_t], edge_at_index, not triangle.is_upper) ## coastal
        edge_bt = landscape_edge(self, [vert_b, vert_t], edge_bt_index, triangle.is_upper) ## coastal
        edge_ct = landscape_edge(self, [vert_c, vert_t], edge_ct_index, far_edge_colour == "red") ## never coastal

        coastal_index = self.coast.index(vert_a)
        self.coastal_edges = self.coastal_edges[:coastal_index] + [edge_at, edge_bt] + self.coastal_edges[coastal_index + 1:]

        if triangle.is_red == triangle.is_upper:
            edge_bc, edge_ca, edge_ab = triangle.edges 
        else:
            edge_ca, edge_ab, edge_bc = triangle.edges

        if tris_ab_are_red == tris_ab_are_upper:
            triangle_a.update_edges( [edge_bc, edge_ct, edge_bt] )
            triangle_b.update_edges( [edge_at, edge_ct, edge_ca] )
        else:
            triangle_a.update_edges( [edge_bt, edge_bc, edge_ct] )
            triangle_b.update_edges( [edge_ca, edge_at, edge_ct] )

        if c_is_red != c_is_upper:
            triangle_c.update_edges( [edge_at, edge_bt, edge_ab] )
        else:
            triangle_c.update_edges( [edge_bt, edge_ab, edge_at] )

        if triangle.is_upper:
            con_tet.set_upper_edge(edge_ct)
            con_tet.set_lower_edge(edge_ab)
        else: 
            con_tet.set_upper_edge(edge_ab)
            con_tet.set_lower_edge(edge_ct)
        con_tet.set_equatorial_edges([edge_at, edge_bt, edge_bc, edge_ca])

        assert edge_ct.links(edge_ab)
        assert not edge_ct.links(edge_at)
        assert not edge_ct.links(edge_bt)
        assert not edge_ct.links(edge_bc)
        assert not edge_ct.links(edge_ca)
        ### could check more...

        ### now for the neighbours

        if tris_ab_are_red == tris_ab_are_upper:
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

        self.num_tetrahedra += 1

        ### now do all the in_fills that we can, coming from the new stuff we added.
        ### This step is not how we do things in [FSS] - there we add tetrahedra only until we cover the triangle we are trying to bury.
        ### Here we must make the continent convex in order to work out where the cusp leaves go as we build.
        self.back_fill(triangle, vert_a, vert_b, vert_t, triangle_a, triangle_b, edge_ab, edge_at, edge_bt, tris_ab_are_upper)

    def back_fill(self, triangle, vert_a, vert_b, vert_t, triangle_a, triangle_b, edge_ab, edge_at, edge_bt, tris_ab_are_upper):    
        """Running this as part of each coastal_fill makes the continent convex"""
        ### at least one of a and b triangles points to the coast, the other we may need to in_fill along

        outer_triangles = [triangle_a, triangle_b]
        while True:
            tri_to_fill = None
            for tri in outer_triangles:
                nb = tri.downriver()
                if tri.is_upper == nb.is_upper:
                    if nb.downriver() == tri:
                        tri_to_fill = tri
                        break
            if tri_to_fill != None:  
                outer_triangles = self.in_fill(tri_to_fill)
            else:
                break 

        ### now the cusp leaves

        ### we may have added lots of faces incident to vert_t. only one of them should have a cusp leaf 
        ### on the same side as we added the filling tetrahedra. That triangle is one of the outer_triangles 

        # print('a', self.vertices.index(vert_a), 'b', self.vertices.index(vert_b), 'c', self.vertices.index(vert_c))
        # print('tris_ab_are_upper', tris_ab_are_upper, 'edge_ab.is_red', edge_ab.is_red)

        last_tri = None
        for tri in outer_triangles:
            if vert_t == tri.vertices[tri.downriver_index()]:
                last_tri, tri_is_coastal = self.flow(tri)
                assert tri_is_coastal
                break
        
        last_edge = last_tri.edges[last_tri.downriver_index()]
        t_leaf = cusp_leaf(self, vert_t, last_edge, tris_ab_are_upper)
        vert_t.cusp_leaves.append(t_leaf)

        ### find where on last_edge.cusp_leaves this cusp_leaf lands
        if tris_ab_are_upper == last_edge.is_red:
            insert_index = 0 # default at start
        else:
            insert_index = None # default at end 
        for i, leaf in enumerate(last_edge.cusp_leaves):
            if leaf.is_upper == t_leaf.is_upper: 
                if are_anticlockwise(t_leaf.cusp.coastal_index(), leaf.cusp.coastal_index(), last_edge.coastal_index() + 0.5):
                    insert_index = i + 1 ## goes after this one, can repeat
                else:
                    if insert_index == None:
                        insert_index = i ## before this one, cannot repeat
                        break
        # print('insert_index', insert_index)
        if insert_index == None:
            last_edge.cusp_leaves.append(t_leaf)
        else:
            last_edge.cusp_leaves.insert(insert_index, t_leaf)

        ### now move forward all the cusp_leaves on the edge that got covered up
        leaves_to_move_forward = edge_ab.cusp_leaves
        edge_ab.cusp_leaves = []  ### clear because no longer coastal

        same_side_leaves = []
        opposite_side_leaves = []
        for leaf in leaves_to_move_forward:
            if leaf.is_upper == triangle.is_upper:
                same_side_leaves.append(leaf)
            else:
                opposite_side_leaves.append(leaf)

        before_cut = []
        after_cut = []
        for leaf in same_side_leaves:
            if are_anticlockwise(vert_t.coastal_index(), last_edge.coastal_index() + 0.5, leaf.cusp.coastal_index()):
                before_cut.append(leaf)
            else:
                after_cut.append(leaf)

        ### add new opposite side cusp leaf
        if edge_ab.is_red == tris_ab_are_upper:
            new_opposite_leaf = cusp_leaf(self, vert_a, edge_bt, not tris_ab_are_upper)
            vert_a.cusp_leaves.append(new_opposite_leaf)
            opposite_side_leaves.insert(0, new_opposite_leaf)
        else:
            new_opposite_leaf = cusp_leaf(self, vert_b, edge_at, not tris_ab_are_upper)
            vert_b.cusp_leaves.insert(0, new_opposite_leaf)
            opposite_side_leaves.append(new_opposite_leaf)

        ### add on opposite_side_leaves
        if tris_ab_are_upper == edge_ab.is_red:
            after_cut = after_cut + opposite_side_leaves
        else:
            before_cut = opposite_side_leaves + before_cut

        edge_at.cusp_leaves = before_cut
        edge_bt.cusp_leaves = after_cut
        for leaf in edge_at.cusp_leaves:
            leaf.coastal_edge = edge_at
        for leaf in edge_bt.cusp_leaves:
            leaf.coastal_edge = edge_bt

    # def make_convex(self): 
    #     ### new triangles are added to the end of the list so this is safe.
    #     ### new sinks created by in-fills have new triangles that point at them, so we get everything.
    #     for tri in self.triangles:
    #         if tri.is_boundary():
    #             downriver_triangle, is_coastal = self.flow(tri)
    #             if not is_coastal:
    #                 self.in_fill(downriver_triangle)

    def is_convex(self):
        for tri in self.triangles:
            if tri.is_boundary():
                downriver_triangle, is_coastal = self.flow(tri)
                if not is_coastal:
                    return False
        return True

    def build_boundary_data(self):
        """Make data used when drawing so only called once"""
        self.upper_landscape_triangles = set([])
        self.lower_landscape_triangles = set([])
        for tri in self.triangles:
            if tri.is_boundary():
                if tri.is_upper:
                    self.upper_landscape_triangles.add(tri)
                else:
                    self.lower_landscape_triangles.add(tri)

        self.upper_landscape_edges = set([])
        self.lower_landscape_edges = set([])

        # for e in self.edges:
        #     e.boundary_triangles = []

        for tri in self.upper_landscape_triangles:
            for e in tri.edges:
                # e.boundary_triangles.append(tri)
                self.upper_landscape_edges.add(e)
        for tri in self.lower_landscape_triangles:
            for e in tri.edges:
                # e.boundary_triangles.append(tri)
                self.lower_landscape_edges.add(e)

    def old_coast(self):
        old_coast = []
        for tri in self.triangles:
            if tri.is_boundary():
                break  ## found an initial unburied tri 
        vert_index = 0  
        initial_vert = tri.vertices[vert_index]
        vert = initial_vert   
        while old_coast == [] or vert != initial_vert:
            if tri.neighbours[(vert_index - 1) % 3].is_upper != tri.is_upper: ## we are coastal
                vert_index = (vert_index + 1) % 3   
                vert = tri.vertices[vert_index]
                old_coast.append(vert)
            else: # walk to the next triangle
                tri = tri.neighbours[(vert_index - 1) % 3]  
                vert_index = tri.vertices.index(vert) 

        #  i+1
        #   *-------* i  
        #    \     /
        #     \   / 
        #      \ /
        #       * i-1 

        # now rotate to put infinity first
        inf_vert_index = old_coast.index( self.infinity )
        old_coast = old_coast[inf_vert_index:] + old_coast[:inf_vert_index]
        return old_coast

    def build_naive(self, max_num_tetrahedra = 50000):  ### just keep building until we hit max tetrahedra
        self.first_non_buried_index = 0
        while self.num_tetrahedra < max_num_tetrahedra:  
            tri = self.triangles[self.first_non_buried_index]  
            self.bury(tri)
            self.first_non_buried_index += 1
            while self.triangles[self.first_non_buried_index].is_buried():
                self.first_non_buried_index += 1
        self.build_boundary_data()  

    def build_triangulation_fundamental_domain(self, max_num_tetrahedra = 50000):
        """Build a continent that contains a fundamental domain as given to us by 
           tree_faces of fundamental_domain.spanning_dual_tree"""
        initial_tet_num = self.tet_face.tet_num
        tree_faces, non_tree_faces, distances_to_root = spanning_dual_tree(self.vt.tri, initial_tet_num = initial_tet_num)
        print(tree_faces)
        initial_continent_tet = self.tetrahedra[0]
        continent_fund_dom_tets = [initial_continent_tet]
        initial_regina_tet = self.vt.tri.tetrahedron(initial_tet_num)
        initial_landscape_triangles = initial_continent_tet.ordered_faces()

        continent_tree_faces_frontier = []
        for i in range(4):
            face_index = initial_regina_tet.face(2, i).index()
            if face_index in tree_faces:
                continent_tree_faces_frontier.append(initial_landscape_triangles[i])
                tree_faces.remove(face_index)

        while len(continent_fund_dom_tets) < self.vt.tri.countTetrahedra():
            triangle = continent_tree_faces_frontier.pop()
            neighbouring_continent_tets = [triangle.upper_tet, triangle.lower_tet] ### could be None
            got_already = [tet in continent_fund_dom_tets for tet in neighbouring_continent_tets]
            assert got_already.count(True) == 1
            new_continent_tet = neighbouring_continent_tets[got_already.index(False)]
            if new_continent_tet == None:
                self.bury(triangle)
            neighbouring_continent_tets = [triangle.upper_tet, triangle.lower_tet]
            new_continent_tet = neighbouring_continent_tets[got_already.index(False)]
            assert new_continent_tet != None
            continent_fund_dom_tets.append(new_continent_tet)
            new_regina_tet = self.vt.tri.tetrahedron(new_continent_tet.index)
            new_landscape_triangles = new_continent_tet.ordered_faces()
            for i in range(4):
                face_index = new_regina_tet.face(2, i).index()
                if face_index in tree_faces:
                    continent_tree_faces_frontier.append(new_landscape_triangles[i])
                    tree_faces.remove(face_index)
        continent_fund_dom_tets.sort(key = lambda t : t.index)
        return continent_fund_dom_tets


    def build_boundary_fundamental_domain_old(self, max_num_tetrahedra = 50000):
        self.first_non_buried_index = 0
        while len(self.desired_vertices) > 0 and self.num_tetrahedra < max_num_tetrahedra:  # will go a little over because we check after each bury, which adds many tetrahedra
            tri = self.triangles[self.first_non_buried_index]  
            self.bury(tri)
            self.first_non_buried_index += 1
            while self.triangles[self.first_non_buried_index].is_buried():
            # while self.triangles[first_non_buried_index].is_buried() or self.triangles[first_non_buried_index].is_upper:
                self.first_non_buried_index += 1
        self.build_boundary_data()  

    ### old version builds lots of things we dont care about, this is much faster.
    def build_boundary_fundamental_domain(self, max_num_tetrahedra = 50000):
        ### fundamental domain for the boundary torus?
        self.first_non_buried_index = 0
        while len(self.desired_vertices) > 0 and self.num_tetrahedra < max_num_tetrahedra:  # will go a little over because we check after each bury, which adds many tetrahedra
            tri = self.triangles[self.first_non_buried_index] 
             ### if this tri is incident to infinity, bury it
            if self.infinity in tri.vertices:
                self.bury(tri)
            self.first_non_buried_index += 1
            while self.triangles[self.first_non_buried_index].is_buried():
            # while self.triangles[first_non_buried_index].is_buried() or self.triangles[first_non_buried_index].is_upper:
                self.first_non_buried_index += 1
        self.build_boundary_data()  



    def build_on_coast(self, max_length = 0.1, max_num_tetrahedra = 50000):  # build until all edges we want to draw are short
        self.max_length = max_length
        # print(('max_length', max_length))

        ## now build

        # while self.num_long_edges > 0 and self.num_tetrahedra < max_num_tetrahedra: 
        while self.num_tetrahedra < max_num_tetrahedra and self.first_non_buried_index < len(self.triangles): 
            tri = self.triangles[self.first_non_buried_index]  

            if tri.ladderpole_descendant_long_coastal_indices() != []:
                self.bury(tri)
            self.first_non_buried_index += 1
            while self.first_non_buried_index < len(self.triangles) and self.triangles[self.first_non_buried_index].is_buried():
                self.first_non_buried_index += 1
        # print(('num_tetrahedra', self.num_tetrahedra))
        hit_max_tetrahedra = self.num_tetrahedra >= max_num_tetrahedra
        # print(('hit max tetrahedra', hit_max_tetrahedra))
        self.build_boundary_data()
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
            while self.first_non_buried_index < len(self.triangles) and self.triangles[self.first_non_buried_index].is_buried():
                self.first_non_buried_index += 1
        # print(('num_tetrahedra', self.num_tetrahedra))
        hit_max_tetrahedra = self.num_tetrahedra >= max_num_tetrahedra
        # print(('hit max tetrahedra', hit_max_tetrahedra))
        self.build_boundary_data()
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
            while self.first_non_buried_index < len(self.triangles) and self.triangles[self.first_non_buried_index].is_buried():
                self.first_non_buried_index += 1
        # print(('num_tetrahedra', self.num_tetrahedra))
        hit_max_tetrahedra = self.num_tetrahedra >= max_num_tetrahedra
        # print(('hit max tetrahedra', hit_max_tetrahedra))
        self.build_boundary_data()
        # print(('num_long_edges_direct_count', self.count_long_edges()))
        # print(('max_coastal_edge_length', self.calculate_max_ladderpole_descendant_coast_edge_length()))
        return hit_max_tetrahedra

    # def build_loxodromics(self, max_length = 0.1, max_num_tetrahedra = 50000):
    #     self.max_length = max_length
    #     print 'max_length', max_length
    #     self.build_boundary_data()

    #     ## now build

    #     while self.first_non_buried_index < len(self.triangles) and self.num_tetrahedra < max_num_tetrahedra: 
    #         tri = self.triangles[self.first_non_buried_index]  

    #         if any( [(edge.is_ladderpole_descendant and tri.dist_from_lox(edge) > self.max_length) for edge in tri.edges] ):
    #             self.bury(tri)
    #         self.first_non_buried_index += 1
    #         while self.first_non_buried_index < len(self.triangles) and self.triangles[self.first_non_buried_index].is_buried():
    #             self.first_non_buried_index += 1
    #     print 'num_tetrahedra', self.num_tetrahedra
    #     self.build_boundary_data()
    #     print 'num_long_edges_direct_count', self.count_long_edges()
    #     print 'max_coastal_edge_length', self.calculate_max_ladderpole_descendant_coast_edge_length()

    def build_long_and_mid(self, max_length = 0.1, max_num_tetrahedra = 50000):  
        self.max_length = max_length
        # print(('max_length', max_length))

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
            while self.first_non_buried_index < len(self.triangles) and self.triangles[self.first_non_buried_index].is_buried():
                self.first_non_buried_index += 1
        # print(('num_tetrahedra', self.num_tetrahedra))
        hit_max_tetrahedra = self.num_tetrahedra >= max_num_tetrahedra
        # print(('hit max tetrahedra', hit_max_tetrahedra))
        self.build_boundary_data()
        # print(('num_long_edges_direct_count', self.count_long_edges()))
        # print(('max_coastal_edge_length', self.calculate_max_ladderpole_descendant_coast_edge_length()))
        return hit_max_tetrahedra

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
            if tri.is_boundary():
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
                                        river_edges = river_edges[:i] + [the_lower_tet.lower_edge] + river_edges[i:]
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
                                    river_edges = river_edges[:i] + [the_lower_tet.lower_edge] + river_edges[i+1:]
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
                                        river_edges = river_edges[:i] + [the_upper_tet.upper_edge] + river_edges[i:]
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
                                    river_edges = river_edges[:i] + [the_upper_tet.upper_edge] + river_edges[i+1:]
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





