from veering.taut import edge_num_to_vert_pair
import continent
import flow_interval

def is_between(i, j, k): ### is j between i and k inclusive
    return (i - j) * (j - k) >= 0

def rect_is_empty(j, i, pj, pi, perm):
    for k in range(j + 1, i):
        pk = perm[k]
        if is_between(pj, pk, pi):
        # if (pj - pk) * (pk - pi) > 0: ### pk is between pi and pj
            return False ### k breaks apart the edge rectangle
    return True

class ordered_rectangle:
    def __init__(self, horiz_ordering, vert_ordering):
        assert len(horiz_ordering) == len(vert_ordering)
        self.horiz_ordering = horiz_ordering
        self.vert_ordering = vert_ordering
        self.edge_rectangle_indices = None
        self.face_rectangle_indices = None
        self.tetrahedron_rectangle_indices = None

    def horiz_to_vert_perm(self):
        return [self.vert_ordering.index(x) for x in self.horiz_ordering]

    def vert_to_horiz_perm(self):
        return [self.horiz_ordering.index(x) for x in self.vert_ordering]

    def __repr__(self):
        return str(self.horiz_to_vert_perm())

    def copy(self):
        other = ordered_rectangle(self.horiz_ordering, self.vert_ordering)
        return other

    def rotate(self):
        self.horiz_ordering.reverse()
        self.vert_ordering.reverse()

    def rotated(self):
        other = self.copy()
        other.rotate()
        return other

    def edge_rectangles(self):
        if self.edge_rectangle_indices != None:
            return self.edge_rectangle_indices
        perm = self.horiz_to_vert_perm()
        out = []
        for i in range(len(self.horiz_ordering)):
            pi = perm[i]
            for j in range(i):
                pj = perm[j]
                if rect_is_empty(j, i, pj, pi, perm):
                    out.append((j, i))
        self.edge_rectangle_indices = out
        return out

    def face_rectangles(self):
        if self.face_rectangle_indices != None:
            return self.face_rectangle_indices
        if self.edge_rectangle_indices == None:
            self.edge_rectangles()
        out = []
        for i in range(len(self.horiz_ordering)):
            for j in range(i):
                if (j,i) in self.edge_rectangle_indices:
                    for k in range(j):
                        if (k,i) in self.edge_rectangle_indices and (k,j) in self.edge_rectangle_indices:
                            out.append((k,j,i))
        self.face_rectangle_indices = out
        return out

    def tetrahedron_rectangles(self):
        if self.tetrahedron_rectangle_indices != None:
            return self.tetrahedron_rectangle_indices
        if self.face_rectangle_indices == None:
            self.face_rectangles()
        out = []
        for i in range(len(self.horiz_ordering)):
            for j in range(i):
                if (j,i) in self.edge_rectangle_indices:
                    for k in range(j):
                        if (k,j,i) in self.face_rectangle_indices:
                            for l in range(k): ### three triangles means the fourth automatically exists
                                if (l,k,j) in self.face_rectangle_indices and (l,k,i) in self.face_rectangle_indices:
                                    out.append((l,k,j,i))
        self.tetrahedron_rectangle_indices = out
        return out

class ordered_edge_rectangle(ordered_rectangle):
    def __init__(self, con, edge_index, horiz_ordering, vert_ordering):
        self.con = con
        self.edge_index = edge_index
        ordered_rectangle.__init__(self, horiz_ordering, vert_ordering)

class ordered_face_rectangle(ordered_rectangle):
    def __init__(self, con, face_index, horiz_ordering, vert_ordering, canonise_orientation = True):
        self.con = con
        self.face_index = face_index
        ordered_rectangle.__init__(self, horiz_ordering, vert_ordering)
        if canonise_orientation:
            if self.corner_location()[0] != 0:
                self.rotate() ### make sure that the corner cusp is on the western side
                assert self.corner_location()[0] == 0

    def corner_location(self):
        W, E = self.horiz_ordering[0], self.horiz_ordering[-1]
        S, N = self.vert_ordering[0], self.vert_ordering[-1]
        if W == S:
            return (0, 0)
        elif E == S:
            return (-1, 0)
        elif E == N: 
            return (-1,-1)
        elif W == N:
            return (0, -1)
        assert False # should not get here

class ordered_tetrahedron_rectangle(ordered_rectangle):
    def __init__(self, con, tet_index, horiz_ordering, vert_ordering):
        self.con = con
        self.tet_index = tet_index
        self.Regina_tet = con.vt.tri.tetrahedron(tet_index)
        ordered_rectangle.__init__(self, horiz_ordering, vert_ordering)
        top_vertices, bottom_vertices = con.vt.tet_vert_posns[tet_index] ### this is currently set in build_continent 
        self.N_ind, self.S_ind = top_vertices ### Regina vertex indices
        self.W_ind, self.E_ind = bottom_vertices
        self.W, self.E = self.horiz_ordering[0], self.horiz_ordering[-1] 
        self.S, self.N = self.vert_ordering[0], self.vert_ordering[-1] 
        self.Regina2cusp = [{self.W_ind: self.W, self.E_ind: self.E, self.S_ind: self.S, self.N_ind: self.N}[i] for i in range(4)]
        self.generate_sub_cells()

    def face(self, ind):
        face_index = self.Regina_tet.triangle(ind).index()
        if ind == self.S_ind:
            southern_ind = min(self.vert_ordering.index(self.W), self.vert_ordering.index(self.E))
            face_vert_ordering = self.vert_ordering[southern_ind:]
            face_horiz_ordering = [x for x in self.horiz_ordering if x in face_vert_ordering]
        elif ind == self.N_ind:
            northern_ind = max(self.vert_ordering.index(self.W), self.vert_ordering.index(self.E))
            face_vert_ordering = self.vert_ordering[:northern_ind + 1]
            face_horiz_ordering = [x for x in self.horiz_ordering if x in face_vert_ordering]
        elif ind == self.W_ind:
            western_ind = min(self.horiz_ordering.index(self.S), self.horiz_ordering.index(self.N))
            face_horiz_ordering = self.horiz_ordering[western_ind:]
            face_vert_ordering = [x for x in self.vert_ordering if x in face_horiz_ordering]
        elif ind == self.E_ind:
            eastern_ind = max(self.horiz_ordering.index(self.S), self.horiz_ordering.index(self.N))
            face_horiz_ordering = self.horiz_ordering[:eastern_ind + 1]
            face_vert_ordering = [x for x in self.vert_ordering if x in face_horiz_ordering]
        face_rect = ordered_face_rectangle(self.con, face_index, face_horiz_ordering, face_vert_ordering)
        return face_rect

    def edge(self, ind):
        Regina_edge = self.Regina_tet.edge(ind)
        edge_index = Regina_edge.index()
        cusps = [self.Regina2cusp[i] for i in edge_num_to_vert_pair[ind]]
        cusp_horiz_indices = [self.horiz_ordering.index(c) for c in cusps]
        cusp_horiz_indices.sort()
        cusp_vert_indices = [self.vert_ordering.index(c) for c in cusps]
        cusp_vert_indices.sort()
        edge_horiz_ordering = self.horiz_ordering[cusp_horiz_indices[0]: cusp_horiz_indices[1] + 1]
        edge_horiz_ordering = [x for x in edge_horiz_ordering if is_between(cusp_vert_indices[0], self.vert_ordering.index(x), cusp_vert_indices[1])]
        edge_vert_ordering = [x for x in self.vert_ordering if x in edge_horiz_ordering]
        edge_rect = ordered_edge_rectangle(self.con, edge_index, edge_horiz_ordering, edge_vert_ordering)

        ### Unlike for faces, which have a way to canonise the orientation of their face rectangles with only local data,
        ### we have to canonise the orientation of the edge rectangle based on the Regina ordering of its vertices.
        perm = self.Regina_tet.edgeMapping(ind)
        if edge_horiz_ordering.index(self.Regina2cusp[perm[0]]) != 0:
            assert edge_horiz_ordering.index(self.Regina2cusp[perm[1]]) == 0
            edge_rect.rotate()
        return edge_rect

    def generate_sub_cells(self):
        self.faces = [self.face(ind) for ind in range(4)]
        self.edges = [self.edge(ind) for ind in range(6)]

    def sub_cell_image(self, sub_cell_dim, sub_cell_ind, new_rectangle_indices):
        """Tries to map new_rectangle_indices into the sub cell (an edge or face of the tet). 
        Returns False if the new rectangle is not contained in the sub cell.
        Otherwise returns the indices for the new rectangle in the sub cell."""
        verts = [self.horiz_ordering[i] for i in new_rectangle_indices]
        if sub_cell_dim == 1:
            sub_cell = self.edges[sub_cell_ind]
        elif sub_cell_dim == 2:
            sub_cell = self.faces[sub_cell_ind]
        else:
            assert sub_cell_dim == 3
            sub_cell = self
        if not all([v in sub_cell.horiz_ordering for v in verts]):
            return False
        else:
            horiz_indices = [sub_cell.horiz_ordering.index(v) for v in verts]
            ### is either ascending or descending list of integers
            vert_indices = [sub_cell.vert_ordering.index(v) for v in verts]
            if horiz_indices[0] > horiz_indices[1]:
                horiz_indices.reverse() ### make increasing
                vert_indices.reverse()
            return (tuple(horiz_indices), tuple(vert_indices))

def build_cardinal_ordering(cusp_order, chunks):
    ordering = [cusp_order[0]]
    for j in range(3):
        ordering.extend(chunks[j])
        ordering.append(cusp_order[j+1])
    return ordering

def build_tetrahedron_rectangle_orderings(con, tetrahedra_cusp_orders, tetrahedra_chunks):
    old_tet_rectangles = []
    for i in range(len(tetrahedra_cusp_orders)):
        horizontal_cusp_order, vertical_cusp_order = tetrahedra_cusp_orders[i]
        we_chunks, sn_chunks = tetrahedra_chunks[i]
        horiz_ordering = build_cardinal_ordering(horizontal_cusp_order, we_chunks)
        vert_ordering = build_cardinal_ordering(vertical_cusp_order, sn_chunks)
        old_tet_rect = ordered_tetrahedron_rectangle(con, i, horiz_ordering, vert_ordering)
        old_tet_rectangles.append(old_tet_rect)
    return old_tet_rectangles

def face_rect_from_face(f, ind, old_tet_rectangles):
    embed = f.embedding(ind)
    tet = embed.simplex()
    tet_num = tet.index()
    tet_face_num = embed.face()
    return old_tet_rectangles[tet_num].face(tet_face_num)

def sanity_check(old_tet_rectangles):
    # print('sanity check')
    r = old_tet_rectangles[0]
    tri = r.con.vt.tri
    for f in tri.triangles():
        # print(f)
        face_rect0 = face_rect_from_face(f, 0, old_tet_rectangles)
        face_rect1 = face_rect_from_face(f, 1, old_tet_rectangles)
        # if face_rect0.corner_location() != face_rect1.corner_location(): 
        #     face_rect1.rotate()    ### unnecessary because we canonise_orientation 
        assert face_rect0.corner_location() == face_rect1.corner_location()
        assert face_rect0.horiz_to_vert_perm() == face_rect1.horiz_to_vert_perm()

def new_rectangle_equiv_class(new_rect, containing_old_tet_rect):
    this_rect_equiv_class = set([])
    ### edges
    for i in range(6):
        im = containing_old_tet_rect.sub_cell_image(1, i, new_rect)
        if im != False: ### the rect fits inside the edge
            Regina_edge = containing_old_tet_rect.Regina_tet.edge(i)
            edge_index = Regina_edge.index()
            new_rect_sig = (1, edge_index, im)  ### 1 means edge
            this_rect_equiv_class.add(new_rect_sig)
    ### faces
    for j in range(4):
        im = containing_old_tet_rect.sub_cell_image(2, j, new_rect)
        if im != False: ### the rect fits inside the face
            Regina_face = containing_old_tet_rect.Regina_tet.triangle(j)
            face_index = Regina_face.index()
            new_rect_sig = (2, face_index, im)  ### 2 means face
            this_rect_equiv_class.add(new_rect_sig)      
    ### the tet itself
    im = containing_old_tet_rect.sub_cell_image(3, 0, new_rect)
    new_rect_sig = (3, containing_old_tet_rect.tet_index, im) 
    this_rect_equiv_class.add(new_rect_sig)  
    return this_rect_equiv_class

# def anticlockwise_ordering(horiz_verts, vertical_verts):
#     W = 0
#     E = 3
#     if vertical_verts[1] < vertical_verts[2]:
#         S = 1
#         N = 2
#     else:
#         S = 2
#         N = 1
#     return [E, N, W, S]



class new_tetrahedron:
    def __init__(self, equiv_class, old_tet_rectangles):
        self.equiv_class = list(equiv_class)
        self.equiv_class.sort()
        self.representative = self.equiv_class[-1]
        self.dim, self.ind, verts = self.representative
        assert self.dim == 3 ### must be in some old tet rect. We want this for the rep for the cusp indices
        self.horiz_verts, self.vertical_verts = verts
        self.anticlockwise_ordering = None
        self.set_anticlockwise_ordering()
        self.index = None
        self.faces = [None, None, None, None]
        self.face_special_corner = [None, None, None, None] ### which of four positions is the special corner of the face
        self.special_vertex_index_in_face = [None, None, None, None]
        self.adjTetFace = [None, None, None, None]
        self.adjGluing = []
        for i in range(4):
            self.adjGluing.append([None, None, None, None])  ### permutation for each gluing
        self.cusp_index = [] ### what is the cusp index in the manifold for this vertex. n for the newly drilled cusp
        
        ### fix
        old_tet_rectangle = old_tet_rectangles[self.ind]
        for i in range(4):
            j = self.anticlockwise_to_horiz[i]
            ind = self.horiz_verts[j]
            if type(old_tet_rectangle.horiz_ordering[ind]) == continent.vertex:
                self.cusp_index.append(old_tet_rectangle.horiz_ordering[ind].Regina_cusp_num)
            else:
                assert type(old_tet_rectangle.horiz_ordering[ind]) == flow_interval.flow_interval
                self.cusp_index.append(old_tet_rectangle.con.vt.tri.countVertices()) ### new cusp
        # print('self.cusp_index', self.cusp_index)

    def set_anticlockwise_ordering(self):
        W, E = 0, 3
        if self.vertical_verts[1] < self.vertical_verts[2]:
            S, N = 1, 2
        else:
            S, N = 2, 1
        self.anticlockwise_to_horiz = [E, N, W, S]
        self.horiz_to_anticlockwise = [self.anticlockwise_to_horiz.index(i) for i in range(4)]  

    def face_rep(self, i):
        horiz_verts = list(self.horiz_verts)
        vertical_verts = list(self.vertical_verts)
        j = self.anticlockwise_to_horiz[i]
        popped_horiz_verts = horiz_verts[:]
        popped_horiz_verts.pop(j)
        popped_vertical_verts = vertical_verts[:]
        popped_vertical_verts.pop(j)
        a, b = min(popped_vertical_verts), max(popped_vertical_verts)
        special_vert_in_tet = None
        if popped_vertical_verts[0] == a:
            self.face_special_corner[i] = (0, 0) ### WS corner
            special_vert_in_tet = self.horiz_to_anticlockwise[vertical_verts.index(a)]
        elif popped_vertical_verts[0] == b:
            self.face_special_corner[i] = (0, -1) ### WN corner
            special_vert_in_tet = self.horiz_to_anticlockwise[vertical_verts.index(b)]
        elif popped_vertical_verts[-1] == a:
            self.face_special_corner[i] = (-1, 0) ### ES corner
            special_vert_in_tet = self.horiz_to_anticlockwise[vertical_verts.index(a)]
        elif popped_vertical_verts[-1] == b:
            self.face_special_corner[i] = (-1, -1) ### EN corner
            special_vert_in_tet = self.horiz_to_anticlockwise[vertical_verts.index(b)]
        assert i != special_vert_in_tet
        if i < special_vert_in_tet:
            self.special_vertex_index_in_face[i] = special_vert_in_tet - 1
        else:
            self.special_vertex_index_in_face[i] = special_vert_in_tet

        return (self.dim, self.ind, (tuple(popped_horiz_verts), tuple(popped_vertical_verts)))

class new_triangle:
    def __init__(self, equiv_class):
        self.equiv_class = equiv_class
        self.tet_face = []

def build_drilled_triangulation_data(old_tet_rectangles):
    sanity_check(old_tet_rectangles)
    new_face_equiv_classes = []
    new_tet_equiv_classes = []
    
    for old_tet_rect in old_tet_rectangles:
        this_tet_new_face_rects = old_tet_rect.face_rectangles()
        for new_face_rect in this_tet_new_face_rects:
            ### find which sub cell(s) it belongs to. Could be in multiple edges.
            this_face_rect_equiv_class = new_rectangle_equiv_class(new_face_rect, old_tet_rect)
            ### now check if we already have this equivalence class
            intersecting_equiv_classes = []
            for equiv_class in new_face_equiv_classes:
                # print('check against equiv_class', equiv_class)
                if not this_face_rect_equiv_class.isdisjoint(equiv_class):
                    intersecting_equiv_classes.append(equiv_class)
            for equiv_class in intersecting_equiv_classes:
                this_face_rect_equiv_class.update(equiv_class) ### union and update
                new_face_equiv_classes.remove(equiv_class)
            new_face_equiv_classes.append(this_face_rect_equiv_class)

        this_tet_new_tet_rects = old_tet_rect.tetrahedron_rectangles()
        for new_tet_rect in this_tet_new_tet_rects:
            ### find which sub cell(s) it belongs to. Could be in multiple edges.
            this_tet_rect_equiv_class = new_rectangle_equiv_class(new_tet_rect, old_tet_rect)
            ### now check if we already have this equivalence class
            intersecting_equiv_classes = []
            for equiv_class in new_tet_equiv_classes:
                if not this_tet_rect_equiv_class.isdisjoint(equiv_class):
                    intersecting_equiv_classes.append(equiv_class)
            for equiv_class in intersecting_equiv_classes:
                this_tet_rect_equiv_class.update(equiv_class) ### union and update
                new_tet_equiv_classes.remove(equiv_class)
            new_tet_equiv_classes.append(this_tet_rect_equiv_class)


    # print('num new faces', len(new_face_equiv_classes))
    # for x in new_face_equiv_classes:
    #     print(x)
    # print('num new tetrahedra', len(new_tet_equiv_classes))
    # for x in new_tet_equiv_classes:
    #     print(x)
    assert len(new_face_equiv_classes) == 2 * len(new_tet_equiv_classes)

    new_tetrahedra = []
    new_faces = []
    for equiv_class in new_tet_equiv_classes:
        new_tetrahedra.append(new_tetrahedron(equiv_class, old_tet_rectangles))
    for equiv_class in new_face_equiv_classes:
        new_faces.append(new_triangle(equiv_class))
    for i, new_tet in enumerate(new_tetrahedra):
        new_tet.index = i
        for j in range(4):
            rep = new_tet.face_rep(j)
            for new_face in new_faces:
                if rep in new_face.equiv_class:
                    new_tet.faces[j] = new_face
                    new_face.tet_face.append((new_tet, j))
                    break

    for new_face in new_faces:
        assert len(new_face.tet_face) == 2
        new_face.other_tet_ind = {new_face.tet_face[0]: new_face.tet_face[1], new_face.tet_face[1]: new_face.tet_face[0]}
    for new_tet in new_tetrahedra:
        # print(new_tet.cusp_index)
        for i in range(4):
            face = new_tet.faces[i]
            new_tet.adjTetFace[i] = face.other_tet_ind[(new_tet, i)]
            other_tet, other_i = new_tet.adjTetFace[i]
            this_tet_face_indices = [0,1,2,3] ### anticlockwise order
            this_tet_face_indices.remove(i)
            other_tet_face_indices = [0,1,2,3]
            other_tet_face_indices.remove(other_i)
            # if new_tet.face_special_corner[i] != other_tet.face_special_corner[other_i]:
            #     other_tet_face_indices.reverse()
            this_special = new_tet.special_vertex_index_in_face[i]
            other_special = other_tet.special_vertex_index_in_face[other_i]
            this_tet_face_indices = [this_tet_face_indices[(k + this_special)%3] for k in range(3)]
            other_tet_face_indices = [other_tet_face_indices[(k + other_special)%3] for k in range(3)]

            perm = [None, None, None, None]
            perm[i] = other_i
            for k in range(3):
                perm[this_tet_face_indices[k]] = other_tet_face_indices[k]
            new_tet.adjGluing[i] = perm

    # for new_tet in new_tetrahedra:
    #     print('new_tet', new_tet.index, [(tf[0].index, tf[1]) for tf in new_tet.adjTetFace], new_tet.adjGluing)
    return (new_tetrahedra, new_faces)














        