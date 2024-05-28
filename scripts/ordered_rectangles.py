


class ordered_rectangle:
    def __init__(self, horiz_ordering, vert_ordering):
        assert len(horiz_ordering) == len(vert_ordering)
        self.horiz_ordering = horiz_ordering
        self.vert_ordering = vert_ordering

    def horiz_to_vert_perm(self):
        return [self.vert_ordering.index(x) for x in self.horiz_ordering]

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

class ordered_face_rectangle(ordered_rectangle):
    def __init__(self, con, face_index, horiz_ordering, vert_ordering):
        self.con = con
        self.face_index = face_index
        ordered_rectangle.__init__(self, horiz_ordering, vert_ordering)

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
        self.N_ind, self.S_ind = top_vertices
        self.W_ind, self.E_ind = bottom_vertices
        self.W, self.E = self.horiz_ordering[0], self.horiz_ordering[-1] 
        self.S, self.N = self.vert_ordering[0], self.vert_ordering[-1] 

    def face(self, ind):
        face_index = self.Regina_tet.face(2, ind).index()
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

def build_cardinal_ordering(cusp_order, chunks):
    ordering = [cusp_order[0]]
    for j in range(3):
        ordering.extend(chunks[j])
        ordering.append(cusp_order[j+1])
    return ordering

def build_tetrahedron_rectangle_orderings(con, tetrahedra_cusp_orders, intervals_inside_tet_rectangles, tetrahedra_chunks):
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
    print('sanity check')
    r = old_tet_rectangles[0]
    tri = r.con.vt.tri
    for f in tri.triangles():
        # print(f)
        face_rect0 = face_rect_from_face(f, 0, old_tet_rectangles)
        face_rect1 = face_rect_from_face(f, 1, old_tet_rectangles)
        if face_rect0.corner_location() != face_rect1.corner_location():
            face_rect1.rotate()
        assert face_rect0.corner_location() == face_rect1.corner_location()
        assert face_rect0.horiz_to_vert_perm() == face_rect1.horiz_to_vert_perm()






        