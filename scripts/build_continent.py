#
# build_continent.py
#

### functions for building continents

from veering.taut import isosig_to_tri_angle, edge_num_to_vert_pair, vert_pair_to_edge_num
from veering.transverse_taut import edge_side_face_collections, get_tet_above_edge, get_tet_top_vert_nums
from veering.veering_tri import veering_triangulation

from continent import continent, continent_tetrahedron
from boundary_triangulation import tet_face

def make_continent_naive(veering_isosig, max_num_tetrahedra = 50):
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle) #, tet_shapes = tet_shapes)
    initial_tet_face = tet_face(vt, 0, 0, verts_pos = [None, None, None, None])
    con = continent( vt, initial_tet_face) #, desired_vertices = desired_vertices )
    con.build_naive(max_num_tetrahedra = max_num_tetrahedra)
    con.build_boundary_data()
    print(len(con.vertices), len(con.edges), len(con.triangles), len(con.tetrahedra))
    return con

def make_continent_fund_dom(veering_isosig, max_num_tetrahedra = 50):
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle) #, tet_shapes = tet_shapes)
    initial_tet_face = tet_face(vt, 0, 0, verts_pos = [None, None, None, None])
    con = continent( vt, initial_tet_face) #, desired_vertices = desired_vertices )
    continent_fund_dom_tets = con.build_triangulation_fundamental_domain(max_num_tetrahedra = max_num_tetrahedra)
    con.build_boundary_data()
    print(len(con.vertices), len(con.edges), len(con.triangles), len(con.tetrahedra))
    return con, continent_fund_dom_tets

# ### this is probably broken, needs to be redone anyway - move around the universal cover building onto the continent as necessary
# def make_continent_drill_dual_cycle(veering_isosig, dual_cycle, num_steps):
#     tri, angle = isosig_to_tri_angle(veering_isosig)
#     vt = veering_triangulation(tri, angle) #, tet_shapes = tet_shapes)
#     ### initial tetrahedron is above face 0 of the dual cycle
#     face0 = vt.tri.triangles()[dual_cycle[0]]
#     embeds = face0.embeddings()
#     tet_0, face_0 = None, None
#     for embed in embeds:
#         tet_num = embed.simplex().index()
#         face_num = embed.face()
#         if vt.coorientations[tet_num][face_num] == -1: ## into the tet
#             tet_0, face_0 = tet_num, face_num
#             break
#     assert tet_0 != None and face_0 != None
#     initial_tet_face = tet_face(vt, tet_0, face_0, verts_pos = [None, None, None, None])
#     con = continent( vt, initial_tet_face)

#     ### identify the triangle corresponding to dual_cycle[1]
#     for triangle in con.triangles:
#         if not triangle.is_upper and triangle.index == dual_cycle[0]:
#             lowest_triangle = triangle
#         if triangle.is_upper and triangle.index == dual_cycle[1]:
#             highest_triangle = triangle
#     triangle_path = [lowest_triangle, highest_triangle]
#     lowest_path_index = 0
#     highest_path_index = 1
#     # for i in range(num_steps):
#     #     last_added_tet = con.bury(highest_triangle)

#     #     ### find next_triangle
#     #     highest_triangle = None
#     #     for triangle in last_added_tet.upper_triangles:
#     #         if triangle.index == dual_cycle[(path_index + 1) % len(dual_cycle)]:
#     #             highest_triangle = triangle
#     #             break
#     #     assert highest_triangle != None
#     #     triangle_path.append(highest_triangle)
#     #     path_index = path_index + 1

#     #     con.build_boundary_data()

#     for i in range(num_steps):
#         if i%2 == 0:
#             last_added_tet = con.bury(triangle_path[-1]) ## go up
#             for triangle in last_added_tet.upper_triangles:
#                 if triangle.index == dual_cycle[(highest_path_index + 1) % len(dual_cycle)]:
#                     triangle_path.append(triangle)
#                     break
#             highest_path_index = highest_path_index + 1
#         else:
#             last_added_tet = con.bury(triangle_path[0]) ## go down
#             for triangle in last_added_tet.lower_triangles:
#                 if triangle.index == dual_cycle[(lowest_path_index - 1) % len(dual_cycle)]:
#                     triangle_path.insert(0, triangle)
#                     break
#             lowest_path_index = lowest_path_index - 1
#     con.build_boundary_data()
#     con.triangle_path = triangle_path
#     return con


def flow_edge_in_continent(con_tet, edge_num):
    tet_vertices = con_tet.vertices
    vert_pair = edge_num_to_vert_pair[edge_num]
    con_vert_pair = [tet_vertices[i] for i in vert_pair]
    upper_tris = con_tet.upper_triangles
    edge_candidates = set(upper_tris[0].edges).union(set(upper_tris[1].edges))
    edge = None
    for e in edge_candidates:
        if set(e.vertices) == set(con_vert_pair): 
            edge = e
            break
    assert edge != None
    return edge

class flow_interval:
    def __init__(self, continent, flow_cycle, init_flow_tetrahedron, init_flow_edge, init_flow_index):
        self.continent = continent
        self.flow_cycle = flow_cycle
        self.tetrahedra = [init_flow_tetrahedron]
        self.edges = [init_flow_edge]
        self.up_index = init_flow_index    ### index in the flow_cycle
        self.down_index = init_flow_index
        self.init_tet = init_flow_tetrahedron
        self.init_edge = init_flow_edge

    def how_far_down(self):
        """how far down have we built the interval from the init tet"""
        ind = self.tetrahedra.index(self.init_tet)
        return ind

    def how_far_up(self):
        """how far up have we built the interval from the init tet"""
        ind = self.tetrahedra.index(self.init_tet)
        return len(self.tetrahedra) - ind - 1

    def position(self, tet):
        """where is tet in the flow interval, where init_flow_tetrahedron is position 0"""
        assert tet in self.tetrahedra
        ind = self.tetrahedra.index(tet)
        return ind - self.tetrahedra.index(self.init_tet)

    def get_tet_at_position(self, i):
        ind = self.tetrahedra.index(self.init_tet)
        return self.tetrahedra[i + ind]

    def get_edge_at_position(self, i):
        ind = self.edges.index(self.init_edge)
        return self.edges[i + ind]

    def go_up(self):
        #flow_cycle, con, flow_edges, flow_tetrahedra, upwards_flow_index):
        edge = self.edges[-1]
        while True:
            upper_boundary_triangles = [t for t in edge.boundary_triangles() if t.is_upper] 
            if len(upper_boundary_triangles) == 0:
                break ## out of while loop
            self.continent.bury(upper_boundary_triangles[0])
        up_tet = edge.upper_tet
        assert up_tet != None
        self.up_index = (self.up_index + 1) % len(self.flow_cycle)
        assert up_tet.index == self.flow_cycle[self.up_index][0]
        self.tetrahedra.append(up_tet)
        self.edges.append(flow_edge_in_continent(up_tet, self.flow_cycle[self.up_index][1]))

    def go_down(self):
        tet = self.tetrahedra[0]
        edge = tet.lower_edge
        self.edges.insert(0, edge) 
        ### now build the continent to get new_lowest_tet
        vt = self.continent.vt
        manifold_edge = vt.tri.edge(edge.index)
        self.down_index = (self.down_index - 1) % len(self.flow_cycle)
        ### is the flow cycle vertical through the tetrahedron? 
        tet_below = get_tet_above_edge(vt.tri, vt.angle, manifold_edge, tet_vert_coorientations = vt.coorientations, get_tet_below_edge = True)
        tet_below_num = tet_below.index()
        top_vert_nums = get_tet_top_vert_nums(vt.coorientations, tet_below_num)
        top_edge_num = vert_pair_to_edge_num[tuple(top_vert_nums)]

        if (tet_below_num, top_edge_num) == self.flow_cycle[self.down_index]: ### flow cycle went straight up
            # print('straight up case')
            if edge.lower_tet == None:
                while True:
                    lower_boundary_triangles = [t for t in edge.boundary_triangles() if not t.is_upper] 
                    if len(lower_boundary_triangles) == 0:
                        break ## out of while loop
                    self.continent.bury(lower_boundary_triangles[0])
            new_tet = edge.lower_tet
        else:
            # print('sideways case')
            ### find which side of the edge our tet is in
            side_tet_collections_at_edge = vt.side_tet_collections[edge.index] ## index in the manifold
            side_face_collections_at_edge = vt.side_face_collections[edge.index]
            downward_path = None
            flow_step = self.flow_cycle[self.down_index] ### pair of (tet index in manifold, index of edge in that tet)
            for i, side_tet_collection in enumerate(side_tet_collections_at_edge):
                if flow_step in side_tet_collection:
                    downward_path = side_tet_collection[:side_tet_collection.index(flow_step) + 1]
                    downward_path_faces = side_face_collections_at_edge[i][:side_tet_collection.index(flow_step) + 1]
            assert downward_path != None
            for j, (tet_num, edge_num) in enumerate(downward_path): 
                lower_boundary_triangles = [t for t in edge.boundary_triangles() if not t.is_upper and t.index == downward_path_faces[j][0]] 
                lower_boundary_triangles = [t for t in edge.boundary_triangles() if not t.is_upper and t.index == downward_path_faces[j][0]] 
                assert len(lower_boundary_triangles) <= 1
                if len(lower_boundary_triangles) == 1:
                    self.continent.bury(lower_boundary_triangles[0])
                for tet in edge.side_tetrahedra:
                    if tet.index == flow_step[0]:
                        if tet.edge(flow_step[1]) == edge:
                            new_tet = tet
                            break
        assert new_tet != None
        self.tetrahedra.insert(0, new_tet)

    def ensure_contains_one_cycle_up(self):
        """Make sure that the interval contains an entire cycle above init_tet"""
        if self.how_far_up() < len(self.flow_cycle):
            for i in range(len(self.flow_cycle) - self.how_far_up()):
                self.go_up()

    def ensure_contains_one_cycle_down(self):
        """Make sure that the interval contains an entire cycle below init_tet"""
        if self.how_far_down() < len(self.flow_cycle):
            for i in range(len(self.flow_cycle) - self.how_far_down()):
                self.go_down()

    def is_inside_edge_rectangle_green_sides(self, con_edge):
        con_edge.ensure_continent_contains_rectangle()
        greens, purples = con_edge.green_purple_rectangle_sides()
        for leaf in greens:
            for v in self.tetrahedra[-1].vertices:
                if leaf.sees_to_its_left(v) == con_edge.is_red:
                    return False
        return True

    def is_inside_edge_rectangle_purple_sides(self, con_edge):
        con_edge.ensure_continent_contains_rectangle()
        greens, purples = con_edge.green_purple_rectangle_sides()                
        for leaf in purples:
            for v in self.tetrahedra[0].vertices:
                if leaf.sees_to_its_left(v) != con_edge.is_red:
                    return False
        return True

    def is_inside_edge_rectangle(self, con_edge):
        return self.is_inside_edge_rectangle_green_sides(con_edge) and self.is_inside_edge_rectangle_purple_sides(con_edge)

    def find_edge_rectangle_we_are_inside(self, tet):
        """Given tet in the fund dom and in our flow cycle, find an edge of tet whose rectangle contains the flow interval's point"""  
        for e in tet.edges():
            if self.is_inside_edge_rectangle(e):
                return e
        assert False

def high_enough(fund_dom_edges, interval):  ### high enough to know where the top tet of flow interval is relative to all the fund dom edges
    for e in fund_dom_edges:
        greens, purples = e.green_purple_rectangle_sides()
        for leaf in greens:
            if not leaf.is_entirely_to_one_side_of(interval.tetrahedra[-1]):
                return False
    return True

def low_enough(fund_dom_edges, interval): ### low enough to know where the bottom tet of flow interval is relative to all the fund dom edges
    for e in fund_dom_edges:
        greens, purples = e.green_purple_rectangle_sides()
        for leaf in purples:
            if not leaf.is_entirely_to_one_side_of(interval.tetrahedra[0]):
                return False
    return True

def make_continent_drill_flow_cycle(veering_isosig, flow_cycle, num_steps = 10):
    ### format for loops: it is a list of tuples, 
    ### each tuple is (tet_index, edge_index within this tet that we exit through)

    con, continent_fund_dom_tets = make_continent_fund_dom(veering_isosig)
    vt = con.vt

    ## install into vt:
    vt.side_face_collections, vt.side_tet_collections = edge_side_face_collections(vt.tri, vt.angle, tet_vert_coorientations = vt.coorientations, return_tets = True, order_bottom_to_top = False)
    # print('sfc', side_face_collections)
    # print('stc', side_tet_collections)
    ### identify the next edge in the cycle 

    # print([t.index for t in continent_fund_dom_tets])
    init_tetrahedron = continent_fund_dom_tets[flow_cycle[0][0]]
    print('init tet index, flow cycle', init_tetrahedron.index, flow_cycle)
    assert init_tetrahedron.index == flow_cycle[0][0]
    init_edge = flow_edge_in_continent(init_tetrahedron, flow_cycle[0][1])

    interval = flow_interval(con, flow_cycle, init_tetrahedron, init_edge, 0)

    fund_dom_edges = set([])
    for tet in continent_fund_dom_tets:
        for e in tet.edges():
            fund_dom_edges.add(e)
    for e in fund_dom_edges:
        e.ensure_continent_contains_rectangle()
    max_steps = 100
    step = 0
    while not high_enough(fund_dom_edges, interval):
        interval.go_up()
        step += 1
        if step > max_steps:
            assert False  ### break runaway while loop
    while not low_enough(fund_dom_edges, interval):
        interval.go_down()
        step += 1
        if step > max_steps:
            assert False  ### break runaway while loop
    ### next find translates of init_tetrahedron along flow interval up one cycle and down one cycle
    interval.ensure_contains_one_cycle_up()
    interval.ensure_contains_one_cycle_down()
    up_translate_tet = interval.get_tet_at_position(len(flow_cycle))
    up_translate_edge = interval.get_edge_at_position(len(flow_cycle))
    down_translate_tet = interval.get_tet_at_position(-len(flow_cycle))



    ### c and gamma^2c are corresponding cusps of up_translate_tet, down_translate_tet

    ### next determine a quadrant that contains the puncture p.
    ### choose a vertex of top edge of up_translate_tet
    up_v = up_translate_tet.upper_edge.vertices[0]
    if up_v not in up_translate_edge.vertices or (up_translate_edge == up_translate_tet.upper_edge):
        ### then quadrant starts at v and contains edge rect for up_translate_tet.upper_edge
        up_quadrant_edge = up_translate_tet.upper_edge
    else:
        ### then quadrant starts at v and contains edge rect for up_translate_edge
        up_quadrant_edge = up_translate_edge

    v_num = up_translate_tet.vertices.index(up_v)
    # other_end_vert_num = up_translate_tet.vertices.index(up_quadrant_edge.other_end(up_v))
    e_num = up_translate_tet.ordered_edges().index(up_quadrant_edge)

    down_v = down_translate_tet.vertices[v_num]
    down_quadrant_edge = down_translate_tet.edge(e_num)

    up_quadrant_edge.ensure_continent_contains_rectangle()
    down_quadrant_edge.ensure_continent_contains_rectangle()

    up_rect_sides = up_quadrant_edge.rectangle_sides()
    up_quadrant_sides = [leaf for leaf in up_rect_sides if leaf.cusp == up_v]
    assert len(up_quadrant_sides) == 2

    down_rect_sides = down_quadrant_edge.rectangle_sides()
    down_quadrant_sides = [leaf for leaf in down_rect_sides if leaf.cusp == down_v]
    assert len(down_quadrant_sides) == 2
    
    # leaves_to_draw = up_quadrant_sides + down_quadrant_sides
    leaves_to_draw = up_quadrant_sides[:]

    for i in range(2):
        interval.go_up()    ### simulate extending the flow interval to get more cusps
        # interval.go_down()

    found_the_first_strip_vertex = False ### first strip vertex is b_NW in the notes
    while not found_the_first_strip_vertex:
        ### get candidate triangles with special vertex (is a corner of the face rectangle) 
        ### at up_v that are also in the quadrant
        candidate_triangles = []
        for f in up_v.triangles_with_special_vertex_here:
            a, b = f.vertices[1:]  ### anticlockwise order as viewed from above is special vertex, then a, then b
            l, m = up_quadrant_sides
            ### The other two vertices should be in the quadrant:
            if l.sees_to_its_left(a) != m.sees_to_its_left(a):
                assert l.sees_to_its_left(b) != m.sees_to_its_left(b)
                candidate_triangles.append(f)

        for f in candidate_triangles:
            # ### One of the other two vertices should be west (same side as up_v) of p (and this will be b_NW) while the other is east (opposite side from up_v) of p

            leaf_dict = f.get_non_special_vertex_cusp_leaves()
            a, b = f.vertices[1:]

            if leaf_dict[a]['internal'].is_upper:
                internal_leaf = leaf_dict[a]['internal']
                boundary_leaves = leaf_dict[b]['boundary']
                candidate_strip_vertex = a
            else:
                internal_leaf = leaf_dict[b]['internal']
                boundary_leaves = leaf_dict[a]['boundary']
                candidate_strip_vertex = b

    ### first_strip_vertex = b_NW is an a or b of some candidate triangle such that:
    ### green leaf is an inner leaf from a (say) that has up_v on one side and p on the other
    ### green leaves on b (say) are outer leaves for the face rect, and they have up_v and p on same side

            if (internal_leaf.weakly_separates([up_v], interval.tetrahedra[-1].vertices) and
                boundary_leaves[0].weakly_separates(interval.tetrahedra[-1].vertices + [up_v], []) and
                boundary_leaves[1].weakly_separates(interval.tetrahedra[-1].vertices + [up_v], [])):
                first_strip_vertex = candidate_strip_vertex
                found_the_first_strip_vertex = True
                leaves_to_draw.append(internal_leaf)
                break
        if not found_the_first_strip_vertex:
            interval.go_up()




    # leaves_to_draw.append(leaf_dict[a][0])
    # leaves_to_draw.append(leaf_dict[b][0])




            

    triangles_to_draw = candidate_triangles

    #### just build up and down some distance
    # for i in range(num_steps):
    #     if i%2 == 0: ### go up
    #         interval.go_up()
    #     else: ### go down
    #         interval.go_down()
    con.build_boundary_data()
    return con, interval, continent_fund_dom_tets, leaves_to_draw, triangles_to_draw

def complete_tetrahedron_rectangles(con, tetrahedra_to_complete):
    """grow the continent so that the given tetrahedra have full tetrahedron rectangles within the continent"""
    k = 0
    for tet in tetrahedra_to_complete:
        for v in tet.vertices:
            # print('tet vert age', con.vertices.index(v))
            # con.build_boundary_data()
            # con.install_thorn_ends()
            sides = tet_rectangle_sides(tet, v)
            for direction in range(2):
                while sides[direction] == None:
                    # print('direction, k', direction, k)
                    e = con.coastal_edges[(v.coastal_index() - direction)%len(con.coast)] 
                    triangles = e.boundary_triangles()  ### grow around this edge
                    if triangles[0].is_upper != (k % 2 == 0): ### xor, alternate which one we add to
                        con.bury(triangles[0])
                    else:
                        con.bury(triangles[1])
                    # con.build_boundary_data()
                    # con.install_thorn_ends()
                    sides = tet_rectangle_sides(tet, v)
                    k += 1
                    if k > 50:
                        print('bail')
                        return None

def get_fund_domain_tetrahedra(con):
    num_tet = con.vt.tri.countTetrahedra()
    fund_dom_tets = list(range(num_tet))
    number_to_find = num_tet
    for con_tet in con.tetrahedra:
        if con_tet.index in fund_dom_tets:
            fund_dom_tets[con_tet.index] = con_tet  ### replace index with the continent tetrahedron
            number_to_find -= 1
        if number_to_find == 0:
            break
    if number_to_find > 0: ## should have found the whole continent
        print('did not find all of the downstairs tetrahedra!')
    update_fund_dom_tet_nums(con, fund_dom_tets)
    return fund_dom_tets

def update_fund_dom_tet_nums(con, fund_dom_tets):
    for v in con.coast:
        v.fund_dom_tet_nums = []
    for con_tet in fund_dom_tets:
        if type(con_tet) == continent_tetrahedron:  ### could be an integer if we didnt find this tet
            for v in con_tet.vertices:
                v.fund_dom_tet_nums.append(con_tet.index)