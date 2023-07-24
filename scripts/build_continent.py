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

### this is probably broken, needs to be redone anyway - move around the universal cover building onto the continent as necessary
def make_continent_drill_dual_cycle(veering_isosig, dual_cycle, num_steps):
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle) #, tet_shapes = tet_shapes)
    ### initial tetrahedron is above face 0 of the dual cycle
    face0 = vt.tri.triangles()[dual_cycle[0]]
    embeds = face0.embeddings()
    tet_0, face_0 = None, None
    for embed in embeds:
        tet_num = embed.simplex().index()
        face_num = embed.face()
        if vt.coorientations[tet_num][face_num] == -1: ## into the tet
            tet_0, face_0 = tet_num, face_num
            break
    assert tet_0 != None and face_0 != None
    initial_tet_face = tet_face(vt, tet_0, face_0, verts_pos = [None, None, None, None])
    con = continent( vt, initial_tet_face)

    ### identify the triangle corresponding to dual_cycle[1]
    for triangle in con.triangles:
        if not triangle.is_upper and triangle.index == dual_cycle[0]:
            lowest_triangle = triangle
        if triangle.is_upper and triangle.index == dual_cycle[1]:
            highest_triangle = triangle
    triangle_path = [lowest_triangle, highest_triangle]
    lowest_path_index = 0
    highest_path_index = 1
    # for i in range(num_steps):
    #     last_added_tet = con.bury(highest_triangle)

    #     ### find next_triangle
    #     highest_triangle = None
    #     for triangle in last_added_tet.upper_triangles:
    #         if triangle.index == dual_cycle[(path_index + 1) % len(dual_cycle)]:
    #             highest_triangle = triangle
    #             break
    #     assert highest_triangle != None
    #     triangle_path.append(highest_triangle)
    #     path_index = path_index + 1

    #     con.build_boundary_data()

    for i in range(num_steps):
        if i%2 == 0:
            last_added_tet = con.bury(triangle_path[-1]) ## go up
            for triangle in last_added_tet.upper_triangles:
                if triangle.index == dual_cycle[(highest_path_index + 1) % len(dual_cycle)]:
                    triangle_path.append(triangle)
                    break
            highest_path_index = highest_path_index + 1
        else:
            last_added_tet = con.bury(triangle_path[0]) ## go down
            for triangle in last_added_tet.lower_triangles:
                if triangle.index == dual_cycle[(lowest_path_index - 1) % len(dual_cycle)]:
                    triangle_path.insert(0, triangle)
                    break
            lowest_path_index = lowest_path_index - 1
    con.build_boundary_data()
    con.triangle_path = triangle_path
    return con

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

def go_up_flow(flow_cycle, con, flow_edges, flow_tetrahedra, upwards_flow_index):
    edge = flow_edges[-1]
    while True:
        con.build_boundary_data()
        upper_boundary_triangles = [t for t in edge.boundary_triangles if t.is_upper] 
        if len(upper_boundary_triangles) == 0:
            break ## out of while loop
        con.bury(upper_boundary_triangles[0])
    last_tet = edge.upper_tet
    assert last_tet != None
    upwards_flow_index = (upwards_flow_index + 1) % len(flow_cycle)
    assert last_tet.index == flow_cycle[upwards_flow_index][0]
    flow_tetrahedra.append(last_tet)
    flow_edges.append(flow_edge_in_continent(last_tet, flow_cycle[upwards_flow_index][1]))
    con.build_boundary_data()
    return upwards_flow_index

def go_down_flow(flow_cycle, con, flow_edges, flow_tetrahedra, downwards_flow_index):
    tet = flow_tetrahedra[0]
    edge = tet.lower_edge()
    flow_edges.insert(0, edge) 
    ### now build the continent to get new_lowest_tet
    vt = con.vt
    manifold_edge = vt.tri.edge(edge.index)
    downwards_flow_index = (downwards_flow_index - 1) % len(flow_cycle)
    ### is the flow cycle vertical through the tetrahedron? 
    tet_below = get_tet_above_edge(vt.tri, vt.angle, manifold_edge, tet_vert_coorientations = vt.coorientations, get_tet_below_edge = True)
    tet_below_num = tet_below.index()
    top_vert_nums = get_tet_top_vert_nums(vt.coorientations, tet_below_num)
    top_edge_num = vert_pair_to_edge_num[tuple(top_vert_nums)]

    if (tet_below_num, top_edge_num) == flow_cycle[downwards_flow_index]: ### flow cycle went straight up
        while True:
            con.build_boundary_data()  
            lower_boundary_triangles = [t for t in edge.boundary_triangles if not t.is_upper] 
            if len(lower_boundary_triangles) == 0:
                break ## out of while loop
            last_tet = con.bury(lower_boundary_triangles[0])
    else:
        ### find which side of the edge our tet is in
        # print('edge index', edge.index)
        side_tet_collections_at_edge = vt.side_tet_collections[edge.index] ## index in the manifold
        side_face_collections_at_edge = vt.side_face_collections[edge.index]
        downward_path = None
        flow_step = flow_cycle[downwards_flow_index]
        for i, side_tet_collection in enumerate(side_tet_collections_at_edge):
            if flow_step in side_tet_collection:
                downward_path = side_tet_collection[:side_tet_collection.index(flow_step) + 1]
                downward_path_faces = side_face_collections_at_edge[i][:side_tet_collection.index(flow_step) + 1]
        assert downward_path != None
        for j, (tet_num, edge_num) in enumerate(downward_path):
            con.build_boundary_data()  
            lower_boundary_triangles = [t for t in edge.boundary_triangles if not t.is_upper and t.index == downward_path_faces[j][0]] 
            assert len(lower_boundary_triangles) == 1
            last_tet = con.bury(lower_boundary_triangles[0])
    assert last_tet != None
    flow_tetrahedra.insert(0, last_tet)
    con.build_boundary_data()
    return downwards_flow_index

def make_continent_drill_flow_cycle(veering_isosig, flow_cycle, num_steps):
    ### format for loops: it is a list of tuples, 
    ### each tuple is (tet_index, edge_index within this tet that we exit through)

    con, continent_fund_dom_tets = make_continent_fund_dom(veering_isosig)
    vt = con.vt

    ## install into vt:
    vt.side_face_collections, vt.side_tet_collections = edge_side_face_collections(vt.tri, vt.angle, tet_vert_coorientations = vt.coorientations, return_tets = True, order_bottom_to_top = False)
    # print('sfc', side_face_collections)
    # print('stc', side_tet_collections)
    ### identify the next edge in the cycle 

    print([t.index for t in continent_fund_dom_tets])
    init_tetrahedron = continent_fund_dom_tets[flow_cycle[0][0]]
    print('init tet index, flow cycle', init_tetrahedron.index, flow_cycle)
    assert init_tetrahedron.index == flow_cycle[0][0]
    init_edge = flow_edge_in_continent(init_tetrahedron, flow_cycle[0][1])

    flow_tetrahedra = [init_tetrahedron]
    flow_edges = [init_edge]
    ### both in the continent

    upwards_flow_index = 0
    downwards_flow_index = 0
    for i in range(num_steps):
        # if i%2 == 0: ### go up
        con.build_boundary_data()
        upwards_flow_index = go_up_flow(flow_cycle, con, flow_edges, flow_tetrahedra, upwards_flow_index)

        # else: ### go down
        #     con.build_boundary_data()
        #     downwards_flow_index = go_down_flow(flow_cycle, con, flow_edges, flow_tetrahedra, downwards_flow_index)
    con.build_boundary_data()
    # con.install_thorn_ends()
    return con, flow_tetrahedra, flow_edges

def complete_tetrahedron_rectangles(con, tetrahedra_to_complete):
    """grow the continent so that the given tetrahedra have full tetrahedron rectangles within the continent"""
    k = 0
    for tet in tetrahedra_to_complete:
        for v in tet.vertices:
            # print('tet vert age', con.vertices.index(v))
            con.build_boundary_data()
            # con.install_thorn_ends()
            sides = tet_rectangle_sides(tet, v)
            for direction in range(2):
                while sides[direction] == None:
                    # print('direction, k', direction, k)
                    e = con.coastal_edges[(v.coastal_index() - direction)%len(con.coast)] 
                    triangles = e.boundary_triangles  ### grow around this edge
                    if triangles[0].is_upper != (k % 2 == 0): ### xor, alternate which one we add to
                        con.bury(triangles[0])
                    else:
                        con.bury(triangles[1])
                    con.build_boundary_data()
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