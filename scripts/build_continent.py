#
# build_continent.py
#

### functions for building continents

from veering.taut import isosig_to_tri_angle
from veering.veering_tri import veering_triangulation
from veering.flow_cycles import is_twisted

from continent import continent, continent_tetrahedron
from boundary_triangulation import tet_face
from flow_interval import flow_interval, uniquify_list_of_flow_intervals, translate_of_interval_from_one_edge_rect_to_another, cmp_to_key

from draw_veering_triangulation_and_mid_annuli import get_consistent_tet_vert_posns

def make_continent_naive(veering_isosig, max_num_tetrahedra = 50):
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle) #, tet_shapes = tet_shapes)
    initial_tet_face = tet_face(vt, 0, 0, verts_pos = [None, None, None, None])
    con = continent( vt, initial_tet_face) #, desired_vertices = desired_vertices )
    con.build_naive(max_num_tetrahedra = max_num_tetrahedra)
    con.build_boundary_data()
    print(len(con.vertices), len(con.edges), len(con.triangles), len(con.tetrahedra))
    return con

def make_continent_fund_dom(veering_isosig, max_num_tetrahedra = 50, for_drawing = False):
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle) #, tet_shapes = tet_shapes)
    initial_tet_face = tet_face(vt, 0, 0, verts_pos = [None, None, None, None])
    con = continent( vt, initial_tet_face, for_drawing = for_drawing) #, desired_vertices = desired_vertices )
    continent_fund_dom_tets = con.build_triangulation_fundamental_domain(max_num_tetrahedra = max_num_tetrahedra)
    con.build_boundary_data()
    # print(len(con.vertices), len(con.edges), len(con.triangles), len(con.tetrahedra))
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


def find_strip_vertex(v, quadrant_sides, interval, is_upper = True):
    ### strip vertex is b_NW in the notes
    while True:
        ### get candidate triangles with special vertex (is a corner of the face rectangle) 
        ### at v that are also in the quadrant
        candidate_triangles = []
        for f in v.triangles_with_special_vertex_here:
            a, b = f.vertices[1:]  ### anticlockwise order as viewed from above is special vertex, then a, then b
            l, m = quadrant_sides
            ### The other two vertices should be in the quadrant:
            if l.sees_to_its_left(a) != m.sees_to_its_left(a):
                assert l.sees_to_its_left(b) != m.sees_to_its_left(b)
                candidate_triangles.append(f)

        for f in candidate_triangles:
            # ### One of the other two vertices should be west (same side as up_v) of p (and this will be b_NW) while the other is east (opposite side from up_v) of p
            leaf_dict = f.get_non_special_vertex_cusp_leaves()
            a, b = f.vertices[1:]

            if leaf_dict[a]['internal'].is_upper == is_upper:
                internal_leaf = leaf_dict[a]['internal']
                boundary_leaves = leaf_dict[b]['boundary']
                other_colour_boundary_leaves = leaf_dict[a]['boundary']
                candidate_strip_vertex = a
                candidate_index = 1
            else:
                internal_leaf = leaf_dict[b]['internal']
                boundary_leaves = leaf_dict[a]['boundary']
                other_colour_boundary_leaves = leaf_dict[b]['boundary']
                candidate_strip_vertex = b
                candidate_index = 2
            assert internal_leaf.is_upper == is_upper

    ### upper_strip_vertex = b_NW is an a or b of some candidate triangle such that:
    ### green leaf is an inner leaf from a (say) that has up_v on one side and p on the other
    ### purple leaves on b (say) are outer leaves for the face rect, and they have up_v and p on same side

            if is_upper:
                tet = interval.tetrahedra[-1]
            else:
                tet = interval.tetrahedra[0]
            if (internal_leaf.weakly_separates([v], tet.vertices) and
                boundary_leaves[0].weakly_separates(tet.vertices + [v], []) and
                boundary_leaves[1].weakly_separates(tet.vertices + [v], [])):
                # print('found strip vertex')

                ### Now candidate_strip_vertex is the strip_vertex. internal leaf is one side of the quadrant at the strip vertex that contains tet

                return (candidate_strip_vertex, internal_leaf, other_colour_boundary_leaves, f, candidate_index)
 
        if is_upper:
            interval.go_up()
        else:
            # print('go down')
            interval.go_down()
    # print('did not find strip vertex')
    # return (None, None, candidate_triangles) ### testing

def make_continent_drill_flow_cycle(veering_isosig, flow_cycle, use_untwisted_speed_up = True, verbose = 0, for_drawing = False):
    ### format for loops: it is a list of tuples, 
    ### each tuple is (tet_index, edge_index within this tet that we exit through)

    con, continent_fund_dom_tets = make_continent_fund_dom(veering_isosig, for_drawing = for_drawing)
    con.for_drawing = for_drawing
    vt = con.vt
    flow_cycle_is_twisted = is_twisted(vt, flow_cycle)
    if verbose > 0:
        print('drill', veering_isosig, flow_cycle, 'is_twisted:', flow_cycle_is_twisted, 'use_untwisted_speed_up:', use_untwisted_speed_up)

    found_parallel = False ### in building the continent, do we find fellow-travelling flow_cycles (if so then drilling naively will not give a hyperbolic manifold)

    ## install into vt:
    # vt.side_face_collections, vt.side_tet_collections = edge_side_face_collections(vt.tri, vt.angle, tet_vert_coorientations = vt.coorientations, return_tets = True, order_bottom_to_top = False)
    # print('sfc', side_face_collections)
    # print('stc', side_tet_collections)
    ### identify the next edge in the cycle 

    # print([t.index for t in continent_fund_dom_tets])
    init_tetrahedron = continent_fund_dom_tets[flow_cycle[0][0]]
    # print('init tet index, flow cycle', init_tetrahedron.index, flow_cycle)
    assert init_tetrahedron.index == flow_cycle[0][0]

    main_interval = flow_interval(con, flow_cycle, init_tetrahedron, 0)

    fund_dom_edges = []
    for tet in continent_fund_dom_tets:
        for e in tet.edges():
            fund_dom_edges.append(e)
    for e in fund_dom_edges:
        e.ensure_continent_contains_rectangle()

    ### next find translates of init_tetrahedron along flow interval up one cycle and down one cycle
    main_interval.ensure_contains_one_cycle_up()
    if flow_cycle_is_twisted or not use_untwisted_speed_up:
        main_interval.ensure_contains_one_cycle_down()  ### does not make much difference experimentally

    if main_interval.is_boundary_parallel():  ### FIX?
        if verbose >= 1:
            print('flow cycle is boundary parallel')
        return None
    up_translate_tet = main_interval.get_tet_at_position(len(flow_cycle))
    up_translate_edge = main_interval.get_edge_at_position(len(flow_cycle))
    if flow_cycle_is_twisted or not use_untwisted_speed_up:
        down_translate_tet = main_interval.get_tet_at_position(-len(flow_cycle))
    else:
        down_translate_tet = main_interval.get_tet_at_position(0)

    ### c and gamma^2c are corresponding cusps of up_translate_tet, down_translate_tet
    ### (or just gamma c if the flow cycle is not twisted)

    ### next determine a quadrant that contains the puncture p.
    ### choose a vertex of top edge of up_translate_tet
    up_v = up_translate_tet.upper_edge.vertices[0] ### this is gamma^2c in the notes
    if up_v not in up_translate_edge.vertices or (up_translate_edge == up_translate_tet.upper_edge):
        ### then quadrant starts at up_v and contains edge rect for up_translate_tet.upper_edge
        up_quadrant_edge = up_translate_tet.upper_edge
    else:
        ### then quadrant starts at up_v and contains edge rect for up_translate_edge
        up_quadrant_edge = up_translate_edge

    v_num = up_translate_tet.vertices.index(up_v)
    # other_end_vert_num = up_translate_tet.vertices.index(up_quadrant_edge.other_end(up_v))
    e_num = up_translate_tet.ordered_edges().index(up_quadrant_edge)

    down_v = down_translate_tet.vertices[v_num] ### this is c in the notes
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
    leaves_to_draw = []
    # leaves_to_draw.extend( up_quadrant_sides[:] + down_quadrant_sides[:] )

    # for i in range(2):
    #     main_interval.go_up()    
### going up a few times before starting to search properly might be more efficient
### because the search below for the upper_strip_vertex involves get_non_special_vertex_cusp_leaves()
### which calls ensure_continent_contains_rectangle, which adds tetrahedra and (it seems) more candidate_triangles
### Maybe there is a more efficient way to check for upper_strip_vertex that doesn't require generating the leaves at 
### the non-special vertices of the candidate triangles? Just look at the position of the cusps on the circle?

        # main_interval.go_down()

    upper_strip_vertex, up_internal_leaf, up_boundary_leaves, up_f, up_candidate_index = find_strip_vertex(up_v, up_quadrant_sides, main_interval, is_upper = True)
    # leaves_to_draw.append(up_internal_leaf)
    # leaves_to_draw.extend(up_boundary_leaves)
    triangles_to_draw = []
    # triangles_to_draw.append(up_f)

    lower_strip_vertex, down_internal_leaf, down_boundary_leaves, down_f, down_candidate_index = find_strip_vertex(down_v, down_quadrant_sides, main_interval, is_upper = False)

    ### Now find the quadrant at upper_strip_vertex that contains the flow line

    ### up_f.edges[0] is the special edge of the face. other_colour_quadrant_sides are the two sides of that edge rect based at upper_strip_vertex

    up_quadrant_sides2 = up_f.edges[0].rectangle_sides_by_vertex()[upper_strip_vertex]

    # leaves_to_draw.extend(up_quadrant_sides2)  # 2 for the edges of opposite slope

    ### Now we have a quadrant based at upper_strip_vertex (b_NW in the notes) that contains p
    up_w = upper_strip_vertex ### rename things to make parallelism clearer


    ### We also need down_w, which is the translate of up_w by \gamma^{-2}

    if up_f.is_upper:
        up_w_tet = up_f.lower_tet
    else:
        up_w_tet = up_f.upper_tet
    f_num = up_w_tet.ordered_faces().index(up_f)

    face_num_path = up_translate_tet.face_num_path_to_other_tet(up_w_tet)
    down_w_tet = down_translate_tet.follow_face_num_path(face_num_path)
    down_f = down_w_tet.ordered_faces()[f_num]
    down_w = down_f.vertices[up_candidate_index]
    down_quadrant_sides2 = down_f.edges[0].rectangle_sides_by_vertex()[down_w]
    # leaves_to_draw.extend(down_quadrant_sides2)

    # triangles_to_draw.append(down_f) 

    upper_strip_vertex2, up_internal_leaf2, up_boundary_leaves2, up_f2, up_candidate_index2 = find_strip_vertex(up_w, up_quadrant_sides2, main_interval, is_upper = True)
    lower_strip_vertex2, down_internal_leaf2, down_boundary_leaves2, down_f2, down_candidate_index2 = find_strip_vertex(down_w, down_quadrant_sides2, main_interval, is_upper = False)

    # triangles_to_draw.append(up_f2) 
    # triangles_to_draw.append(down_f2) 
    # leaves_to_draw.append(up_internal_leaf2)
    # leaves_to_draw.append(down_internal_leaf2)

       

    # print('up_v', up_v, 'down_v', down_v, 'up_w', up_w, 'down_w', down_w, 'lower_strip_vertex', lower_strip_vertex, 'upper_strip_vertex2', upper_strip_vertex2, 'lower_strip_vertex2', lower_strip_vertex2)

    ### Continent is now big enough to contain a translate of any edge whose rectangle contains the puncture.

    candidate_drilled_edges = con.edges[:] # shallow copy
    drilled_continent_edges = []
            
    ### find edges in continent that have a puncture
    for e in candidate_drilled_edges:
        if main_interval.is_inside_edge_rectangle(e):
            drilled_continent_edges.append(e)
            # print('main interval found in edge with index', e.index)
            # print('edge is red?', e.is_red)
    # print('continent_size', len(con.tetrahedra))

    ### translate punctures to fund dom edges (those edges below fund dom tets)
    flow_intervals = []
    for e in drilled_continent_edges:
        if e.upper_tet == None:
            e.ensure_continent_contains_tet_above()
        e2 = continent_fund_dom_tets[e.upper_tet.index].lower_edge
        new_interval = translate_of_interval_from_one_edge_rect_to_another(e, e2, main_interval)
        flow_intervals.append( new_interval )

    for interval in flow_intervals:
        interval.extend_within_continent()

    flow_intervals, found_parallel_this_time = uniquify_list_of_flow_intervals(flow_intervals)
    if not found_parallel and found_parallel_this_time:
        found_parallel = True
    if verbose >= 2:
        print('num unique flow_intervals', len(flow_intervals))

    fund_dom_unique_edges = [t.lower_edge for t in continent_fund_dom_tets]
    fund_dom_unique_edge_indices = [e.index for e in fund_dom_unique_edges]

    ### figure out which edge rectangles contain which unique punctures
    intervals_inside_edge_rectangles = []
    for e in fund_dom_unique_edges:
        flow_intervals_in_e = []
        for interval in flow_intervals:
            if interval.is_inside_edge_rectangle(e):
                flow_intervals_in_e.append(interval)
        intervals_inside_edge_rectangles.append(flow_intervals_in_e)
    for i, intervals in enumerate(intervals_inside_edge_rectangles):
        e = fund_dom_unique_edges[i]
        if verbose >= 2:
            print('edge below tet', i, 'is index', e.index, 'is red', e.is_red, 'contains intervals', len(intervals))

    ### copy punctures around to get all punctures going through a tet rect of fund dom

    # print(len(continent_fund_dom_tets))
    intervals_inside_tet_rectangles = []
    for i, t in enumerate(continent_fund_dom_tets):
        flow_intervals_in_t = []
        for interval in intervals_inside_edge_rectangles[i]: ### start with intervals in the lower edge of t
            flow_intervals_in_t.append(interval.copy())  ### we don't want the same interval showing up in multiple tetrahedra - naming them gets messy
        edges_to_add = t.equatorial_edges
        # edges_to_add.append(t.upper_edge)  ### for testing, this should already be covered by other tet rects
        for e_eq in edges_to_add:   ### then add intervals in equatorial edges
            j = fund_dom_unique_edge_indices.index(e_eq.index)
            e = fund_dom_unique_edges[j]  
            flow_intervals_in_e = intervals_inside_edge_rectangles[j]
            for interval in flow_intervals_in_e:
                new_interval = translate_of_interval_from_one_edge_rect_to_another(e, e_eq, interval) ### these are copies as well
                flow_intervals_in_t.append(new_interval)
        for interval in flow_intervals_in_t:
            interval.extend_within_continent()
        flow_intervals_in_t, found_parallel_this_time = uniquify_list_of_flow_intervals(flow_intervals_in_t)
        if not found_parallel and found_parallel_this_time:
            found_parallel = True
        intervals_inside_tet_rectangles.append(flow_intervals_in_t)

    # for i, intervals in enumerate(intervals_inside_tet_rectangles):
    #     print('tet', i, 'contains intervals', len(intervals))   

    ### Now start sorting the intervals inside the tet rectangles.
    ### Do horizontal order first. Start with positions of intervals relative to the S, N cusps

    tet_vert_posns, _, _, _ = get_consistent_tet_vert_posns(vt.tri, vt.angle, vt.tet_types, vt.coorientations)
    vt.tet_vert_posns = tet_vert_posns

    ###           top[0]
    ###          /   |   \
    ### bottom[0]--- | ---bottom[1]
    ###          \   |   /
    ###           top[1]
    
    tetrahedra_cusp_orders = []
    tetrahedra_chunks = []
    for i, t in enumerate(continent_fund_dom_tets):
        top_vertices, bottom_vertices = tet_vert_posns[i]
        N_ind, S_ind = top_vertices
        W_ind, E_ind = bottom_vertices

        S_leaf, N_leaf = t.upper_edge.green_purple_rectangle_sides()[0] ## green sides are the cusp_leaves coming from S, N cusps
        S, N = S_leaf.cusp, N_leaf.cusp
        if t.vertices[N_ind] != N:
            assert t.vertices[N_ind] == S
            S_leaf, N_leaf = N_leaf, S_leaf
            S, N = N, S

        W_leaf, E_leaf = t.lower_edge.green_purple_rectangle_sides()[1] 
        W, E = W_leaf.cusp, E_leaf.cusp
        if S_leaf.sees_to_its_left(E):
            W, E = E, W
            W_leaf, E_leaf = E_leaf, W_leaf
        assert S_leaf.sees_to_its_left(W) and not S_leaf.sees_to_its_left(E)
        assert t.vertices[S_ind] == S
        assert t.vertices[N_ind] == N
        assert t.vertices[W_ind] == W
        assert t.vertices[E_ind] == E

        if S_leaf.sees_to_its_left(N):
            assert N_leaf.sees_to_its_left(S)
            horizontal_cusp_order = [W, N, S, E]
        else:
            assert not N_leaf.sees_to_its_left(S)
            horizontal_cusp_order = [W, S, N, E]
        if W_leaf.sees_to_its_left(E):
            assert E_leaf.sees_to_its_left(W)
            vertical_cusp_order = [S, W, E, N]
        else:
            assert not E_leaf.sees_to_its_left(W)
            vertical_cusp_order = [S, E, W, N]

        ### Sanity check:
        horiz_to_vert_perm = [vertical_cusp_order.index(horizontal_cusp_order[i]) for i in range(4)]
        if t.upper_edge.is_red:
            if t.lower_edge.is_red:
                assert horiz_to_vert_perm == [1,0,3,2]
            else:
                assert horiz_to_vert_perm == [2,0,3,1]
        else:
            if t.lower_edge.is_red:
                assert horiz_to_vert_perm == [1,3,0,2]
            else:
                assert horiz_to_vert_perm == [2,3,0,1]

        tetrahedra_cusp_orders.append([horizontal_cusp_order, vertical_cusp_order])

        # print('tet', i, 'cusp orders', horizontal_cusp_order, vertical_cusp_order)
        east_edge = None
        north_edge = None
        west_edge = None
        south_edge = None
        for e in t.equatorial_edges:
            # print('e.vertices', e.vertices)
            if horizontal_cusp_order[0] in e.vertices and horizontal_cusp_order[1] in e.vertices:
                west_edge = e
            if horizontal_cusp_order[2] in e.vertices and horizontal_cusp_order[3] in e.vertices:
                east_edge = e
            if vertical_cusp_order[2] in e.vertices and vertical_cusp_order[3] in e.vertices:
                north_edge = e 
            if vertical_cusp_order[0] in e.vertices and vertical_cusp_order[1] in e.vertices:
                south_edge = e
        # print(east_edge, north_edge, west_edge, south_edge)
        assert east_edge != None and west_edge != None and north_edge != None and south_edge != None

        # print('tet', i, 'horizontal_cusp_order', horizontal_cusp_order, 'vertical_cusp_order', vertical_cusp_order)
        flow_intervals_in_t = intervals_inside_tet_rectangles[i]
        west_intervals = []
        we_middle_intervals = []
        east_intervals = []

        south_intervals = []
        sn_middle_intervals = []
        north_intervals = []

        for interval in flow_intervals_in_t:
            if interval.is_inside_edge_rectangle_green_sides(t.upper_edge):
                we_middle_intervals.append(interval)
                interval.we_chunk_num = 1
            elif interval.is_inside_edge_rectangle_green_sides(west_edge):
                west_intervals.append(interval)
                interval.we_chunk_num = 0
            else:
                assert interval.is_inside_edge_rectangle_green_sides(east_edge)
                east_intervals.append(interval)
                interval.we_chunk_num = 2

            if interval.is_inside_edge_rectangle_purple_sides(t.lower_edge):
                sn_middle_intervals.append(interval)
                interval.sn_chunk_num = 1
            elif interval.is_inside_edge_rectangle_purple_sides(south_edge):
                south_intervals.append(interval)
                interval.sn_chunk_num = 0
            else:
                assert interval.is_inside_edge_rectangle_purple_sides(north_edge)
                north_intervals.append(interval)
                interval.sn_chunk_num = 2            

        # print('tet', i, 'contains west, middle, east intervals', len(west_intervals), len(we_middle_intervals), len(east_intervals)) 
        # print('tet', i, 'contains south, middle, north intervals', len(south_intervals), len(sn_middle_intervals), len(north_intervals))
        
        for interval in west_intervals:
            assert interval.we_chunk_num == 0
        for interval in we_middle_intervals:
            assert interval.we_chunk_num == 1
        for interval in east_intervals:
            assert interval.we_chunk_num == 2
        for interval in south_intervals:
            assert interval.sn_chunk_num == 0
        for interval in sn_middle_intervals:
            assert interval.sn_chunk_num == 1
        for interval in north_intervals:
            assert interval.sn_chunk_num == 2

        # for k, interval in enumerate(west_intervals):
        #     print('west interval ' + str(k), interval)
        #     for tet in interval.tetrahedra:
        #         print(tet)
        we_chunks = [west_intervals, we_middle_intervals, east_intervals]
        for j, we_chunk in enumerate(we_chunks):
            # print('we_chunk', j, we_chunk)
            if len(we_chunk) > 1: ### sort within the chunk
                chunk_W_cusp, chunk_E_cusp = horizontal_cusp_order[j], horizontal_cusp_order[j+1]
                we_chunk.sort(key=cmp_to_key(chunk_W_cusp, chunk_E_cusp, is_horizontal = True))
            # for k, interval in enumerate(we_chunk):
            #     interval.name = str(i) + str(j) + str(k)
        # print('tet', i, 'we_chunks', we_chunks[0], we_chunks[1], we_chunks[2])
            # print('sorted_we_chunk', j, we_chunk)

        sn_chunks = [south_intervals, sn_middle_intervals, north_intervals]
        for j, sn_chunk in enumerate(sn_chunks):
            # print(j, sn_chunk)
            if len(sn_chunk) > 1: ### sort within the chunk
                chunk_S_cusp, chunk_N_cusp = vertical_cusp_order[j], vertical_cusp_order[j+1]
                sn_chunk.sort(key=cmp_to_key(chunk_S_cusp, chunk_N_cusp, is_horizontal = False))
        # print('tet', i, 'sn_chunks', sn_chunks[0], sn_chunks[1], sn_chunks[2])

        tetrahedra_chunks.append([we_chunks, sn_chunks])

        for interval in flow_intervals_in_t:
            interval.we_in_chunk_index = we_chunks[interval.we_chunk_num].index(interval)
            interval.sn_in_chunk_index = sn_chunks[interval.sn_chunk_num].index(interval)
            # print('tet', i, 'wechn', interval.we_chunk_num, 'weind', interval.we_in_chunk_index, 'snchn', interval.sn_chunk_num, 'snind', interval.sn_in_chunk_index)
            interval.owning_tet_index = i
            interval.name = str(i) + str(interval.we_chunk_num) + str(interval.we_in_chunk_index) + str(interval.sn_chunk_num) + str(interval.sn_in_chunk_index) 

    all_intervals = []
    for ints in intervals_inside_tet_rectangles:
        all_intervals.extend(ints)
    for interval in all_intervals:
        interval.extend_within_continent()

    con.build_boundary_data()
    if verbose >= 1:
        print('continent num tetrahedra:', len(con.tetrahedra))
    return con, tetrahedra_cusp_orders, tetrahedra_chunks, intervals_inside_tet_rectangles, all_intervals, continent_fund_dom_tets, leaves_to_draw, triangles_to_draw, found_parallel

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