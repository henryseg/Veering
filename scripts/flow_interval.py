from veering.taut import edge_num_to_vert_pair, vert_pair_to_edge_num
from veering.transverse_taut import get_tet_above_edge, get_tet_top_vert_nums

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
    def __init__(self, continent, flow_cycle, init_flow_tetrahedron, init_flow_index):
        self.continent = continent
        self.flow_cycle = flow_cycle
        self.tetrahedra = [init_flow_tetrahedron]
        init_flow_edge = flow_edge_in_continent(init_flow_tetrahedron, flow_cycle[0][1])
        self.edges = [init_flow_edge]
        self.up_index = init_flow_index    ### index in the flow_cycle
        self.down_index = init_flow_index
        self.init_flow_index = init_flow_index  ### never alter after initialization
        self.init_tet = init_flow_tetrahedron  ### never alter after initialization
        self.init_edge = init_flow_edge  ### never alter after initialization
        self.name = 'int_' + self.init_tet.__repr__()

    def __repr__(self):
        return self.name
    
    def __str__(self):
        return self.__repr__()

    def copy(self):
        other = flow_interval(self.continent, self.flow_cycle, self.init_tet, self.init_flow_index)
        other.extend_within_continent()
        return other

    def equals(self, other):
        """assuming both intervals have been extended within the continent, checks if they are the same"""
        if self.flow_cycle != other.flow_cycle:
            return False
        self.extend_within_continent()
        other.extend_within_continent()
        out = (self.tetrahedra[-1] == other.tetrahedra[-1]) and (self.up_index == other.up_index)
        if out:
            # print('interval lengths', len(self.tetrahedra), len(other.tetrahedra))
            # if self.tetrahedra[0] != other.tetrahedra[0]:
            #     print('self init tet index', self.tetrahedra.index(self.init_tet))
            #     print('other init tet index', other.tetrahedra.index(other.init_tet))
            assert self.tetrahedra[0] == other.tetrahedra[0]
            assert self.down_index == other.down_index
        return out

    def ensure_contains_one_cycle(self):
        while len(self.tetrahedra) <= len(self.flow_cycle):
            self.go_up()
            self.go_down()  ### make longer, we don't really care how.

    def fellow_travels(self, other):
        """Do the two flow cycles bound an annulus?"""
        ### This will detect flow intervals based on different flow cycles that bound an annulus. 
        ### But what if some multiple of one flow cycle cobounds an annulus with some multiple of another?
        ### If this happens then the code will run forever, trying to separate the flow intervals
        self.ensure_contains_one_cycle()
        other.ensure_contains_one_cycle()
        this_lowest = self.tetrahedra[0]
        this_one_cycle_up = self.tetrahedra[len(self.flow_cycle)]
        other_lowest = other.tetrahedra[0]
        other_one_cycle_up = other.tetrahedra[len(other.flow_cycle)]
        path = this_lowest.face_num_path_to_other_tet(other_lowest)
        return this_one_cycle_up.follow_face_num_path(path) == other_one_cycle_up

    def is_in_list(self, intervals_list):
        for other in intervals_list:
            if self.equals(other):
                return True, False  ### first is the result, second is whether or not we found a parallel interval
            if self.fellow_travels(other):
                # print('found fellow travelling flow intervals')
                return True, True
        return False, False

    def is_boundary_parallel(self):  ### FIX we should prove that this works... or better check it before building continents
    #     ### use def tri_loop_is_boundary_parallel(tri_loop, tri)?? from veering.flow_cycles
    #     ### if you go straight up through a tetrahedron you do so infinitely many times so you are not boundary parallel
    #     ### if you are boundary parallel then you are trapped inside one ladder of the boundary triangulation, so you 
    #     ### follow the ladder pole slope.
    #     self.ensure_contains_one_cycle_up()
    #     self.ensure_contains_one_cycle_down() ### two cycles should be enough to separate the vertices if not boundary parallel. 
    #     ### One might not be enough if there is a rotation by pi as we go up
    #     lowest = self.tetrahedra[0]
    #     two_cycles_up = self.tetrahedra[2 * len(self.flow_cycle)]
    #     # print(lowest.vertices, one_cycle_up.vertices)
    #     # print(not set(lowest.vertices).isdisjoint(set(one_cycle_up.vertices)))
    #     return not set(lowest.vertices).isdisjoint(set(two_cycles_up.vertices))

        self.ensure_contains_one_cycle()
        lowest = self.tetrahedra[0]
        one_cycle_up = self.tetrahedra[len(self.flow_cycle)]
        for i in range(4):
            if lowest.vertices[i] == one_cycle_up.vertices[i]:
                return True
        return False

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

    def go_up(self, expand_continent = True):
        edge = self.edges[-1]
        if expand_continent:
            edge.ensure_continent_contains_tet_above()
        if edge.upper_tet != None:
            up_tet = edge.upper_tet
            assert up_tet != None
            self.up_index = (self.up_index + 1) % len(self.flow_cycle)
            assert up_tet.index == self.flow_cycle[self.up_index][0]
            self.tetrahedra.append(up_tet)
            self.edges.append(flow_edge_in_continent(up_tet, self.flow_cycle[self.up_index][1]))
            return True
        else:
            return False ### we ran out of continent to extend the flow interval into

    def go_down(self, expand_continent = True):
        # print('going down')
        tet = self.tetrahedra[0]
        the_edge = tet.lower_edge
        ### now build the continent to get new_lowest_tet
        vt = self.continent.vt
        manifold_edge = vt.tri.edge(the_edge.index)
        next_down_index = (self.down_index - 1) % len(self.flow_cycle)
        ### is the flow cycle vertical through the tetrahedron? 
        tet_below = get_tet_above_edge(vt.tri, vt.angle, manifold_edge, tet_vert_coorientations = vt.coorientations, get_tet_below_edge = True)
        tet_below_num = tet_below.index()
        top_vert_nums = get_tet_top_vert_nums(vt.coorientations, tet_below_num)
        top_edge_num = vert_pair_to_edge_num[tuple(top_vert_nums)]

        if (tet_below_num, top_edge_num) == self.flow_cycle[next_down_index]: ### flow cycle went straight up
            # print('straight down case')
            if the_edge.lower_tet == None:
                if expand_continent:
                    while True:
                        lower_boundary_triangles = [t for t in the_edge.boundary_triangles() if not t.is_upper] 
                        if len(lower_boundary_triangles) == 0:
                            break ## out of while loop
                        self.continent.bury(lower_boundary_triangles[0])
                else:
                    return False ### we ran out of continent to extend the flow interval into
            new_tet = the_edge.lower_tet
        else:
            # print('sideways case')
            #### this was too complicated... let's be less efficient and just make it work
            # ### find which side of the edge our tet is in
            # side_tet_collections_at_edge = vt.side_tet_collections[the_edge.index] ## index in the manifold
            # side_face_collections_at_edge = vt.side_face_collections[the_edge.index] ## here these are ordered top to bottom
            # downward_path = None
            # flow_step = self.flow_cycle[next_down_index] ### pair of (tet index in manifold, index of edge in that tet)
            # for i, side_tet_collection in enumerate(side_tet_collections_at_edge):
            #     if flow_step in side_tet_collection:
            #         downward_path = side_tet_collection[:side_tet_collection.index(flow_step) + 1]
            #         downward_path_faces = side_face_collections_at_edge[i][:side_tet_collection.index(flow_step) + 1]
            # assert downward_path != None
            # new_tet = None
            # for j, (tet_num, edge_num) in enumerate(downward_path): 
            #     lower_boundary_triangles = [t for t in the_edge.boundary_triangles() if not t.is_upper and t.index == downward_path_faces[j][0]] 
            #     assert len(lower_boundary_triangles) <= 1
            #     for tet in the_edge.side_tetrahedra:  ### check to see if we already have the new_tet
            #         if tet.index == flow_step[0]:
            #             if tet.edge(flow_step[1]) == the_edge:
            #                 new_tet = tet
            #                 break
            #     if expand_continent:  
            #         if len(lower_boundary_triangles) == 1:
            #             self.continent.bury(lower_boundary_triangles[0])
            #     else:  ### did not find new_tet without expanding the continent
            #         return False
            if expand_continent:  ### this edge is on the bottom of some tet so it has either 2 or 0 lower boundary triangles incident to it, bury them until we get 0.
                while True:
                    lower_boundary_triangles = [t for t in the_edge.boundary_triangles() if not t.is_upper] 
                    if len(lower_boundary_triangles) == 0:
                        break ## out of while loop
                    self.continent.bury(lower_boundary_triangles[0])
            new_tet = None
            flow_step = self.flow_cycle[next_down_index] ### pair of (tet index in manifold, index of edge in that tet)
            for tet in the_edge.side_tetrahedra:
                if tet.index == flow_step[0]:
                    if tet.edge(flow_step[1]) == the_edge:
                        new_tet = tet
                        break
            if new_tet == None: ### we did not find the tet
                return False

        assert new_tet != None
        self.tetrahedra.insert(0, new_tet)
        self.edges.insert(0, the_edge) 
        self.down_index = next_down_index
        return True

    def extend_within_continent(self): 
        """extend this flow interval as much as possible without growing the continent"""
        while True:
            succeeded = self.go_up(expand_continent = False)
            if not succeeded:
                break ### out of while loop
        while True:
            succeeded = self.go_down(expand_continent = False)
            if not succeeded:
                break ### out of while loop

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

    def crossing_leaves(self): ### find green leaves of tetrahedra[-1] that cross purple leaves of tetrahedra[0]
        green_boundary = self.tetrahedra[-1].get_boundary_cusp_leaves()[0]
        purple_boundary = self.tetrahedra[0].get_boundary_cusp_leaves()[1]
        out_green = set([])
        out_purple = set([])
        for g in green_boundary:
            for p in purple_boundary:
                if g.links(p):
                    out_green.add(g)
                    out_purple.add(p)
        return (list(out_green), list(out_purple))

    def is_inside_edge_rectangle_green_sides(self, con_edge):
        con_edge.ensure_continent_contains_rectangle()
        greens, purples = con_edge.green_purple_rectangle_sides()

        for leaf in greens:
            while not leaf.is_entirely_to_one_side_of(self.tetrahedra[-1]):
                self.go_up()
            v_is_outside = [leaf.sees_to_its_left(v) == con_edge.is_red for v in self.tetrahedra[-1].vertices]
            assert all(v_is_outside) or not any(v_is_outside)
            if v_is_outside[0]:
                return False
        return True

    def is_inside_edge_rectangle_purple_sides(self, con_edge):
        con_edge.ensure_continent_contains_rectangle()
        greens, purples = con_edge.green_purple_rectangle_sides()                
        for leaf in purples:
            while not leaf.is_entirely_to_one_side_of(self.tetrahedra[0]):
                self.go_down()
            v_is_outside = [leaf.sees_to_its_left(v) != con_edge.is_red for v in self.tetrahedra[0].vertices]
            assert all(v_is_outside) or not any(v_is_outside)
            if v_is_outside[0]:
                return False
        return True

    def is_inside_edge_rectangle(self, con_edge):
        return self.is_inside_edge_rectangle_green_sides(con_edge) and self.is_inside_edge_rectangle_purple_sides(con_edge)

    # def find_edge_rectangle_we_are_inside(self, tet):
    #     """Given tet in the fund dom and in our flow cycle, find an edge of tet whose rectangle contains the flow interval's point"""  
    #     for e in tet.edges():
    #         if self.is_inside_edge_rectangle(e):
    #             return e
    #     assert False

    def is_green_comparable_with(self, other): ### check if upper tets do not overlap horizontally
        self.extend_within_continent()
        other.extend_within_continent()
        green_boundary = self.tetrahedra[-1].get_boundary_cusp_leaves()[0]
        if all([leaf.is_entirely_to_one_side_of(other.tetrahedra[-1]) for leaf in green_boundary]):
            green_boundary_other = other.tetrahedra[-1].get_boundary_cusp_leaves()[0]
            if all([leaf.is_entirely_to_one_side_of(self.tetrahedra[-1]) for leaf in green_boundary_other]): 
                return True ### both checks necessary. If one tet spans the other then each leaf of the outer tet sees the inner tet verts all to one side
        else:
            return False

    def is_purple_comparable_with(self, other): ### check if upper tets do not overlap vertically
        self.extend_within_continent()
        other.extend_within_continent()
        purple_boundary = self.tetrahedra[0].get_boundary_cusp_leaves()[1]
        if all([leaf.is_entirely_to_one_side_of(other.tetrahedra[0]) for leaf in purple_boundary]):
            purple_boundary_other = other.tetrahedra[0].get_boundary_cusp_leaves()[1]
            if all([leaf.is_entirely_to_one_side_of(self.tetrahedra[0]) for leaf in purple_boundary_other]):
                return True ### both checks are necessary - for an example see drilling_flow_cycle.drill_flow_cycle('gLLAQbecdfffhhnkqnc_120012', [(0, 4), (4, 5), (2, 4), (1, 2), (5, 1), (0, 4), (4, 0), (0, 4), (4, 5), (2, 4), (1, 2), (5, 1)]
        else:
            return False

    def make_green_comparable_with(self, other): ### make sure that upper tets do not overlap horizontally
        assert not self.equals(other)
        while not self.is_green_comparable_with(other):
            self.go_up()
            other.go_up()

    def make_purple_comparable_with(self, other): ### make sure that lower tets do not overlap vertically
        assert not self.equals(other)
        while not self.is_purple_comparable_with(other):
            self.go_down()
            other.go_down()

    def order_relative_to_W_E(self, other, W, E): ### are we closer to W or to E than other?
        """Assuming that W and E are cusps and self, other are between W and E, which is closer to which"""
        if self.equals(other):
            return 0  ## for use in sorting, we want answers of -1, 0, and 1.
        # print('making green comparable with')
        self.make_green_comparable_with(other)
        # print('done making green comparable with')
        ### find a leaf of self that separates W from E
        green_boundary = self.tetrahedra[-1].get_boundary_cusp_leaves()[0]
        test_leaf = None
        for leaf in green_boundary:
            if leaf.sees_to_its_left(W) != leaf.sees_to_its_left(E):
                test_leaf = leaf
                break
        assert test_leaf != None
        if test_leaf.separates([W], other.tetrahedra[-1].vertices):
            return -1  ## for use in sorting, we want answers of -1, 0, and 1.
        else:
            assert test_leaf.separates([E], other.tetrahedra[-1].vertices)
            return +1  ## for use in sorting, we want answers of -1, 0, and 1.

    def order_relative_to_S_N(self, other, S, N): ### are we closer to S or to N than other?
        """Assuming that S and N are cusps and self, other are between S and N, which is closer to which"""
        if self.equals(other):
            return 0  ## for use in sorting, we want answers of -1, 0, and 1.
        self.make_purple_comparable_with(other)
        ### find a leaf of self that separates S from N
        purple_boundary = self.tetrahedra[0].get_boundary_cusp_leaves()[1]
        test_leaf = None
        for leaf in purple_boundary:
            if leaf.sees_to_its_left(S) != leaf.sees_to_its_left(N):
                test_leaf = leaf
                break
        assert test_leaf != None
        if test_leaf.separates([S], other.tetrahedra[0].vertices):
            return -1  ## for use in sorting, we want answers of -1, 0, and 1.
        else:
            assert test_leaf.separates([N], other.tetrahedra[0].vertices)
            return +1  ## for use in sorting, we want answers of -1, 0, and 1.

def cmp_to_key(small_cusp, big_cusp, is_horizontal = True):  ### see https://stackoverflow.com/questions/32752739/how-does-the-functools-cmp-to-key-function-work
    '''Convert a cmp= function into a key= function''' ### small is either W or S, big is either E or N
    if is_horizontal:
        class K(object):
            def __init__(self, obj, *args):
                # print('obj created with ',obj)
                self.obj = obj
            def __lt__(self, other):
                # print('comparing less than ',self.obj)
                return self.obj.order_relative_to_W_E(other.obj, small_cusp, big_cusp) < 0
            def __gt__(self, other):
                # print('comparing greater than ',self.obj)
                return self.obj.order_relative_to_W_E(other.obj, small_cusp, big_cusp) > 0
            def __eq__(self, other):
                # print('comparing equal to ',self.obj)
                return self.obj.equals(other.obj)
            def __le__(self, other):
                 # print('comparing less than equal ',self.obj)
                return self.obj.order_relative_to_W_E(other.obj, small_cusp, big_cusp) <= 0
            def __ge__(self, other):
                # print('comparing greater than equal',self.obj)
               return self.obj.order_relative_to_W_E(other.obj, small_cusp, big_cusp) >= 0
            def __ne__(self, other):
                # print('comparing not equal ',self.obj)
                return not self.obj.equals(other.obj)
        return K
    else:
        class K(object):
            def __init__(self, obj, *args):
                # print('obj created with ',obj)
                self.obj = obj
            def __lt__(self, other):
                # print('comparing less than ',self.obj)
                return self.obj.order_relative_to_S_N(other.obj, small_cusp, big_cusp) < 0
            def __gt__(self, other):
                # print('comparing greater than ',self.obj)
                return self.obj.order_relative_to_S_N(other.obj, small_cusp, big_cusp) > 0
            def __eq__(self, other):
                # print('comparing equal to ',self.obj)
                return self.obj.equals(other.obj)
            def __le__(self, other):
                 # print('comparing less than equal ',self.obj)
                return self.obj.order_relative_to_S_N(other.obj, small_cusp, big_cusp) <= 0
            def __ge__(self, other):
                # print('comparing greater than equal',self.obj)
               return self.obj.order_relative_to_S_N(other.obj, small_cusp, big_cusp) >= 0
            def __ne__(self, other):
                # print('comparing not equal ',self.obj)
                return not self.obj.equals(other.obj)
        return K

def uniquify_list_of_flow_intervals(flow_intervals):
    unique_flow_intervals = []
    found_parallel = False
    for interval in flow_intervals:
        result, found_parallel_here = interval.is_in_list(unique_flow_intervals)
        if not found_parallel and found_parallel_here:
            found_parallel = True
        if not result:
            unique_flow_intervals.append(interval)
    return unique_flow_intervals, found_parallel

def translate_of_interval_from_one_edge_rect_to_another(e1, e2, interval):
    assert interval.is_inside_edge_rectangle(e1)
    assert e1.index == e2.index
    if e1.upper_tet == None:
        e1.ensure_continent_contains_tet_above() 
    if e2.upper_tet == None:
        e2.ensure_continent_contains_tet_above()
    t1 = e1.upper_tet
    t2 = e2.upper_tet
    assert t1.index == t2.index
    path = t1.face_num_path_to_other_tet(interval.continent.init_tet) + interval.continent.init_tet.face_num_path_to_other_tet(interval.init_tet)
    translated_interval_init_tet = t2.follow_face_num_path(path)
    new_interval = flow_interval(interval.continent, interval.flow_cycle, translated_interval_init_tet, 0)
    assert new_interval.is_inside_edge_rectangle(e2)
    return new_interval


