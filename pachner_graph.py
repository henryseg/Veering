
import regina
from taut import isosig_to_tri_angle, isosig_from_tri_angle
from pachner import twoThreeMove, threeTwoMove
from branched_surface import upper_branched_surface, isosig_from_tri_angle_branch, isosig_to_tri_angle_branch, has_non_sing_semiflow
from veering import is_veering
from flow_cycles import flow_cycle_to_triangle_loop, find_flow_cycles, tri_loop_is_boundary_parallel
from drill import drill

class pachner_node():
    def __init__(self, isoSig, tri, angle, branch, drilled_cusp_index = None, ceiling = 9999, floor = 0):
        self.isoSig = isoSig
        self.neighbour_moves_up_faces = {}
        self.neighbour_moves_down_edges = {}
        self.neighbour_moves_up_tri_angle_branch = {}
        self.neighbour_moves_down_tri_angle_branch = {}
        self.neighbour_moves_tri_angle_branch = None
        self.tri, self.angle, self.branch = tri, angle, branch
        # self.is_frontier = False #is a node on the boundary of what we have explored
        if drilled_cusp_index != None:
            self.drilled_cusp_index = drilled_cusp_index
        self.generate_neighbour_isoSigs(ceiling, floor)

    def all_neighbour_isoSigs(self):
        # return self.neighbour_moves_up_faces | self.neighbour_moves_down_edges  #union
        return set(self.neighbour_moves_up_faces.keys()) | set(self.neighbour_moves_down_edges.keys())

    def generate_neighbour_isoSigs(self, ceiling, floor):
        # tri, angle = self.getTriang()
        tri, angle, branch = self.tri, self.angle, self.branch
        num_tetrahedra = tri.countTetrahedra()
        assert num_tetrahedra <= ceiling
        assert num_tetrahedra >= floor

        if num_tetrahedra < ceiling:
            for face_index in range(tri.countTriangles()):
                ### only go this way if it expands the drilled cusp
                face = tri.triangles()[face_index]
                embed0 = face.embedding(0)
                tet0 = embed0.simplex()
                tet_0_face_num = embed0.face()
                embed1 = face.embedding(1)
                tet1 = embed1.simplex()
                tet_1_face_num = embed1.face()
                if (tet0.vertex(tet_0_face_num).index() == self.drilled_cusp_index) or (tet1.vertex(tet_1_face_num).index() == self.drilled_cusp_index):
                    tri_copy = regina.Triangulation3(tri)
                    angle_copy = angle[:]
                    branch_copy = branch[:]
                    output = twoThreeMove(tri_copy, face_index, angle = angle_copy, branch = branch_copy, return_vertex_perm = True)
                    # print('2-3', face_index, output)
                    if output != False:
                        tri_new, angle_new, branches_new, vertex_perm = output
                        drilled_cusp_index_new = vertex_perm[self.drilled_cusp_index]
                        for branch_new in branches_new:
                            sig_new = isosig_from_tri_angle_branch(tri_new, angle_new, branch_new)
                            self.neighbour_moves_up_faces[sig_new] = 'f' + str(face_index)
                            self.neighbour_moves_up_tri_angle_branch[sig_new] = (tri_new, angle_new, branch_new, drilled_cusp_index_new)

        if num_tetrahedra > floor:
            for edge_index in range(tri.countEdges()):
                ### only go this way if it expands the drilled cusp
                edge = tri.edges()[edge_index]
                if (edge.vertex(0).index() != self.drilled_cusp_index) and (edge.vertex(1).index() != self.drilled_cusp_index):
                    tri_copy = regina.Triangulation3(tri)
                    angle_copy = angle[:]
                    branch_copy = branch[:]
                    output = threeTwoMove(tri_copy, edge_index, angle = angle_copy, branch = branch_copy, return_vertex_perm = True)
                    if output != False:
                        tri_new, angle_new, branch_new, vertex_perm = output
                        drilled_cusp_index_new = vertex_perm[self.drilled_cusp_index]
                        sig_new = isosig_from_tri_angle_branch(tri_new, angle_new, branch_new)
                        self.neighbour_moves_down_edges[sig_new] = 'e' + str(edge_index)
                        self.neighbour_moves_down_tri_angle_branch[sig_new] = (tri_new, angle_new, branch_new, drilled_cusp_index_new)
        self.neighbour_moves_tri_angle_branch = {**self.neighbour_moves_up_tri_angle_branch, **self.neighbour_moves_down_tri_angle_branch} ## union of the dictionaries

def print_path(target_isoSig, big_dict_of_nodes):
    ### Given the target isoSig and the big dict of nodes, each node of which 
    ### knows where it came from in the search, build the path from target back
    ### to the start. Print out the list of nodes and save them as regina files.
    target_node = big_dict_of_nodes[target_isoSig]
    path = [target_node]
    current_node = target_node
    while current_node.came_from != None:
        current_node = big_dict_of_nodes[current_node.came_from]
        path.append(current_node)
    path.reverse()  ### now goes from start to target
    for i, node in enumerate(path):
        print(node.isoSig)
        # filename = node.isoSig + '|' + "".join([str(num) for num in node.angle]) + '.rga'
        ### First part is the tri-angle isoSig, second part is the angle structure as it comes
        ### in the sequence of triangulations without relabelling vertices as in the canonical
        ### isoSig triangulation.
        # print(filename)
        # node.tri.save('Output/' + filename)
        # if i < len(path) - 1:  ## not the last one
        #     next_node = path[i+1]
        #     if next_node.isoSig in node.neighbour_moves_up_faces:
        #         print(node.neighbour_moves_up_faces[next_node.isoSig])
        #     else:
        #         assert next_node.isoSig in node.neighbour_moves_down_edges
        #         print(node.neighbour_moves_down_edges[next_node.isoSig])

# def search_Pachner_graph_for_shortest_path(start_isoSig, name=None, search_depth = 3, ceiling = 5, check_property = False, property = None, save_dir = None):
#     tri, angle, branch = isosig_to_tri_angle_branch(start_isoSig) 
def search_Pachner_graph_for_shortest_path(start_isoSig, tri, angle, branch, name=None, search_depth = 3, ceiling = 5, drilled_cusp_index = None, check_property = False, property = None, save_dir = None):
    start_node = pachner_node( start_isoSig, tri, angle, branch, drilled_cusp_index = drilled_cusp_index, ceiling = ceiling )
    start_node.came_from = None

    big_dict_of_nodes = {start_isoSig : start_node}
    frontier_isoSigs = set([start_isoSig])
    # print(len(big_dict_of_nodes), len(frontier_isoSigs))
    for counter in range(search_depth):
        if len(frontier_isoSigs) == 0: #we are done...
            print('done')
            break
        new_frontier_isoSigs = set([]) 
        # for each element in the frontier check to see if it appears on the big_list_of_sigs if not we add it to the big list 
        for cur_isoSig in frontier_isoSigs:
            current_node = big_dict_of_nodes[cur_isoSig]
            neighbour_isoSigs = current_node.all_neighbour_isoSigs()
            for nb_isoSig in neighbour_isoSigs: 
                if not nb_isoSig in big_dict_of_nodes:
                    nb_tri, nb_angle, nb_branch, nb_drilled_cusp_index = current_node.neighbour_moves_tri_angle_branch[nb_isoSig]
                    #print('nb drilled cusp', nb_drilled_cusp_index)
                    new_node = pachner_node(nb_isoSig, nb_tri, nb_angle, nb_branch, drilled_cusp_index = nb_drilled_cusp_index, ceiling = ceiling)
                    new_node.came_from = cur_isoSig
                    if counter == search_depth - 1: #last layer
                        new_node.is_frontier = True

                    new_frontier_isoSigs.add(nb_isoSig)
                    big_dict_of_nodes[nb_isoSig] = new_node

                    if is_veering(nb_tri, nb_angle):
                        print('veering!', isosig_from_tri_angle_branch(nb_tri, nb_angle, nb_branch))
                        upper_branch = upper_branched_surface(nb_tri, nb_angle, return_lower = False)
                        lower_branch = upper_branched_surface(nb_tri, nb_angle, return_lower = True)
                        print('upper:', isosig_from_tri_angle_branch(nb_tri, nb_angle, upper_branch))
                        print('lower:', isosig_from_tri_angle_branch(nb_tri, nb_angle, lower_branch))
                        print_path(nb_isoSig, big_dict_of_nodes)
                        ## break
                        return None 

        frontier_isoSigs = new_frontier_isoSigs
        print(len(big_dict_of_nodes), len(frontier_isoSigs))

    print('did not find veering')
    return None

def main():
    depth = 100
    ceiling = 10

    print('depth', depth)
    print('ceiling', ceiling)

    # sig = 'cPcbbbdxm_10'
    # sig = 'dLQacccjsnk_200'
    sig = 'dLQbccchhfo_122'
    tri, angle = isosig_to_tri_angle(sig) 
    branch = upper_branched_surface(tri, angle)
    # tl = flow_cycle_to_triangle_loop(tri, branch, [(0, 2)]) 
    # drilled_cusp_index = drill(tri, tl, angle = angle, branch = branch) 
    # print('angle', angle, 'branch', branch, 'drilled_cusp_index', drilled_cusp_index)
    # # assert has_non_sing_semiflow(tri, branch)

    # # start veering: cPcbbbdxm_10_dl loop [(0, 2)] tri_loop [(2, 102)]
    # # drill: eLMkbbddddhapu_2100_fjek

    # start_isoSig = isosig_from_tri_angle_branch(tri, angle, branch)
    # assert start_isoSig == 'eLMkbbddddhapu_2100_fjek'  
    # tri, angle, branch = isosig_to_tri_angle_branch(start_isoSig)     ### fails?? 
    

    # graph = search_Pachner_graph_for_shortest_path(start_isoSig, tri, angle, branch,  name=None, search_depth = depth, ceiling = ceiling, drilled_cusp_index = drilled_cusp_index, check_property = False, property = None, save_dir = None)


    loops = find_flow_cycles(tri, branch)
    tri_loops = [flow_cycle_to_triangle_loop(tri, branch, loop) for loop in loops]
    
    for tri_loop in tri_loops:
        if tri_loop != False: # False means that tri_loop goes more than once  along the same triangle - not currently implemented
            tri, angle = isosig_to_tri_angle(sig)
            if tri_loop_is_boundary_parallel(tri_loop, tri) == False: # if a loop is boundary parallel then we don't drill
                branch = upper_branched_surface(tri, angle)
                drilled_cusp_index = drill(tri, tri_loop, angle = angle, branch = branch)
                start_isoSig = isosig_from_tri_angle_branch(tri, angle, branch)
                print(start_isoSig, 'angle', angle, 'branch', branch, 'drilled_cusp_index', drilled_cusp_index)
                graph = search_Pachner_graph_for_shortest_path(start_isoSig, tri, angle, branch,  name=None, search_depth = depth, ceiling = ceiling, drilled_cusp_index = drilled_cusp_index, check_property = False, property = None, save_dir = None)










