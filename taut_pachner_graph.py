
import regina
from taut import isosig_to_tri_angle, isosig_from_tri_angle
from taut_pachner import twoThreeMove, threeTwoMove

class taut_pachner_node():
    def __init__(self, isoSig, ceiling = 9999, floor = 0):
        self.isoSig = isoSig
        self.neighbour_isoSigs_down = {}
        self.neighbour_isoSigs_up = {}
        # self.is_frontier = False #is a node on the boundary of what we have explored
        self.generate_neighbour_isoSigs(ceiling, floor)

    def get_neighbour_isoSigs(self):
        # return self.neighbour_isoSigs_up | self.neighbour_isoSigs_down  #union
        return set(self.neighbour_isoSigs_up.keys()) | set(self.neighbour_isoSigs_down.keys())

    def getTriang(self):
        # return regina.NTriangulation(self.isoSig)
        return isosig_to_tri_angle(self.isoSig)

    def generate_neighbour_isoSigs(self, ceiling, floor):
        tri, angle = self.getTriang()
        num_tetrahedra = tri.countTetrahedra()
        assert num_tetrahedra <= ceiling
        assert num_tetrahedra >= floor

        if num_tetrahedra < ceiling:
            for face_index in range(tri.countTriangles()):
                tri_copy = regina.NTriangulation(tri)
                angle_copy = angle[:]
                output = twoThreeMove(tri_copy, angle_copy, face_index)
                if output != False:
                    tri_new, angle_new = output
                    self.neighbour_isoSigs_up[isosig_from_tri_angle(tri_new, angle_new)] = 'f' + str(face_index)
        # else:
        #     self.is_frontier = True  #frontier either because we are at ceiling # tet, or because this is last layer of search

        if num_tetrahedra > floor:
            for edge_index in range(tri.countEdges()):
                tri_copy = regina.NTriangulation(tri)
                angle_copy = angle[:]
                output = threeTwoMove(tri_copy, angle_copy, edge_index)
                if output != False:
                    tri_new, angle_new = output
                    self.neighbour_isoSigs_down[isosig_from_tri_angle(tri_new, angle_new)] = 'e' + str(edge_index)

def print_path_home(target_isoSig, big_dict_of_nodes):
    current_node = big_dict_of_nodes[target_isoSig]
    while True:
        print(current_node.isoSig)
        if current_node.came_from == None:
            break
        else:
            if current_node.came_from in current_node.neighbour_isoSigs_down:
                print(current_node.neighbour_isoSigs_down[current_node.came_from])
            else:
                print(current_node.neighbour_isoSigs_up[current_node.came_from])
            current_node = big_dict_of_nodes[current_node.came_from]

def search_Pachner_graph_for_shortest_path(start_isoSig, target_isoSig, name=None, search_depth = 3, ceiling = 5, check_property = False, property = None, save_dir = None):
    start_node = taut_pachner_node( start_isoSig, ceiling = ceiling )
    start_node.came_from = None

    big_dict_of_nodes = {start_isoSig : start_node}
    frontier_isoSigs = set([start_isoSig])
    print(len(big_dict_of_nodes), len(frontier_isoSigs))
    for counter in range(search_depth):
        if len(frontier_isoSigs) == 0: #we are done...
            break
        new_frontier_isoSigs = set([]) 
        # for each element in the frontier check to see if it appears on the big_list_of_sigs if not we add it to the big list 
        for cur_isoSig in frontier_isoSigs:
            neighbour_isoSigs = big_dict_of_nodes[cur_isoSig].get_neighbour_isoSigs()
            for nb_isoSig in neighbour_isoSigs: 

                # if nb_isoSig == target_isoSig:
                #         print('found', target_isoSig)
                #         print_path_home(cur_isoSig, big_dict_of_nodes)
                #         print('---')

                if not nb_isoSig in big_dict_of_nodes:
                    new_node = taut_pachner_node(nb_isoSig, ceiling = ceiling)
                    new_node.came_from = cur_isoSig
                    # if check_property:
                    #     if property(start_node, name, save_dir):
                    #         print start_node.isoSig
                    if counter == search_depth - 1: #last layer
                        new_node.is_frontier = True
                    # update_neighbours(new_node, frontier_nodes) 
                    # # neighbours of new_node cannot be in new_frontier_nodes by parity of num tet, must be in frontier_nodes
                    # # cannot be in rest of nodes, since their neighbours have all been found 
                    new_frontier_isoSigs.add(nb_isoSig)
                    big_dict_of_nodes[nb_isoSig] = new_node

                    if nb_isoSig.split('_')[0] == target_isoSig:
                        print_path_home(nb_isoSig, big_dict_of_nodes)
                        break

                    

        frontier_isoSigs = new_frontier_isoSigs
        print(len(big_dict_of_nodes), len(frontier_isoSigs))

    return None

def main():
    depth = 1000
    ceiling = 8

    print('depth', depth)
    print('ceiling', ceiling)

    target_isoSig = 'gLLPQceeffefiiaellu'  ### drilled
    start_isoSig = 'gLLPQccdfeffhggaagb_201022'  ### veering
    graph = search_Pachner_graph_for_shortest_path(start_isoSig, target_isoSig, name=None, search_depth = depth, ceiling = ceiling, check_property = False, property = None, save_dir = None)




