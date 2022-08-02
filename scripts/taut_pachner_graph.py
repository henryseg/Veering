
import regina

from veering.taut import isosig_to_tri_angle, isosig_from_tri_angle
from veering.taut_pachner import twoThreeMove, threeTwoMove
from veering.veering_tri import is_veering

class taut_pachner_node():
    def __init__(self, isoSig, tri = None, angle = None, ceiling = 9999, floor = 0):
        self.isoSig = isoSig
        self.neighbour_moves_up_faces = {}
        self.neighbour_moves_down_edges = {}
        self.neighbour_moves_up_tri_angles = {}
        self.neighbour_moves_down_tri_angles = {}
        self.neighbour_moves_tri_angles = None
        if tri == None:
            print('gen new triang')
            self.tri, self.angle = self.getTriang()
        else:
            self.tri, self.angle = tri, angle
        # self.is_frontier = False #is a node on the boundary of what we have explored
        self.generate_neighbour_isoSigs(ceiling, floor)

    def all_neighbour_isoSigs(self):
        # return self.neighbour_moves_up_faces | self.neighbour_moves_down_edges  #union
        return set(self.neighbour_moves_up_faces.keys()) | set(self.neighbour_moves_down_edges.keys())

    def getTriang(self):
        # return regina.Triangulation3(self.isoSig)
        return isosig_to_tri_angle(self.isoSig)

    def generate_neighbour_isoSigs(self, ceiling, floor):
        # tri, angle = self.getTriang()
        tri, angle = self.tri, self.angle
        num_tetrahedra = tri.countTetrahedra()
        assert num_tetrahedra <= ceiling
        assert num_tetrahedra >= floor

        if num_tetrahedra < ceiling:
            for face_index in range(tri.countTriangles()):
                tri_copy = regina.Triangulation3(tri)
                angle_copy = angle[:]
                output = twoThreeMove(tri_copy, angle_copy, face_index)
                if output != False:
                    tri_new, angle_new = output
                    self.neighbour_moves_up_faces[isosig_from_tri_angle(tri_new, angle_new)] = 'f' + str(face_index)
                    self.neighbour_moves_up_tri_angles[isosig_from_tri_angle(tri_new, angle_new)] = (tri_new, angle_new)

        if num_tetrahedra > floor:
            for edge_index in range(tri.countEdges()):
                tri_copy = regina.Triangulation3(tri)
                angle_copy = angle[:]
                output = threeTwoMove(tri_copy, angle_copy, edge_index)
                if output != False:
                    tri_new, angle_new = output
                    self.neighbour_moves_down_edges[isosig_from_tri_angle(tri_new, angle_new)] = 'e' + str(edge_index)
                    self.neighbour_moves_down_tri_angles[isosig_from_tri_angle(tri_new, angle_new)] = (tri_new, angle_new)
        self.neighbour_moves_tri_angles = {**self.neighbour_moves_up_tri_angles, **self.neighbour_moves_down_tri_angles} ## union of the dictionaries

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
        filename = node.isoSig + '|' + "".join([str(num) for num in node.angle]) + '.rga'
        ### First part is the tri-angle isoSig, second part is the angle structure as it comes
        ### in the sequence of triangulations without relabelling vertices as in the canonical
        ### isoSig triangulation.
        print(filename)
        node.tri.save('Output/' + filename)
        if i < len(path) - 1:  ## not the last one
            next_node = path[i+1]
            if next_node.isoSig in node.neighbour_moves_up_faces:
                print(node.neighbour_moves_up_faces[next_node.isoSig])
            else:
                assert next_node.isoSig in node.neighbour_moves_down_edges
                print(node.neighbour_moves_down_edges[next_node.isoSig])

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
            current_node = big_dict_of_nodes[cur_isoSig]
            neighbour_isoSigs = current_node.all_neighbour_isoSigs()
            for nb_isoSig in neighbour_isoSigs: 
                if not nb_isoSig in big_dict_of_nodes:
                    nb_tri, nb_angle = current_node.neighbour_moves_tri_angles[nb_isoSig]
                    new_node = taut_pachner_node(nb_isoSig, tri = nb_tri, angle = nb_angle, ceiling = ceiling)
                    new_node.came_from = cur_isoSig
                    if counter == search_depth - 1: #last layer
                        new_node.is_frontier = True

                    new_frontier_isoSigs.add(nb_isoSig)
                    big_dict_of_nodes[nb_isoSig] = new_node

                    if nb_isoSig == target_isoSig:
                        print_path(nb_isoSig, big_dict_of_nodes)
                        break 

        frontier_isoSigs = new_frontier_isoSigs
        print(len(big_dict_of_nodes), len(frontier_isoSigs))

    return None

def find_veering_from_drilled(start_isoSig, name=None, search_depth = 100, ceiling = 8, check_property = False, property = None, save_dir = None):
    #first check if the first one is not veering
    tri, angle = isosig_to_tri_angle(start_isoSig)
    if is_veering(tri, angle):
        print('start_isoSig is veering')
    else:
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
                current_node = big_dict_of_nodes[cur_isoSig]
                neighbour_isoSigs = current_node.all_neighbour_isoSigs()
                for nb_isoSig in neighbour_isoSigs: 
                    if not nb_isoSig in big_dict_of_nodes:
                        nb_tri, nb_angle = current_node.neighbour_moves_tri_angles[nb_isoSig]
                        new_node = taut_pachner_node(nb_isoSig, tri = nb_tri, angle = nb_angle, ceiling = ceiling)
                        new_node.came_from = cur_isoSig
                        if counter == search_depth - 1: #last layer
                            new_node.is_frontier = True

                        new_frontier_isoSigs.add(nb_isoSig)
                        big_dict_of_nodes[nb_isoSig] = new_node

                        if is_veering(nb_tri, nb_angle):
                            print_path(nb_isoSig, big_dict_of_nodes)
                            print('veering:', isosig_from_tri_angle(nb_tri, nb_angle))
                            break 

            frontier_isoSigs = new_frontier_isoSigs
            print(len(big_dict_of_nodes), len(frontier_isoSigs))

        return None
