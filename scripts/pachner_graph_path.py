### find shortest path in the pachner graph from a start to a target triangulation
### No taut or branch structures, so we use plain Regina isosigs

import regina
from veering.pachner import twoThreeMove, threeTwoMove

class pachner_node():
    def __init__(self, isoSig, tri, ceiling = 9999, floor = 0):
        self.isoSig = isoSig
        self.neighbour_moves_up_faces = {}
        self.neighbour_moves_down_edges = {}
        self.neighbour_moves_up_tri = {}
        self.neighbour_moves_down_tri = {}
        self.neighbour_moves_tri = None
        self.tri = tri 
        # self.is_frontier = False #is a node on the boundary of what we have explored
        self.generate_neighbour_isoSigs(ceiling, floor)

    def all_neighbour_isoSigs(self):
        return set(self.neighbour_moves_up_faces.keys()) | set(self.neighbour_moves_down_edges.keys())

    def generate_neighbour_isoSigs(self, ceiling, floor):
        tri = self.tri
        num_tetrahedra = tri.countTetrahedra()
        assert num_tetrahedra <= ceiling
        assert num_tetrahedra >= floor

        if num_tetrahedra < ceiling:
            for face_index in range(tri.countTriangles()):
                face = tri.triangles()[face_index]
                embed0 = face.embedding(0)
                tet0 = embed0.simplex()
                tet_0_face_num = embed0.face()
                embed1 = face.embedding(1)
                tet1 = embed1.simplex()
                tet_1_face_num = embed1.face()

                tri_copy = regina.Triangulation3(tri) 
                output = twoThreeMove(tri_copy, face_index)
   
                if output != False:
                    tri_new = output[0]
                    sig_new = tri_new.isoSig()  
                    self.neighbour_moves_up_faces[sig_new] = 'f' + str(face_index)
                    self.neighbour_moves_up_tri[sig_new] = tri_new

        if num_tetrahedra > floor:
            for edge_index in range(tri.countEdges()):
                edge = tri.edges()[edge_index]
                tri_copy = regina.Triangulation3(tri)
                output = threeTwoMove(tri_copy, edge_index)
                if output != False:
                    tri_new = output[0]
                    sig_new = tri_new.isoSig()  
                    self.neighbour_moves_down_edges[sig_new] = 'e' + str(edge_index)
                    self.neighbour_moves_down_tri[sig_new] = tri_new
        self.neighbour_moves_tri = {**self.neighbour_moves_up_tri, **self.neighbour_moves_down_tri} ## union of the dictionaries

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
        # filename = node.isoSig + '.rga'
        # print(filename)
        # node.tri.save('Output/' + filename)
        if i < len(path) - 1:  ## not the last one
            next_node = path[i+1]
            if next_node.isoSig in node.neighbour_moves_up_faces:
                print(node.neighbour_moves_up_faces[next_node.isoSig])
            else:
                assert next_node.isoSig in node.neighbour_moves_down_edges
                print(node.neighbour_moves_down_edges[next_node.isoSig])

def search_Pachner_graph_for_shortest_path(start_isoSig, target_isoSig, name = None, search_depth = 3, ceiling = 5, save_dir = None):
    start_tri = regina.Triangulation3.fromIsoSig(start_isoSig)
    start_tri.orient()
    start_node = pachner_node( start_isoSig, start_tri, ceiling = ceiling )
    start_node.came_from = None

    big_dict_of_nodes = {start_isoSig : start_node}
    frontier_isoSigs = set([start_isoSig])
    # print(len(big_dict_of_nodes), len(frontier_isoSigs))
    for counter in range(search_depth):
        if len(frontier_isoSigs) == 0: #we are done...
            # print('done')
            break
        new_frontier_isoSigs = set([]) 
        # for each element in the frontier check to see if it appears on the big_list_of_sigs if not we add it to the big list 
        for cur_isoSig in frontier_isoSigs:
            current_node = big_dict_of_nodes[cur_isoSig]
            neighbour_isoSigs = current_node.all_neighbour_isoSigs()
            for nb_isoSig in neighbour_isoSigs: 
                if not nb_isoSig in big_dict_of_nodes:
                    nb_tri = current_node.neighbour_moves_tri[nb_isoSig]
                    new_node = pachner_node(nb_isoSig, nb_tri, ceiling = ceiling)
                    new_node.came_from = cur_isoSig
                    if counter == search_depth - 1: #last layer
                        new_node.is_frontier = True

                    new_frontier_isoSigs.add(nb_isoSig)
                    big_dict_of_nodes[nb_isoSig] = new_node

                    if nb_isoSig == target_isoSig:
                        print('found path')
                        print_path(nb_isoSig, big_dict_of_nodes)
                        ## break
                        return None 

        frontier_isoSigs = new_frontier_isoSigs
        print(len(big_dict_of_nodes), len(frontier_isoSigs))

    print('did not find path')
    return None