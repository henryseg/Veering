
import regina
from taut import isosig_from_tri_angle, isosig_to_tri_angle
from branched_surface import isosig_to_tri_angle_branch, isosig_from_tri_angle_branch, upper_branched_surface
from taut_polytope import is_layered
import taut_pachner
import branched_pachner
import string

class taut_branched_pachner_node():
    def __init__(self, isoSig, tri = None, angle = None, branch = None, ceiling = 9999, floor = 0):
        self.isoSig = isoSig
        self.neighbour_moves_up_faces = {}
        self.neighbour_moves_down_edges = {}
        self.neighbour_moves_up_tri_angle_branches = {}
        self.neighbour_moves_down_tri_angle_branches = {}
        self.neighbour_moves_tri_angle_branches = None
        if tri == None:
            print('gen new triang')
            self.tri, self.angle, self.branch = self.getTriang()
        else:
            self.tri, self.angle, self.branch = tri, angle, branch
        # self.is_frontier = False #is a node on the boundary of what we have explored
        self.generate_neighbour_isoSigs(ceiling, floor)

    def all_neighbour_isoSigs(self):
        # return self.neighbour_moves_up_faces | self.neighbour_moves_down_edges  #union
        return set(self.neighbour_moves_up_faces.keys()) | set(self.neighbour_moves_down_edges.keys())

    def getTriang(self):
        return isosig_to_tri_angle_branch(self.isoSig)  ### broken

    def generate_neighbour_isoSigs(self, ceiling, floor):
        # tri, angle = self.getTriang()
        tri, angle, branch = self.tri, self.angle, self.branch
        num_tetrahedra = tri.countTetrahedra()
        assert num_tetrahedra <= ceiling
        assert num_tetrahedra >= floor

        if num_tetrahedra < ceiling:
            for face_index in range(tri.countTriangles()):
                # print('face_index', face_index)
                tri_copy = regina.Triangulation3(tri)
                angle_copy = angle[:]
                output_taut = taut_pachner.twoThreeMove(tri_copy, angle_copy, face_index)
                # print('output_taut', output_taut)
                tri_copy2 = regina.Triangulation3(tri)
                branch_copy = branch[:]
                output_branch = branched_pachner.twoThreeMove(tri_copy2, branch_copy, face_index)
                # print('output_branch', output_branch)
                if output_taut != False and output_branch != False:
                    tri_new, angle_new = output_taut
                    _, branch_new_list = output_branch
                    for branch_new in branch_new_list:
                        new_isosig = isosig_from_tri_angle_branch(tri_new, angle_new, branch_new)
                        # print('new_isosig', new_isosig)
                        self.neighbour_moves_up_faces[new_isosig] = 'f' + str(face_index)
                        self.neighbour_moves_up_tri_angle_branches[new_isosig] = (tri_new, angle_new, branch_new)

        if num_tetrahedra > floor:
            # print('now do threeTwoMove')
            for edge_index in range(tri.countEdges()):
                # print('edge_index', edge_index)
                tri_copy = regina.Triangulation3(tri)
                angle_copy = angle[:]
                output_taut = taut_pachner.threeTwoMove(tri_copy, angle_copy, edge_index)
                # print('output_taut', output_taut)
                tri_copy2 = regina.Triangulation3(tri)
                branch_copy = branch[:]
                output_branch = branched_pachner.threeTwoMove(tri_copy2, branch_copy, edge_index)
                # print('output_branch', output_branch)
                if output_taut != False and output_branch != False:
                    tri_new, angle_new = output_taut
                    _, branch_new = output_branch
                    new_isosig = isosig_from_tri_angle_branch(tri_new, angle_new, branch_new)
                    # print('new_isosig', new_isosig)
                    self.neighbour_moves_down_edges[new_isosig] = 'e' + str(edge_index)
                    self.neighbour_moves_down_tri_angle_branches[new_isosig] = (tri_new, angle_new, branch_new)
        self.neighbour_moves_tri_angle_branches = {**self.neighbour_moves_up_tri_angle_branches, **self.neighbour_moves_down_tri_angle_branches} ## union of the dictionaries

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
        filename = node.isoSig + '|' + "".join([str(num) for num in node.angle]) + '|' + "".join([string.ascii_lowercase[b] for b in node.branch]) + '.rga'
        ### First part is the tri-angle-branch isoSig, second part is the angle structure as it comes
        ### in the sequence of triangulations without relabelling vertices as in the canonical
        ### isoSig triangulation.
        ### Third part is the branched surface 
        print(filename)
        # node.tri.save('Output/' + filename)
        print(is_layered(node.tri, node.angle))
        if i < len(path) - 1:  ## not the last one
            next_node = path[i+1]
            if next_node.isoSig in node.neighbour_moves_up_faces:
                print(node.neighbour_moves_up_faces[next_node.isoSig])
            else:
                assert next_node.isoSig in node.neighbour_moves_down_edges
                print(node.neighbour_moves_down_edges[next_node.isoSig])

def search_Pachner_graph_for_shortest_path(start_isoSig, start_tri, start_angle, start_branch, target_isoSig, name=None, search_depth = 3, ceiling = 5, check_property = False, property = None, save_dir = None):
    start_node = taut_branched_pachner_node( start_isoSig, tri = start_tri, angle = start_angle, branch = start_branch , ceiling = ceiling )
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
                    nb_tri, nb_angle, nb_branch = current_node.neighbour_moves_tri_angle_branches[nb_isoSig]
                    new_node = taut_branched_pachner_node(nb_isoSig, tri = nb_tri, angle = nb_angle, branch = nb_branch, ceiling = ceiling)
                    new_node.came_from = cur_isoSig
                    if counter == search_depth - 1: #last layer
                        new_node.is_frontier = True

                    new_frontier_isoSigs.add(nb_isoSig)
                    big_dict_of_nodes[nb_isoSig] = new_node

                    if nb_isoSig.split('_')[:2] == target_isoSig.split('_'):
                        print_path(nb_isoSig, big_dict_of_nodes)
                        return 'done'

        frontier_isoSigs = new_frontier_isoSigs
        print(len(big_dict_of_nodes), len(frontier_isoSigs))

    return None
