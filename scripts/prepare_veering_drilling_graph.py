import regina
import pickle

from veering.file_io import output_to_pickle, read_from_pickle
from drill_veering_census_analyse import parse_line

class veering_node():
    def __init__(self, isoSig):
        self.isoSig = isoSig
        self.neighbour_isoSigs_drill_dict = {} ### keys are the flow_cycle we drill along
        self.neighbour_isoSigs_fill_set = set([]) ### pairs of (undrilled_sig, flow_cycle)

        tri = self.getReginaTriang()
        self.num_tetrahedra = tri.countTetrahedra()
        self.num_cusps = tri.countVertices()

    # def get_neighbour_isoSigs(self):
    #     return self.neighbour_isoSigs_up | self.neighbour_isoSigs_down  #union

    def getReginaTriang(self):
        regina_isosig = self.isoSig.split('_')[0]
        return regina.Triangulation3(regina_isosig)

def make_graph(max_count = 100):
    
    filename = '/Users/segerman/Library/CloudStorage/Dropbox/Data/drillings_census_4_ladders_max_cycle_len_5.txt'

    graph = dict()

    with open(filename, 'r') as file:
        count = 0
        while True:
            count += 1
            if count % 500 == 0:
                print(count)
            line = file.readline()
            if len(line) == 0 or count == max_count:
                break ### EOF

            drill_data_str, drill_data, undrilled_sig = parse_line(line, return_set = False)
            if undrilled_sig in graph:
                node = graph[undrilled_sig]
            else:
                node = veering_node(undrilled_sig)
                graph[undrilled_sig] = node
            for pair in drill_data:
                drilled_sig, flow_cycle = pair
                if not drilled_sig in graph:
                    neighbour = veering_node(drilled_sig)
                    graph[drilled_sig] = neighbour
                else:
                    neighbour = graph[drilled_sig]
                node.neighbour_isoSigs_drill_dict[flow_cycle] = drilled_sig
                neighbour.neighbour_isoSigs_fill_set.add((undrilled_sig, flow_cycle))

    output_to_pickle(graph, 'data/veering_drilling_graph.pkl')
    return graph
