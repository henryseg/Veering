import regina
import pickle

from veering.file_io import output_to_pickle, read_from_pickle
from drill_veering_census_analyse import parse_line

class veering_node():
    def __init__(self, isoSig):
        self.isoSig = isoSig
        self.neighbour_isoSigs_drill = set([])
        # self.neighbour_isoSigs_up = set([])
        # self.is_frontier = False #is a node on the boundary of what we have explored
        # self.generate_neighbour_isoSigs(ceiling, floor)

        tri = self.getReginaTriang()
        self.num_tetrahedra = tri.countTetrahedra()
        self.num_cusps = tri.countVertices()

    # def get_neighbour_isoSigs(self):
    #     return self.neighbour_isoSigs_up | self.neighbour_isoSigs_down  #union

    def getReginaTriang(self):
        regina_isosig = self.isoSig.split('_')[0]
        return regina.Triangulation3(regina_isosig)


def main(max_count = 100):
    
    filename = '/Users/segerman/Library/CloudStorage/Dropbox/Data/drillings_census_4_ladders_max_cycle_len_5.txt'

    graph = dict()

    with open(filename, 'r') as file:
        count = 0
        while True:
            count += 1
            line = file.readline()
            if len(line) == 0 or count == max_count:
                break ### EOF

            drilled_sigs, undrilled_sig = parse_line(line)
            if undrilled_sig in graph:
                node = graph[undrilled_sig]
            else:
                node = veering_node(undrilled_sig)
                graph[undrilled_sig] = node
            for sig in drilled_sigs:
                if not sig in graph:
                    neighbour = veering_node(sig)
                    node.neighbour_isoSigs_drill.add(sig)
                    graph[sig] = neighbour

    output_to_pickle(graph, 'data/veering_drilling_graph.pkl')
    return graph

