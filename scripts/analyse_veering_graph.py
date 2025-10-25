import regina
import pickle
import veering
import drilling_flow_cycle

def read_from_pickle(filename):
    f = open(filename, 'rb')
    data = pickle.load(f)
    f.close()
    return data

def prepare_table(graph, max_tetrahedra = 5, max_cusps = 3):
    print('make table')
    table = []
    nodes_to_draw = []
    for i in range(max_cusps + 1):
        row = []
        for j in range(max_tetrahedra + 1):
            row.append([])
        table.append(row)

    for node in graph.values():
        if node.num_tetrahedra <= max_tetrahedra and node.num_cusps <= max_cusps:
            table[node.num_cusps][node.num_tetrahedra].append(node)
            nodes_to_draw.append(node)
    for i in range(max_cusps + 1):
        for j in range(max_tetrahedra + 1):
            table[i][j].sort(key = lambda n: n.isoSig)
    return table, nodes_to_draw

def main(max_tetrahedra = 16, max_cusps = 5):
    census = veering.veering_census()
    filename = 'data/veering_drilling_graph.pkl'
    graph = read_from_pickle(filename)
    
    
    # most_tetrahedra = 0
    # best_sig = None
    largest_drillings = []
    for i in range(17):
        largest_drillings.append([0, None])
    for node in graph.values():
        if node.num_tetrahedra <= 16:
            assert node.isoSig in census

            for flow_cycle in node.neighbour_isoSigs_drill_dict.keys():
                neighbour_sig = node.neighbour_isoSigs_drill_dict[flow_cycle]
                neighbour_num_tet = graph[neighbour_sig].num_tetrahedra
                if largest_drillings[node.num_tetrahedra][0] < neighbour_num_tet:
                    largest_drillings[node.num_tetrahedra] = [neighbour_num_tet, [(node.isoSig, flow_cycle)]]
                elif largest_drillings[node.num_tetrahedra][0] == neighbour_num_tet:
                    largest_drillings[node.num_tetrahedra][1].append((node.isoSig, flow_cycle))
    for line in largest_drillings:
        print(line)
        to_drill = line[1]
        if to_drill != None:
            for pair in to_drill:
                sig, fc = pair
                drilling_flow_cycle.drill_flow_cycle(sig, fc, draw_rectangles = True)

    #     if node.num_tetrahedra >= most_tetrahedra:
    #         most_tetrahedra = node.num_tetrahedra
    #         best_sig = node.isoSig
    # print(most_tetrahedra, best_sig)





