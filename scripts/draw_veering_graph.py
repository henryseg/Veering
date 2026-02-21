import regina
import pickle
from pyx import path, trafo, canvas, style, text, color, deco
from math import sqrt, ceil, floor, cos, sin, pi

def read_from_pickle(filename):
    f = open(filename, 'rb')
    data = pickle.load(f)
    f.close()
    return data

def make_line(a, b):
    p = path.path( path.moveto(a.real, a.imag) )
    p.append( path.lineto(b.real, b.imag) )
    return p

def make_curve(a, at, bt, b):
    a1 = a + at
    b1 = b - bt
    p = path.path( path.moveto(a.real, a.imag) )
    p.append( path.curveto(a1.real, a1.imag, b1.real, b1.imag, b.real, b.imag) )
    return p

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

def transp(target_dot_rad, max_dot_rad, max_transp = 0.9, min_transp = 0.6):
    param = target_dot_rad/max_dot_rad
    # param = param**1.2
    return (1.0 - param)*max_transp + param*min_transp

def draw(graph, name = 'foo', max_tetrahedra = 5, max_cusps = 3, max_transp = 0.6, min_transp = 0.5):

    canv = canvas.canvas()

    max_dot_rad = 0.02

    red = color.rgb(0.9, 0.1, 0)
    blue = color.rgb(0, 0.3, 0.9)
    green = color.rgb(0, 0.7, 0.0)
    purple = color.rgb(0.5, 0.0, 0.5)
    yellow = color.rgb(0.7, 0.7, 0.0)
    black = color.rgb(0.0, 0.0, 0.0)
    white = color.rgb(0.9, 0.9, 0.9)
    # default_col = white
    # background_col = black
    default_col = black
    # background_col = white

    global_scale_up = 5.0
    scl = trafo.trafo(matrix=((global_scale_up, 0), (0, global_scale_up)), vector=(0, 0))
    vertical_spacing_scale = 1.5

    r = 0.2
    theta = 0.2*pi
    tangent = r*complex(cos(theta), sin(theta)) # tangent for this kind of edge
    # phi = 0.4*pi
    # t14 = 2*r*complex(cos(phi), sin(phi))
    # psi = -0.1*pi
    # t02 = 2*r*complex(cos(psi), sin(psi))

    ### background rectangle
    # p = path.rect(1.5,0.5,max_tetrahedra + 1, max_vertices + 1)
    # p = p.transformed(scl)
    # canv.stroke(p, [deco.filled(), style.linewidth(0), background_col])

    text_pos = global_scale_up * complex(2.35, 0.75 - (1 - vertical_spacing_scale))
    canv.text(text_pos.real, text_pos.imag, "tetrahedra", textattrs=[text.size(3), default_col, text.halign.center, text.valign.middle])
    # canv.text(text_pos.real, text_pos.imag, "vertices", textattrs=[text.size(3), default_col, text.halign.center, text.valign.middle])

    text_pos = global_scale_up * complex(1.75, 1.35 - (1 - vertical_spacing_scale))
    canv.text(text_pos.real, text_pos.imag, "cusps", textattrs=[text.size(3), default_col, text.halign.center, text.valign.middle, trafo.rotate(90)])
    # canv.text(text_pos.real, text_pos.imag, "bubbles", textattrs=[text.size(3), default_col, text.halign.center, text.valign.middle, trafo.rotate(90)])


    for i in range(2, max_tetrahedra + 1):
        # print('i', i)
        text_pos = global_scale_up * complex(i + 0.35, 0.9 - (1 - vertical_spacing_scale))
        canv.text(text_pos.real, text_pos.imag, "$"+str(i)+"$", textattrs=[text.size(4), default_col, text.halign.center, text.valign.middle])
    for j in range(1, max_cusps + 1):
        # print('j', j)
        text_pos = global_scale_up * complex(1.9, vertical_spacing_scale*j + 0.35)
        canv.text(text_pos.real, text_pos.imag, "$"+str(j)+"$", textattrs=[text.size(4), default_col, text.halign.center, text.valign.middle])

    table, nodes_to_draw = prepare_table(graph, max_tetrahedra = max_tetrahedra, max_cusps = max_cusps)

    for row in table:
        print([len(entry) for entry in row])

    print('calculate coordinates and radii')
    for j, row in enumerate(table):
        for i, entry in enumerate(row):    
            num_dots = len(entry)
            side_length = ceil(sqrt(num_dots))
            # print(len(entry), side_length)

            square_offset = complex(i, vertical_spacing_scale*j)
            for k, node in enumerate(entry):
                x = floor(k/side_length) + 0.5
                y = k % side_length + 0.5
                node.dot_coords = complex(x, y)/(side_length * 1.5) + square_offset
                node.dot_rad = min(max_dot_rad, 0.15/side_length)

    print('drawing connections')
    for i, node in enumerate(nodes_to_draw):
        if i % 100 == 0:
            print(i/len(nodes_to_draw))
        for neighbour_sig in node.neighbour_isoSigs_drill_dict.values():
            if neighbour_sig in graph.keys():
                neighbour = graph[neighbour_sig]
                if neighbour in nodes_to_draw:
                    p = make_curve(node.dot_coords, tangent, tangent, neighbour.dot_coords)
                    p = p.transformed(scl)
                    canv.stroke(p, [style.linewidth(global_scale_up * 1.0 * neighbour.dot_rad), style.linecap.round, red, color.transparency(1.0 * transp(neighbour.dot_rad, max_dot_rad, max_transp = max_transp, min_transp = min_transp))])

    print('drawing dots')
    for node in nodes_to_draw:
        col = default_col
        p = path.circle(node.dot_coords.real, node.dot_coords.imag, node.dot_rad)
        p = p.transformed(scl)
        canv.stroke(p, [deco.filled(), style.linewidth(0), col])

    output_filename = 'Images/Graphs/' + name + '.pdf'
    print(output_filename)
    canv.writePDFfile(output_filename)

def draw_from_filename(filename, max_tetrahedra = 5, max_cusps = 2):
    name = filename.split('.')[0] + '_draw_max_' + str(max_tetrahedra) + '_' + str(max_cusps)
    name = name.split('/')[1]

    graph = read_from_pickle(filename)
    max_transp = 0.99999
    min_transp = 0.8
    draw(graph, name = name, max_tetrahedra = max_tetrahedra, max_cusps = max_cusps, max_transp = max_transp, min_transp = min_transp)


def main(max_tetrahedra = 16, max_cusps = 5):
    filename = 'data/veering_drilling_graph.pkl'
    draw_from_filename(filename, max_tetrahedra = max_tetrahedra, max_cusps = max_cusps)
    





