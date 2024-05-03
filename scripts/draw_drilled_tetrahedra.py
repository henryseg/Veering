
import pyx ### vector graphics 
from pyx import path, trafo, canvas, style, text, color, deco

def make_line(a, b):
    p = path.path( path.moveto(a.real, a.imag) )
    p.append( path.lineto(b.real, b.imag) )
    return p

def draw_square(i, canv, global_scale_up, scl, edge_thickness, green, purple):
    offset = complex(4*i,0)
    for j in range(4):
        p = make_line(complex(0,j) + offset, complex(3,j) + offset)
        p = p.transformed(scl)
        canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, purple])
    for j in range(4):
        p = make_line(complex(j,0) + offset, complex(j,3) + offset)
        p = p.transformed(scl)
        canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, green])
    text_pos = global_scale_up * (complex(1.5, 3.5) + offset)
    canv.text(text_pos.real, text_pos.imag, "$"+str(i)+"$", textattrs=[text.size(3), text.halign.center, text.valign.middle])

def draw_drilled_tetrahedra(con, name = "", tetrahedra_cusp_orders = None, intervals_inside_tet_rectangles = None, tetrahedra_chunks = None,
    edge_thickness = 0.02,
    leaf_thickness = 0.03,
    transparency = 0.9):

    global_scale_up = 1.0
    red = color.rgb(0.9,0.3,0)
    blue = color.rgb(0,0.3,0.9)
    edge_colours = {True: red, False: blue}
    green = color.rgb(0.0,0.5,0.0)
    purple = color.rgb(0.5,0.0,0.5)
    leaf_colours = {True: green, False: purple}
    
    scl = trafo.trafo(matrix=((global_scale_up, 0), (0, global_scale_up)), vector=(0, 0))
    canv = canvas.canvas()
    # canv.stroke(path.circle(0,0,global_scale_up), [style.linewidth(edge_thickness)])
    
    for i in range(len(tetrahedra_cusp_orders)):
        draw_square(i, canv, global_scale_up, scl, edge_thickness, green, purple)

        horizontal_cusp_order, vertical_cusp_order = tetrahedra_cusp_orders[i]
        offset = complex(4*i,0)
        W_coords = complex(0, vertical_cusp_order.index(horizontal_cusp_order[0])) 
        E_coords = complex(3, vertical_cusp_order.index(horizontal_cusp_order[3])) 
        S_coords = complex(horizontal_cusp_order.index(vertical_cusp_order[0]), 0) 
        N_coords = complex(horizontal_cusp_order.index(vertical_cusp_order[3]), 3)
        top_vertices, bottom_vertices = con.vt.tet_vert_posns[i] ### this is currently set in build_continent 
        N_ind, S_ind = top_vertices
        W_ind, E_ind = bottom_vertices

        for coords, ind in [(W_coords, W_ind), (E_coords, E_ind), (S_coords, S_ind), (N_coords, N_ind)]:
            
            dot_coords = coords + offset
            symbol_coords = 1.13*(coords - complex(1.5, 1.5)) + complex(1.5, 1.5) + offset
            canv.stroke(path.circle(dot_coords.real, dot_coords.imag, 0.05), [deco.filled(), style.linewidth(0)])
            canv.text(symbol_coords.real, symbol_coords.imag, "$"+str(ind)+"$", textattrs=[text.size(-2), text.halign.center, text.valign.middle])

        we_chunks, sn_chunks = tetrahedra_chunks[i]
        widths = [len(chunk) + 1 for chunk in we_chunks]
        heights = [len(chunk) + 1 for chunk in sn_chunks]
        # print('tet', i, 'widths', widths, 'heights', heights)
        tet_intervals = intervals_inside_tet_rectangles[i]
        for interval in tet_intervals:
            # print('we chunk num', interval.we_chunk_num, 'we_in_chunk_index', interval.we_in_chunk_index, 'sn chunk num', interval.sn_chunk_num, 'sn_in_chunk_index', interval.sn_in_chunk_index, 'interval', interval)
            x = interval.we_chunk_num + (interval.we_in_chunk_index + 1)/widths[interval.we_chunk_num]
            y = interval.sn_chunk_num + (interval.sn_in_chunk_index + 1)/heights[interval.sn_chunk_num]
            coords = offset + complex(x, y)
            canv.stroke(path.circle(coords.real, coords.imag, 0.035), [deco.filled(), style.linewidth(0)])



    output_filename = 'Images/DrilledTetrahedra/' + name + '.pdf'
    print(output_filename)
    canv.writePDFfile(output_filename)
        