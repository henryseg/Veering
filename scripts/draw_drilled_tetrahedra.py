
import pyx ### vector graphics 
from pyx import path, trafo, canvas, style, text, color, deco

def make_line(a, b):
    p = path.path( path.moveto(a.real, a.imag) )
    p.append( path.lineto(b.real, b.imag) )
    return p

def draw_square(i, canv, global_scale_up, scl, edge_thickness, green, purple, draw_inner_lines = True):
    offset = complex(4*i,0)
    if draw_inner_lines:
        indices = [0,1,2,3]
    else:
        indices = [0,3]
    for j in indices:
        p = make_line(complex(0,j) + offset, complex(3,j) + offset)
        p = p.transformed(scl)
        canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, purple, color.transparency(0.75)])
    for j in indices:
        p = make_line(complex(j,0) + offset, complex(j,3) + offset)
        p = p.transformed(scl)
        canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, green, color.transparency(0.75)])
    text_pos = global_scale_up * (complex(1.5, 3.8) + offset)
    canv.text(text_pos.real, text_pos.imag, "$"+str(i)+"$", textattrs=[text.size(3), text.halign.center, text.valign.middle])

def draw_edge_rectangle(canv, corners, scl, edge_thickness, green, purple, rect_scale_in = 0.9):
    center = 0.5*(corners[0] + corners[1])
    corners = [rect_scale_in*(corner - center) + center for corner in corners]
    x0, x1 = sorted([corner.real for corner in corners])
    y0, y1 = sorted([corner.imag for corner in corners])
    p = make_line(complex(x0, y0), complex(x1, y0))
    p = p.transformed(scl)
    canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, purple])
    p = make_line(complex(x1, y0), complex(x1, y1))
    p = p.transformed(scl)
    canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, green])
    p = make_line(complex(x1, y1), complex(x0, y1))
    p = p.transformed(scl)
    canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, purple])
    p = make_line(complex(x0, y1), complex(x0, y0))
    p = p.transformed(scl)
    canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, green])

def draw_rectangle(canv, x_interval, y_interval, scl, edge_thickness, green, purple, rect_scale_in = 0.9):
    corners = [complex(x_interval[i], y_interval[i]) for i in range(2)]
    draw_edge_rectangle(canv, corners, scl, edge_thickness, green, purple, rect_scale_in = rect_scale_in)

def draw_drilled_tetrahedra(con, name = "", tetrahedra_cusp_orders = None, 
    intervals_inside_tet_rectangles = None, 
    tetrahedra_chunks = None,
    old_tet_rectangles = [],
    draw_edge_rectangles = False,
    draw_face_rectangles = False, 
    draw_tetrahedron_rectangles = True,
    draw_square_inner_lines = True,
    draw_vertex_numbers = False, 
    edge_thickness = 0.02,
    transparency = 0.9):

    global_scale_up = 1.0
    red = color.rgb(0.9,0.3,0)
    blue = color.rgb(0,0.3,0.9)
    edge_colours = {True: red, False: blue}
    green = color.rgb(0.0,0.5,0.0)
    purple = color.rgb(0.5,0.0,0.5)
    leaf_colours = {True: green, False: purple}
    grey = color.rgb(0.7,0.7,0.7)
    
    scl = trafo.trafo(matrix=((global_scale_up, 0), (0, global_scale_up)), vector=(0, 0))
    canv = canvas.canvas()
    
    for i in range(len(tetrahedra_cusp_orders)):
        draw_square(i, canv, global_scale_up, scl, edge_thickness, green, purple, draw_inner_lines = draw_square_inner_lines)
        old_tet_rect = old_tet_rectangles[i]
        tet_rect_coords = [None] * len(old_tet_rect.horiz_ordering)

        horizontal_cusp_order, vertical_cusp_order = tetrahedra_cusp_orders[i]
        offset = complex(4*i,0)
        W_coords = complex(0, vertical_cusp_order.index(horizontal_cusp_order[0])) 
        W_mid_coords = complex(0, 1.5) 
        E_coords = complex(3, vertical_cusp_order.index(horizontal_cusp_order[3])) 
        E_mid_coords = complex(3, 1.5) 
        S_coords = complex(horizontal_cusp_order.index(vertical_cusp_order[0]), 0) 
        S_mid_coords = complex(1.5, 0) 
        N_coords = complex(horizontal_cusp_order.index(vertical_cusp_order[3]), 3)
        N_mid_coords = complex(1.5, 3)
        top_vertices, bottom_vertices = con.vt.tet_vert_posns[i] ### this is currently set in build_continent 
        N_ind, S_ind = top_vertices
        W_ind, E_ind = bottom_vertices

        tet_rect_coords[0] = global_scale_up*(W_coords + offset)
        tet_rect_coords[-1] = global_scale_up*(E_coords + offset)
        v2hperm = old_tet_rect.vert_to_horiz_perm()
        tet_rect_coords[v2hperm[0]] = global_scale_up*(S_coords + offset)
        tet_rect_coords[v2hperm[-1]] = global_scale_up*(N_coords + offset)

        ### tetrahedron vertex indices

        for coords, ind in [(W_coords, W_ind), (E_coords, E_ind), (S_coords, S_ind), (N_coords, N_ind)]:
            dot_coords = coords + offset
            canv.stroke(path.circle(dot_coords.real, dot_coords.imag, 0.05), [deco.filled(), style.linewidth(0)])
            if draw_vertex_numbers:
                symbol_coords = global_scale_up*(1.13*(coords - complex(1.5, 1.5)) + complex(1.5, 1.5) + offset)
                canv.text(symbol_coords.real, symbol_coords.imag, "$"+str(ind)+"$", textattrs=[text.size(-4), grey, text.halign.center, text.valign.middle])

        Regina_tet = con.vt.tri.tetrahedron(i)
        for coords, ind in [(W_mid_coords, E_ind), (E_mid_coords, W_ind), (S_mid_coords, N_ind), (N_mid_coords, S_ind)]:
            symbol_coords = global_scale_up*(1.13*(coords - complex(1.5, 1.5)) + complex(1.5, 1.5) + offset)
            tri_num = Regina_tet.triangle(ind).index()
            canv.text(symbol_coords.real, symbol_coords.imag, "$"+str(tri_num)+"$", textattrs=[text.size(-2), text.halign.center, text.valign.middle])


        we_chunks, sn_chunks = tetrahedra_chunks[i]
        widths = [len(chunk) + 1 for chunk in we_chunks]
        heights = [len(chunk) + 1 for chunk in sn_chunks]
        # print('tet', i, 'widths', widths, 'heights', heights)
        tet_intervals = intervals_inside_tet_rectangles[i]
        for interval in tet_intervals:
            # print('we chunk num', interval.we_chunk_num, 'we_in_chunk_index', interval.we_in_chunk_index, 'sn chunk num', interval.sn_chunk_num, 'sn_in_chunk_index', interval.sn_in_chunk_index, 'interval', interval)
            x = interval.we_chunk_num + (interval.we_in_chunk_index + 1)/widths[interval.we_chunk_num]
            y = interval.sn_chunk_num + (interval.sn_in_chunk_index + 1)/heights[interval.sn_chunk_num]
            coords = global_scale_up*(offset + complex(x, y))
            tet_rect_coords[old_tet_rect.horiz_ordering.index(interval)] = coords
            canv.stroke(path.circle(coords.real, coords.imag, 0.035), [deco.filled(), style.linewidth(0)])

        if draw_edge_rectangles:
            edge_rects = old_tet_rect.edge_rectangles()
            print('tet', i, 'num edge rects', len(edge_rects))
            for edge_rect in edge_rects:
                corners = [tet_rect_coords[j] for j in edge_rect]
                draw_edge_rectangle(canv, corners, scl, 0.4*edge_thickness, green, purple)

        if draw_face_rectangles:
            face_rects = old_tet_rect.face_rectangles()
            print('tet', i, 'num face rects', len(face_rects))
            for rect in face_rects:
                verts = [tet_rect_coords[j] for j in rect]
                x_coords = [v.real for v in verts]
                y_coords = [v.imag for v in verts]
                x_interval = (min(x_coords), max(x_coords))
                y_interval = (min(y_coords), max(y_coords))
                draw_rectangle(canv, x_interval, y_interval, scl, 0.4*edge_thickness, green, purple)

        if draw_tetrahedron_rectangles:
            tet_rects = old_tet_rect.tetrahedron_rectangles()
            print('tet', i, 'num tet rects', len(tet_rects))
            for rect in tet_rects:
                verts = [tet_rect_coords[j] for j in rect]
                x_coords = [v.real for v in verts]
                y_coords = [v.imag for v in verts]
                x_interval = (min(x_coords), max(x_coords))
                y_interval = (min(y_coords), max(y_coords))
                draw_rectangle(canv, x_interval, y_interval, scl, 0.4*edge_thickness, green, purple)


    output_filename = 'Images/DrilledTetrahedra/' + name + '_tet_rects' + '.pdf'
    print(output_filename)
    canv.writePDFfile(output_filename)
        