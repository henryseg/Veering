
import pyx ### vector graphics 
from pyx import path, trafo, canvas, style, text, color, deco
from math import sqrt

def HSLtoRGB(Hk,S,L): ### Hue, Saturation, Luminance between 0.0 and 1.0
    if (S == 0.0):
        return [L,L,L]
    Q = 0.0
    if L < 0.5:
        Q = L*(1.0+S)
    else:
        Q = L+S-(L*S)
    P = 2.0 * L - Q
    r = CalcHSLRGB(Hk + (1.0/3.0), P, Q)
    g = CalcHSLRGB(Hk, P, Q)
    b = CalcHSLRGB(Hk - (1.0/3.0), P, Q)
    
    return color.rgb(r,g,b)  

# This function is needed by HSLtoRGB().
def CalcHSLRGB(T, P, Q):
    if (T < 0.0):
        T += 1.0
    elif (T > 1.0):
        T -= 1.0    
    c = 0
    if (T < (1.0/6.0)):
        c = P + ((Q-P)*6.0*T)
    elif (T < (0.5)):
        c = Q
    elif (T < (2.0/3.0)):
        c = P + ((Q-P)*((2.0/3.0)-T)*6.0)
    else:
        c = P
    return c     

def make_line(a, b):
    p = path.path( path.moveto(a.real, a.imag) )
    p.append( path.lineto(b.real, b.imag) )
    return p

def draw_square(i, canv, global_scale_up, global_scale_matrix, edge_thickness, green, purple, draw_inner_lines = True):
    offset = complex(4*i,0)
    if draw_inner_lines:
        indices = [0,1,2,3]
    else:
        indices = [0,3]
    for j in indices:
        p = make_line(complex(0,j) + offset, complex(3,j) + offset)
        p = p.transformed(global_scale_matrix)
        canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, purple, color.transparency(0.75)])
    for j in indices:
        p = make_line(complex(j,0) + offset, complex(j,3) + offset)
        p = p.transformed(global_scale_matrix)
        canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, green, color.transparency(0.75)])
    text_pos = global_scale_up * (complex(1.5, 3.8) + offset)
    canv.text(text_pos.real, text_pos.imag, "$"+str(i)+"$", textattrs=[text.size(3), text.halign.center, text.valign.middle])

def draw_edge_rectangle(canv, corners, global_scale_matrix, edge_thickness, green, purple, global_delta = None):
    center = 0.5*(corners[0] + corners[1])
    x0, x1 = sorted([corner.real for corner in corners])
    y0, y1 = sorted([corner.imag for corner in corners])

    x0, x1 = x0 + global_delta, x1 - global_delta
    y0, y1 = y0 + global_delta, y1 - global_delta
    
    p = make_line(complex(x0, y0), complex(x1, y0))
    p = p.transformed(global_scale_matrix)
    canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, purple])
    
    p = make_line(complex(x1, y1), complex(x0, y1))
    p = p.transformed(global_scale_matrix)
    canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, purple])

    p = make_line(complex(x1, y0), complex(x1, y1))
    p = p.transformed(global_scale_matrix)
    canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, green])
    
    p = make_line(complex(x0, y1), complex(x0, y0))
    p = p.transformed(global_scale_matrix)
    canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, green])

def draw_rectangle(canv, x_interval, y_interval, global_scale_matrix, edge_thickness, green, purple, global_delta = None):
    corners = [complex(x_interval[i], y_interval[i]) for i in range(2)]
    draw_edge_rectangle(canv, corners, global_scale_matrix, edge_thickness, green, purple, global_delta = global_delta)

class rect_wrapper:
    def __init__(self, rect):
        self.rect = rect

def draw_drilled_tetrahedra(con, name = "", tetrahedra_cusp_orders = None, 
    intervals_inside_tet_rectangles = None, 
    tetrahedra_chunks = None,
    old_tet_rectangles = [],
    draw_drilled = True, ### if false, just draw the old tetrahedron rectangles
    draw_edge_rectangles = False,
    draw_face_rectangles = False, 
    draw_tetrahedron_rectangles = True,
    draw_square_inner_lines = True,
    draw_vertex_numbers = False, 
    edge_thickness = 0.02,
    transparency = 0.9,
    global_scale_up = 1.0,
    scale_drawn_elements = 1.0):

    global_delta = 0.04      ### distance to move new rectangles in from dots
    rectangle_spacing = 0.02 ### additional distance for multiple rectangles to one side of a dot
    big_dot_rad = 0.03
    small_dot_rad = 0.02

    red = color.rgb(0.9,0.3,0)
    blue = color.rgb(0,0.3,0.9)
    edge_colours = {True: red, False: blue}
    green = color.rgb(0.0,0.5,0.0)
    purple = color.rgb(0.5,0.0,0.5)
    leaf_colours = {True: green, False: purple}
    grey = color.rgb(0.5,0.5,0.5)
    
    global_scale_matrix = trafo.trafo(matrix=((global_scale_up, 0), (0, global_scale_up)), vector=(0, 0))
    canv = canvas.canvas()
    
    ### find all flow cycles involved
    flow_cycles = set([])
    for tet_intervals in intervals_inside_tet_rectangles:
        for interval in tet_intervals:
            flow_cycles.add(interval.flow_cycle)
    flow_cycles = list(flow_cycles)
    flow_cycles.sort()

    n = len(flow_cycles)
    colours = [HSLtoRGB(float(i)/float(n), 1.0, 0.33) for i in range(n)]
    for i in range(n):
        coords = global_scale_up * (complex(0.0, -0.6 - 0.25*float(i)))
        text_coords = global_scale_up * (complex(0.075, -0.6 - 0.25*float(i)))
        canv.stroke(path.circle(coords.real, coords.imag, big_dot_rad), [deco.filled(), colours[i], style.linewidth(0)])
        canv.text(text_coords.real, text_coords.imag, "$"+str(flow_cycles[i])+"$", textattrs=[text.size(-4), color.rgb(0.4,0.4,0.4), text.halign.left, text.valign.middle])

    ### now that we have drawn the dots for the key, scale everything down.
    global_delta *= scale_drawn_elements
    rectangle_spacing *= scale_drawn_elements
    big_dot_rad *= scale_drawn_elements
    small_dot_rad *= scale_drawn_elements
    edge_thickness *= scale_drawn_elements

    for i in range(len(tetrahedra_cusp_orders)):
        draw_square(i, canv, global_scale_up, global_scale_matrix, edge_thickness, green, purple, draw_inner_lines = draw_square_inner_lines)
        old_tet_rect = old_tet_rectangles[i]
        tet_rect_coords = [None] * len(old_tet_rect.horiz_ordering)

        if draw_drilled:
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

            tet_rect_coords[0] = (W_coords + offset)
            tet_rect_coords[-1] = (E_coords + offset)
            v2hperm = old_tet_rect.vert_to_horiz_perm()
            tet_rect_coords[v2hperm[0]] = (S_coords + offset)
            tet_rect_coords[v2hperm[-1]] = (N_coords + offset)

            Regina_tet = con.vt.tri.tetrahedron(i)
            ### Draw triangle numbers
            for coords, ind in [(W_mid_coords, E_ind), (E_mid_coords, W_ind), (S_mid_coords, N_ind), (N_mid_coords, S_ind)]:
                symbol_coords = global_scale_up*(1.13*(coords - complex(1.5, 1.5)) + complex(1.5, 1.5) + offset)
                tri_num = Regina_tet.triangle(ind).index()
                canv.text(symbol_coords.real, symbol_coords.imag, "$"+str(tri_num)+"$", textattrs=[text.size(-1), text.halign.center, text.valign.middle])

            we_chunks, sn_chunks = tetrahedra_chunks[i]
            widths = [len(chunk) + 1 for chunk in we_chunks]
            heights = [len(chunk) + 1 for chunk in sn_chunks]
            # print('tet', i, 'widths', widths, 'heights', heights)
            tet_intervals = intervals_inside_tet_rectangles[i]
            for interval in tet_intervals:
                # print('we chunk num', interval.we_chunk_num, 'we_in_chunk_index', interval.we_in_chunk_index, 'sn chunk num', interval.sn_chunk_num, 'sn_in_chunk_index', interval.sn_in_chunk_index, 'interval', interval)
                x = interval.we_chunk_num + (interval.we_in_chunk_index + 1)/widths[interval.we_chunk_num]
                y = interval.sn_chunk_num + (interval.sn_in_chunk_index + 1)/heights[interval.sn_chunk_num]
                # coords = global_scale_up*(offset + complex(x, y))
                coords = (offset + complex(x, y))
                tet_rect_coords[old_tet_rect.horiz_ordering.index(interval)] = coords
                # canv.stroke(path.circle(coords.real, coords.imag, 0.035), [deco.filled(), style.linewidth(0)])

            if draw_edge_rectangles:
                edge_rects = old_tet_rect.edge_rectangles()
                # print('tet', i, 'num edge rects', len(edge_rects))
                for edge_rect in edge_rects:
                    corners = [tet_rect_coords[j] for j in edge_rect]
                    draw_edge_rectangle(canv, corners, global_scale_matrix, 0.4*edge_thickness, green, purple)

            if draw_face_rectangles:
                face_rects = old_tet_rect.face_rectangles()
                # print('tet', i, 'num face rects', len(face_rects))
                for rect in face_rects:
                    verts = [tet_rect_coords[j] for j in rect]
                    x_coords = [v.real for v in verts]
                    y_coords = [v.imag for v in verts]
                    x_interval = (min(x_coords), max(x_coords))
                    y_interval = (min(y_coords), max(y_coords))
                    draw_rectangle(canv, x_interval, y_interval, global_scale_matrix, 0.4*edge_thickness, green, purple)

            if draw_tetrahedron_rectangles:

                rectws = [rect_wrapper(rect) for rect in old_tet_rect.tetrahedron_rectangles()]

                ### work out which rectangles hit which side of which dots

                ### rectangles incident to this dot on the dot's E, W, N, S:
                dot_sides = [] 
                for i in range(len(old_tet_rect.horiz_ordering)):
                    dot_sides.append({'E':[], 'W':[], 'N':[], 'S':[]})
                # print('dot sides', dot_sides)

                for rectw in rectws:
                    verts = [tet_rect_coords[j] for j in rectw.rect]
                    rectw.x_coords = [v.real for v in verts]
                    rectw.y_coords = [v.imag for v in verts]

                    rectw.W_ind = min(range(len(rectw.x_coords)), key=rectw.x_coords.__getitem__)
                    rectw.W_dot = rectw.rect[rectw.W_ind]
                    dot_sides[rectw.W_dot]['E'].append(rectw)  ### west side of rect is east side of dot
                    rectw.E_ind = max(range(len(rectw.x_coords)), key=rectw.x_coords.__getitem__)
                    rectw.E_dot = rectw.rect[rectw.E_ind]
                    dot_sides[rectw.E_dot]['W'].append(rectw)
                    rectw.S_ind = min(range(len(rectw.y_coords)), key=rectw.y_coords.__getitem__)
                    rectw.S_dot = rectw.rect[rectw.S_ind]
                    dot_sides[rectw.S_dot]['N'].append(rectw)
                    rectw.N_ind = max(range(len(rectw.y_coords)), key=rectw.y_coords.__getitem__)
                    rectw.N_dot = rectw.rect[rectw.N_ind]
                    dot_sides[rectw.N_dot]['S'].append(rectw)
                # print('dot sides', dot_sides)

                for dot in dot_sides:
                    dot['E'].sort(key = lambda rectw: rectw.x_coords[rectw.E_ind]) ### sort rects to east of dot by where their eastern sides are
                    dot['W'].sort(key = lambda rectw: -rectw.x_coords[rectw.W_ind]) ### - to reverse order
                    dot['N'].sort(key = lambda rectw: rectw.y_coords[rectw.N_ind]) ### sort rects to east of dot by where their eastern sides are
                    dot['S'].sort(key = lambda rectw: -rectw.y_coords[rectw.S_ind]) ### - to reverse order                

                for rectw in rectws:
                    W_place = dot_sides[rectw.W_dot]['E'].index(rectw)
                    E_place = dot_sides[rectw.E_dot]['W'].index(rectw)
                    S_place = dot_sides[rectw.S_dot]['N'].index(rectw)
                    N_place = dot_sides[rectw.N_dot]['S'].index(rectw)
                    x_interval = [min(rectw.x_coords) + rectangle_spacing*W_place, max(rectw.x_coords) - rectangle_spacing*E_place]
                    y_interval = [min(rectw.y_coords) + rectangle_spacing*S_place, max(rectw.y_coords) - rectangle_spacing*N_place]
                    draw_rectangle(canv, x_interval, y_interval, global_scale_matrix, 0.4*edge_thickness, green, purple, global_delta = global_delta)

            ### drilling dots
            for interval in tet_intervals:
                # print('we chunk num', interval.we_chunk_num, 'we_in_chunk_index', interval.we_in_chunk_index, 'sn chunk num', interval.sn_chunk_num, 'sn_in_chunk_index', interval.sn_in_chunk_index, 'interval', interval)
                x = interval.we_chunk_num + (interval.we_in_chunk_index + 1)/widths[interval.we_chunk_num]
                y = interval.sn_chunk_num + (interval.sn_in_chunk_index + 1)/heights[interval.sn_chunk_num]
                coords = global_scale_up*(offset + complex(x, y))
                ### draw drilled dots
                flow_cycle_ind = flow_cycles.index(interval.flow_cycle)
                canv.stroke(path.circle(coords.real, coords.imag, small_dot_rad), [deco.filled(), colours[flow_cycle_ind], style.linewidth(0)])
                # canv.text(coords.real, coords.imag, flow_cycle_ind, textattrs=[text.size(-3), grey, text.halign.center, text.valign.middle])

        ### Draw dots and indices for the vertices of the original tetrahedra
        for coords, ind in [(W_coords, W_ind), (E_coords, E_ind), (S_coords, S_ind), (N_coords, N_ind)]:
            dot_coords = global_scale_up*(coords + offset)
            canv.stroke(path.circle(dot_coords.real, dot_coords.imag, big_dot_rad), [deco.filled(), style.linewidth(0)])
            if draw_vertex_numbers:
                symbol_coords = global_scale_up*(1.11*(coords - complex(1.5, 1.5)) + complex(1.5, 1.5) + offset)
                canv.text(symbol_coords.real, symbol_coords.imag, "$"+str(ind)+"$", textattrs=[text.size(-4), grey, text.halign.center, text.valign.middle])

    output_filename = 'Images/DrilledTetrahedra/' + name + '_tet_rects' + '.pdf'
    print(output_filename)
    canv.writePDFfile(output_filename)
        