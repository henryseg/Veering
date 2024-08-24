
import pyx ### vector graphics 
from pyx import path, trafo, canvas, style, text, color, deco
from math import pi, cos, sin, acos, tan, atan2, sqrt

from continent import continent_tetrahedron
from circular_order import are_anticlockwise, are_linking

def unitize(a):
    return a/abs(a)

def to_complex(a):
    return complex(a[0], a[1])

def make_line(a, b):
    p = path.path( path.moveto(a.real, a.imag) )
    p.append( path.lineto(b.real, b.imag) )
    return p

def dot(p, q):
    return p.real*q.real + p.imag*q.imag

def arc_center_rad(a, b):
    c = b - a
    d = b * dot(a, a)
    e = a * dot(b, b)
    det = 4*(a.real*b.imag - a.imag*b.real)
    center = (2.0/det) * complex(c.imag + d.imag - e.imag, -c.real - d.real + e.real)
    r = abs(a - center)
    return center, r

def reflect_line(a, b, c):
    """reflect c in the line a b"""
    normal = complex(0,1) * (b-a)
    offset = dot(normal, c - a)
    return c - normal * (2*offset)/dot(normal, normal)

def reflect_arc(a, b, c):
    """reflect c in the arc from a to b"""
    center, r = arc_center_rad(a, b)
    cp = c - center
    R_in = dot(cp, cp)
    return ( cp * r * r / (R_in) ) + center

def sign(x):
    if x<0: return -1
    else: return 1

def circle_circle_isects(c1, r1, c2, r2): #worked out in some mathematica notebook "circle_circle_intersect.nb", you'll never find it again
    """intersection points of two circles"""
    c2 = c2 - c1
    x2, y2 = c2.real, c2.imag
    foox = sqrt(-(y2**2)*(-(r1-r2)**2 + x2**2 + y2**2)*(-(r1+r2)**2 + x2**2 + y2**2) )
    fooy = sign(y2)*sqrt(-(-(r1-r2)**2 + x2**2 + y2**2)*(-(r1+r2)**2 + x2**2 + y2**2) )

    u1 = (x2*(r1**2 - r2**2 + x2**2 + y2**2) - foox)/(2*(x2**2 + y2**2))
    u2 = (x2*(r1**2 - r2**2 + x2**2 + y2**2) + foox)/(2*(x2**2 + y2**2))
    v1 = (y2*(r1**2 - r2**2 + x2**2 + y2**2) + x2*fooy)/(2*(x2**2 + y2**2))
    v2 = (y2*(r1**2 - r2**2 + x2**2 + y2**2) - x2*fooy)/(2*(x2**2 + y2**2))
    w1 = complex(u1,v1) + c1
    w2 = complex(u2,v2) + c1
    return w1, w2

def geodesic_isect(a1, b1, a2, b2):
    """find intersection point between arcs in the Poincare disk from a1 to b1 and from a2 to b2"""
    c1, r1 = arc_center_rad(a1, b1)
    c2, r2 = arc_center_rad(a2, b2)
    if r1 + r2 < abs(c1 - c2):
        return None
    else:
        w1, w2 = circle_circle_isects(c1, r1, c2, r2)
        if abs(w1) < 1.0:
            return w1
        else:
            assert abs(w2) < 1.0
            return w2

def incenter(a, b, c):
    det_ab = (a.real*b.imag - a.imag*b.real)
    if abs(det_ab) < 0.0001: 
        cp = reflect_line(a, b, c)
    else:
        cp = reflect_arc(a, b, c)
    det_bc = (b.real*c.imag - b.imag*c.real)
    if abs(det_bc) < 0.0001: 
        ap = reflect_line(b, c, a)
    else:
        ap = reflect_arc(b, c, a)
    pc = make_arc(c, cp)
    pa = make_arc(a, ap)
    isect = pc.intersect(pa)
    x, y = pc.at(isect[0][0])
    # print(x.true())
    # print(float(x), float(y))
    # return complex(x, y)
    return pc.at(isect[0][0])  ### does not return a complex number, it's some horrible pyx thing...

def make_arc(a, b, return_midpt = False):
    det = 4*(a.real*b.imag - a.imag*b.real)
    if abs(det) < 0.0001: 
        if not return_midpt:
            return make_line(a, b)
        else:
            return ( make_line(a, b), 0.5*(a+b) )
    else:
        # theta = 0.5 * acos( a[0]*b[0] + a[1]*b[1] )
        # r = tan(theta)
        # p = path.path( path.moveto(a[0], a[1]) )
        # p.append( path.arct(0.0,0.0,b[0], b[1], r))

        center, r = arc_center_rad(a, b)
        A2 = a - center
        B2 = b - center
        r = abs(A2)
        dot_prod = (A2.real*B2.real + A2.imag*B2.imag)
        theta = 0.5 * acos( dot_prod/(abs(A2)*abs(B2)) )
        direction = (0.5*(a + b) - center)
        direction = (1.0/abs(direction)) * direction
        tangent_point = center + (abs(A2)/cos(theta)) * direction
        
        p = path.path( path.moveto(a.real, a.imag) )
        p.append( path.arct(tangent_point.real, tangent_point.imag, b.real, b.imag, r))

        # p = path.path( path.moveto(a[0], a[1]) )
        # p.append( path.lineto(center.real, center.imag) )
        # p.append( path.lineto(b[0], b[1]) )
        if not return_midpt:
            return p
        else:
            midpt = center + r * direction
            return (p, midpt)

def end_pos(e2, e1, offset = 0.0):
    ### position of end e2 in edge e1
    t = float(e1.ends.index(e2) + 1 + offset)/float(len(e1.ends) + 1)
    u1, v1 = e1.vertices
    u1, v1 = u1.circle_pos, v1.circle_pos     
    return unitize((1-t) * u1 + t * v1)

def end_pos2(thorn_end):
    ### position of ind along edge
    edge, ind = thorn_end
    t = float(ind + 1)/float(len(edge.ends) + 1)
    u1, v1 = edge.vertices
    u1, v1 = u1.circle_pos, v1.circle_pos     
    return unitize((1-t) * u1 + t * v1)

def circle_position(t, n):
    s = 2*pi*t/n
    return complex(cos(s), sin(s))

def tet_rectangle_sides(tet, v):
    ### a side is a pair (v, thorn_end), where v is the cusp and thorn_end is (coastal_arc, index_on_that_arc)
    ### we are looking for at most two cusp leaves to form the sides

    if v in tet.upper_edge.vertices:
        b = tet.upper_edge.other_end[v]
        c, d = tet.lower_edge.vertices      ### FIX, this is old...
        thorn_ends = v.green_thorn_ends
        other_colour_thorn_ends = v.purple_thorn_ends
    else:
        b = tet.lower_edge.other_end[v]
        c, d = tet.upper_edge.vertices
        thorn_ends = v.purple_thorn_ends
        other_colour_thorn_ends = v.green_thorn_ends
    
    for thorn_end in thorn_ends:
        thorn_end_location = thorn_end[0].coastal_index() + 0.5
        # print(v.coastal_index, thorn_end_location, c.coastal_index, d.coastal_index)
        if are_linking(v.coastal_index(), thorn_end_location, c.coastal_index(), d.coastal_index()):
            ### we found the cusp leaf that goes through the tetrahedron
            if len(other_colour_thorn_ends) == 0:
                # print('nn')
                return [ None, None ]
            for i, other_colour_thorn_end in enumerate(other_colour_thorn_ends):
                # print('i', i, v.coastal_index, thorn_end_location, other_colour_thorn_end[0].coastal_index + 0.5, are_anticlockwise(v.coastal_index, thorn_end_location, other_colour_thorn_ends[0][0].coastal_index + 0.5))
                if are_anticlockwise(v.coastal_index(), thorn_end_location, other_colour_thorn_end[0].coastal_index() + 0.5):
                    if i == 0:  ### 0th other colour thorn is after thorn_end
                        # print('ns')
                        return [ None, (v, other_colour_thorn_ends[0]) ]
                    else:
                        # print('ss')
                        return [ (v, other_colour_thorn_ends[i-1]), (v, other_colour_thorn_ends[i]) ]
            ### didn't find an other colour thorn after thorn_end
            # print('sn')
            return [ (v, other_colour_thorn_ends[-1]), None ]

def complete_tetrahedron_rectangles(con, tetrahedra_to_complete):
    """grow the continent so that the given tetrahedra have full tetrahedron rectangles within the continent"""
    k = 0
    for tet in tetrahedra_to_complete:
        for v in tet.vertices:
            # print('tet vert age', con.vertices.index(v))
            # con.build_boundary_data()
            # con.install_thorn_ends()
            sides = tet_rectangle_sides(tet, v)
            for direction in range(2):
                while sides[direction] == None:
                    # print('direction, k', direction, k)
                    e = con.coastal_edges[(v.coastal_index() - direction)%len(con.coast)] 
                    triangles = e.boundary_triangles()  ### grow around this edge
                    if triangles[0].is_upper != (k % 2 == 0): ### xor, alternate which one we add to
                        con.bury(triangles[0])
                    else:
                        con.bury(triangles[1])
                    # con.build_boundary_data()
                    # con.install_thorn_ends()
                    sides = tet_rectangle_sides(tet, v)
                    k += 1
                    if k > 50:
                        print('bail')
                        return None

def tet_purple_rectangle_sides(tet, actually_do_green = False):
    out = []
    if not actually_do_green:
        verts = tet.upper_edge().vertices
    else:
        verts = tet.lower_edge().vertices
    for v in verts:
        out.append(tet_rectangle_sides(tet, v))
    return out

def draw_edge(e, edge_thickness, edge_colours, scl, canv):
    col = edge_colours[e.is_red]
    start, end = e.drawing_end_positions()
    start = circle_position(start, len(e.continent.coast))
    end = circle_position(end, len(e.continent.coast))
    p = make_arc(start, end)
    p = p.transformed(scl)
    canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, col])

def draw_triangle(f, edge_thickness, edge_colours, scl, canv):
    for e in f.edges:
        draw_edge(e, edge_thickness, edge_colours, scl, canv)

def draw_tetrahedron(t, edge_thickness, edge_colours, scl, canv):
    draw_edge(t.lower_edge, edge_thickness, edge_colours, scl, canv)
    for e in t.equatorial_edges:
        draw_edge(e, edge_thickness, edge_colours, scl, canv)
    draw_edge(t.upper_edge, edge_thickness, edge_colours, scl, canv)

def draw_leaf(leaf, leaf_thickness, leaf_colours, scl, canv):
    col = leaf_colours[leaf.is_upper]
    start, end = leaf.drawing_end_positions()
    start = circle_position(start, len(leaf.continent.coast))
    end = circle_position(end, len(leaf.continent.coast))
    p = make_arc(start, end)
    p = p.transformed(scl)
    canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, col])

def draw_edge_rectangle_half(l, m, leaf_thickness, leaf_colours, scl, canv):
    if l == None and m != None:
        draw_leaf(m, leaf_thickness, leaf_colours, scl, canv)
        ### doesnt return anything, can make other things break...
    elif m == None and l != None:
        draw_leaf(l, leaf_thickness, leaf_colours, scl, canv)
    elif l != None and m != None:
        l_start, l_end = l.drawing_end_positions()
        l_start = circle_position(l_start, len(l.continent.coast))
        l_end = circle_position(l_end, len(l.continent.coast))
        m_start, m_end = m.drawing_end_positions()
        m_start = circle_position(m_start, len(m.continent.coast))
        m_end = circle_position(m_end, len(m.continent.coast))

        intersection = geodesic_isect(l_start, l_end, m_start, m_end)
        if intersection == None:
            print('no intersection')
            draw_leaf(m, 2*leaf_thickness, leaf_colours, scl, canv)
            draw_leaf(l, 2*leaf_thickness, leaf_colours, scl, canv)
            return None
        # assert intersection != None

        p = make_arc(l_start, intersection)
        p = p.transformed(scl)
        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, leaf_colours[l.is_upper]])
        q = make_arc(intersection, m_start)
        q = q.transformed(scl)
        canv.stroke(q, [style.linewidth(leaf_thickness), style.linecap.round, leaf_colours[m.is_upper]])
        return [p,q]

def draw_continent_circle(con, name = "", draw_labels = True, draw_upper_landscape = False, draw_lower_landscape = False, 
    draw_coastal_edges = True, draw_all_edges = False,
    draw_cusp_leaves = True,
    shade_triangles = False, draw_fund_domain = False, fund_dom_tets = None, draw_fund_domain_edges = False,
    edge_rectangles_to_draw = [],
    tetrahedron_rectangles_to_draw = [],
    tetrahedron_rectangles_to_shade = [],
    triangles_to_draw = [],
    tetrahedra_to_draw = [],
    draw_edges_for_edge_rectangles = False,
    leaves_to_draw = [],
    intervals_to_draw = [],
    edge_thickness = 0.02,
    leaf_thickness = 0.03,
    transparency = 0.9,
    text_size = -4):
    
    # ### check edge rectangle
    # check_edge = con.edges[240]
    # print(len(con.tetrahedra))
    # check_edge.ensure_continent_contains_rectangle()
    # print(len(con.tetrahedra))

    global_scale_up = 25.0
    red = color.rgb(0.9,0.3,0)
    blue = color.rgb(0,0.3,0.9)
    edge_colours = {True: red, False: blue}
    green = color.rgb(0.0,0.5,0.0)
    purple = color.rgb(0.5,0.0,0.5)
    leaf_colours = {True: green, False: purple}
    
    scl = trafo.trafo(matrix=((global_scale_up, 0), (0, global_scale_up)), vector=(0, 0))
    canv = canvas.canvas()
    canv.stroke(path.circle(0,0,global_scale_up), [style.linewidth(edge_thickness)])

    con.build_boundary_data()

    n = len(con.coast)
    for v in con.coast:
        i = v.coastal_index()
        t = 2*pi*float(i)/float(n)
        v.circle_pos = complex(cos(t), sin(t))
        vert_pos = v.circle_pos * 1.01 * global_scale_up
        if draw_labels:
            # label = "$"+str(con.vertices.index(v))+"$"
            label = "$"+str(v.Regina_cusp_num)+"$"
            canv.text(vert_pos.real, vert_pos.imag, label, textattrs=[text.size(text_size), text.halign.left, text.valign.middle,trafo.rotate((180/pi) * atan2(vert_pos.imag, vert_pos.real))])

        # vert_pos2 = v.circle_pos * 1.2 * global_scale_up
        # p = path.path(path.moveto(vert_pos.real, vert_pos.imag), path.lineto(vert_pos2.real, vert_pos2.imag))
        # canv.stroke(p, [deco.curvedtext("$"+str(con.vertices.index(v))+"$")])

    for e in con.edges:  
    ### would be nice to get edges drawn in the correct order. upper bdry edges are on top, lower bdry edges are below, coast doesnt matter
    ### internal edges (we think) have their edge rectangles in cusp leaves. Then check who spans who.
    ### alternatively, draw lower landscape first then flip up through tetrahedra.
        if draw_all_edges or (draw_coastal_edges and e.is_coastal()):
            draw_edge(e, edge_thickness, edge_colours, scl, canv)

    ### highlight vertices of tetrahedra in a fundamental domain
    if draw_fund_domain:
        for con_tet in fund_dom_tets:
            if type(con_tet) == continent_tetrahedron:  ### could be an integer if we didnt find this tet
                if draw_fund_domain_edges:
                    for e in con_tet.edges():
                        draw_edge(e, 1.5*edge_thickness, edge_colours, scl, canv)

        for v in con.coast:
            v.fund_dom_tet_nums = []
        for con_tet in fund_dom_tets:
            for v in con_tet.vertices:
                v.fund_dom_tet_nums.append(con_tet.index)
        for v in [v for v in con.coast if v.fund_dom_tet_nums != []]:
            vert_pos = v.circle_pos * 1.03 * global_scale_up
            if draw_labels:
                canv.text(vert_pos.real, vert_pos.imag, "$"+str(v.fund_dom_tet_nums)+"$", textattrs=[text.size(text_size), text.halign.left, text.valign.middle,trafo.rotate((180/pi) * atan2(vert_pos.imag, vert_pos.real))])


    # lower_colours = {True: color.rgb(0.5,0.3,0), False: color.rgb(0,0.3,0.5)}
    # upper_colours = {True: color.rgb(0.9,0.3,0), False: color.rgb(0,0.3,0.9)}


    landscape_edges = [con.lower_landscape_edges, con.upper_landscape_edges]

    # colours = [lower_colours, upper_colours]

    upper_tris = con.upper_landscape_triangles
    lower_tris = con.lower_landscape_triangles
    boundary_tris = [lower_tris, upper_tris]

    if shade_triangles:
        u,v,w = con.triangle_path[0].vertices
        p = make_arc(u.circle_pos, v.circle_pos)
        q = make_arc(v.circle_pos, w.circle_pos)
        r = make_arc(w.circle_pos, u.circle_pos)
        p.append(q[1])
        p.append(r[1])  ### remove extraneous moveto commands
        p = p.transformed(scl)
        canv.stroke(p, [deco.filled([color.transparency(transparency)]), style.linewidth(0)])
        u,v,w = con.triangle_path[-1].vertices
        p = make_arc(u.circle_pos, v.circle_pos)
        q = make_arc(v.circle_pos, w.circle_pos)
        r = make_arc(w.circle_pos, u.circle_pos)
        p.append(q[1])
        p.append(r[1])  ### remove extraneous moveto commands
        p = p.transformed(scl)
        canv.stroke(p, [deco.filled([color.transparency(transparency)]), style.linewidth(0)])
        # for triangle in con.triangle_path:
        #     u,v,w = triangle.vertices
        #     p = make_arc(u.circle_pos, v.circle_pos)
        #     q = make_arc(v.circle_pos, w.circle_pos)
        #     r = make_arc(w.circle_pos, u.circle_pos)
        #     p.append(q[1])
        #     p.append(r[1])  ### remove extraneous moveto commands
        #     p = p.transformed(scl)
        #     canv.stroke(p, [deco.filled([color.transparency(0.8)]), style.linewidth(0)])

    to_do = []
    if draw_lower_landscape:
        to_do.append(0)
    if draw_upper_landscape:
        to_do.append(1)
    for i in to_do:
        for e in landscape_edges[i]:
            col = edge_colours[e.is_red]
            transp = []
            if i == 0:
                transp = [color.transparency(0.75)]
            u, v = e.vertices
            p = make_arc(u.circle_pos, v.circle_pos)
            p = p.transformed(scl)
            canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, col] + transp)
        for tri in boundary_tris[i]:
            center = incenter(tri.vertices[0].circle_pos, tri.vertices[1].circle_pos, tri.vertices[2].circle_pos)
            # canv.fill(path.circle(global_scale_up*center[0], global_scale_up*center[1], 0.1)) 
            tri_name = str(tri.index) ## in quotient manifold
            tri_name = str(tri.continent.triangles.index(tri)) ## in continent
            if draw_labels:
                canv.text(global_scale_up*center[0], global_scale_up*center[1], "$"+tri_name+"$", textattrs=[text.size(-2), text.halign.center, text.valign.middle] + transp)

    if draw_cusp_leaves:
        leaves = { green: [], purple: []}
        for i, c in enumerate(con.coast):
            for leaf in c.cusp_leaves:
                draw_leaf(leaf, leaf_thickness, leaf_colours, scl, canv)
        #         if leaf.is_upper:
        #             leaves[green].append(leaf)
        #         else:
        #             leaves[purple].append(leaf)
        # for col in [purple, green]:
        #     for leaf in leaves[col]:
        #         start, end = leaf.drawing_end_positions()
        #         start = circle_position(start, len(con.coast))
        #         end = circle_position(end, len(con.coast))
        #         p = make_arc(start, end)
        #         p = p.transformed(scl)
        #         canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, col])

    for leaf in leaves_to_draw:
        draw_leaf(leaf, 3*leaf_thickness, leaf_colours, scl, canv)

    for e in edge_rectangles_to_draw:
        if draw_edges_for_edge_rectangles:
            draw_edge(e, 3*edge_thickness, edge_colours, scl, canv)

        #     -3-v
        #    |  /|
        #    0 e 2
        #    |/  |
        #    u-1-
        a, b, c, d = e.rectangle_sides()
        draw_edge_rectangle_half(a, d, leaf_thickness, leaf_colours, scl, canv)
        draw_edge_rectangle_half(b, c, leaf_thickness, leaf_colours, scl, canv)

    for f in triangles_to_draw:
        draw_triangle(f, edge_thickness, edge_colours, scl, canv)

    for t in tetrahedra_to_draw:
        draw_tetrahedron(t, edge_thickness, edge_colours, scl, canv)

    for t in tetrahedron_rectangles_to_draw:
        eq_edges = t.equatorial_edges
        t_rect_arcs = []
        for i, e in enumerate(eq_edges):
            before_e_verts = eq_edges[(i-1) % 4].vertices
            after_e_verts = eq_edges[(i+1) % 4].vertices
            u, v = e.vertices
            if u in before_e_verts:
                assert v in after_e_verts
                p = e.rectangle_sides()[1]
                q = e.rectangle_sides()[2]
            else:
                assert u in after_e_verts
                assert v in before_e_verts
                p = e.rectangle_sides()[3]
                q = e.rectangle_sides()[0]
            arcs = draw_edge_rectangle_half(p, q, leaf_thickness, leaf_colours, scl, canv)
            t_rect_arcs.extend(arcs)

        assert len(t_rect_arcs) == 8
        if t in tetrahedron_rectangles_to_shade:
            polygon = t_rect_arcs[0] << t_rect_arcs[1] << t_rect_arcs[2] << t_rect_arcs[3] << t_rect_arcs[4] << t_rect_arcs[5] << t_rect_arcs[6] << t_rect_arcs[7]
            canv.stroke(polygon, [deco.filled([color.transparency(transparency)]), style.linewidth(0)])

    for interval in intervals_to_draw:
        g, p = interval.crossing_leaves() ### careful, this can grow the continent. So we run it first outside the draw function, in scripts.py
        assert len(g) == 2 and len(p) == 2

        # for leaf in g:
        #     draw_leaf(leaf, leaf_thickness, leaf_colours, scl, canv)
        # for leaf in p:
        #     draw_leaf(leaf, leaf_thickness, leaf_colours, scl, canv)

        intersections2 = []
        for l in g:
            intersections1 = []
            l_start, l_end = l.drawing_end_positions()
            l_start = circle_position(l_start, len(l.continent.coast))
            l_end = circle_position(l_end, len(l.continent.coast))
            for m in p:
                m_start, m_end = m.drawing_end_positions()
                m_start = circle_position(m_start, len(m.continent.coast))
                m_end = circle_position(m_end, len(m.continent.coast))

                intersection = geodesic_isect(l_start, l_end, m_start, m_end)
                intersections1.append(intersection)
            intersections2.append(intersections1)
        arcs = []
        p = make_arc(intersections2[0][0], intersections2[0][1])
        p = p.transformed(scl)
        arcs.append(p)
        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, green])

        q = make_arc(intersections2[0][1], intersections2[1][1])
        q = q.transformed(scl)
        arcs.append(q)
        canv.stroke(q, [style.linewidth(leaf_thickness), style.linecap.round, purple])

        p = make_arc(intersections2[1][1], intersections2[1][0])
        p = p.transformed(scl)
        arcs.append(p)
        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, green])

        q = make_arc(intersections2[1][0], intersections2[0][0])
        q = q.transformed(scl)
        arcs.append(q)
        canv.stroke(q, [style.linewidth(leaf_thickness), style.linecap.round, purple])

        polygon = arcs[0] << arcs[1] << arcs[2] << arcs[3]
        canv.stroke(polygon, [deco.filled([color.transparency(transparency)]), style.linewidth(0)])

        text_pos = intersections2[0][0]

        text_offset = 0.25
        text_pos = global_scale_up * text_pos + complex(0, text_offset*interval.owning_tet_index) ### move labels off of each other 
        canv.text(text_pos.real, text_pos.imag, "$"+interval.name+"$", textattrs=[text.size(text_size), text.halign.left, text.valign.middle])





        ### draw full leaves
        # for leaf in e.rectangle_sides():
        #     if leaf != None:
        #         draw_leaf(leaf, leaf_thickness, leaf_colours, scl, canv)






                # if leaf.is_upper:
                #     col = green
                # else:
                #     col = purple
                # start, end = leaf.drawing_end_positions()
                # start = circle_position(start, len(con.coast))
                # end = circle_position(end, len(con.coast))
                # p = make_arc(start, end)
                # p = p.transformed(scl)
                # canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, col])


    # ### draw edge rectangle
    # col = edge_colours[check_edge.is_red]
    # u, v = check_edge.vertices
    # p = make_arc(u.circle_pos, v.circle_pos)
    # p = p.transformed(scl)
    # canv.stroke(p, [style.linewidth(2*edge_thickness), style.linecap.round, col])
    # cusp_leaves = check_edge.rectangle_sides()
    # for leaf in cusp_leaves:
    #     if leaf != None:
    #         if leaf.is_upper:
    #             col = green
    #         else:
    #             col = purple
    #         start, end = leaf.drawing_end_positions()
    #         start = circle_position(start, len(con.coast))
    #         end = circle_position(end, len(con.coast))
    #         p = make_arc(start, end)
    #         p = p.transformed(scl)
    #         canv.stroke(p, [style.linewidth(2*leaf_thickness), style.linecap.round, col])



    #         for tet in draw_tetrahedron_rectangles:
    #             a, c = tet.upper_edge().vertices
    #             b, d = tet.lower_edge().vertices
    #             if not are_anticlockwise(a.coastal_index, b.coastal_index, c.coastal_index):
    #                 b, d = d, b   ### now a, b, c, d are anticlockwise
    #             # print(con.vertices.index(a), con.vertices.index(b), con.vertices.index(c), con.vertices.index(d))
    #             all_sides_discrete = [tet_rectangle_sides(tet, v) for v in [a,b,c,d]]
    #             all_sides_geometry = []
    #             for vertex_sides_discrete in all_sides_discrete:
    #                 vertex_sides_geometry = []
    #                 for side_discrete in vertex_sides_discrete:
    #                     a, b = side_discrete
    #                     if side_discrete != None:
    #                         v, thorn_end = side_discrete
    #                         # print(thorn_end)
    #                         # print((v.circle_pos, end_pos2(thorn_end)))
    #                         vertex_sides_geometry.append( (v.circle_pos, end_pos2(thorn_end)) )
    #                     else:
    #                         vertex_sides_geometry.append( None )
    #                 all_sides_geometry.append(vertex_sides_geometry)
    #             for i in range(4):
    #                 if all_sides_geometry[i][0] != None and all_sides_geometry[(i+1)%4][1] != None:
    #                     v1, t1 = all_sides_geometry[i][0]
    #                     v2, t2 = all_sides_geometry[(i+1)%4][1]
    #                     intersection = geodesic_isect(v1, t1, v2, t2)
    #                     assert intersection != None
    #                     all_sides_geometry[i][0] = (v1, intersection)
    #                     all_sides_geometry[(i+1)%4][1] = (v2, intersection)

    #             for i, vertex_sides_geometry in enumerate(all_sides_geometry):
    #                 if i%2 == 0:
    #                     col = purple
    #                 else:
    #                     col = green
    #                 for side_geometry in vertex_sides_geometry: 
    #                     if side_geometry != None:
    #                         # print(side_geometry)
    #                         v_pos, t_pos = side_geometry 
    #                         p = make_arc(v_pos, t_pos)
    #                         p = p.transformed(scl)
    #                         canv.stroke(p, [style.linewidth(3*leaf_thickness), style.linecap.round, col])

    output_filename = 'Images/CircleContinent/' + name + '.pdf'
    print(output_filename)
    canv.writePDFfile(output_filename)



