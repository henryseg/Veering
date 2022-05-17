
import pyx ### vector graphics 
from pyx import path, trafo, canvas, style, text, color, deco
from math import pi, cos, sin, acos, tan, atan2, sqrt
from taut import isosig_to_tri_angle, edge_num_to_vert_pair, vert_pair_to_edge_num
from transverse_taut import edge_side_face_collections, get_tet_above_edge, get_tet_top_vert_nums
from veering import veering_triangulation
from continent import continent, continent_tetrahedron
from boundary_triangulation import tet_face

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
            return ( make_line(a, b), (0.5*(a[0]+b[0]), 0.5*(a[1]+b[1])) )
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

def are_anticlockwise(a, b, c): ### given indices in this order, are they anticlockwise on the circle?
    return (b - a)*(c - b)*(a - c) < 0  ### either one or two are negative

def are_linking(a1, a2, b1, b2): ### given indices on a circle, do the a's link the b's?
    return (a1 - b1)*(a1 - b2)*(a2 - b1)*(a2 - b2) < 0 ### same as the sign of the cross-ratio

def install_thorn_ends(con):
    """For each cusp c, install c.purple_thorn_ends = [] ### [coastal arc, position along that arc] 
    and same for green_thorn_ends"""
    purple_train_routes = []  ### pairs of coastal edges corresponding to a train route
    green_train_routes = []
  
    for edge in con.lower_landscape_edges:
        leaf_end_edges = []
        if edge.is_coastal():
            if not edge.is_coastal_sink(upper = False):
                leaf_end_edges.append(edge)
                for tri in edge.boundary_triangles:
                    if not tri.is_upper:
                        last_tri = con.flow(tri)[0]
                        last_edge = last_tri.edges[last_tri.downriver_index()]
                        leaf_end_edges.append(last_edge)
        else:
            if edge.is_watershed():
                for tri in edge.boundary_triangles:
                    last_tri = con.flow(tri)[0]
                    last_edge = last_tri.edges[last_tri.downriver_index()]
                    leaf_end_edges.append(last_edge)
        if len(leaf_end_edges) == 2:
            purple_train_routes.append(leaf_end_edges)

        
    for edge in con.upper_landscape_edges:
        leaf_end_edges = []
        if edge.is_coastal():
            if not edge.is_coastal_sink(upper = True):
                leaf_end_edges.append(edge)
                for tri in edge.boundary_triangles:
                    if tri.is_upper:
                        last_tri = con.flow(tri)[0]
                        last_edge = last_tri.edges[last_tri.downriver_index()]
                        leaf_end_edges.append(last_edge)
        else:
            if edge.is_watershed():
                for tri in edge.boundary_triangles:
                    last_tri = con.flow(tri)[0]
                    last_edge = last_tri.edges[last_tri.downriver_index()]
                    leaf_end_edges.append(last_edge)
        if len(leaf_end_edges) == 2:
            green_train_routes.append(leaf_end_edges)
                
    for e in con.coastal_edges:
        e.purple_ends = []  ### of train routes
        e.green_ends = []
    for e1, e2 in purple_train_routes:
        assert e1.is_coastal()
        assert e2.is_coastal()
        e1.purple_ends.append(e2)
        e2.purple_ends.append(e1)
    for e1, e2 in green_train_routes:
        e1.green_ends.append(e2)
        e2.green_ends.append(e1)
    for i, e in enumerate(con.coastal_edges):
        rotated_coastal_edges = con.coastal_edges[i:] + con.coastal_edges[:i]
        e.purple_ends.sort(key = lambda e_other:rotated_coastal_edges.index(e_other), reverse = True)
        e.green_ends.sort(key = lambda e_other:rotated_coastal_edges.index(e_other), reverse = True)
        if e.is_red:
            e.ends = e.green_ends + e.purple_ends
        else:
            e.ends = e.purple_ends + e.green_ends
    
    for i, c in enumerate(con.coast):
        c.purple_thorn_ends = [] ### [coastal arc, position along that arc]
        e = con.coastal_edges[i]
        e1 = e.purple_ends[0]
        while True:
            index = e1.purple_ends.index(e)
            if index == len(e1.purple_ends) - 1:
                assert con.coast[ (con.coastal_edges.index(e1) + 1) % len(con.coast) ] == c
                break
            else:
                c.purple_thorn_ends.append( (e1, e1.ends.index(e)) )
                e, e1 = e1, e1.purple_ends[index + 1]
        
    for i, c in enumerate(con.coast):
        c.green_thorn_ends = [] ### [coastal arc, position along that arc]
        e = con.coastal_edges[i]  ### immediately after the cusp
        e1 = e.green_ends[0] ### of train routes, ordered counterclockwise along the edge
        while True: # go around the crown counterclockwise, meaning that the thorn_ends are ordered clockwise around c
            index = e1.green_ends.index(e)
            if index == len(e1.green_ends) - 1:
                assert con.coast[ (con.coastal_edges.index(e1) + 1) % len(con.coast) ] == c
                break
            else:
                c.green_thorn_ends.append( (e1, e1.ends.index(e)) )
                e, e1 = e1, e1.green_ends[index + 1]

def tet_rectangle_sides(tet, v):
    ### a side is a pair (v, thorn_end), where v is the cusp and thorn_end is (coastal_arc, index_on_that_arc)
    ### we are looking for at most two cusp leaves to form the sides

    if v in tet.upper_edge().vertices:
        b = tet.upper_edge().other_end[v]
        c, d = tet.lower_edge().vertices
        thorn_ends = v.green_thorn_ends
        other_colour_thorn_ends = v.purple_thorn_ends
    else:
        b = tet.lower_edge().other_end[v]
        c, d = tet.upper_edge().vertices
        thorn_ends = v.purple_thorn_ends
        other_colour_thorn_ends = v.green_thorn_ends
    
    for thorn_end in thorn_ends:
        thorn_end_location = thorn_end[0].coastal_index + 0.5
        # print(v.coastal_index, thorn_end_location, c.coastal_index, d.coastal_index)
        if are_linking(v.coastal_index, thorn_end_location, c.coastal_index, d.coastal_index):
            ### we found the cusp leaf that goes through the tetrahedron
            if len(other_colour_thorn_ends) == 0:
                # print('nn')
                return [ None, None ]
            for i, other_colour_thorn_end in enumerate(other_colour_thorn_ends):
                # print('i', i, v.coastal_index, thorn_end_location, other_colour_thorn_end[0].coastal_index + 0.5, are_anticlockwise(v.coastal_index, thorn_end_location, other_colour_thorn_ends[0][0].coastal_index + 0.5))
                if are_anticlockwise(v.coastal_index, thorn_end_location, other_colour_thorn_end[0].coastal_index + 0.5):
                    if i == 0:  ### 0th other colour thorn is after thorn_end
                        # print('ns')
                        return [ None, (v, other_colour_thorn_ends[0]) ]
                    else:
                        # print('ss')
                        return [ (v, other_colour_thorn_ends[i-1]), (v, other_colour_thorn_ends[i]) ]
            ### didn't find an other colour thorn after thorn_end
            # print('sn')
            return [ (v, other_colour_thorn_ends[-1]), None ]

def tet_purple_rectangle_sides(tet, actually_do_green = False):
    out = []
    if not actually_do_green:
        verts = tet.upper_edge().vertices
    else:
        verts = tet.lower_edge().vertices
    for v in verts:
        out.append(tet_rectangle_sides(tet, v))
    return out

def draw_continent_circle(con, name = "", draw_upper_landscape = True, draw_lower_landscape = False, 
    draw_upper_green = True, draw_lower_purple = False, draw_train_tracks = False, 
    draw_foliation = True, foliation_style_old = False, foliation_style_split = False, 
    foliation_style_cusp_leaves = True, foliation_style_boundary_leaves = True,
    shade_triangles = False, draw_fund_domain = False, fund_dom_tets = None, draw_fund_domain_edges = False,
    draw_tetrahedron_rectangles = []):
    
    global_scale_up = 10.0
    edge_thickness = 0.02
    track_thickness = 0.02
    leaf_thickness = 0.03
    edge_colours = {True: color.rgb(0.9,0.3,0), False: color.rgb(0,0.3,0.9)}
    green = color.rgb(0.0,0.5,0.0)
    purple = color.rgb(0.5,0.0,0.5)

    scl = trafo.trafo(matrix=((global_scale_up, 0), (0, global_scale_up)), vector=(0, 0))
    canv = canvas.canvas()
    canv.stroke(path.circle(0,0,global_scale_up), [style.linewidth(0.02)])

    n = len(con.coast)
    for v in con.coast:
        i = v.coastal_index
        t = 2*pi*float(i)/float(n)
        v.circle_pos = complex(cos(t), sin(t))
        vert_pos = v.circle_pos * 1.01 * global_scale_up
        canv.text(vert_pos.real, vert_pos.imag, "$"+str(con.vertices.index(v))+"$", textattrs=[text.size(-4), text.halign.left, text.valign.middle,trafo.rotate((180/pi) * atan2(vert_pos.imag, vert_pos.real))])

        # vert_pos2 = v.circle_pos * 1.2 * global_scale_up
        # p = path.path(path.moveto(vert_pos.real, vert_pos.imag), path.lineto(vert_pos2.real, vert_pos2.imag))
        # canv.stroke(p, [deco.curvedtext("$"+str(con.vertices.index(v))+"$")])

    ### highlight vertices of tetrahedra in a fundamental domain
    if draw_fund_domain:
        if fund_dom_tets == None:
            fund_dom_tets = get_fund_domain_tetrahedra(con)
        for con_tet in fund_dom_tets:
            if type(con_tet) == continent_tetrahedron:  ### could be an integer if we didnt find this tet
                if draw_fund_domain_edges:
                    for e in con_tet.edges():
                        col = edge_colours[e.is_red]
                        u, v = e.vertices
                        p = make_arc(u.circle_pos, v.circle_pos)
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(edge_thickness), style.linecap.round, col])

        update_fund_dom_tet_nums(con, fund_dom_tets)
        for v in [v for v in con.coast if v.fund_dom_tet_nums != []]:
            vert_pos = v.circle_pos * 1.03 * global_scale_up
            canv.text(vert_pos.real, vert_pos.imag, "$"+str(v.fund_dom_tet_nums)+"$", textattrs=[text.size(-4), text.halign.left, text.valign.middle,trafo.rotate((180/pi) * atan2(vert_pos.imag, vert_pos.real))])


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
        canv.stroke(p, [deco.filled([color.transparency(0.8)]), style.linewidth(0)])
        u,v,w = con.triangle_path[-1].vertices
        p = make_arc(u.circle_pos, v.circle_pos)
        q = make_arc(v.circle_pos, w.circle_pos)
        r = make_arc(w.circle_pos, u.circle_pos)
        p.append(q[1])
        p.append(r[1])  ### remove extraneous moveto commands
        p = p.transformed(scl)
        canv.stroke(p, [deco.filled([color.transparency(0.8)]), style.linewidth(0)])
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
            canv.text(global_scale_up*center[0], global_scale_up*center[1], "$"+str(tri.index)+"$", textattrs=[text.size(-2), text.halign.center, text.valign.middle] + transp)

    ### train tracks...

    purple_train_routes = []  ### pairs of coastal edges corresponding to a train route
    green_train_routes = []
    if draw_lower_purple:
        if draw_train_tracks:
            for tri in lower_tris:
                midpts = []
                is_reds = []
                for e in tri.edges:
                    is_reds.append(e.is_red)
                    u, v = e.vertices
                    p, midpt = make_arc(u.circle_pos, v.circle_pos, return_midpt = True)
                    midpts.append(midpt)
                for i in range(3):
                    if (is_reds[i] == is_reds[(i+1)%3]) or (not is_reds[i] and is_reds[(i+1)%3]):
                        p = make_arc(midpts[i], midpts[(i+1)%3])
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(track_thickness), style.linecap.round, purple])    
        if draw_foliation:   
            for edge in con.lower_landscape_edges:
                leaf_end_edges = []
                if edge.is_coastal():
                    if not edge.is_coastal_sink(upper = False):
                        leaf_end_edges.append(edge)
                        for tri in edge.boundary_triangles:
                            if not tri.is_upper:
                                last_tri = con.flow(tri)[0]
                                last_edge = last_tri.edges[last_tri.downriver_index()]
                                leaf_end_edges.append(last_edge)
                else:
                    if edge.is_watershed():
                        for tri in edge.boundary_triangles:
                            last_tri = con.flow(tri)[0]
                            last_edge = last_tri.edges[last_tri.downriver_index()]
                            leaf_end_edges.append(last_edge)
                if len(leaf_end_edges) == 2:
                    purple_train_routes.append(leaf_end_edges)
                    if foliation_style_old:
                        leaf_ends = []
                        for e in leaf_end_edges:
                            endpts = e.vertices
                            _, midpt = make_arc(endpts[0].circle_pos, endpts[1].circle_pos, return_midpt = True)
                            leaf_ends.append(midpt)
                        p = make_arc(leaf_ends[0], leaf_ends[1])
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, purple])

    if draw_upper_green:
        if draw_train_tracks:
            for tri in upper_tris:
                midpts = []
                is_reds = []
                for e in tri.edges:
                    is_reds.append(e.is_red)
                    u, v = e.vertices
                    p, midpt = make_arc(u.circle_pos, v.circle_pos, return_midpt = True)
                    midpts.append(midpt)
                for i in range(3):
                    if (is_reds[i] == is_reds[(i+1)%3]) or (is_reds[i] and not is_reds[(i+1)%3]):
                        p = make_arc(midpts[i], midpts[(i+1)%3])
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(track_thickness), style.linecap.round, green])
        if draw_foliation:
            for edge in con.upper_landscape_edges:
                leaf_end_edges = []
                if edge.is_coastal():
                    if not edge.is_coastal_sink(upper = True):
                        leaf_end_edges.append(edge)
                        for tri in edge.boundary_triangles:
                            if tri.is_upper:
                                last_tri = con.flow(tri)[0]
                                last_edge = last_tri.edges[last_tri.downriver_index()]
                                leaf_end_edges.append(last_edge)
                else:
                    if edge.is_watershed():
                        for tri in edge.boundary_triangles:
                            last_tri = con.flow(tri)[0]
                            last_edge = last_tri.edges[last_tri.downriver_index()]
                            leaf_end_edges.append(last_edge)
                if len(leaf_end_edges) == 2:
                    green_train_routes.append(leaf_end_edges)
                    if foliation_style_old:
                        leaf_ends = []
                        for e in leaf_end_edges:
                            endpts = e.vertices
                            _, midpt = make_arc(endpts[0].circle_pos, endpts[1].circle_pos, return_midpt = True)
                            leaf_ends.append(midpt)
                        p = make_arc(leaf_ends[0], leaf_ends[1])
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, green])
    
    if draw_foliation and (foliation_style_split or foliation_style_cusp_leaves or foliation_style_boundary_leaves):
        for e in con.coastal_edges:
            e.purple_ends = []
            e.green_ends = []
        for e1, e2 in purple_train_routes:
            e1.purple_ends.append(e2)
            e2.purple_ends.append(e1)
        for e1, e2 in green_train_routes:
            e1.green_ends.append(e2)
            e2.green_ends.append(e1)
        for i, e in enumerate(con.coastal_edges):
            rotated_coastal_edges = con.coastal_edges[i:] + con.coastal_edges[:i]
            e.purple_ends.sort(key = lambda e_other:rotated_coastal_edges.index(e_other), reverse = True)
            e.green_ends.sort(key = lambda e_other:rotated_coastal_edges.index(e_other), reverse = True)
            if e.is_red:
                e.ends = e.green_ends + e.purple_ends
            else:
                e.ends = e.purple_ends + e.green_ends
        if foliation_style_split:
            for e1, e2 in purple_train_routes:
                p1 = end_pos(e2, e1)
                p2 = end_pos(e1, e2)
                p = make_arc(p1, p2)
                p = p.transformed(scl)
                canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, purple])
            for e1, e2 in green_train_routes:
                p1 = end_pos(e2, e1)
                p2 = end_pos(e1, e2)
                p = make_arc(p1, p2)
                p = p.transformed(scl)
                canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, green])
        if foliation_style_cusp_leaves or foliation_style_boundary_leaves:
            for i, c in enumerate(con.coast):
                c.purple_thorn_end_positions = [] ### complex numbers
                c.purple_thorn_ends = [] ### [coastal arc, position along that arc]
                e = con.coastal_edges[i]
                e1 = e.purple_ends[0]
                while True:
                    index = e1.purple_ends.index(e)
                    if index == len(e1.purple_ends) - 1:
                        break
                    else:
                        c.purple_thorn_end_positions.append( end_pos(e, e1, offset = 0.5) )
                        c.purple_thorn_ends.append( (e1, e1.ends.index(e)) )
                        e, e1 = e1, e1.purple_ends[index + 1]

                if foliation_style_boundary_leaves:
                    e_before = con.coastal_edges[(i-1)%len(con.coast)]
                    e_after = con.coastal_edges[i]
                    first_pos = end_pos(e_after.purple_ends[0], e_after, offset = -0.25)
                    last_pos = end_pos(e_before.purple_ends[-1], e_before, offset = 0.25)
                    c.purple_thorn_end_positions = [first_pos] + c.purple_thorn_end_positions + [last_pos]
                    arcs = []
                    for i in range(len(c.purple_thorn_end_positions) - 1):
                        arcs.append(make_arc(c.purple_thorn_end_positions[i], c.purple_thorn_end_positions[i+1]))
                    for p in arcs:
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, purple])

                if foliation_style_cusp_leaves:
                    for thorn_end in c.purple_thorn_ends:
                        thorn_end_pos = end_pos2(thorn_end)
                        p = make_arc(c.circle_pos, thorn_end_pos)
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, purple])
                
            for i, c in enumerate(con.coast):
                c.green_thorn_end_positions = [] ### complex numbers
                c.green_thorn_ends = [] ### [coastal arc, position along that arc]
                e = con.coastal_edges[i]
                e1 = e.green_ends[0]
                while True:
                    index = e1.green_ends.index(e)
                    if index == len(e1.green_ends) - 1:
                        break
                    else:
                        c.green_thorn_end_positions.append( end_pos(e, e1, offset = 0.5) )
                        c.green_thorn_ends.append( (e1, e1.ends.index(e)) )
                        e, e1 = e1, e1.green_ends[index + 1]
                if foliation_style_boundary_leaves:
                    e_before = con.coastal_edges[(i-1)%len(con.coast)]
                    e_after = con.coastal_edges[i]
                    first_pos = end_pos(e_after.green_ends[0], e_after, offset = -0.25)
                    last_pos = end_pos(e_before.green_ends[-1], e_before, offset = 0.25)
                    c.green_thorn_end_positions = [first_pos] + green_thorn_end_positions + [last_pos]
                    arcs = []
                    for i in range(len(c.green_thorn_end_positions) - 1):
                        arcs.append(make_arc(c.green_thorn_end_positions[i], c.green_thorn_end_positions[i+1]))
                    for p in arcs:
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, green])

                if foliation_style_cusp_leaves:
                    for thorn_end in c.green_thorn_ends:
                        thorn_end_pos = end_pos2(thorn_end)
                        p = make_arc(c.circle_pos, thorn_end_pos)
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, green])

            for tet in draw_tetrahedron_rectangles:
                a, c = tet.upper_edge().vertices
                b, d = tet.lower_edge().vertices
                if not are_anticlockwise(a.coastal_index, b.coastal_index, c.coastal_index):
                    b, d = d, b   ### now a, b, c, d are anticlockwise
                # print(con.vertices.index(a), con.vertices.index(b), con.vertices.index(c), con.vertices.index(d))
                all_sides_discrete = [tet_rectangle_sides(tet, v) for v in [a,b,c,d]]
                all_sides_geometry = []
                for vertex_sides_discrete in all_sides_discrete:
                    vertex_sides_geometry = []
                    for side_discrete in vertex_sides_discrete:
                        a, b = side_discrete
                        if side_discrete != None:
                            v, thorn_end = side_discrete
                            # print(thorn_end)
                            # print((v.circle_pos, end_pos2(thorn_end)))
                            vertex_sides_geometry.append( (v.circle_pos, end_pos2(thorn_end)) )
                        else:
                            vertex_sides_geometry.append( None )
                    all_sides_geometry.append(vertex_sides_geometry)
                for i in range(4):
                    if all_sides_geometry[i][0] != None and all_sides_geometry[(i+1)%4][1] != None:
                        v1, t1 = all_sides_geometry[i][0]
                        v2, t2 = all_sides_geometry[(i+1)%4][1]
                        intersection = geodesic_isect(v1, t1, v2, t2)
                        assert intersection != None
                        all_sides_geometry[i][0] = (v1, intersection)
                        all_sides_geometry[(i+1)%4][1] = (v2, intersection)

                for i, vertex_sides_geometry in enumerate(all_sides_geometry):
                    if i%2 == 0:
                        col = purple
                    else:
                        col = green
                    for side_geometry in vertex_sides_geometry: 
                        if side_geometry != None:
                            # print(side_geometry)
                            v_pos, t_pos = side_geometry 
                            p = make_arc(v_pos, t_pos)
                            p = p.transformed(scl)
                            canv.stroke(p, [style.linewidth(3*leaf_thickness), style.linecap.round, col])

    output_filename = 'Images/CircleContinent/' + name + '.pdf'
    canv.writePDFfile(output_filename)

def make_continent_naive(veering_isosig, max_num_tetrahedra = 50):
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle) #, tet_shapes = tet_shapes)
    initial_tet_face = tet_face(vt, 0, 0, verts_pos = [None, None, None, None])
    con = continent( vt, initial_tet_face) #, desired_vertices = desired_vertices )
    con.build_naive(max_num_tetrahedra = max_num_tetrahedra)
    con.make_convex()
    print(len(con.triangles), len(con.vertices), len(con.tetrahedra))
    return con

def make_continent_drill_dual_cycle(veering_isosig, dual_cycle, num_steps):
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle) #, tet_shapes = tet_shapes)
    ### initial tetrahedron is above face 0 of the dual cycle
    face0 = vt.tri.triangles()[dual_cycle[0]]
    embeds = face0.embeddings()
    tet_0, face_0 = None, None
    for embed in embeds:
        tet_num = embed.simplex().index()
        face_num = embed.face()
        if vt.coorientations[tet_num][face_num] == -1: ## into the tet
            tet_0, face_0 = tet_num, face_num
            break
    assert tet_0 != None and face_0 != None
    initial_tet_face = tet_face(vt, tet_0, face_0, verts_pos = [None, None, None, None])
    con = continent( vt, initial_tet_face)

    ### identify the triangle corresponding to dual_cycle[1]
    for triangle in con.triangles:
        if not triangle.is_upper and triangle.index == dual_cycle[0]:
            lowest_triangle = triangle
        if triangle.is_upper and triangle.index == dual_cycle[1]:
            highest_triangle = triangle
    triangle_path = [lowest_triangle, highest_triangle]
    lowest_path_index = 0
    highest_path_index = 1
    # for i in range(num_steps):
    #     last_added_tet = con.bury(highest_triangle)

    #     ### find next_triangle
    #     highest_triangle = None
    #     for triangle in last_added_tet.upper_triangles:
    #         if triangle.index == dual_cycle[(path_index + 1) % len(dual_cycle)]:
    #             highest_triangle = triangle
    #             break
    #     assert highest_triangle != None
    #     triangle_path.append(highest_triangle)
    #     path_index = path_index + 1

    #     con.make_convex() ### could this take us further up dual_cycle?

    for i in range(num_steps):
        if i%2 == 0:
            last_added_tet = con.bury(triangle_path[-1]) ## go up
            for triangle in last_added_tet.upper_triangles:
                if triangle.index == dual_cycle[(highest_path_index + 1) % len(dual_cycle)]:
                    triangle_path.append(triangle)
                    break
            highest_path_index = highest_path_index + 1
        else:
            last_added_tet = con.bury(triangle_path[0]) ## go down
            for triangle in last_added_tet.lower_triangles:
                if triangle.index == dual_cycle[(lowest_path_index - 1) % len(dual_cycle)]:
                    triangle_path.insert(0, triangle)
                    break
            lowest_path_index = lowest_path_index - 1

        con.make_convex() ### could this take us further up dual_cycle?

    con.update_boundary()
    con.triangle_path = triangle_path
    return con

def flow_edge_in_continent(con_tet, edge_num):
    tet_vertices = con_tet.ordered_vertices()
    vert_pair = edge_num_to_vert_pair[edge_num]
    vertex_pair = [tet_vertices[i] for i in vert_pair]
    upper_tris = con_tet.upper_triangles
    edge_candidates = set(upper_tris[0].edges).union(set(upper_tris[1].edges))
    edge = None
    for e in edge_candidates:
        if set(e.vertices) == set(vertex_pair): 
            edge = e
            break
    assert edge != None
    return edge

def make_continent_drill_flow_cycle(veering_isosig, flow_cycle, num_steps):
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle) #, tet_shapes = tet_shapes)
    ### format for loops: it is a list of tuples, 
    ### each tuple is (tet_index, edge_index within this tet that we exit through)
    tet_0 = flow_cycle[0][0]
    face_0 = 0 ### only used for initial numbering of vertices of the continent, so who cares
    initial_tet_face = tet_face(vt, tet_0, face_0, verts_pos = [None, None, None, None])
    con = continent( vt, initial_tet_face)

    side_face_collections, side_tet_collections = edge_side_face_collections(vt.tri, vt.angle, tet_vert_coorientations = vt.coorientations, return_tets = True, order_bottom_to_top = False)
    # print('sfc', side_face_collections)
    # print('stc', side_tet_collections)
    ### identify the next edge in the cycle 

    init_tetrahedron = con.tetrahedra[0]
    init_edge = flow_edge_in_continent(init_tetrahedron, flow_cycle[0][1])

    flow_tetrahedra = [init_tetrahedron]
    flow_edges = [init_edge]
    ### both in the continent

    upwards_flow_index = 0
    downwards_flow_index = 0
    for i in range(num_steps):
        # print(i)
        if i%2 == 0: ### go up
            edge = flow_edges[-1]
            last_tet = None
            while True:
                con.update_boundary()  
                upper_boundary_triangles = [t for t in edge.boundary_triangles if t.is_upper] 
                if len(upper_boundary_triangles) == 0:
                    break ## out of while loop
                last_tet = con.bury(upper_boundary_triangles[0])
            assert last_tet != None
            upwards_flow_index = (upwards_flow_index + 1) % len(flow_cycle)
            assert last_tet.index == flow_cycle[upwards_flow_index][0]
            flow_tetrahedra.append(last_tet)
            flow_edges.append(flow_edge_in_continent(last_tet, flow_cycle[upwards_flow_index][1]))
            con.make_convex() ### could this take us further up dual_cycle? We don't think so
        else: ### go down
            tet = flow_tetrahedra[0]
            edge = tet.lower_edge()
            flow_edges = [edge] + flow_edges
            ### now build the continent to get new_lowest_tet
            manifold_edge = tri.edge(edge.index)
            downwards_flow_index = (downwards_flow_index - 1) % len(flow_cycle)
            ### is the flow cycle vertical through the tetrahedron? 
            tet_below = get_tet_above_edge(vt.tri, vt.angle, manifold_edge, tet_vert_coorientations = vt.coorientations, get_tet_below_edge = True)
            tet_below_num = tet_below.index()
            top_vert_nums = get_tet_top_vert_nums(vt.coorientations, tet_below_num)
            top_edge_num = vert_pair_to_edge_num[tuple(top_vert_nums)]

            if (tet_below_num, top_edge_num) == flow_cycle[downwards_flow_index]: ### flow cycle went straight up
                while True:
                    con.update_boundary()  
                    lower_boundary_triangles = [t for t in edge.boundary_triangles if not t.is_upper] 
                    if len(lower_boundary_triangles) == 0:
                        break ## out of while loop
                    last_tet = con.bury(lower_boundary_triangles[0])
            else:
                ### find which side of the edge our tet is in
                # print('edge index', edge.index)
                side_tet_collections_at_edge = side_tet_collections[edge.index] ## index in the manifold
                side_face_collections_at_edge = side_face_collections[edge.index]
                downward_path = None
                flow_step = flow_cycle[downwards_flow_index]
                for i, side_tet_collection in enumerate(side_tet_collections_at_edge):
                    if flow_step in side_tet_collection:
                        downward_path = side_tet_collection[:side_tet_collection.index(flow_step) + 1]
                        downward_path_faces = side_face_collections_at_edge[i][:side_tet_collection.index(flow_step) + 1]
                assert downward_path != None
                for j, (tet_num, edge_num) in enumerate(downward_path):
                    con.update_boundary()  
                    lower_boundary_triangles = [t for t in edge.boundary_triangles if not t.is_upper and t.index == downward_path_faces[j][0]] 
                    assert len(lower_boundary_triangles) == 1
                    last_tet = con.bury(lower_boundary_triangles[0])
            assert last_tet != None
            flow_tetrahedra = [last_tet] + flow_tetrahedra
            con.make_convex() ### could this take us further down dual_cycle? We don't think so

    con.update_boundary()
    return con, flow_tetrahedra, flow_edges

def complete_tetrahedron_rectangles(con, tetrahedra_to_complete):
    """grow the continent so that the given tetrahedra have full tetrahedron rectangles within the continent"""
    k = 0
    for tet in tetrahedra_to_complete:
        for v in tet.vertices():
            # print('tet vert age', con.vertices.index(v))
            con.update_boundary()
            install_thorn_ends(con)
            sides = tet_rectangle_sides(tet, v)
            for direction in range(2):
                while sides[direction] == None:
                    # print('direction, k', direction, k)
                    e = con.coastal_edges[(v.coastal_index - direction)%len(con.coast)] 
                    triangles = e.boundary_triangles  ### grow around this edge
                    if triangles[0].is_upper != (k % 2 == 0): ### xor, alternate which one we add to
                        con.bury(triangles[0])
                    else:
                        con.bury(triangles[1])
                    con.make_convex()
                    con.update_boundary()
                    install_thorn_ends(con)
                    sides = tet_rectangle_sides(tet, v)
                    k += 1
                    if k > 50:
                        print('bail')
                        return None

def get_fund_domain_tetrahedra(con):
    num_tet = con.vt.tri.countTetrahedra()
    fund_dom_tets = list(range(num_tet))
    number_to_find = num_tet
    for con_tet in con.tetrahedra:
        if con_tet.index in fund_dom_tets:
            fund_dom_tets[con_tet.index] = con_tet  ### replace index with the continent tetrahedron
            number_to_find -= 1
        if number_to_find == 0:
            break
    if number_to_find > 0: ## should have found the whole continent
        print('did not find all of the downstairs tetrahedra!')
    update_fund_dom_tet_nums(con, fund_dom_tets)
    return fund_dom_tets

def update_fund_dom_tet_nums(con, fund_dom_tets):
    for v in con.coast:
        v.fund_dom_tet_nums = []
    for con_tet in fund_dom_tets:
        if type(con_tet) == continent_tetrahedron:  ### could be an integer if we didnt find this tet
            for v in con_tet.vertices():
                v.fund_dom_tet_nums.append(con_tet.index)

def main():
    veering_isosig = 'cPcbbbdxm_10' 
    flow_cycle = [(0, 2)]

    # veering_isosig = 'eLAkaccddjsnak_2001'
    # flow_cycle = [(1, 0), (2, 5)]

    # for num_steps in range(10):
    num_steps = 5
    con, flow_tetrahedra, flow_edges = make_continent_drill_flow_cycle(veering_isosig, flow_cycle, num_steps)
    fund_dom_tets = get_fund_domain_tetrahedra(con)
    complete_tetrahedron_rectangles(con, fund_dom_tets)
    print(len(flow_tetrahedra))
    name = veering_isosig + '_' + str(flow_cycle) + '_' + str(num_steps) + '_cusp_leaves'
    # tets_to_draw = [flow_tetrahedra[0], flow_tetrahedra[-1]]
    tets_to_draw = fund_dom_tets
    draw_continent_circle(con, name = name,
        draw_upper_landscape = False, draw_lower_landscape = False, 
        draw_upper_green = True, draw_lower_purple = True,
        draw_train_tracks = False, draw_foliation = True, 
        foliation_style_old = False, foliation_style_split = False, 
        foliation_style_cusp_leaves = True, foliation_style_boundary_leaves = False,
        shade_triangles = False, draw_fund_domain = True, fund_dom_tets = fund_dom_tets,
        draw_fund_domain_edges = True, draw_tetrahedron_rectangles = tets_to_draw)



    # veering_isosig = 'cPcbbbiht_12'
    # veering_isosig = 'dLQacccjsnk_200'
    # max_num_tetrahedra = 50
    # con = make_continent_naive(veering_isosig, max_num_tetrahedra = max_num_tetrahedra)
    # name = veering_isosig + '_' + str(max_num_tetrahedra) + '_cusp_leaves'
    # draw_continent_circle(con, name = name,
    #     draw_upper_landscape = False, draw_lower_landscape = False, 
    #     draw_upper_green = True, draw_lower_purple = True,
    #     draw_train_tracks = False, draw_foliation = True,
    #     foliation_style_old = False, foliation_style_split = False, 
    #     foliation_style_cusp_leaves = True, foliation_style_boundary_leaves = False,
    #     draw_fund_domain = True)

#    veering_isosig = 'dLQacccjsnk_200'
#    dual_cycle = [4,5]
#    for num_steps in range(20):
#    # num_steps = 7
#        con = make_continent_drill(veering_isosig, dual_cycle, num_steps)
#        name = veering_isosig + '_' + str(dual_cycle) + '_' + str(num_steps) + '_cusp_leaves'
#        draw_continent_circle(con, name = name,
#            draw_upper_landscape = False, draw_lower_landscape = False, 
#            draw_upper_green = True, draw_lower_purple = True,
#            draw_train_tracks = False, draw_foliation = True, 
#            foliation_style_old = False, foliation_style_split = False, 
#            foliation_style_cusp_leaves = True, foliation_style_boundary_leaves = False,
#            shade_triangles = True, draw_fund_domain = True,
#            draw_fund_domain_edges = True)
        
    # veering_isosig = 'cPcbbbdxm_10' # example where after drilling we get an edge between the new cusp and itself
    # dual_cycle = [1,2]
    # for num_steps in range(20):
    # # num_steps = 7
    #     con = make_continent_drill_dual_cycle(veering_isosig, dual_cycle, num_steps)
    #     name = veering_isosig + '_' + str(dual_cycle) + '_' + str(num_steps) + '_cusp_leaves'
    #     draw_continent_circle(con, name = name,
    #         draw_upper_landscape = False, draw_lower_landscape = False, 
    #         draw_upper_green = True, draw_lower_purple = True,
    #         draw_train_tracks = False, draw_foliation = True, 
    #         foliation_style_old = False, foliation_style_split = False, 
    #         foliation_style_cusp_leaves = True, foliation_style_boundary_leaves = False,
    #         shade_triangles = True, draw_fund_domain = True,
    #         draw_fund_domain_edges = True)




