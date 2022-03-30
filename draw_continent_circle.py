
import pyx ### vector graphics 
from pyx import path, trafo, canvas, style, text, color, deco
from math import pi, cos, sin, acos, tan, atan2, sqrt
from taut import isosig_to_tri_angle
from veering import veering_triangulation
from continent import continent
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

def draw_continent_circle(con, name = "", draw_upper_landscape = True, draw_lower_landscape = False, 
    draw_upper_green = True, draw_lower_purple = False, draw_train_tracks = False, 
    draw_foliation = True, foliation_style_old = False, foliation_style_split = False, 
    foliation_style_cusp_leaves = True, foliation_style_boundary_leaves = True,
    shade_triangles = False):
    
    global_scale_up = 10.0
    edge_thickness = 0.01
    track_thickness = 0.02
    leaf_thickness = 0.03

    scl = trafo.trafo(matrix=((global_scale_up, 0), (0, global_scale_up)), vector=(0, 0))
    canv = canvas.canvas()
    canv.stroke(path.circle(0,0,global_scale_up), [style.linewidth(0.02)])

    n = len(con.coast)
    for i,v in enumerate(con.coast):
        v.coastal_index = i
        t = 2*pi*float(i)/float(n)
        v.circle_pos = complex(cos(t), sin(t))
        vert_pos = v.circle_pos * 1.01 * global_scale_up
        canv.text(vert_pos.real, vert_pos.imag, "$"+str(con.vertices.index(v))+"$", textattrs=[text.size(-4), text.halign.left, text.valign.middle,trafo.rotate((180/pi) * atan2(vert_pos.imag, vert_pos.real))])

        # vert_pos2 = v.circle_pos * 1.2 * global_scale_up
        # p = path.path(path.moveto(vert_pos.real, vert_pos.imag), path.lineto(vert_pos2.real, vert_pos2.imag))
        # canv.stroke(p, [deco.curvedtext("$"+str(con.vertices.index(v))+"$")])

    # lower_colours = {True: color.rgb(0.5,0.3,0), False: color.rgb(0,0.3,0.5)}
    # upper_colours = {True: color.rgb(0.9,0.3,0), False: color.rgb(0,0.3,0.9)}
    edge_colours = {True: color.rgb(0.9,0.3,0), False: color.rgb(0,0.3,0.9)}
    green = color.rgb(0.0,0.5,0.0)
    purple = color.rgb(0.5,0.0,0.5)

    landscape_edges = [con.lower_landscape_edges, con.upper_landscape_edges]

    # colours = [lower_colours, upper_colours]

    upper_tris = con.upper_landscape_triangles
    lower_tris = con.lower_landscape_triangles
    boundary_tris = [lower_tris, upper_tris]

    if shade_triangles:
        u,v,w = con.lowest_triangle.vertices
        p = make_arc(u.circle_pos, v.circle_pos)
        q = make_arc(v.circle_pos, w.circle_pos)
        r = make_arc(w.circle_pos, u.circle_pos)
        p.append(q[1])
        p.append(r[1])  ### remove extraneous moveto commands
        p = p.transformed(scl)
        canv.stroke(p, [deco.filled([color.transparency(0.8)]), style.linewidth(0)])
        u,v,w = con.highest_triangle.vertices
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

    purple_leaves = []
    green_leaves = []
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
                    purple_leaves.append(leaf_end_edges)
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
                    green_leaves.append(leaf_end_edges)
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
        for e1, e2 in purple_leaves:
            e1.purple_ends.append(e2)
            e2.purple_ends.append(e1)
        for e1, e2 in green_leaves:
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
            for e1, e2 in purple_leaves:
                p1 = end_pos(e2, e1)
                p2 = end_pos(e1, e2)
                p = make_arc(p1, p2)
                p = p.transformed(scl)
                canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, purple])
            for e1, e2 in green_leaves:
                p1 = end_pos(e2, e1)
                p2 = end_pos(e1, e2)
                p = make_arc(p1, p2)
                p = p.transformed(scl)
                canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, green])
        if foliation_style_cusp_leaves or foliation_style_boundary_leaves:
            for i, c in enumerate(con.coast):
                purple_thorn_ends = []
                e = con.coastal_edges[i]
                e1 = e.purple_ends[0]
                while True:
                    index = e1.purple_ends.index(e)
                    if index == len(e1.purple_ends) - 1:
                        break
                    else:
                        purple_thorn_ends.append( end_pos(e, e1, offset = 0.5) )
                        e, e1 = e1, e1.purple_ends[index + 1]
                if foliation_style_cusp_leaves:
                    for thorn_end in purple_thorn_ends:
                        p = make_arc(c.circle_pos, thorn_end)
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, purple])
                if foliation_style_boundary_leaves:
                    e_before = con.coastal_edges[(i-1)%len(con.coast)]
                    e_after = con.coastal_edges[i]
                    first_pos = end_pos(e_after.purple_ends[0], e_after, offset = -0.25)
                    last_pos = end_pos(e_before.purple_ends[-1], e_before, offset = 0.25)
                    purple_thorn_ends = [first_pos] + purple_thorn_ends + [last_pos]
                    arcs = []
                    for i in range(len(purple_thorn_ends) - 1):
                        arcs.append(make_arc(purple_thorn_ends[i], purple_thorn_ends[i+1]))
                    for p in arcs:
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, purple])

            for i, c in enumerate(con.coast):
                green_thorn_ends = []
                e = con.coastal_edges[i]
                e1 = e.green_ends[0]
                while True:
                    index = e1.green_ends.index(e)
                    if index == len(e1.green_ends) - 1:
                        break
                    else:
                        green_thorn_ends.append( end_pos(e, e1, offset = 0.5) )
                        e, e1 = e1, e1.green_ends[index + 1]
                if foliation_style_cusp_leaves:
                    for thorn_end in green_thorn_ends:
                        p = make_arc(c.circle_pos, thorn_end)
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, green])
                if foliation_style_boundary_leaves:
                    e_before = con.coastal_edges[(i-1)%len(con.coast)]
                    e_after = con.coastal_edges[i]
                    first_pos = end_pos(e_after.green_ends[0], e_after, offset = -0.25)
                    last_pos = end_pos(e_before.green_ends[-1], e_before, offset = 0.25)
                    green_thorn_ends = [first_pos] + green_thorn_ends + [last_pos]
                    arcs = []
                    for i in range(len(green_thorn_ends) - 1):
                        arcs.append(make_arc(green_thorn_ends[i], green_thorn_ends[i+1]))
                    for p in arcs:
                        p = p.transformed(scl)
                        canv.stroke(p, [style.linewidth(leaf_thickness), style.linecap.round, green])


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

def make_continent_drill(veering_isosig, dual_cycle, num_steps):
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
            next_triangle = triangle
    next_triangles = []
    path_index = 1
    for i in range(num_steps):
        last_added_tet = con.bury(next_triangle)

        ### find next_triangle
        next_triangle = None
        for triangle in last_added_tet.upper_triangles:
            if triangle.index == dual_cycle[(path_index + 1) % len(dual_cycle)]:
                next_triangle = triangle
                break
        assert next_triangle != None
        next_triangles.append(next_triangle)
        path_index = path_index + 1

        con.make_convex() ### could this take us further up dual_cycle?

    con.update_boundary()
    con.lowest_triangle = lowest_triangle
    con.triangle_path = next_triangles
    con.highest_triangle = next_triangle
    return con


def main():
    # veering_isosig = 'cPcbbbiht_12'
    veering_isosig = 'dLQacccjsnk_200'
    max_num_tetrahedra = 50
    con = make_continent_naive(veering_isosig, max_num_tetrahedra = max_num_tetrahedra)
    name = veering_isosig + '_' + str(max_num_tetrahedra) + '_boundary_leaves'
    draw_continent_circle(con, name = name,
        draw_upper_landscape = False, draw_lower_landscape = False, 
        draw_upper_green = True, draw_lower_purple = True,
        draw_train_tracks = False, draw_foliation = True,
        foliation_style_old = False, foliation_style_split = False, 
        foliation_style_cusp_leaves = False, foliation_style_boundary_leaves = True)

    # veering_isosig = 'dLQacccjsnk_200'
    # dual_cycle = [4,5]
    # for num_steps in range(20):
    # # num_steps = 7
    #     con = make_continent_drill(veering_isosig, dual_cycle, num_steps)
    #     name = veering_isosig + '_' + str(dual_cycle) + '_' + str(num_steps) + '_cusp_leaves'
    #     draw_continent_circle(con, name = name,
    #         draw_upper_landscape = True, draw_lower_landscape = True, 
    #         draw_upper_green = True, draw_lower_purple = True,
    #         draw_train_tracks = False, draw_foliation = True, foliation_style_old = False,
    #         foliation_style_split = False, foliation_style_cusp_leaves = True,
    #         shade_triangles = True)




