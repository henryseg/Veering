
import pyx ### vector graphics 
from math import pi, cos, sin, acos, tan, atan2, sqrt
from taut import isosig_to_tri_angle
from veering import veering_triangulation
from continent import continent
from boundary_triangulation import tet_face

def to_complex(a):
    return complex(a[0], a[1])

def make_line(a, b):
    p = pyx.path.path( pyx.path.moveto(a.real, a.imag) )
    p.append( pyx.path.lineto(b.real, b.imag) )
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
        # p = pyx.path.path( pyx.path.moveto(a[0], a[1]) )
        # p.append( pyx.path.arct(0.0,0.0,b[0], b[1], r))

        center, r = arc_center_rad(a, b)
        A2 = a - center
        B2 = b - center
        r = abs(A2)
        dot_prod = (A2.real*B2.real + A2.imag*B2.imag)
        theta = 0.5 * acos( dot_prod/(abs(A2)*abs(B2)) )
        direction = (0.5*(a + b) - center)
        direction = (1.0/abs(direction)) * direction
        tangent_point = center + (abs(A2)/cos(theta)) * direction
        
        p = pyx.path.path( pyx.path.moveto(a.real, a.imag) )
        p.append( pyx.path.arct(tangent_point.real, tangent_point.imag, b.real, b.imag, r))

        # p = pyx.path.path( pyx.path.moveto(a[0], a[1]) )
        # p.append( pyx.path.lineto(center.real, center.imag) )
        # p.append( pyx.path.lineto(b[0], b[1]) )
        if not return_midpt:
            return p
        else:
            midpt = center + r * direction
            return (p, midpt)



def draw_continent_circle(veering_isosig, max_num_tetrahedra = 50, draw_upper_landscape = True, draw_lower_landscape = False, draw_upper_green = True, draw_lower_purple = False, draw_train_tracks = False, draw_foliation = True):
    global_scale_up = 10.0
    edge_thickness = 0.01
    track_thickness = 0.02
    leaf_thickness = 0.03

    scl = pyx.trafo.trafo(matrix=((global_scale_up, 0), (0, global_scale_up)), vector=(0, 0))
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle) #, tet_shapes = tet_shapes)

    canv = pyx.canvas.canvas()
    canv.stroke(pyx.path.circle(0,0,global_scale_up), [pyx.style.linewidth(0.02)])

    initial_tet_face = tet_face(vt, 0, 0, verts_pos = [None, None, None, None])
    con = continent( vt, initial_tet_face) #, desired_vertices = desired_vertices )
    con.build_naive(max_num_tetrahedra = max_num_tetrahedra)
    con.make_convex()
    print(len(con.triangles), len(con.vertices), len(con.tetrahedra))

    n = len(con.coast)
    for i,v in enumerate(con.coast):
        v.coastal_index = i
        t = 2*pi*float(i)/float(n)
        v.circle_pos = complex(cos(t), sin(t))

    lower_colours = {True: pyx.color.rgb(0.5,0.3,0), False: pyx.color.rgb(0,0.3,0.5)}
    upper_colours = {True: pyx.color.rgb(0.9,0.3,0), False: pyx.color.rgb(0,0.3,0.9)}
    green = pyx.color.rgb(0.0,0.5,0.0)
    purple = pyx.color.rgb(0.5,0.0,0.5)

    landscape_edges = [con.lower_landscape_edges, con.upper_landscape_edges]

    colours = [lower_colours, upper_colours]

    upper_tris = con.upper_landscape_triangles
    lower_tris = con.lower_landscape_triangles
    boundary_tris = [lower_tris, upper_tris]

    to_do = []
    if draw_lower_landscape:
        to_do.append(0)
    if draw_upper_landscape:
        to_do.append(1)
    for i in to_do:
        for e in landscape_edges[i]:
            col = colours[i][e.is_red]
            u, v = e.vertices
            p = make_arc(u.circle_pos, v.circle_pos)
            p = p.transformed(scl)
            canv.stroke(p, [pyx.style.linewidth(edge_thickness), pyx.style.linecap.round, col])
        for tri in boundary_tris[i]:
            center = incenter(tri.vertices[0].circle_pos, tri.vertices[1].circle_pos, tri.vertices[2].circle_pos)
            # canv.fill(pyx.path.circle(global_scale_up*center[0], global_scale_up*center[1], 0.1)) ### text has bottom left corner here
            canv.text(global_scale_up*center[0], global_scale_up*center[1], "$"+str(tri.index)+"$", textattrs=[pyx.text.size(0)])


    ### train tracks...

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
                        canv.stroke(p, [pyx.style.linewidth(track_thickness), pyx.style.linecap.round, purple])    
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
                    leaf_ends = []
                    for e in leaf_end_edges:
                        endpts = e.vertices
                        _, midpt = make_arc(endpts[0].circle_pos, endpts[1].circle_pos, return_midpt = True)
                        leaf_ends.append(midpt)
                    p = make_arc(leaf_ends[0], leaf_ends[1])
                    p = p.transformed(scl)
                    canv.stroke(p, [pyx.style.linewidth(leaf_thickness), pyx.style.linecap.round, purple])

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
                        canv.stroke(p, [pyx.style.linewidth(track_thickness), pyx.style.linecap.round, green])
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
                    leaf_ends = []
                    for e in leaf_end_edges:
                        endpts = e.vertices
                        _, midpt = make_arc(endpts[0].circle_pos, endpts[1].circle_pos, return_midpt = True)
                        leaf_ends.append(midpt)
                    p = make_arc(leaf_ends[0], leaf_ends[1])
                    p = p.transformed(scl)
                    canv.stroke(p, [pyx.style.linewidth(leaf_thickness), pyx.style.linecap.round, green])

    output_filename = 'Images/CircleContinent/' + veering_isosig + '.pdf'
    canv.writePDFfile(output_filename)


def main():
    # veering_isosig = 'cPcbbbiht_12'
    veering_isosig = 'dLQacccjsnk_200'
    draw_continent_circle(veering_isosig, max_num_tetrahedra = 500, 
        draw_upper_landscape = False, draw_lower_landscape = False, 
        draw_upper_green = True, draw_lower_purple = True,
        draw_train_tracks = False, draw_foliation = True)




