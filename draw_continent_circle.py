
import pyx ### vector graphics 
from math import pi, cos, sin, acos, tan, atan2
from taut import isosig_to_tri_angle
from veering import veering_triangulation
from continent import continent
from boundary_triangulation import tet_face

def make_line(a, b):
    p = pyx.path.path( pyx.path.moveto(a[0], a[1]) )
    p.append( pyx.path.lineto(b[0], b[1]) )
    return p

def make_arc(a, b, return_midpt = False):
    det = 4*(a[0]*b[1] - a[1]*b[0])
    if abs(det) < 0.0001: 
        if not return_midpt:
            return make_line(a, b)
        else:
            return (make_line(a, b), (0.0,0.0))
    else:
        # theta = 0.5 * acos( a[0]*b[0] + a[1]*b[1] )
        # r = tan(theta)
        # p = pyx.path.path( pyx.path.moveto(a[0], a[1]) )
        # p.append( pyx.path.arct(0.0,0.0,b[0], b[1], r))

        A = complex(a[0], a[1])
        B = complex(b[0], b[1])
        C = B - A
        D = B * (a[0]*a[0] + a[1]*a[1])
        E = A * (b[0]*b[0] + b[1]*b[1])
        center = (2.0/det) * complex(C.imag + D.imag - E.imag, -C.real - D.real + E.real)
        A2 = A - center
        B2 = B - center
        r = abs(A2)
        dot_prod = (A2.real*B2.real + A2.imag*B2.imag)
        theta = 0.5 * acos( dot_prod/(abs(A2)*abs(B2)) )
        direction = (0.5*(A + B) - center)
        direction = (1.0/abs(direction)) * direction
        tangent_point = center + (abs(A2)/cos(theta)) * direction
        
        p = pyx.path.path( pyx.path.moveto(a[0], a[1]) )
        p.append( pyx.path.arct(tangent_point.real, tangent_point.imag, b[0], b[1], r))

        # p = pyx.path.path( pyx.path.moveto(a[0], a[1]) )
        # p.append( pyx.path.lineto(center.real, center.imag) )
        # p.append( pyx.path.lineto(b[0], b[1]) )
        if not return_midpt:
            return p
        else:
            midpt = center + r * direction
            return (p, (midpt.real, midpt.imag))

def draw_continent_circle(veering_isosig, max_num_tetrahedra = 50, draw_upper_landscape = True, draw_lower_landscape = False, draw_upper_green = True, draw_lower_purple = False):
    global_scale_up = 10.0
    scl = pyx.trafo.trafo(matrix=((global_scale_up, 0), (0, global_scale_up)), vector=(0, 0))
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle) #, tet_shapes = tet_shapes)

    canv = pyx.canvas.canvas()
    canv.stroke(pyx.path.circle(0,0,global_scale_up), [pyx.style.linewidth(0.02)])

    initial_tet_face = tet_face(vt, 0, 0, verts_pos = [None, None, None, None])
    con = continent( vt, initial_tet_face) #, desired_vertices = desired_vertices )
    con.build_naive(max_num_tetrahedra = max_num_tetrahedra)
    print(len(con.triangles), len(con.vertices), len(con.tetrahedra))

    n = len(con.coast)
    for i,v in enumerate(con.coast):
        v.coastal_index = i
        t = 2*pi*float(i)/float(n)
        v.circle_pos = (cos(t), sin(t))

    lower_colours = {True: pyx.color.rgb(0.5,0.3,0), False: pyx.color.rgb(0,0.3,0.5)}
    upper_colours = {True: pyx.color.rgb(0.9,0.3,0), False: pyx.color.rgb(0,0.3,0.9)}
    green = pyx.color.rgb(0.0,0.5,0.0)
    purple = pyx.color.rgb(0.5,0.0,0.5)

    landscape_edges = con.boundary_landscape_edges()

    colours = [lower_colours, upper_colours]

    upper_tris = []
    lower_tris = []
    for tri in con.triangles:
        if not tri.is_buried:
            if tri.is_upper:
                upper_tris.append(tri)
            else:
                lower_tris.append(tri)
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
            canv.stroke(p, [pyx.style.linewidth(0.001), col])
        for tri in boundary_tris[i]:
            midpts = []
            for e in tri.edges:
                u, v = e.vertices
                p, midpt = make_arc(u.circle_pos, v.circle_pos, return_midpt = True)
                midpts.append(midpt)
            center = ( (1.0/3.0)*(midpts[0][0] + midpts[1][0] + midpts[2][0]), (1.0/3.0)*(midpts[0][1] + midpts[1][1] + midpts[2][1]) )
            canv.text(global_scale_up*center[0], global_scale_up*center[1], "$"+str(tri.index)+"$", textattrs=[pyx.text.halign.center, pyx.text.vshift.middlezero, pyx.text.size(0)])


    ### train tracks...

    if draw_lower_purple:
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
                    canv.stroke(p, [pyx.style.linewidth(0.004), purple])       
    if draw_upper_green:
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
                    canv.stroke(p, [pyx.style.linewidth(0.004), green])


    #     midpts.append(midpt)
    #     for i in range(3):
    #         draw_line(canv, midpts[i], midpts[(i+1)%3], [pyx.style.linewidth(0.004)])


    output_filename = 'Images/CircleContinent/' + veering_isosig + '.pdf'
    canv.writePDFfile(output_filename)


def main():
    # veering_isosig = 'cPcbbbiht_12'
    veering_isosig = 'dLQacccjsnk_200'
    draw_continent_circle(veering_isosig, max_num_tetrahedra = 20, 
        draw_upper_landscape = True, draw_lower_landscape = False, 
        draw_upper_green = True, draw_lower_purple = False)




