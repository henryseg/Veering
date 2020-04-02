
import pyx ### vector graphics 

from file_io import parse_data_file, read_from_pickle
from taut import isosig_to_tri_angle
from veering import veering_triangulation
from continent import continent
from boundary_triangulation import boundary_triangulation


# def pre_draw_transformation( z, ladder_holonomy ):
    # return z/ladder_holonomy

def draw_path(canv, path_C, draw_options):
    p = pyx.path.path( pyx.path.moveto(path_C[0].real, path_C[0].imag) )
    for coord in path_C[1:]:  ### this is how path drawing works...
        p.append( pyx.path.lineto(coord.real, coord.imag) )
    canv.stroke(p, draw_options)

def draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, output_filename = None, draw_args = None, lw = 0.005 ):

    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle, tet_shapes = tet_shapes)
    B = boundary_triangulation(vt)
    B.generate_canvases(args = draw_args)

    # print B.torus_triangulation_list[0].ladder_holonomy
    T = B.torus_triangulation_list[0]

    ### make initial_tet_face be in the lower left of the fundamental domain
    # initial_tet_face = T.ladder_list[0].ladder_unit_list[0]
    
    ### make initial_tet_face be in the middle of the fundamental domain
    num_ladders = len(T.ladder_list)
    L = T.ladder_list[int(num_ladders/2 - 1)]  ## -1 because we split the last ladder between the left and right
    num_ladder_units = len(L.ladder_unit_list)
    initial_tet_face = L.ladder_unit_list[int(num_ladder_units/2)]
    print 'initial_tet_face', initial_tet_face

    ladderpoles_vertices = T.left_ladder_pole_vertices() 
    desired_vertices = [v for L in ladderpoles_vertices for v in L]

    con = continent( vt, initial_tet_face, desired_vertices = desired_vertices )
    
    # con.build(until_have_desired_vertices = False, max_num_tetrahedra = max_num_tetrahedra)
    con.build_fundamental_domain()  ## expand the continent until we have all vertices of the boundary triangulation fundamental domain

    print 'unfound desired_vertices', con.desired_vertices
    print 'num_tetrahedra', con.num_tetrahedra

    # now replace ladderpoles_vertices with the continent's corresponding vertices 
    epsilon = 0.001
    for ladderpole_vertices in ladderpoles_vertices:
        for i, w in enumerate(ladderpole_vertices):
            for v in con.boundary_triangulation_vertices:
                if abs(v.pos.complex() - w) < epsilon:
                    ladderpole_vertices[i] = v    
                    break     

    interesting_segments = []
    for ladderpole_vertices in ladderpoles_vertices:
        segment = [con.coast.index(ladderpole_vertices[0]), con.coast.index(ladderpole_vertices[-1])]
        segment.sort()
        interesting_segments.append( segment )

    print interesting_segments

    con.mark_interesting_segments(interesting_segments)

    # con.build(max_interesting_edge_length = 0.2, max_num_tetrahedra = max_num_tetrahedra)

    # eq = con.segment_between( ladderpoles_vertices[0][0], ladderpoles_vertices[0][1] )   ## segment under one edge of ladderpole
    # eq = con.segment_between( ladderpoles_vertices[0][0], ladderpoles_vertices[0][-1] )   ## 0th ladderpole

    grad = pyx.color.gradient.Hue
    # colours = {'L':pyx.color.rgb.blue, 'R':pyx.color.rgb.red}  

    # draw_options = [pyx.style.linewidth(lw), colour])  ## needs to know which colour if we do this
    draw_options = [pyx.style.linewidth(lw), pyx.deco.colorgradient(grad)]

    canv = T.canv

    for ladderpole_vertices in ladderpoles_vertices:
        for i in range(len(ladderpole_vertices) - 1):  # fenceposts
            path = con.segment_between(ladderpole_vertices[i], ladderpole_vertices[i+1])  
            path_C = [ T.drawing_scale * v.pos.complex() for v in path ]
            draw_path(canv, path_C, draw_options)  

    # adj_verts, adj_edges = con.vertices_and_edges_adjacent_to_infinity()  

    ### continent drawing the boundary triangulation
    # lines for triangles meeting infinity
    
    # for endpoints, veering_colour in adj_edges:
    #     z, w = [T.drawing_scale * v.pos.complex() for v in endpoints]
    #     pyx_stroke_col = pyx.deco.stroked([colours[veering_colour]])
    #     canv.stroke(pyx.path.line( z.real, z.imag, w.real, w.imag),  [pyx.style.linewidth(lw * 5), pyx_stroke_col]  )

    # # dots for edges from infinity
    # for v,veering_colour in adj_verts:
    #     z = v.pos.complex()
    #     pyx_fill_col = pyx.deco.filled([colours[veering_colour]])
    #     canv.fill(pyx.path.circle(z.real, z.imag, 0.02), [pyx_fill_col])

    ### continent drawing the left_ladder_pole_vertices
    # T = B.torus_triangulation_list[0]
    # for L in T.ladder_list:
    #     for v in L.left_ladder_pole_vertices():
    #         v *= T.drawing_scale
    #         if L.is_even:
    #             col = colours['L']
    #         else:
    #             col = colours['R']
    #         canv.stroke(pyx.path.circle(v.real, v.imag, 0.15), [pyx.style.linewidth(lw * 3), col])

    ### circles around found vertices
    # for v in con.found_vertices:
    #     pos = v.pos.complex()
    #     pos *= T.drawing_scale
    #     canv.stroke(pyx.path.circle(pos.real, pos.imag, 0.2), [pyx.style.linewidth(lw * 3)])

    canv.writePDFfile(output_filename)

def draw_cannon_thurston_from_veering_isosigs_file(veering_isosigs_filename, output_dirname, max_num_tetrahedra = 500, num_to_draw = None):
    veering_isosigs_list = parse_data_file(veering_isosigs_filename)
    if num_to_draw != None:
        to_draw = veering_isosigs_list[:num_to_draw]
    else:
        to_draw = veering_isosigs_list

    shapes_data = read_from_pickle('Data/veering_shapes_up_to_ten_tetrahedra.pkl')
    for veering_isosig in to_draw:
        print veering_isosig
        tet_shapes = shapes_data[veering_isosig]
        filename = output_dirname + '/' + veering_isosig + '_' + str(max_num_tetrahedra) + '.pdf'
        draw_continent(veering_isosig, tet_shapes, max_num_tetrahedra, output_filename = filename)


if __name__ == '__main__':
    args = {'draw_boundary_triangulation':False, 'global_drawing_scale': 1.5, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps
    # veering_isosig = 'cPcbbbiht_12'
    veering_isosig = 'dLQacccjsnk_200'
    # veering_isosig = 'iLLLAQccdffgfhhhqgdatgqdm_21012210' ## no symmetry - helps us spot errors
    shapes_data = read_from_pickle('Data/veering_shapes_up_to_ten_tetrahedra.pkl')
    tet_shapes = shapes_data[veering_isosig]

    max_num_tetrahedra = 50000
    filename = 'Images/Cannon-Thurston/' + veering_isosig + '_' + str(max_num_tetrahedra) + '.pdf'
    draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, output_filename = filename, draw_args = args )

    # draw_cannon_thurston_from_veering_isosigs_file('Data/veering_census.txt', 'Images/Cannon-Thurston', max_num_tetrahedra = 10000, num_to_draw = 10)
    


