
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

def draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, draw_CT_curve = True, draw_lightning_curve = False, draw_landscapes = False, max_length = 0.1, output_filename = None, draw_args = None, build_type = None ):

    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle, tet_shapes = tet_shapes)
    B = boundary_triangulation(vt)
    B.generate_canvases(args = draw_args)

    for i,T in enumerate(B.torus_triangulation_list):
        print(('cusp', i))
        ### make initial_tet_face be in the lower left of the fundamental domain
        # initial_tet_face = T.ladder_list[0].ladder_unit_list[0]
        
        ### make initial_tet_face be in the middle of the fundamental domain
        num_ladders = len(T.ladder_list)
        L = T.ladder_list[int(num_ladders/2 - 1)]  ## -1 because we split the last ladder between the left and right
        num_ladder_units = len(L.ladder_unit_list)
        initial_tet_face = L.ladder_unit_list[int(num_ladder_units/2)]
        print(('initial_tet_face', initial_tet_face))

        ladderpoles_vertices = T.left_ladder_pole_vertices() 
        desired_vertices = [v for L in ladderpoles_vertices for v in L]

        con = continent( vt, initial_tet_face, desired_vertices = desired_vertices )
        
        con.build_fundamental_domain()  ## expand the continent until we have all vertices of the boundary triangulation fundamental domain

        print(('unfound desired_vertices', con.desired_vertices))
        print(('num_tetrahedra', con.num_tetrahedra))

        # now replace ladderpoles_vertices with the continent's corresponding vertices 
        epsilon = 0.001
        for ladderpole_vertices in ladderpoles_vertices:
            for i, w in enumerate(ladderpole_vertices):
                for v in con.boundary_triangulation_vertices:
                    if abs(v.pos.complex() - w) < epsilon:
                        ladderpole_vertices[i] = v    
                        break     

        ### the following is the list with correctly replaced vertices
        all_ladderpole_vertices = [v for L in ladderpoles_vertices for v in L]

        ladderpole_descendant_segments = []
        for ladderpole_vertices in ladderpoles_vertices:
            segment = [con.coast.index(ladderpole_vertices[0]), con.coast.index(ladderpole_vertices[-1])]
            segment.sort()
            ladderpole_descendant_segments.append( segment )

        print(ladderpole_descendant_segments)

        con.mark_ladderpole_descendants(ladderpole_descendant_segments)

##### experimenting
        # # print 'orig verts', original_desired_vertices
        # important_vertices = []
        # con.build_long_and_mid(max_length = max_length, max_num_tetrahedra = 1000)
        # con.update_coast()
        # for i, v in enumerate(con.coast):
        #     if len(v.ladderpole_ancestors) > 0:
        #         w = con.coast[(i+1)]
        #         if len(v.ladderpole_ancestors.intersection(w.ladderpole_ancestors)) > 0:
        #             if not v in all_ladderpole_vertices and not w in all_ladderpole_vertices:
        #                 ## we found a ladderpole descendant edge that does not touch one of the original vertices
        #                 for u in con.coast:
        #                     u.ladderpole_ancestors = set()
        #                 v.ladderpole_ancestors.add(0)
        #                 w.ladderpole_ancestors.add(0)
        #                 important_vertices = [v,w]
        #                 break

        # con.build_long_and_mid(max_length = max_length, max_num_tetrahedra = max_num_tetrahedra)

        # print 'important verts', important_vertices
        # for v in important_vertices:
        #     z = T.drawing_scale * v.pos.complex()
        #     pyx_fill_col = pyx.deco.filled([pyx.color.rgb.black])
        #     T.canv.fill(pyx.path.circle(z.real, z.imag, 0.02), [pyx_fill_col])



#### end experimenting

        hit_max_tetrahedra = False ### default assumption is that we had enough tetrahedra to get the max_length we want.
        print(build_type)
        if build_type == 'build_naive':
            con.build_naive(max_num_tetrahedra = max_num_tetrahedra)
        elif build_type == 'build_on_coast':
            hit_max_tetrahedra = con.build_on_coast(max_length = max_length, max_num_tetrahedra = max_num_tetrahedra)
        elif build_type == 'build_make_long_descendant_edges_internal':
            hit_max_tetrahedra = con.build_make_long_descendant_edges_internal(max_length = max_length, max_num_tetrahedra = max_num_tetrahedra)
        elif build_type == 'build_explore_prongs':
            hit_max_tetrahedra = con.build_explore_prongs(max_length = max_length, max_num_tetrahedra = max_num_tetrahedra)
        elif build_type == 'build_long_and_mid':
            hit_max_tetrahedra = con.build_long_and_mid(max_length = max_length, max_num_tetrahedra = max_num_tetrahedra)
        
        if hit_max_tetrahedra:
            output_filename = output_filename[:-4] + '_hitmax.pdf'

        #######

        # eq = con.segment_between( ladderpoles_vertices[0][0], ladderpoles_vertices[0][1] )   ## segment under one edge of ladderpole
        # eq = con.segment_between( ladderpoles_vertices[0][0], ladderpoles_vertices[0][-1] )   ## 0th ladderpole

        grad = pyx.color.gradient.Hue
        # colours = {"blue":pyx.color.rgb.blue, "red":pyx.color.rgb.red}  
        colours = {"blue":pyx.color.rgb(0,0,0.5), "red":pyx.color.rgb(0.5,0,0)}

        ct_lw = draw_args['ct_lw']

        draw_options = [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round, pyx.deco.colorgradient(grad)] ## this may get overwritten with colour information for the ladder

        ##### draw the Cannon-Thurston curve

        if draw_CT_curve:
            if draw_args['only_draw_ladderpoles']:
                for j, ladderpole_vertices in enumerate(ladderpoles_vertices):
                    # if j%2 == 0:
                    #     col = colours["red"]
                    # else:
                    #     col = colours["blue"]
                    # draw_options = [pyx.style.linewidth(ct_lw), col]
                    for i in range(len(ladderpole_vertices) - 1):  # fenceposts
                        path = con.segment_between(ladderpole_vertices[i], ladderpole_vertices[i+1])  
                        for v in path[:-1]:
                            assert v.is_ladderpole_descendant()
                        path_C = [ T.drawing_scale * v.pos.complex() for v in path ]
                        draw_path(T.canv, path_C, draw_options)  
            else:
                path = con.coast
                path.remove(con.infinity)
                path_C = [ T.drawing_scale * v.pos.complex() for v in path ]
                draw_path(T.canv, path_C, draw_options)  


        ##############
        # adj_verts, adj_edges = con.vertices_and_edges_adjacent_to_infinity()  

        ### continent drawing the boundary triangulation (we don't generally do this because we have boundary_triangulation to do this)
        # lines for triangles meeting infinity

        # for endpoints, veering_colour in adj_edges:
        #     z, w = [T.drawing_scale * v.pos.complex() for v in endpoints]
        #     pyx_stroke_col = pyx.deco.stroked([colours[veering_colour]])
        #     T.canv.stroke(pyx.path.line( z.real, z.imag, w.real, w.imag),  [pyx.style.linewidth(lw * 5), pyx_stroke_col]  )

        # # dots for edges from infinity
        # for v,veering_colour in adj_verts:
        #     z = v.pos.complex()
        #     pyx_fill_col = pyx.deco.filled([colours[veering_colour]])
        #     T.canv.fill(pyx.path.circle(z.real, z.imag, 0.02), [pyx_fill_col])

        ### continent drawing the left_ladder_pole_vertices
        # T = B.torus_triangulation_list[0]
        # for L in T.ladder_list:
        #     for v in L.left_ladder_pole_vertices():
        #         v *= T.drawing_scale
        #         if L.is_even:
        #             col = colours["blue"]
        #         else:
        #             col = colours["red"]
        #         T.canv.stroke(pyx.path.circle(v.real, v.imag, 0.15), [pyx.style.linewidth(lw * 3), col])

        ### circles around found vertices
        # for v in con.found_vertices:
        #     pos = v.pos.complex()
        #     pos *= T.drawing_scale
        #     T.canv.stroke(pyx.path.circle(pos.real, pos.imag, 0.2), [pyx.style.linewidth(lw * 3)])



        #### draw upper and lower landscapes for the continent

        if draw_landscapes:
            lower_colours = {True: pyx.color.rgb(0.5,0,0), False: pyx.color.rgb(0,0,0.5)}
            upper_colours = {True: pyx.color.rgb(0.9,0,0), False: pyx.color.rgb(0,0,0.9)}

            landscape_edges = con.boundary_landscape_edges()

            colours = [lower_colours, upper_colours]
            # for i in range(2):
            i = 1
            for e in landscape_edges[i]:
                col = colours[i][e.is_red]
                u, v = e.vertices
                if u == con.infinity or v == con.infinity:
                    if u == con.infinity:
                        z = T.drawing_scale * v.pos.complex()
                    else:
                        z = T.drawing_scale * u.pos.complex()
                    T.canv.fill(pyx.path.circle(z.real, z.imag, 0.05), [col])
                else:
                    draw_path(T.canv, [T.drawing_scale * u.pos.complex(), T.drawing_scale * v.pos.complex()], [pyx.style.linewidth(0.5 * ct_lw), col])

        #### draw lightning curves

        if draw_lightning_curve:
            lightning_colours = [pyx.color.rgb(0,0.5,0), pyx.color.rgb(0.5,0,0.5)]
            lightning_curves = con.lightning_curves([])  ## only lightning curves for infinity
            # lightning_curves = con.lightning_curves(all_ladderpole_vertices)  ## add in more lightning curves
            for i in range(2):
                for crv in lightning_curves[i]:
                    ## remove any infinities
                    crv = [c for c in crv if c != con.infinity]

                    crv = [ T.drawing_scale * c.pos.complex() for c in crv ]
                    draw_path(T.canv, crv, [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round, lightning_colours[i]])
                    # ## trim to ladder poles
                    # ladderpole_vertex_indices = []
                    # for i, v in enumerate(crv):
                    #     if v in all_ladderpole_vertices:
                    #         ladderpole_vertex_indices.append(i)
                    # if len(ladderpole_vertex_indices) > 0:
                    #     crv = crv[ladderpole_vertex_indices[0]: ladderpole_vertex_indices[-1] + 1]
                    #     # for e in crv:
                    #         # pts = [T.drawing_scale * v.pos.complex() for v in e.vertices]
                    #         # draw_path(T.canv, pts, [pyx.style.linewidth(ct_lw)])  
                    #     crv = [ T.drawing_scale * c.pos.complex() for c in crv ]
                    #     draw_path(T.canv, crv, [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round]) 




    out_canvas = pyx.canvas.canvas()
    height_offset = 0.0
    canvases = [T.canv for T in B.torus_triangulation_list]
    for i,c in enumerate(canvases):
        out_canvas.insert(c, attrs=[pyx.trafo.translate(-c.bbox().left(), height_offset - c.bbox().bottom())])
        height_offset += c.bbox().height() + 0.05 ### add a tiny bit to stop crashes due to line width
    out_canvas.writePDFfile(output_filename)

def draw_cannon_thurston_from_veering_isosigs_file(veering_isosigs_filename, output_dirname, max_num_tetrahedra = 500, max_length = 0.1, num_to_draw = None, draw_args = None, build_type = None):
    veering_isosigs_list = parse_data_file(veering_isosigs_filename)
    if num_to_draw != None:
        to_draw = veering_isosigs_list[:num_to_draw]
    else:
        to_draw = veering_isosigs_list

    shapes_data = read_from_pickle('Data/veering_shapes_up_to_ten_tetrahedra.pkl')
    for veering_isosig in to_draw:
        print(veering_isosig)
        tet_shapes = shapes_data[veering_isosig]
        filename = output_dirname + '/' + veering_isosig + '_' + str(max_num_tetrahedra) + '_' + str(max_length) + '_' + build_type + '.pdf'
        draw_continent(veering_isosig, tet_shapes, max_num_tetrahedra, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )


if __name__ == '__main__':
    # draw_args = {'draw_boundary_triangulation':True, 'only_draw_ladderpoles': True, 'ct_lw': 0.002, 'global_drawing_scale': 4, 'draw_labels': False, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps
    draw_args = {'draw_boundary_triangulation':True, 'only_draw_ladderpoles': True, 'ct_lw': 0.02, 'global_drawing_scale': 4, 'draw_labels': False, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps

    
    # max_num_tetrahedra = 50000
    # max_num_tetrahedra = 100000
    # max_num_tetrahedra = 400000
    max_num_tetrahedra = 2000000
    # max_length = 0.15
    # max_length = 0.1
    # max_length = 0.07
    # max_length = 0.06
    # max_length = 0.02
    max_length = 0.0

    draw_args['ct_lw'] = 0.2 * max_length 

    # build_type = 'build_naive'
    # build_type = 'build_on_coast'
    # build_type = 'build_make_long_descendant_edges_internal'
    # build_type = 'build_explore_prongs'
    build_type = 'build_long_and_mid'

    veering_isosig = 'cPcbbbiht_12'
    # # # # veering_isosig = 'cPcbbbdxm_10'
    # # # # veering_isosig = 'dLQacccjsnk_200'
    # # # veering_isosig = 'eLMkbcddddedde_2100'
    # # # # veering_isosig = 'eLAkaccddjsnak_2001'
    # # # veering_isosig = 'gLAMPbbcdeffdhwqqqj_210202'
    # veering_isosig = 'gLLAQbecdfffhhnkqnc_120012'
    # # # # veering_isosig = 'iLLLAQccdffgfhhhqgdatgqdm_21012210' ## no symmetry - helps us spot errors
    # # # veering_isosig = 'iLLPwQcccdfehghhhggaahhbg_20102211'

    shapes_data = read_from_pickle('Data/veering_shapes_up_to_ten_tetrahedra.pkl')
    tet_shapes = shapes_data[veering_isosig]
    filename = 'Images/Cannon-Thurston/' + veering_isosig + '_' + str(max_num_tetrahedra) + '_' + str(max_length) + '_' + build_type + '.pdf'
    draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, draw_CT_curve = True, draw_lightning_curve = False, draw_landscapes = False, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )
    # draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, draw_CT_curve = False, draw_lightning_curve = True, draw_landscapes = False, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )
    

    ### draw many:

    # num_to_draw = 30
    # draw_cannon_thurston_from_veering_isosigs_file('Data/veering_census.txt', 'Images/Cannon-Thurston', max_num_tetrahedra = max_num_tetrahedra, max_length = max_length, num_to_draw = num_to_draw, draw_args = draw_args, build_type = build_type)
    


