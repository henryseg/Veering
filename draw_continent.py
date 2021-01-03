
import pyx ### vector graphics 
import cmath

from file_io import parse_data_file, read_from_pickle, output_to_pickle
from taut import isosig_to_tri_angle
from veering import veering_triangulation
from continent import continent
from boundary_triangulation import boundary_triangulation


# def pre_draw_transformation( z, ladder_holonomy ):
    # return z/ladder_holonomy

def draw_path(canv, path_C, draw_options, fill = False):
    p = pyx.path.path( pyx.path.moveto(path_C[0].real, path_C[0].imag) )
    for coord in path_C[1:]:  ### this is how path drawing works...
        p.append( pyx.path.lineto(coord.real, coord.imag) )
    if fill:
        canv.fill(p, draw_options)
    else:
        canv.stroke(p, draw_options)

def dividers_to_lightning_curves_spiky(dividers):
    # print(('dividers', len(dividers[0]), len(dividers[1])))
    ### now find lightning curve from the edges of the divider
    out = []
    for i in range(2):  ## upper and lower
        curves = []
        for divider in dividers[i]:
            if len(divider) <= 2:
                ### add nothing
                continue
            a, b = divider[:2]
            x = a.shared_vertex(b)
            lightning_curve = [x]
            for c in divider[2:]:
                y = b.shared_vertex(c)
                if x == y:
                    a, b = a, c
                else:
                    lightning_curve.append(y)
                    a, b = b, c
                    x = y
            curves.append(lightning_curve)
        out.append(curves)    

    return out

def special_vertex_in_and_out(divider_list, v, v_inds):
    ### v_inds are indices of edges in divider_list that are incident to v
    ### answers how we get to v and back out when we are doing the non-spiky lightning curves
    if len(v_inds) % 2 == 1:
        # print 'odd'
        ind = (v_inds[0] + v_inds[-1])//2
        return ( [v.pos.complex()], ind, ind ) ### last midpoint is at ind, then we do the extra point in the list, then start midpoints up again with ind
    else:
        # print 'even'
        ind = (v_inds[0] + v_inds[-1] + 1)//2
        mid_step = (divider_list[ind-1].midpoint() + divider_list[ind].midpoint()) * 0.5
        return ( [mid_step, v.pos.complex(), mid_step], ind - 1, ind ) ### more complicated since we have an extra point on the midpoints curve

def lightning_curve_from_dividers(dividers, a, b, special_vertices = [], spiky = True):
    #### give it only upper dividers, or only lower dividers
    curve = []
    for divider_list in dividers:
        ### check if both s and e are endpoints of this divider list
        a_inds = []
        b_inds = []
        for i, div in enumerate(divider_list):
            if a in div.vertices:
                a_inds.append(i) 
            if b in div.vertices:
                b_inds.append(i) 
        if len(a_inds) > 0 and len(b_inds) > 0:  # we found the correct divider_list
            divider_list = divider_list[:] ### make a copy
            a_to_b = a_inds[0] < b_inds[0]

            if spiky:
                if a_to_b:
                    divider_list = divider_list[a_inds[-1]:b_inds[0]]
                else:
                    divider_list = divider_list[b_inds[-1]:a_inds[0]]
                p, q = divider_list[:2]
                x = p.shared_vertex(q)
                lightning_curve = [x.pos.complex()]
                for r in divider_list[2:]:
                    y = q.shared_vertex(r)
                    if x == y:
                        p, q = p, r
                    else:
                        lightning_curve.append(y.pos.complex())
                        p, q = q, r
                        x = y
                if a_to_b:
                    lightning_curve = [a.pos.complex()] + lightning_curve + [b.pos.complex()]
                else:
                    lightning_curve = [b.pos.complex()] + lightning_curve + [a.pos.complex()]
                return lightning_curve

            else:     ### have to hit a, b, and any special vertices along the way
                ### find all special vertices along the path from a to b (or b to a), go in and out of each in turn...
                if a_to_b:
                    s_inds, e_inds = a_inds, b_inds
                    s, e = a, b
                else:
                    s_inds, e_inds = b_inds, a_inds
                    s, e = b, a
                special_verts = special_vertices[:]
                if a in special_vertices:
                    special_verts.remove(a)
                if b in special_vertices:
                    special_verts.remove(b)
                visited_special_vertices = [s]
                visited_special_vertex_inds = [s_inds]
                for j in range(s_inds[-1] + 1, e_inds[0]):
                    edge = divider_list[j]
                    p, q = edge.vertices
                    v = None
                    if p in special_vertices:
                        v = p
                    elif q in special_vertices:
                        v = q
                    if v != None:
                        if v == visited_special_vertices[-1]:
                            visited_special_vertex_inds[-1].append(j)
                        else:
                            visited_special_vertices.append(v)
                            visited_special_vertex_inds.append([j])
                visited_special_vertices.append(e)
                visited_special_vertex_inds.append(e_inds)
                all_extra_steps = []
                all_in_indices = []
                all_out_indices = []
                for v, v_inds in zip(visited_special_vertices, visited_special_vertex_inds):
                    extra_steps, in_ind, out_ind = special_vertex_in_and_out(divider_list, v, v_inds)
                    all_extra_steps.append(extra_steps)
                    all_in_indices.append(in_ind)
                    all_out_indices.append(out_ind)
                lightning_curve = []
                for j in range(len(all_extra_steps) - 1):
                    lightning_curve.extend(all_extra_steps[j])
                    lightning_curve.extend([edge.midpoint() for edge in divider_list[all_out_indices[j] : all_in_indices[j+1] + 1]])
                lightning_curve.extend(all_extra_steps[-1])
                if len(s_inds) % 2 == 0:
                    lightning_curve = lightning_curve[1:]
                if len(e_inds) % 2 == 0:
                    lightning_curve = lightning_curve[:-1]
                return lightning_curve




def draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, draw_CT_curve = True, draw_lightning_curve = False, draw_jordan_curve = False, draw_landscapes = False, draw_box_for_cohom_frac = False, draw_alignment_dots = False, max_length = 0.1, output_filename = None, draw_args = None, build_type = None, more = False ):

    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle, tet_shapes = tet_shapes)
    B = boundary_triangulation(vt)
    B.generate_canvases(args = draw_args)

    out_data = []

    for i,T in enumerate(B.torus_triangulation_list):
        # print(('cusp', i))
        ### make initial_tet_face be in the lower left of the fundamental domain
        # initial_tet_face = T.ladder_list[0].ladder_unit_list[0]
        
        ### make initial_tet_face be in the middle of the fundamental domain
        num_ladders = len(T.ladder_list)
        L = T.ladder_list[int(num_ladders/2 - 1)]  ## -1 because we split the last ladder between the left and right
        num_ladder_units = len(L.ladder_unit_list)
        initial_tet_face = L.ladder_unit_list[int(num_ladder_units/2)]

        # print(('initial_tet_face', initial_tet_face)) 
        # print(('origin_in_C', initial_tet_face.origin_in_C))
        # print(('verts_pos', initial_tet_face.verts_pos))
        ### want to draw a box which will contain a fund dom, will be what we render as a cohom frac

        ladderpoles_vertices = T.left_ladder_pole_vertices() 
        left_extra_ladderpole = [z - T.sideways_once_holonomy for z in ladderpoles_vertices[-1]]
        right_extra_ladderpole = [z + T.sideways_once_holonomy for z in ladderpoles_vertices[0]]

        ### doing the following crashes develop_ideal_hyperbolic_tetrahedra every time :(

        # print T.sideways_once_holonomy.imag
        # if abs(T.sideways_once_holonomy.imag) > 0.001:   
        #     if T.sideways_once_holonomy.imag > 0.0: ### add a copy below to the right side, add a copy above to the left side 
        #     ### (this probably isn entirely correct, but maybe its good enough...)
        #         right_extra_ladderpole = [z - T.ladder_holonomy for z in right_extra_ladderpole[:-1]] + right_extra_ladderpole
        #         left_extra_ladderpole = left_extra_ladderpole + [z + T.ladder_holonomy for z in left_extra_ladderpole[1:]]
        #     else: ### add a copy above to the right side, add a copy below to the left side 
        #         right_extra_ladderpole = right_extra_ladderpole + [z + T.ladder_holonomy for z in right_extra_ladderpole[1:]]  
        #         left_extra_ladderpole = [z - T.ladder_holonomy for z in left_extra_ladderpole[:-1]] + left_extra_ladderpole
   
        if more:
            more_ladderpoles_vertices = [left_extra_ladderpole] + ladderpoles_vertices[:] + [right_extra_ladderpole]
            ### extend all ladders by one ^b^b^b two... above and below
            more_ladderpoles_vertices = [[L[-3] - T.ladder_holonomy, L[-2] - T.ladder_holonomy] + L + [L[1] + T.ladder_holonomy, L[2] + T.ladder_holonomy] for L in more_ladderpoles_vertices]


        if more:
            desired_vertices = [v for L in more_ladderpoles_vertices for v in L]
        else:
            desired_vertices = [v for L in ladderpoles_vertices for v in L]

        con = continent( vt, initial_tet_face, desired_vertices = desired_vertices )
        
        con.build_fundamental_domain()  ## expand the continent until we have all vertices of the boundary triangulation fundamental domain

        # print(('unfound desired_vertices', con.desired_vertices))
        # assert con.desired_vertices == [] ### found all the desired vertices
        if con.desired_vertices != []:
            print(veering_isosig, 'did not find all torus triangulation vertices')
            return False

        # now replace ladderpoles_vertices with the continent's corresponding vertices 
        epsilon = 0.001
        for ladderpole_vertices in ladderpoles_vertices:
            for i, w in enumerate(ladderpole_vertices):
                for v in con.boundary_triangulation_vertices:
                    if abs(v.pos.complex() - w) < epsilon:
                        ladderpole_vertices[i] = v    
                        break     

        if more:
            for ladderpole_vertices in more_ladderpoles_vertices:
                for i, w in enumerate(ladderpole_vertices):
                    for v in con.boundary_triangulation_vertices:
                        if abs(v.pos.complex() - w) < epsilon:
                            ladderpole_vertices[i] = v    
                            break     

        ### the following is the list with correctly replaced vertices
        all_ladderpole_vertices = [v for L in ladderpoles_vertices for v in L] ### dont need this to be bigger when doing "more" 


        ladderpole_descendant_segments = []

        if more:
            foo = more_ladderpoles_vertices
        else:
            foo = ladderpoles_vertices
        for ladderpole_vertices in foo:
            segment = [con.coast.index(ladderpole_vertices[0]), con.coast.index(ladderpole_vertices[-1])]
            segment.sort()
            ladderpole_descendant_segments.append( segment )



        # print(ladderpole_descendant_segments)

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
        # print(build_type)
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

        draw_options = [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round, pyx.deco.colorgradient(grad, steps = 256)] ## this may get overwritten with colour information for the ladder

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

        #### draw upper and lower landscapes for the continent

        if draw_landscapes:
            lower_colours = {True: pyx.color.rgb(0.5,0.3,0), False: pyx.color.rgb(0,0.3,0.5)}
            upper_colours = {True: pyx.color.rgb(0.9,0.3,0), False: pyx.color.rgb(0,0.3,0.9)}

            landscape_edges = con.boundary_landscape_edges()

            colours = [lower_colours, upper_colours]
            for i in range(2):
            # i = 1
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
            # lightning_colours = [pyx.color.rgb(0,0.5,0), pyx.color.rgb(0.5,0,0.5)]  ### green, purple
            # lightning_colours = [pyx.color.rgb(1,0,0), pyx.color.rgb(1,0,0)]  ### red, red
            lightning_colours = [pyx.color.rgb(0,0,0), pyx.color.rgb(0,0,0)]  ### black, black       

            dividers = con.lightning_dividers([])  ## only lightning curves for infinity
            # dividers = con.lightning_dividers(all_ladderpole_vertices)  ## add in more lightning curves

            for ladder in ladderpoles_vertices:
                for j in range(len(ladder) - 1):
                    s, e = ladder[j], ladder[j+1]
                    for i in range(2):
                        lightning_curve = lightning_curve_from_dividers(dividers[i], s, e, special_vertices = all_ladderpole_vertices, spiky = False)
                        lightning_curve = [ T.drawing_scale * c for c in lightning_curve ]
                        draw_path(T.canv, lightning_curve, [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round, lightning_colours[i]])


            # for i in range(2):
            #     for crv in lightning_curves[i]:
            #         ## remove any infinities
            #         crv = [c for c in crv if c != con.infinity]
            #         # if draw_args['only_draw_ladderpoles']:
            #         #     ladderpole_indices = []
            #         #     for j, c in enumerate(crv):
            #         #         if c in all_ladderpole_vertices:
            #         #             ladderpole_indices.append(j)
            #         #     if len(ladderpole_indices) <= 1:
            #         #         crv = []
            #         #     else:
            #         #         crv = crv[ladderpole_indices[0]:ladderpole_indices[-1] + 1]
            #         crv = [ T.drawing_scale * c.pos.complex() for c in crv ]
            #         if len(crv) > 0:
            #             draw_path(T.canv, crv, [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round, lightning_colours[i]])
            #         # ## trim to ladder poles
            #         # ladderpole_vertex_indices = []
            #         # for i, v in enumerate(crv):
            #         #     if v in all_ladderpole_vertices:
            #         #         ladderpole_vertex_indices.append(i)
            #         # if len(ladderpole_vertex_indices) > 0:
            #         #     crv = crv[ladderpole_vertex_indices[0]: ladderpole_vertex_indices[-1] + 1]
            #         #     # for e in crv:
            #         #         # pts = [T.drawing_scale * v.pos.complex() for v in e.vertices]
            #         #         # draw_path(T.canv, pts, [pyx.style.linewidth(ct_lw)])  
            #         #     crv = [ T.drawing_scale * c.pos.complex() for c in crv ]
            #         #     draw_path(T.canv, crv, [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round]) 

        if draw_jordan_curve:
            jordan_colours = [pyx.color.rgb(0,0,0), pyx.color.rgb(1,1,1)]  ### black, white   
            dividers = con.lightning_dividers([])  ## only lightning curves for infinity

            layer1 = T.canv.layer("layer1")
            layer2 = T.canv.layer("layer2", below = "layer1")

            for ladder in ladderpoles_vertices:
                for j in range(len(ladder) - 1):
                    s, e = ladder[j], ladder[j+1]
                    CT_ladder = con.segment_between(s, e)
                    if CT_ladder[0] != s:
                        CT_ladder.reverse()
                    assert CT_ladder[0] == s
                    assert CT_ladder[-1] == e
                    for i in range(2):
                        lightning_curve = lightning_curve_from_dividers(dividers[i], s, e, special_vertices = all_ladderpole_vertices, spiky = False)
                        CT_ladder_pos = [v.pos.complex() for v in CT_ladder ]
                        if CT_ladder_pos[0] == lightning_curve[0]:
                            lightning_curve.reverse() 
                        assert CT_ladder_pos[0] == lightning_curve[-1] and CT_ladder_pos[-1] == lightning_curve[0]
                        lightning_curve = lightning_curve[1:-1]
                        loop = CT_ladder_pos + lightning_curve
                        path_C = [ T.drawing_scale * c for c in loop ]
                        # draw_path(T.canv, path_C, [jordan_colours[i]], fill = True) 
                        draw_path(layer1, path_C, [jordan_colours[i]], fill = True)     
            box = T.canv.bbox()
            layer2.fill(box.rect(), [pyx.color.rgb(0.5,0.5,0.5)])

            for divider_list in dividers[0]:
                for edge in divider_list:
                    edgeC = [T.drawing_scale * v.pos.complex() for v in edge.vertices]
                    draw_path(layer1, edgeC, [pyx.color.rgb(1,0,0), pyx.style.linewidth(0.005)])
            for divider_list in dividers[1]:
                for edge in divider_list:
                    edgeC = [T.drawing_scale * v.pos.complex() for v in edge.vertices]
                    draw_path(layer1, edgeC, [pyx.color.rgb(0,0,1), pyx.style.linewidth(0.005)])


        if draw_box_for_cohom_frac:
            box = T.canv.bbox()
            diam = pyx.unit.tocm( max([box.right() - box.left(), box.top() - box.bottom()]) )
            ### pyx stores lengths in a weird format, we have to convert to cm to get a float out
            box_center = complex(pyx.unit.tocm( 0.5*(box.right() + box.left()) ), pyx.unit.tocm( 0.5*(box.top() + box.bottom()) ))
            box_right_midpoint = box_center + complex(0.5*diam, 0)

            # T.canv.stroke(box.rect())
            T.canv.stroke(pyx.path.rect(box_center.real - 0.5*diam, box_center.imag - 0.5*diam, diam, diam))  # lower left corner coords, width, height


            inf_vert = initial_tet_face.face
            zero_vert = inf_vert - ((inf_vert%2)*2 - 1)  ## swaps 0 with 1, 2 with 3
            one_vert = inf_vert - 2*(((inf_vert/2) % 2)*2 - 1) ## swaps 0 with 2, 1 with 3

            zero_vert_pos = T.drawing_scale * initial_tet_face.verts_pos[zero_vert].complex()
            one_vert_pos = T.drawing_scale * initial_tet_face.verts_pos[one_vert].complex()

            half_pos = 0.5 * (zero_vert_pos + one_vert_pos)  
            ### we must rotate, translate, and scale the cohom fractal picture to fit in the box

            if draw_alignment_dots:
                T.canv.fill(pyx.path.circle(zero_vert_pos.real, zero_vert_pos.imag, 0.2))
                T.canv.fill(pyx.path.circle(half_pos.real, half_pos.imag, 0.15))
                T.canv.fill(pyx.path.circle(box_center.real, box_center.imag, 0.1))
                T.canv.fill(pyx.path.circle(box_right_midpoint.real, box_right_midpoint.imag, 0.1))

            ### need to send zero_vert_pos to box_center, and half_pos to box_right_midpoint



            # print veering_isosig
            out_data.append(veering_isosig)
            # print 'tet face', initial_tet_face
            out_data.append((initial_tet_face.tet_num, initial_tet_face.face))

            picture_unit = one_vert_pos - zero_vert_pos
            translation = (box_center - zero_vert_pos)/picture_unit ## need to write this in coord system of one_vert_pos and zero_vert_pos

            # print 'translation', [translation.real, translation.imag]
            out_data.append((translation.real, translation.imag))

            
            
            ### first rotate and scale, then do parabolicBy2DVector(v)
            complex_scale = (box_right_midpoint - box_center)/(half_pos - zero_vert_pos)

            # print 'scale, rotate', cmath.polar(complex_scale)
            out_data.append(cmath.polar(complex_scale))

            # pic_center = zero_vert_pos
            # print pyx.unit.tocm(box.right())
            # print type(pyx.unit.tocm(box.right()))
            # rad = pyx.unit.tocm( max([box.right() - pic_center.real, pic_center.real - box.left(), box.top() - pic_center.imag, pic_center.imag - box.bottom()]) )
            ### pyx stores lengths in a weird format, we have to convert to cm to get a float out
        

            # T.canv.stroke(pyx.path.rect(pic_center.real - rad, pic_center.imag - rad, 2*rad, 2*rad))  # lower left corner coords, width, height

            # right_midpoint = pic_center + complex(rad,0)
            # complex_scale = (right_midpoint - pic_center) / (half_pos - pic_center)
            # print "complex_scale", complex_scale

            ### to do: 
            ### 1. remove bits of lightning curve beyond the fund domain DONE
            ### 2. how to get the data into the cohom frac code DONE

        ## circles around found vertices
        for pos in desired_vertices:
            # pos = v.pos.complex()
            pos *= T.drawing_scale
            T.canv.stroke(pyx.path.circle(pos.real, pos.imag, 0.2), [pyx.style.linewidth(0.1), pyx.color.rgb.black])




    out_canvas = pyx.canvas.canvas()
    height_offset = 0.0
    canvases = [T.canv for T in B.torus_triangulation_list]
    for i,c in enumerate(canvases):
        out_canvas.insert(c, attrs=[pyx.trafo.translate(-c.bbox().left(), height_offset - c.bbox().bottom())])
        height_offset += c.bbox().height() + 0.05 ### add a tiny bit to stop crashes due to line width
    out_canvas.writePDFfile(output_filename)
    return out_data

def draw_cannon_thurston_from_veering_isosigs_file(veering_isosigs_filename, output_dirname, max_num_tetrahedra = 500, max_length = 0.1, interval_to_draw = None, draw_args = None, build_type = None):
    veering_isosigs_list = parse_data_file(veering_isosigs_filename)
    if interval_to_draw != None:
        to_draw = veering_isosigs_list[interval_to_draw[0]:interval_to_draw[1]]
    else:
        to_draw = veering_isosigs_list

    shapes_data = read_from_pickle('Data/veering_shapes.pkl')
    for veering_isosig in to_draw:
        print(veering_isosig)
        tet_shapes = shapes_data[veering_isosig]
        filename = output_dirname + '/' + veering_isosig + '_' + str(max_num_tetrahedra) + '_' + str(max_length) + '_' + build_type + '.pdf'
        # draw_continent(veering_isosig, tet_shapes, max_num_tetrahedra, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )
        try:
            draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, draw_CT_curve = False, draw_lightning_curve = False, draw_jordan_curve = True, draw_landscapes = False, draw_box_for_cohom_frac = False, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )
        except: 
            print('failed to draw ' + veering_isosig)

def draw_jigsaw_from_veering_isosigs_file(veering_isosigs_filename, output_dirname, jigsaw_data_out_filename = "jigsaw_data.pkl", max_num_tetrahedra = 2000000, max_length = 0.2, interval_to_draw = None):
    veering_isosigs_list = parse_data_file(veering_isosigs_filename)
    if interval_to_draw != None:
        to_draw = veering_isosigs_list[interval_to_draw[0]:interval_to_draw[1]]
    else:
        to_draw = veering_isosigs_list

    build_type = 'build_long_and_mid'
    # draw_args = {'draw_boundary_triangulation':False, 'only_draw_ladderpoles': True, 'ct_lw': 0.2 * max_length, 'global_drawing_scale': 4, 'draw_labels': False, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps
    draw_args = {'draw_boundary_triangulation':True, 'draw_labels': True, 'only_draw_ladderpoles': True, 'ct_lw': 0.2 * max_length, 'global_drawing_scale': 4, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps
    
    # shapes_data = read_from_pickle('Data/veering_shapes.pkl')
    shapes_data = read_from_pickle('Data/shapes_jig_no_symm.pkl')

    data_for_cohom_fracs = {}
    for i, veering_isosig in enumerate(to_draw):
        # print(veering_isosig)
        if i%50 == 0:
            print(i)
        tet_shapes = shapes_data[veering_isosig]
        # print 'tet_shapes', tet_shapes


        filename = output_dirname + '/' + veering_isosig + '_' + str(max_num_tetrahedra) + '_' + str(max_length) + '_' + build_type + '.pdf'
        # draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, draw_CT_curve = True, draw_lightning_curve = False, draw_landscapes = False, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )
        out = draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, draw_CT_curve = False, draw_lightning_curve = True, draw_landscapes = False, draw_box_for_cohom_frac = True, draw_alignment_dots = True, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )
        if out != False:
            data_for_cohom_fracs[out[0]] = out[1:]
    output_to_pickle(data_for_cohom_fracs, jigsaw_data_out_filename)

if __name__ == '__main__':
    # draw_args = {'draw_boundary_triangulation':True, 'only_draw_ladderpoles': True, 'ct_lw': 0.002, 'global_drawing_scale': 4, 'draw_labels': False, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps
    # draw_args = {'draw_boundary_triangulation':True, 'only_draw_ladderpoles': True, 'ct_lw': 0.02, 'global_drawing_scale': 4, 'draw_labels': True, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps
    draw_args = {'draw_boundary_triangulation':False, 'only_draw_ladderpoles': True, 'ct_lw': 0.02, 'global_drawing_scale': 4, 'draw_labels': False, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps

    
    # max_num_tetrahedra = 5000
    # max_num_tetrahedra = 50000
    # max_num_tetrahedra = 100000
    # max_num_tetrahedra = 400000
    max_num_tetrahedra = 2000000
    max_length = 0.4
    # max_length = 0.2
    # max_length = 0.15
    # max_length = 0.1
    # max_length = 0.07
    # max_length = 0.06
    # max_length = 0.02
    # max_length = 0.01

    draw_args['ct_lw'] = 0.2 * max_length 

    # build_type = 'build_naive'
    # build_type = 'build_on_coast'
    # build_type = 'build_make_long_descendant_edges_internal'
    # build_type = 'build_explore_prongs'
    build_type = 'build_long_and_mid'

    # veering_isosig = 'cPcbbbiht_12'
    # # # # veering_isosig = 'cPcbbbdxm_10'
    # veering_isosig = 'dLQacccjsnk_200'
    # veering_isosig = 'eLMkbcddddedde_2100'
    # veering_isosig = 'eLAkaccddjsnak_2001'
    # veering_isosig = 'gLAMPbbcdeffdhwqqqj_210202'
    # veering_isosig = 'gLLAQbecdfffhhnkqnc_120012'
    # # # # # veering_isosig = 'iLLLAQccdffgfhhhqgdatgqdm_21012210' ## no symmetry - helps us spot errors
    # # # # veering_isosig = 'iLLPwQcccdfehghhhggaahhbg_20102211'
    # # veering_isosig = 'jLAwwAQbcbdfghihihhwhnaaxrn_211211021' ## first non geometric
    # # veering_isosig = 'nLLwMLPMMkbeefeihjkjlmlmhhaaaektxnaqrs_0111000011220'  ### quite big negative shape
    # # veering_isosig = 'qLvPvvMQQLQkccgkgjkmlknpooppoqjaajqqhhqqaqxhhh_0222110112222211'
    # veering_isosig = 'fLLQcbeddeehhnkhh_21112'
    veering_isosig = 'eLAkbbcdddhwqj_2102'

    # # # veering_isosig = 'mLvLLLQQQbegikhjiilkllhiardrnnkxeif_120000112222'
    # veering_isosig = 'mLvLLMMQQcehfhjlklkjlktilbbjumhtfai_011220220111'
    # # veering_isosig = 'mLLvLQLQQbeffjglhlkjklxxxjsfqjhhoqo_102210101022'
    # # # veering_isosig ='kLLLAPPkcdgfehhjijjhfhaqiphffj_2010222001'

    shapes_data = read_from_pickle('Data/veering_shapes_up_to_twelve_tetrahedra.pkl')
    # shapes_data = read_from_pickle('Data/veering_shapes.pkl')
    # shapes_data = read_from_pickle('Data/shapes_jig_no_symm.pkl')
    tet_shapes = shapes_data[veering_isosig]
    # print tet_shapes
    filename = 'Images/Cannon-Thurston/' + veering_isosig + '_' + str(max_num_tetrahedra) + '_' + str(max_length) + '_' + build_type + '.pdf'
    
    more = True ### generate bigger continent to get lightning curves right
    # # draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, draw_CT_curve = True, draw_lightning_curve = False, draw_landscapes = False, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )
    # draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, draw_CT_curve = True, draw_lightning_curve = True, draw_landscapes = False, draw_box_for_cohom_frac = True, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )
    # draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, draw_CT_curve = False, draw_lightning_curve = True, draw_landscapes = False, draw_box_for_cohom_frac = False, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type, more = more )
    
    draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, draw_CT_curve = False, draw_lightning_curve = False, draw_jordan_curve = True, draw_landscapes = True, draw_box_for_cohom_frac = False, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type, more = more )
    


    ## draw many:
    # start_num = 50
    # end_num = 1000
    # draw_cannon_thurston_from_veering_isosigs_file('Data/veering_census.txt', 'Images/Jordan_curve', max_num_tetrahedra = max_num_tetrahedra, max_length = max_length, interval_to_draw = (start_num, end_num), draw_args = draw_args, build_type = build_type)
    
    ### jigsaws

    # draw_jigsaw_from_veering_isosigs_file('Data/veering_for_jigsaws.txt', 'Images/Jigsaw', num_to_draw = 2)
    # draw_jigsaw_from_veering_isosigs_file('Data/veering_for_jigsaws.txt', 'Images/Jigsaw')
    # draw_jigsaw_from_veering_isosigs_file('Data/layered_one_cusp.txt', 'Images/Jigsaw', jigsaw_data_out_filename = "jigsaw_data_layered_one_cusp.pkl", num_to_draw = 500)
    # draw_jigsaw_from_veering_isosigs_file('Data/sigs_for_jigs_no_symm.txt', 'Images/Jigsaw', jigsaw_data_out_filename = "jigsaw_data_no_symm.pkl", max_length = max_length, num_to_draw = 876)  # all up through n's is 876. The 281th has trouble developing


