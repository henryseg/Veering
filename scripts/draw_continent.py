
import pyx ### vector graphics 
import cmath

from file_io import parse_data_file, read_from_pickle, output_to_pickle
from taut import isosig_to_tri_angle
from veering import veering_triangulation
from continent import continent
from boundary_triangulation import boundary_triangulation, tet_face
from snappy_util import shapes


# def pre_draw_transformation( z, ladder_holonomy ):
    # return z/ladder_holonomy

def draw_path(canv, path_C, draw_options, fill = False, closed = False):
    p = pyx.path.path( pyx.path.moveto(path_C[0].real, path_C[0].imag) )
    for coord in path_C[1:]:  ### this is how path drawing works...
        p.append( pyx.path.lineto(coord.real, coord.imag) )
    if closed:
        p.closepath()
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
    #### result is oriented not in terms of order of a and b, but in order of the dividers
    print(len(dividers))
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
    assert False, 'no correct divider list?'

### function to get lightning curve to the left or right of a ladderpole edge. look for tet_faces you need near the edge. Calculate any needed offsets then and there.


def explore_side(vt, tet, face, pivot_vert, leading_vert, trailing_vert, ladderpole_colour):
    tet_faces_pivots = []
    while True:
        tet_faces_pivots.append( (tet_face(vt, tet.index(), face), pivot_vert) )
        if vt.get_edge_between_verts_colour(tet.index(), (pivot_vert, leading_vert)) != ladderpole_colour:
            break
        gluing = tet.adjacentGluing(trailing_vert)
        tet = tet.adjacentTetrahedron(trailing_vert)
        face = gluing[face]
        pivot_vert = gluing[pivot_vert]
        leading_vert, trailing_vert = gluing[trailing_vert], gluing[leading_vert]
    tet_faces_pivots.reverse() ## seems to play nice with continent coast order
    return tet_faces_pivots

def lightning_curve_for_ladder_unit(dividers, lu, offset):
    non_inf_verts = [0,1,2,3]
    non_inf_verts.remove(lu.face)
    edge_colours = []
    for i in range(3):
        edge_colours.append( lu.ladder.vt.get_edge_between_verts_colour(lu.tet_num, (non_inf_verts[(i+1)%3], non_inf_verts[(i+2)%3]) ) )
    for i in range(3):
        if edge_colours[i] == edge_colours[(i+1)%3]:
            odd_one_out = (i+2)%3
    s, e = lu.con_verts[non_inf_verts[(odd_one_out + 1)%3]], lu.con_verts[non_inf_verts[(odd_one_out + 2)%3]] 

    T = lu.ladder.torus_triangulation
    ladder_parity = T.ladder_list.index(lu.ladder) % 2
    lightning_curve = lightning_curve_from_dividers(dividers[ladder_parity], s, e, special_vertices = [], spiky = False)  ### no need for extra ladderpole vertices - there are none between these two
    print('len lightning_curve', len(lightning_curve))
    lightning_curve = [c + offset for c in lightning_curve]
    return lightning_curve
    # lightning_curve_scaled = [ T.drawing_scale * (c + offset) for c in lightning_curve ]
    # draw_path(T.canv, lightning_curve_scaled, [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round, lightning_colours[ladder_parity%2]])

def lightning_curves_around_ladder_unit_ladder_pole_edge(dividers, lu, T):
    ### describe a ladderpole edge by a neighbouring ladder unit. Draw lightning curves around it
    vt = T.vt
    tet_num, face = lu.tet_num, lu.face
    if lu.is_on_left():
        opp_vert = lu.right_vertices[0]
        ladderpole_verts = lu.left_vertices
    else:
        opp_vert = lu.left_vertices[0]
        ladderpole_verts = lu.right_vertices

    ladderpole_verts_complex = [lu.verts_pos[i].complex() for i in ladderpole_verts]

    ladderpole_colour = vt.get_edge_between_verts_colour(tet_num, ladderpole_verts)
    if (ladderpole_colour == "red" and lu.is_on_left()) or (ladderpole_colour == "blue" and lu.is_on_right()):
        pivot_vert, trailing_vert = ladderpole_verts
    else:
        trailing_vert, pivot_vert = ladderpole_verts 
        ladderpole_verts_complex.reverse()
        
    leading_vert = opp_vert  

    tet = T.vt.tri.tetrahedron(tet_num)
    gluing = tet.adjacentGluing(opp_vert)
    nb_tet = tet.adjacentTetrahedron(opp_vert)
    nb_face = gluing[face]
    nb_pivot_vert, nb_trailing_vert = gluing[trailing_vert], gluing[pivot_vert] ### swap
    nb_leading_vert = gluing[leading_vert]

    ### if the edge is blue, rotate up to the right to find relevant tet_faces to its right, and down to the left.
    ### vice versa for red ladderpole edge

    our_side_tet_faces_pivots = explore_side(vt, tet, face, pivot_vert, leading_vert, trailing_vert, ladderpole_colour)
    nb_side_tet_faces_pivots = explore_side(vt, nb_tet, nb_face, nb_pivot_vert, nb_leading_vert, nb_trailing_vert, ladderpole_colour)

    ### for each of these tet_faces, use the edge for the corresponding ladder unit in T, translated as necessary using the offset
    ### from either end of the ladderpole edge.

    if (ladderpole_colour == "red" and lu.is_on_left()) or (ladderpole_colour == "blue" and lu.is_on_right()):
        all_tet_faces_pivots = [our_side_tet_faces_pivots, nb_side_tet_faces_pivots]
    else:
        all_tet_faces_pivots = [nb_side_tet_faces_pivots, our_side_tet_faces_pivots]

    out = []
    for k in range(2):
        full_crv = []
        for j, (tf, pivot_vert) in enumerate(all_tet_faces_pivots[k]):
            lu = T.bt.get_ladder_unit_from_tet_face(tf)
            offset = ladderpole_verts_complex[k] - lu.verts_pos[pivot_vert].complex()
            crv = lightning_curve_for_ladder_unit(dividers, lu, offset)
            if j > 0:
                crv = crv[1:]
            full_crv.extend(crv)
        out.append(full_crv)
    return out

def uniquify_list(dup_list, subtract = [], epsilon = 0.001):
    unique_list = []
    for z in dup_list:
        keep = True
        for w in unique_list + subtract:
            if abs(z-w)<0.001:
                keep = False
                break
        if keep:
            unique_list.append(z)
    return unique_list

def replace_with_continent_vertices(v_list, con, epsilon = 0.001):
    for i, w in enumerate(v_list):
        for v in con.boundary_triangulation_vertices:
            # print(type(v), type(w))
            if abs(v.pos.complex() - w) < epsilon:
                v_list[i] = v    
                break    

def assign_continent_vertices_to_tet_faces(T, con):
    for L in T.ladder_list:
        for lu in L.ladder_unit_list:
            lu.con_verts = [v.complex() for v in lu.verts_pos]
            replace_with_continent_vertices(lu.con_verts, con)

def draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, max_length, output_filename = None, draw_args = None, build_type = None ):
    draw_CT_curve, draw_lightning_curve, draw_jordan_curve = draw_args['draw_CT_curve'], draw_args['draw_lightning_curve'], draw_args['draw_jordan_curve']
    draw_dividers, draw_landscapes, draw_box_for_cohom_frac = draw_args['draw_dividers'], draw_args['draw_landscapes'], draw_args['draw_box_for_cohom_frac']
    draw_alignment_dots, draw_desired_vertices, expand_fund_dom = draw_args['draw_alignment_dots'], draw_args['draw_desired_vertices'], draw_args['expand_fund_dom']

    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle, tet_shapes = tet_shapes)
    B = boundary_triangulation(vt)
    B.generate_canvases(args = draw_args)

    out_data = []

    for i,T in enumerate(B.torus_triangulation_list):
        print(T.sideways_index_shift)
        print(T.sideways_once_holonomy)
        print(T.sideways_holonomy)
        print(T.ladder_holonomy)
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

        ladderpoles_vertices = T.left_ladder_pole_vertices() ### everything is on left of the ladders...

        if expand_fund_dom:
            left_ladder = T.ladder_list[0]
            right_ladder = T.ladder_list[-1]

            left_ladder_nbd = []
            for lu in left_ladder.ladder_unit_list:
                if lu.is_on_left():
                    # print(lu, lu.left_vertices[0])
                    left_ladder_nbd.extend( lu.vert_positions_around_corner( lu.left_vertices[0] ) )
            last = left_ladder.ladder_unit_list[-1]
            # print('extra', last, last.left_vertices[-1])
            left_ladder_nbd.extend( last.vert_positions_around_corner( last.left_vertices[-1] ) )
            # print('right')
            right_ladder_nbd = []
            for lu in right_ladder.ladder_unit_list:
                if lu.is_on_left():
                    # print(lu, lu.left_vertices[0])
                    right_ladder_nbd.extend( lu.vert_positions_around_corner( lu.left_vertices[0] ) )
            last = right_ladder.ladder_unit_list[-1]
            # print('extra', last, last.left_vertices[-1])
            right_ladder_nbd.extend( last.vert_positions_around_corner( last.left_vertices[-1] ) )
            
            left_ladder_nbd = uniquify_list(left_ladder_nbd)
            right_ladder_nbd = uniquify_list(right_ladder_nbd)

            bottom_nbd = []
            top_nbd = []
            for i, L in enumerate(T.ladder_list):
                if i%2 == 0:
                    lu_bottom = L.ladder_unit_list[0]
                    lu_top = L.ladder_unit_list[-1]

                    bottom_nbd.extend( lu_bottom.vert_positions_around_corner(lu_bottom.left_vertices[0]) )
                    bottom_nbd.extend( lu_bottom.vert_positions_around_corner(lu_bottom.right_vertices[0]) )
                    top_nbd.extend( lu_top.vert_positions_around_corner(lu_top.left_vertices[-1]) )
                    top_nbd.extend( lu_top.vert_positions_around_corner(lu_top.right_vertices[-1]) )
            bottom_nbd = uniquify_list(bottom_nbd)
            top_nbd = uniquify_list(top_nbd)

            all_ladderpole_vertices = [v for L in ladderpoles_vertices for v in L]
            left_ladder_nbd = uniquify_list(left_ladder_nbd, subtract = all_ladderpole_vertices)
            right_ladder_nbd = uniquify_list(right_ladder_nbd, subtract = all_ladderpole_vertices + left_ladder_nbd)  ### right and left can overlap
            bottom_nbd = uniquify_list(bottom_nbd, subtract = all_ladderpole_vertices + left_ladder_nbd + right_ladder_nbd) ### top and bottom cannot
            top_nbd = uniquify_list(top_nbd, subtract = all_ladderpole_vertices + left_ladder_nbd + right_ladder_nbd)

            desired_vertices = all_ladderpole_vertices + left_ladder_nbd + right_ladder_nbd + bottom_nbd + top_nbd
            desired_vertices = uniquify_list(desired_vertices) ## just in case...
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
        for ladderpole_vertices in ladderpoles_vertices:
            replace_with_continent_vertices(ladderpole_vertices, con)    
        if expand_fund_dom:
            replace_with_continent_vertices(left_ladder_nbd, con)
            replace_with_continent_vertices(right_ladder_nbd, con)
            replace_with_continent_vertices(bottom_nbd, con)
            replace_with_continent_vertices(top_nbd, con)
            nbd = left_ladder_nbd + bottom_nbd + top_nbd + right_ladder_nbd
            nbd.sort(key = lambda v: con.coast.index(v))

        ### the following is the list with correctly replaced vertices
        all_ladderpole_vertices = [v for L in ladderpoles_vertices for v in L] ### dont need this to be bigger when doing "expand_fund_dom" 

        ladderpole_descendant_segments = []
        if expand_fund_dom:  
            ### we must sort the extra vertices of nbd into the correct ladderpoles... just do by colour of the edge connecting to infty
            ladderpole_is_red = nbd[0].edge_between(con.infinity).is_red
            segment_start = con.coast.index(nbd[0])
            for i,v in enumerate(nbd):
                if v.edge_between(con.infinity).is_red != ladderpole_is_red:
                    ladderpole_descendant_segments.append( [segment_start, con.coast.index(nbd[i-1])] )
                    segment_start = con.coast.index(nbd[i])
                    ladderpole_is_red = not ladderpole_is_red
            ladderpole_descendant_segments.append( [segment_start, con.coast.index(nbd[-1])] )
        else:
            for ladderpole_vertices in ladderpoles_vertices:
                segment = [con.coast.index(ladderpole_vertices[0]), con.coast.index(ladderpole_vertices[-1])]
                segment.sort()
                ladderpole_descendant_segments.append( segment )

        con.mark_ladderpole_descendants(ladderpole_descendant_segments)

        # print 'important verts', important_vertices
        # for v in important_vertices:
        #     z = T.drawing_scale * v.pos.complex()
        #     pyx_fill_col = pyx.deco.filled([pyx.color.rgb.black])
        #     T.canv.fill(pyx.path.circle(z.real, z.imag, 0.02), [pyx_fill_col])

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

        assign_continent_vertices_to_tet_faces(T, con)
        dividers = con.lightning_dividers([])  ## only lightning curves for infinity
        layer1 = T.canv.layer("layer1")
        layer2 = T.canv.layer("layer2", below = "layer1")

        if draw_jordan_curve:
            jordan_colours = [pyx.color.rgb(0,0,0), pyx.color.rgb(1,1,1)]  ### black, white   


            for k, L in enumerate(T.ladder_list):
                if k%2 == 0:
                    for lu in L.ladder_unit_list:
                        if lu.is_on_left():
                            s, e = lu.con_verts[lu.left_vertices[0]], lu.con_verts[lu.left_vertices[1]]
                        else:
                            s, e = lu.con_verts[lu.right_vertices[0]], lu.con_verts[lu.right_vertices[1]]
                        CT_ladderpole = con.segment_between(s, e)
                        print(lu, T)
                        curves = lightning_curves_around_ladder_unit_ladder_pole_edge(dividers, lu, T)
                        if CT_ladderpole[0] != e:
                            CT_ladderpole.reverse()
                        assert CT_ladderpole[0] == e
                        assert CT_ladderpole[-1] == s

                        CT_ladderpole_pos = [v.pos.complex() for v in CT_ladderpole ]
                        assert abs(CT_ladderpole_pos[0] - curves[0][-1]) < 0.0001
                        assert abs(CT_ladderpole_pos[-1] - curves[0][0]) < 0.0001
                        loop0 = curves[0] + CT_ladderpole_pos[1:-1]
                        path_C = [ T.drawing_scale * c for c in loop0 ]
                        draw_path(layer1, path_C, [jordan_colours[0]], fill = True)  
                        assert abs(CT_ladderpole_pos[0] - curves[1][0]) < 0.0001
                        assert abs(CT_ladderpole_pos[-1] - curves[1][-1]) < 0.0001
                        CT_ladderpole_pos.reverse()
                        loop1 = curves[1] + CT_ladderpole_pos[1:-1]
                        path_C = [ T.drawing_scale * c for c in loop1 ]
                        draw_path(layer1, path_C, [jordan_colours[1]], fill = True)  
                        

            # for k, ladderpole in enumerate(ladderpoles_vertices):
            #     for j in range(len(ladderpole) - 1):
            #         s, e = ladderpole[j], ladderpole[j+1]
            #         CT_ladderpole = con.segment_between(s, e)
            #         if CT_ladderpole[0] != s:
            #             CT_ladderpole.reverse()
            #         assert CT_ladderpole[0] == s
            #         assert CT_ladderpole[-1] == e
            #         for i in range(2):
            #             lightning_curve = lightning_curve_from_dividers(dividers[i], s, e, special_vertices = all_ladderpole_vertices, spiky = False)
            #             CT_ladderpole_pos = [v.pos.complex() for v in CT_ladderpole ]
            #             if CT_ladderpole_pos[0] == lightning_curve[0]:
            #                 lightning_curve.reverse() 
            #             assert CT_ladderpole_pos[0] == lightning_curve[-1] and CT_ladderpole_pos[-1] == lightning_curve[0]
            #             lightning_curve = lightning_curve[1:-1]
            #             loop = CT_ladderpole_pos + lightning_curve
            #             path_C = [ T.drawing_scale * c for c in loop ]
            #             draw_path(layer1, path_C, [jordan_colours[i]], fill = True)  
                        # if k == 0:
                        #     L = T.ladder_list[0]
                        #     offset = L.left_ladderpole_index_to_ladder_unit(j).gluing_offset
                        #     path_C = [ T.drawing_scale * (c + offset) for c in loop ]
                        #     draw_path(layer1, path_C, [jordan_colours[i]], fill = True) 

            box = T.canv.bbox()
            layer2.fill(box.enlarged(1.0).rect(), [pyx.color.rgb(0.5,0.5,0.5)])

        # lightning_colours = [pyx.color.rgb(0,0.5,0), pyx.color.rgb(0.5,0,0.5)]  ### green, purple
        # lightning_colours = [pyx.color.rgb(1,0,0), pyx.color.rgb(1,0,0)]  ### red, red
        lightning_colours = [pyx.color.rgb(0,0,0), pyx.color.rgb(0,0,0)]  ### black, black   

        if draw_dividers:
            for k in range(2):
                for divider_list in dividers[k]:
                    for edge in divider_list:
                        edgeC = [T.drawing_scale * v.pos.complex() for v in edge.vertices]
                        draw_path(layer1, edgeC, [lightning_colours[k], pyx.style.linewidth(0.005)])

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

            landscape_edges = [con.lower_landscape_edges, con.upper_landscape_edges]

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

        lightning_colours = [pyx.color.rgb(0,0.5,0), pyx.color.rgb(0.5,0,0.5)]  ### green, purple
        # lightning_colours = [pyx.color.rgb(1,0,0), pyx.color.rgb(1,0,0)]  ### red, red
        # lightning_colours = [pyx.color.rgb(0,0,0), pyx.color.rgb(0,0,0)]  ### black, black   

        if draw_lightning_curve:    
            dividers = con.lightning_dividers([])  ## only lightning curves for infinity
            # dividers = con.lightning_dividers(all_ladderpole_vertices)  ## add in more lightning curves

            # k = 0
            for k, L in enumerate(T.ladder_list):
                for lu in L.ladder_unit_list:
                    if k%2 == 0:
                        curves = lightning_curves_around_ladder_unit_ladder_pole_edge(dividers, lu, T)
                        for j, crv in enumerate(curves):
                            lightning_curve_scaled = [ T.drawing_scale * c for c in crv ]
                            draw_path(T.canv, lightning_curve_scaled, [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round, lightning_colours[j]])



                    # non_inf_verts = [0,1,2,3]
                    # non_inf_verts.remove(lu.face)
                    # edge_colours = []
                    # for i in range(3):
                    #     edge_colours.append( T.vt.get_edge_between_verts_colour(lu.tet_num, (non_inf_verts[(i+1)%3], non_inf_verts[(i+2)%3]) ) )
                    # for i in range(3):
                    #     if edge_colours[i] == edge_colours[(i+1)%3]:
                    #         odd_one_out = (i+2)%3
                    # s, e = lu.con_verts[non_inf_verts[(odd_one_out + 1)%3]], lu.con_verts[non_inf_verts[(odd_one_out + 2)%3]] 

                    # lightning_curve = lightning_curve_from_dividers(dividers[k%2], s, e, special_vertices = all_ladderpole_vertices, spiky = False)
                    # lightning_curve_scaled = [ T.drawing_scale * c for c in lightning_curve ]
                    # draw_path(T.canv, lightning_curve_scaled, [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round, lightning_colours[k%2]])
                    # if k == len(T.ladder_list) - 1:
                    #     lightning_curve_shifted = [c + lu.gluing_offset for c in lightning_curve]
                    #     lightning_curve_shifted_scaled = [ T.drawing_scale * c for c in lightning_curve_shifted ]
                    #     draw_path(T.canv, lightning_curve_shifted_scaled, [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round, lightning_colours[k%2]])


               

            # for k, ladderpole in enumerate(ladderpoles_vertices):
            #     for j in range(len(ladderpole) - 1):
            #         s, e = ladderpole[j], ladderpole[j+1]
            #         for i in range(2):
            #             lightning_curve = lightning_curve_from_dividers(dividers[i], s, e, special_vertices = all_ladderpole_vertices, spiky = False)
            #             lightning_curve_scaled = [ T.drawing_scale * c for c in lightning_curve ]
            #             draw_path(T.canv, lightning_curve_scaled, [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round, lightning_colours[i]])

                        # if k == 0:
                        #     L = T.ladder_list[0]
                        #     offset = L.left_ladderpole_index_to_ladder_unit(j).gluing_offset
                        #     lightning_curve_scaled_shifted = [ T.drawing_scale * (c + offset) for c in lightning_curve ]
                        #     draw_path(T.canv, lightning_curve_scaled_shifted, [pyx.style.linewidth(ct_lw), pyx.style.linejoin.round, pyx.color.rgb(0,1,1)]) 

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
            one_vert = inf_vert - 2*(((inf_vert//2) % 2)*2 - 1) ## swaps 0 with 2, 1 with 3

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
        if draw_desired_vertices:
            for pos in desired_vertices:
                # pos = v.pos.complex()
                pos *= T.drawing_scale
                T.canv.stroke(pyx.path.circle(pos.real, pos.imag, 0.3), [pyx.style.linewidth(0.1), pyx.color.rgb.green])

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

    # shapes_data = read_from_pickle('Data/veering_shapes.pkl')
    for veering_isosig in to_draw:
        print(veering_isosig)
        # tet_shapes = shapes_data[veering_isosig]
        tri, angle = isosig_to_tri_angle(veering_sig)
        tet_shapes = shapes(tri)
        filename = output_dirname + '/' + veering_isosig + '_' + str(max_num_tetrahedra) + '_' + str(max_length) + '_' + build_type + '.pdf'
        try:
            draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )
        except: 
            print('failed to draw ' + veering_isosig)

def draw_jigsaw_from_veering_isosigs_list(veering_isosigs_list, output_dirname, jigsaw_data_out_filename = "jigsaw_data.pkl", max_num_tetrahedra = 2000000, max_length = 0.2, draw_args = None):
    build_type = 'build_long_and_mid'
    # shapes_data = read_from_pickle('Data/shapes_jig_no_symm.pkl')

    data_for_cohom_fracs = {}
    for i, veering_isosig in enumerate(veering_isosigs_list):
        # print(veering_isosig)
        if i%25 == 0:
            print(i)
        # tet_shapes = shapes_data[veering_isosig]
        tri, angle = isosig_to_tri_angle(veering_sig)
        tet_shapes = shapes(tri)
        # print 'tet_shapes', tet_shapes

        filename = output_dirname + '/' + veering_isosig + '_' + str(max_num_tetrahedra) + '_' + str(max_length) + '_' + build_type + '.pdf'
        expand_fund_dom = True
        out = draw_continent( veering_isosig, tet_shapes, max_num_tetrahedra, max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )

        if out != False:
            data_for_cohom_fracs[out[0]] = out[1:]
    output_to_pickle(data_for_cohom_fracs, jigsaw_data_out_filename)


def draw_jigsaw_from_veering_isosigs_file(veering_isosigs_filename, output_dirname, jigsaw_data_out_filename = "jigsaw_data.pkl", max_num_tetrahedra = 2000000, max_length = 0.2, interval_to_draw = None, draw_args = None):
    veering_isosigs_list = parse_data_file(veering_isosigs_filename)
    if interval_to_draw != None:
        to_draw = veering_isosigs_list[interval_to_draw[0]:interval_to_draw[1]]
    else:
        to_draw = veering_isosigs_list
    draw_jigsaw_from_veering_isosigs_list(to_draw, output_dirname, jigsaw_data_out_filename = jigsaw_data_out_filename, max_num_tetrahedra = max_num_tetrahedra, max_length = max_length, draw_args = draw_args)

