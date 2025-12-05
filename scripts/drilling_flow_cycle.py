from build_continent import make_continent_drill_flow_cycles
from draw_drilled_tetrahedra import draw_drilled_tetrahedra
from ordered_rectangles import build_tetrahedron_rectangle_orderings, sanity_check, build_drilled_triangulation_data

from veering.taut import is_taut, isosig_from_tri_angle
from veering.veering_tri import is_veering
from veering.flow_cycles import generate_flow_cycles, flow_cycle_to_dual_edge_loop
import regina 

def triangulation_data_to_tri_angle(new_tetrahedra, new_faces):
    tri = regina.Triangulation3() ### empty
    Regina_tets = []
    for new_tet in new_tetrahedra:
        Regina_tets.append(tri.newTetrahedron())

    for i, new_tet in enumerate(new_tetrahedra):
        Regina_tet = Regina_tets[i]
        for j in range(4):
            if Regina_tet.adjacentTetrahedron(j) == None: # don't glue if already glued
                other_tet, _ = new_tet.adjTetFace[j]
                other_Regina_tet = Regina_tets[other_tet.index]
                Regina_perm = regina.Perm4(new_tet.adjGluing[j])
                Regina_tet.join(j, other_Regina_tet, Regina_perm)                
    assert tri.isValid()
    assert tri.isIdeal()
    assert tri.isOriented()
    angle = [1] * len(Regina_tets)
    # print('angle struct', angle)
    # print('countBoundaryComponents', tri.countBoundaryComponents())
    assert is_veering(tri, angle)
    cusp_mapping = []
    # print('tri.countVertices()', tri.countVertices())
    for i in range(tri.countVertices()):
        # print('vertex', i)
        v = tri.vertex(i)
        cusp_ind_measurements = []
        for embed in v.embeddings():
            tet_index = embed.simplex().index()
            vert_num = embed.face()
            # print(tet_index, vert_num)
            cusp_ind_measurements.append(new_tetrahedra[tet_index].cusp_index[vert_num])
        cusp_ind = cusp_ind_measurements[0]
        # print(cusp_ind_measurements)
        assert all([cusp_ind == cusp_ind_m for cusp_ind_m in cusp_ind_measurements]) 
        cusp_mapping.append(cusp_ind)

    return tri, angle, cusp_mapping

def drill_flow_cycles(veering_isosig, flow_cycles, save_filename = None, return_isosig_tri_angle = False, return_built_tri_angle = False, return_triangulation_data = False, draw_rectangles = False, draw_tetrahedron_rectangles = False, global_scale_up = 1.0, scale_drawn_elements = 1.0, return_found_parallel = False, return_cusp_mapping = False, use_untwisted_speed_up = True, verbose = 0):
    out = make_continent_drill_flow_cycles(veering_isosig, flow_cycles, use_untwisted_speed_up = use_untwisted_speed_up, verbose = verbose)
    if out == None:  ### flow cycle is boundary parallel
        # print('flow_cycle is boundary parallel')
        return None
    con, tetrahedra_cusp_orders, tetrahedra_chunks, intervals_inside_tet_rectangles, _, _, _, _, found_parallel = out
    old_tet_rectangles = build_tetrahedron_rectangle_orderings(con, tetrahedra_cusp_orders, tetrahedra_chunks)
    new_tetrahedra, new_faces = build_drilled_triangulation_data(old_tet_rectangles)
    
    if draw_rectangles:
        if save_filename == None:
            save_filename = veering_isosig + '' + str(flow_cycles) + '_drill'
        draw_drilled_tetrahedra(con, name = save_filename, draw_vertex_numbers = True, 
        tetrahedra_cusp_orders = tetrahedra_cusp_orders, 
        intervals_inside_tet_rectangles = intervals_inside_tet_rectangles, 
        tetrahedra_chunks = tetrahedra_chunks,
        old_tet_rectangles = old_tet_rectangles,
        draw_tetrahedron_rectangles = draw_tetrahedron_rectangles,
        global_scale_up = global_scale_up,
        scale_drawn_elements = scale_drawn_elements)
    
    built_tri, built_angle, built_to_original_cusp_mapping = triangulation_data_to_tri_angle(new_tetrahedra, new_faces)

    new_sig, isom, isosig_tri, isosig_angle = isosig_from_tri_angle(built_tri, built_angle, return_isom = True, return_Regina_tri = True, return_isosig_angle = True)

    if return_cusp_mapping:
        isosig_to_built_cusp_mapping = [None] * isosig_tri.countVertices()
        for i in range(built_tri.countVertices()):
            v = built_tri.vertex(i)
            mapped_cusp_indices = []
            for embed in v.embeddings():
                tet_index = embed.simplex().index()
                vert_num = embed.face()
                mapped_tet_index = isom.tetImage(tet_index)
                mapped_vert_num = isom.facetPerm(tet_index)[vert_num]
                mapped_v = isosig_tri.tetrahedron(mapped_tet_index).vertex(mapped_vert_num)
                mapped_cusp_indices.append(mapped_v.index())
            cusp_ind = mapped_cusp_indices[0]
            assert all([cusp_ind == cusp_ind_m for cusp_ind_m in mapped_cusp_indices])
            isosig_to_built_cusp_mapping[cusp_ind] = i

        isosig_to_original_cusp_mapping = [built_to_original_cusp_mapping[ind] for ind in isosig_to_built_cusp_mapping]

    ### We could also ask which flow cycle corresponds to which cusp of the drilled manifold.
    ### This is not currently implemented.

    if verbose >= 1:
        print('drilled sig:', new_sig)

    out = [new_sig]
    if return_isosig_tri_angle:
        out.extend([isosig_tri, isosig_angle])
    if return_built_tri_angle:
        out.extend([built_tri, built_angle])
    if return_triangulation_data:
        out.extend([new_tetrahedra, new_faces])
    if return_found_parallel:
        out.append(found_parallel)
    if return_cusp_mapping:
        out.append(isosig_to_original_cusp_mapping)
    return out

def drill_flow_cycle_list(veering_isosig, max_length = 2, monochromatic_only = False, min_length = None):
    cycles = generate_flow_cycles(veering_isosig, max_length = max_length, monochromatic_only = monochromatic_only, max_length_only = max_length_only)
    for cycle, num_steps_up in cycles:
        print(veering_isosig, num_steps_up, cycle, drill_flow_cycles(veering_isosig, [cycle])) 

def drill_all_flow_cycles(veering_isosig, min_length, max_length, draw_tetrahedron_rectangles = True, global_scale_up = 1.0, scale_drawn_elements = 1.0, verbose = 0):
    cycles = generate_flow_cycles(veering_isosig, max_length = max_length, min_length = min_length)
    print('flow_cycles', cycles)
    save_filename = veering_isosig + '_drill_all_' + str(min_length) + '_' + str(max_length)
    drilled = drill_flow_cycles(veering_isosig, cycles, save_filename = save_filename, 
        draw_rectangles = True, draw_tetrahedron_rectangles = draw_tetrahedron_rectangles, 
        global_scale_up = global_scale_up, scale_drawn_elements = scale_drawn_elements, verbose = verbose)
    return drilled

def main():
    ### examples where drilling the flow cycle and drilling the geodesic give different answers:
    ### drilling cPcbbbiht_12 along ((0, 0), (0, 0), (0, 5), (0, 0), (0, 5)) gives different results jLLwQLQbeefgehiiixxxaaxxxcv [o9_40888(0,0)(0,0)] nLvALzAAQkbeffhhikjlkmmmhaihggfhujcvcf [L14n33639(0,0)(0,0)]
    ### snappy word: 'bbCabCa'

    ### drilling dLQacccjsnk_200 along ((0, 4), (2, 2), (2, 5), (1, 1)) gives different results pLLPwvAPPAQccdfejhmjklnmnooqffaakvachckcvhw [o9_43267(0,0)(0,0)] oLLzMLLzQQcaceefiljkmnlnmnjkxccnabqqarggr [o9_41941(0,0)(0,0)]

    import snappy
    from snappy_drill_homotopic import tet_and_face_indices_to_word, drill_tet_and_face_indices
    from veering.taut import isosig_to_tri_angle

    sig = 'cPcbbbiht_12'
    flow_cycle = ((0, 0), (0, 0), (0, 5), (0, 0), (0, 5))
    tri, angle = isosig_to_tri_angle(sig)
    tet_and_face_indices = flow_cycle_to_dual_edge_loop(tri, angle, flow_cycle)
    mfd = snappy.Manifold(tri)
    word = tet_and_face_indices_to_word(mfd, tet_and_face_indices)
    print(sig, flow_cycle, word)
    drilled = drill_tet_and_face_indices(mfd, tet_and_face_indices)
    print(drilled.identify())

    # output: cPcbbbiht_12 ((0, 0), (0, 0), (0, 5), (0, 0), (0, 5)) bbCabCa
    #         [L14n33639(0,0)(0,0)]





