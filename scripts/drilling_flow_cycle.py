from build_continent import make_continent_drill_flow_cycle
from draw_drilled_tetrahedra import draw_drilled_tetrahedra
from ordered_rectangles import build_tetrahedron_rectangle_orderings, sanity_check, build_drilled_triangulation_data

from veering.taut import is_taut, isosig_from_tri_angle
from veering.veering_tri import is_veering
from veering.flow_cycles import generate_flow_cycles
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

def drill_flow_cycle(veering_isosig, flow_cycle, return_tri_angle = False, draw_rectangles = False, return_found_parallel = False, return_cusp_mapping = False, use_untwisted_speed_up = True, verbose = 0):
    out = make_continent_drill_flow_cycle(veering_isosig, flow_cycle, use_untwisted_speed_up = use_untwisted_speed_up, verbose = verbose)
    if out == None:  ### flow cycle is boundary parallel
        # print('flow_cycle is boundary parallel')
        return None
    con, tetrahedra_cusp_orders, tetrahedra_chunks, intervals_inside_tet_rectangles, _, _, _, _, found_parallel = out
    old_tet_rectangles = build_tetrahedron_rectangle_orderings(con, tetrahedra_cusp_orders, tetrahedra_chunks)
    new_tetrahedra, new_faces = build_drilled_triangulation_data(old_tet_rectangles)
    
    if draw_rectangles:
        name = veering_isosig + '' + str(flow_cycle) + '_drill'
        draw_drilled_tetrahedra(con, name = name, draw_vertex_numbers = False, 
        tetrahedra_cusp_orders = tetrahedra_cusp_orders, 
        intervals_inside_tet_rectangles = intervals_inside_tet_rectangles, 
        tetrahedra_chunks = tetrahedra_chunks,
        old_tet_rectangles = old_tet_rectangles)
    
    tri, angle, cusp_mapping = triangulation_data_to_tri_angle(new_tetrahedra, new_faces)

    new_sig, isom, isosig_tri = isosig_from_tri_angle(tri, angle, return_isom = True, return_Regina_tri = True)

    if return_cusp_mapping:
        cusp_mapping_isosig = []
        for i in range(tri.countVertices()):
            v = tri.vertex(i)
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
            cusp_mapping_isosig.append(cusp_ind)

    if verbose >= 1:
        print('drilled sig:', new_sig)

    out = [new_sig]
    if return_tri_angle:
        out.extend([tri, angle])
    if return_found_parallel:
        out.append(found_parallel)
    if return_cusp_mapping:
        out.append(cusp_mapping_isosig)
    return out

def drill_flow_cycles(veering_isosig, max_length = 2, monochromatic_only = False, max_length_only = False):
    cycles = generate_flow_cycles(veering_isosig, max_length = max_length, monochromatic_only = monochromatic_only, max_length_only = max_length_only)
    for cycle, num_steps_up in cycles:
        print(veering_isosig, num_steps_up, cycle, drill_flow_cycle(veering_isosig, cycle)) 






