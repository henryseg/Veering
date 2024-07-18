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
    return tri, angle

def drill_flow_cycle(veering_isosig, flow_cycle, return_tri_angle = False, draw_rectangles = False, return_found_parallel = False, verbose = 0):
    out = make_continent_drill_flow_cycle(veering_isosig, flow_cycle, verbose = verbose)
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
    
    tri, angle = triangulation_data_to_tri_angle(new_tetrahedra, new_faces)

    if return_tri_angle:
        if return_found_parallel:
            return (tri, angle, found_parallel)
        else:
            return (tri, angle)
    else:
        if return_found_parallel:
            return (isosig_from_tri_angle(tri, angle), found_parallel)
        else:
            return isosig_from_tri_angle(tri, angle)
    

def drill_flow_cycles(veering_isosig, max_length = 2, monochromatic_only = False, max_length_only = False):
    cycles = generate_flow_cycles(veering_isosig, max_length = max_length, monochromatic_only = monochromatic_only, max_length_only = max_length_only)
    for cycle, num_steps_up in cycles:
        print(veering_isosig, num_steps_up, cycle, drill_flow_cycle(veering_isosig, cycle)) 
