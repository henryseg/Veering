# from draw_continent_circle import draw_continent_circle
from build_continent import make_continent_drill_flow_cycle
# from draw_drilled_tetrahedra import draw_drilled_tetrahedra
from ordered_rectangles import build_tetrahedron_rectangle_orderings, sanity_check, build_drilled_triangulation_data

from veering.taut import is_taut, isosig_from_tri_angle
from veering.veering_tri import is_veering
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
    # 
    # print(sig)

def drill_flow_cycle(veering_isosig, flow_cycle):
    con, tetrahedra_cusp_orders, tetrahedra_chunks, _, _, _, _, _ = make_continent_drill_flow_cycle(veering_isosig, flow_cycle, verbose = 0)
    new_tetrahedra, new_faces = build_drilled_triangulation_data(con, tetrahedra_cusp_orders, tetrahedra_chunks)
    tri, angle = triangulation_data_to_tri_angle(new_tetrahedra, new_faces)
    sig = isosig_from_tri_angle(tri, angle)
    return sig


