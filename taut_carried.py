#
# taut_carried.py
#

# functions for working with surfaces carried by the two-skeleton of a 
# transverse taut triangulation.

import regina # needed inside of imported files
from taut import liberal
from transverse_taut import is_transverse_taut, convert_tetrahedron_coorientations_to_faces
from taut_homology import edge_side_face_collections

@liberal
def boundary_cycles_from_surface(tri, angle, surface, tet_vert_coorientations = None):
    """ 
    Takes a carried surface. For each cusp of tri, look at the boundary curve
    of the surface on the boundary torus for that cusp. Push it up slightly, record 
    which faces of tri it goes through.
    """

    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    # set up output vectors
    out = []
    for vertex in tri.faces(0):  # 0 is the dimension of the face, so this is cusps
        out.append([0] * tri.countFaces(2))

    edge_sides = edge_side_face_collections(tri, angle, tet_vert_coorientations = tet_vert_coorientations)

    face_coorientations = convert_tetrahedron_coorientations_to_faces(tri, tet_vert_coorientations)

    for f in tri.faces(2):
        weight = int(surface[f.index()])
        for vert_num in range(3):
            cusp_index = f.face(0, vert_num).index()

            sgn = face_coorientations[f.index()]
            trailing_vert_num = (vert_num + sgn) % 3
            trailing_edge = f.face(1, trailing_vert_num) 
            leading_vert_num = (vert_num - sgn) % 3
            leading_edge = f.face(1, leading_vert_num) 
            ## these are with respect to right hand rule on surface as viewed from above
            ## the faces above us on the trailing edge get -weight, 
            ## the faces above us on the leading edge get +weight.

            left, right = edge_sides[trailing_edge.index()]
            pair = (f.index(), trailing_vert_num)
            if pair in left:
                side = left
            else:
                assert pair in right
                side = right
            above = side[side.index(pair) + 1:]
            for (f_ind, vert) in above:
                out[cusp_index][f_ind] -= weight

            left, right = edge_sides[leading_edge.index()]
            pair = (f.index(), leading_vert_num)
            if pair in left:
                side = left
            else:
                assert pair in right
                side = right
            above = side[side.index(pair) + 1:]
            for (f_ind, vert) in above:
                out[cusp_index][f_ind] += weight
    return out



