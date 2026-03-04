#
# ordered_tri.py
#

import regina
from .taut import isosig_to_tri_angle
from .veering_tri import veering_triangulation
from .edge_orientability import is_edge_orientable
from .transverse_taut import is_transverse_taut, get_top_and_bottom_vert_nums

@liberal
def ordered_tri(tri, angle):
    """
    Given a transverse veering triangulation, create a new triangulation with vertices labelled so that the 
    triangulation is ordered in a way consistent with the upper branched surface.

    """
    vt = veering_triangulation(tri, angle)

    out = regina.Triangulation3.fromIsoSig(isosig)
    frontier_tet_nums = [0]
    done_tet_nums = []



    while len(frontier_tet_nums) > 0:
        tet_num = frontier_tet_nums.pop()
        tri_tet = vt.tri.tetrahedron(tet_num)
        out_tet = out.tetrahedron(tet_num)

        (t0,t1), (b0,b1) = get_top_and_bottom_vert_nums(vt.coorientations, tet_num) 
        ### the ordering within these pairs is arbitrary
        (ot0,ot1), (ob0,ob1) = (1,2), (0,3)
        

        ###      1
        ###      *      top edge is blue
        ###     /|`.
        ###    /  | `.
        ### 0 *---|---* 3
        ###    `. |  /
        ###      `.|/
        ###        *
        ###        2

        ###        2
        ###        *    top edge is red
        ###      ,'|\
        ###    ,' |  \
        ### 0 *---|---* 3
        ###    \  | ,'
        ###     \|,'
        ###      *
        ###      1


