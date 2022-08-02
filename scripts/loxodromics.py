#
# loxodromics.py
#

from veering.basic_math import CP1, move_in_PSL

from develop_ideal_hyperbolic_tetrahedra import develop_verts_pos


def mob_tsfm_from_tets(init_verts_pos, final_verts_pos):
    a,b,c,d = init_verts_pos
    p,q,r,s = final_verts_pos
    M = move_in_PSL(a,b,c,p,q,r)
    assert M(d).is_close_to(s)
    return M

def get_tet_above(vt, tet_num, verts_pos):
    """
    given a veering triang with shapes, and a tet with its vert posns,
    go up from the tet once, and return the upper tet and its vert
    posns
    """
    coorientations = vt.coorientations
    tet_shapes = vt.tet_shapes

    coors = coorientations[tet_num]
    top_verts = []
    bottom_verts = []
    for i in range(4):
        if coors[tet_num][i] == 1:
            bottom_verts.append(i)
        else:
            top_verts.append(i)

    trailing_vert = bottom_verts[0]
    leading_vert = bottom_verts[1]
    current_tet = vt.tri.tetrahedron(tet_num)

    while True:
        next_tet = current_tet.adjacentTetrahedron(trailing_vert)
        gluing = current_tet.adjacentGluing(trailing_vert)
        verts_pos = develop_verts_pos(verts_pos, gluing, trailing_vert, tet_shapes[next_tet.index()])

        trailing_vert, leading_vert = gluing[leading_vert], gluing[trailing_vert]

        current_tet = next_tet

        assert coorientations[current_tet.index()][leading_vert] == -1
        if coorientations[current_tet.index()][trailing_vert] == -1: ### we are now at the top of the edge
            return (current_tet.index(), verts_pos)

def loxodromic_from_tet(vt, init_tet_num, init_verts_pos = None):
    """
    given a veering triang with shapes, and a tet, go up from the tet
    until you find it again, return the Mob tsfm
    """
    current_tet_num = init_tet_num
    if init_verts_pos == None:
        init_verts_pos = [CP1((1,0)), CP1((0,1)), CP1((1,1)), CP1((vt.tet_shapes[tet_num],1))]
    current_verts_pos = init_verts_pos
    while True:
        current_tet_num, current_verts_pos = get_tet_above(vt, current_tet_num, current_verts_pos)
        if current_tet_num == init_tet_num:
            return mob_tsfm_from_tets(init_verts_pos, current_verts_pos)

def loxodromic_from_flag(vt, tet_num, face_vert, edge_vert, init_verts_pos = None):
    """
    given a veering triang with shapes and a flag in a tet, walk
    around it until we can call the other function
    """
    coorientations = vt.coorientations
    tet_shapes = vt.tet_shapes
    if init_verts_pos == None:
        verts_pos = [CP1((1,0)), CP1((0,1)), CP1((1,1)), CP1((vt.tet_shapes[tet_num],1))]
    else:
        verts_pos = init_verts_pos
    leading_vert = face_vert
    trailing_vert = edge_vert
    current_tet = vt.tri.tetrahedron(tet_num)

    while coorientations[current_tet.index()][leading_vert] != coorientations[current_tet.index()][trailing_vert]:
        next_tet = current_tet.adjacentTetrahedron(trailing_vert)
        gluing = current_tet.adjacentGluing(trailing_vert)
        verts_pos = develop_verts_pos(verts_pos, gluing, trailing_vert, tet_shapes[next_tet.index()])
        trailing_vert, leading_vert = gluing[leading_vert], gluing[trailing_vert]
        current_tet = next_tet

    return loxodromic_from_tet(vt, current_tet.index(), init_verts_pos = verts_pos)
