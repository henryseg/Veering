#
# ordered_tri.py
#

import regina
from veering.taut import isosig_to_tri_angle, liberal
from veering.veering_tri import veering_triangulation
from veering.edge_orientability import is_edge_orientable
from veering.transverse_taut import is_transverse_taut, get_top_and_bottom_vert_nums

def is_increasing(L):
    a, b, c = L
    return a < b and b < c

def preserves_order(face_num, gluing): # does gluing along face_num with gluing permutation gluing preserve the ordering?
    vert_nums = [0,1,2,3]
    vert_nums.remove(face_num)
    new_vert_nums = [gluing[v] for v in vert_nums]
    return is_increasing(new_vert_nums)

def print_gluings(tri):
    for i in range(tri.countTetrahedra()):
        tet = tri.tetrahedron(i)
        for j in range(4):
            if tet.face(2, j).isBoundary():
                print(i,j,'boundary')
            else:
                gluing = tet.adjacentGluing(j)
                next_tet = tet.adjacentTetrahedron(j)
                print(i, j, next_tet.index(), gluing, preserves_order(j, gluing))

def get_new_top_bottom_verts(vt, tet_num):
    (t0,t1), (b0,b1) = get_top_and_bottom_vert_nums(vt.coorientations, tet_num) 
    ### the ordering within these pairs is increasing vertex number
    # print('tet num', tet_num, 'old tb', (t0, t1), (b0, b1))
    assert t0 < t1 and b0 < b1
    (ot0,ot1), (ob0,ob1) = (1,2), (0,3)
    tri_to_out_perm = regina.Perm4(t0, ot0, t1, ot1, b0, ob0, b1, ob1)
    ### The veering triangulations we generate are always orientable. These ordered triangulations reverse orientation depending on top edge colour.
    ### if the perm is orientation preserving we want the top edge to be red. 
    if (tri_to_out_perm.sign() == 1) != (vt.get_edge_between_verts_colour(tet_num, (t0, t1)) == 'red'):
        (ot0,ot1), (ob0,ob1) = (2,1), (0,3)
    tri_to_out_perm = regina.Perm4(t0, ot0, t1, ot1, b0, ob0, b1, ob1)
    assert (tri_to_out_perm.sign() == 1) == (vt.get_edge_between_verts_colour(tet_num, (t0, t1)) == 'red')
    return ((ot0,ot1), (ob0,ob1))

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

@liberal
def ordered_tri(tri, angle, verbose = 0):
    """
    Given a transverse veering triangulation, create a new triangulation with vertices labelled so that the 
    triangulation is ordered in a way consistent with the upper branched surface.

    """
    vt = veering_triangulation(tri, angle)
    assert is_edge_orientable(vt.tri, angle)
    N = tri.countTetrahedra()

    if verbose > 0:
        print_gluings(tri)
        print('top, bottom verts')
        for i in range(N):
            print(get_top_and_bottom_vert_nums(vt.coorientations, i))
        print('-')

    out = regina.Triangulation3()
    for i in range(N):
        out.newSimplex()
    
    tri_to_out_perms = ["-"] * N
    tet_num = 0
    tri_tet = vt.tri.tetrahedron(tet_num)
    out_tet = out.tetrahedron(tet_num)

    (t0,t1), (b0,b1) = get_top_and_bottom_vert_nums(vt.coorientations, tet_num) 
    ((ot0,ot1), (ob0,ob1)) = get_new_top_bottom_verts(vt, tet_num)
    tri_to_out_perms[0] = regina.Perm4(t0, ot0, t1, ot1, b0, ob0, b1, ob1)

    if verbose > 0:
        print('top, bottom', (ot0,ot1), (ob0,ob1))
        print('tri_to_out_perms[0]', tri_to_out_perms[0], tri_to_out_perms[0].sign())

    frontier_tet_face_pairs = [(0,b0), (0,b1)]  ### initially, going from tet 0 upwards
    done_tet_nums = [0]

    while len(frontier_tet_face_pairs) > 0:
        if verbose > 1:
            print('frontier_tet_face_pairs', frontier_tet_face_pairs)
        tet_num, face_num = frontier_tet_face_pairs.pop()
        if verbose > 1:
            print('chosen tet face', (tet_num, face_num))
        tri_tet = vt.tri.tetrahedron(tet_num)
        tri_gluing = tri_tet.adjacentGluing(face_num)
        next_tri_tet = tri_tet.adjacentTetrahedron(face_num)
        next_tet_num = next_tri_tet.index()
        if verbose > 1:
            print('next_tet_num, tri_gluing', next_tet_num, tri_gluing)

        out_tet = out.tetrahedron(tet_num)
        next_out_tet = out.tetrahedron(next_tet_num)

        tri_to_out_perm = tri_to_out_perms[tet_num]
        if verbose > 1:
            print('tri_to_out_perm', tri_to_out_perm)

        if not next_tet_num in done_tet_nums:
            (t0,t1), (b0,b1) = get_top_and_bottom_vert_nums(vt.coorientations, next_tet_num) 
            # ### the ordering within these pairs is arbitrary
            ((ot0,ot1), (ob0,ob1)) = get_new_top_bottom_verts(vt, next_tet_num)
            next_tri_to_out_perm = regina.Perm4(t0, ot0, t1, ot1, b0, ob0, b1, ob1)

            if verbose > 1:
                print('top, bottom', (ot0,ot1), (ob0,ob1))
            out_gluing_perm = next_tri_to_out_perm * tri_gluing * tri_to_out_perm.inverse()  # rightmost happens first
            face_vert_nums = [0, 1, 2, 3]
            face_vert_nums.remove(tri_to_out_perm[face_num])
            next_face_vert_nums = [out_gluing_perm[f] for f in face_vert_nums]
            if verbose > 1:
                print('next_tri_to_out_perm', next_tri_to_out_perm)
                print('next_face_vert_nums', next_face_vert_nums)
            if not is_increasing(next_face_vert_nums): ### then we need to rotate the next tet labels
                (ot0,ot1), (ob0,ob1) = (ot1,ot0), (ob1,ob0)
                next_tri_to_out_perm = regina.Perm4(t0, ot0, t1, ot1, b0, ob0, b1, ob1)

                if verbose > 1:
                    out_gluing_perm = next_tri_to_out_perm * tri_gluing * tri_to_out_perm.inverse()  # rightmost happens first
                    face_vert_nums = [0, 1, 2, 3]
                    face_vert_nums.remove(tri_to_out_perm[face_num])
                    next_face_vert_nums = [out_gluing_perm[f] for f in face_vert_nums]
                    print('next_tri_to_out_perm', next_tri_to_out_perm)
                    print('next_face_vert_nums', next_face_vert_nums)
                    # assert is_increasing(next_face_vert_nums)

            tri_to_out_perms[next_tet_num] = next_tri_to_out_perm
            frontier_tet_face_pairs.extend([(next_tet_num, b0), (next_tet_num, b1)])
            done_tet_nums.append(next_tet_num)

        next_tri_to_out_perm = tri_to_out_perms[next_tet_num]
        out_gluing_perm = next_tri_to_out_perm * tri_gluing * tri_to_out_perm.inverse()  # rightmost happens first
        if verbose > 1:
            print('tri_to_out_perms', tri_to_out_perms) #, [p.sign() for p in tri_to_out_perms])
            print('to glue: tet', out_tet.index(), 'along face', tri_to_out_perm[face_num], 'to tet', next_out_tet.index(), 'with perm', out_gluing_perm)
        out_tet.join(tri_to_out_perm[face_num], next_out_tet, out_gluing_perm)
        if verbose > 1:
            print_gluings(out)

    if verbose > 0:
        print('tri_to_out_perms', tri_to_out_perms, [p.sign() for p in tri_to_out_perms])
    # assert out.isIsomorphicTo(vt.tri)  ### slow!
    assert out.isOrdered()
    return out





