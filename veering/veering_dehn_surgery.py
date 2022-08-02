#
# veering_dehn_surgery.py
#

import regina

from .veering_tri import is_veering
from .taut import is_taut, isosig_from_tri_angle, isosig_to_tri_angle, vert_pair_to_edge_num


def get_mobius_strip_indices(triangulation):
    out = []
    for face in triangulation.triangles():
        if face.isMobiusBand():
            out.append(face.index())
    return out


def shares_pi_with(ang, vert):
    # Within a single tetrahedron with angle structure ang (a member
    # of {0, 1, 2}), given a vertex vert (a member of {0, 1, 2, 3}),
    # find the vertex at the other end of the pi edge from vert.
    other_part = [1, 2, 3]
    other_part.remove(ang + 1)
    partition = [[0, ang + 1], other_part]
    for part in partition:
        if vert in part:
            part.remove(vert)
            out = part.pop()
    return out


def veering_mobius_dehn_surgery(triangulation, angle_struct, face_num):
    tri = regina.Triangulation3(triangulation)  # make a copy
    angle = list(angle_struct)  # make a copy
    face = tri.triangle(face_num)
    assert face.isMobiusBand()
    # Note that dunce caps cannot appear in a veering triangulation

    # Find which vertex is on both copies of the identified edge of the face
    edges = [face.edge(i) for i in range(3)]  # edge i is opposite vertex i, i in [0, 1, 2]
    for j in range(3):
        if edges[j]==edges[(j+1)%3]:
            B = (j+2)%3
            break

    embed0 = face.embedding(0)
    embed1 = face.embedding(1)
    tet0 = embed0.tetrahedron()
    tet1 = embed1.tetrahedron()
    embed0_verts = embed0.vertices()
    embed1_verts = embed1.vertices()

    # In tet0: B gives "b".  Let "c" be the edge sharing a
    # pi with "b".  Let "d" be the vertex not meeting the given face.
    # Let "a" be the remaining vertex.

    b = embed0.vertices()[B]
    c = shares_pi_with(angle[tet0.index()], b)
    d = embed0.vertices()[3]  # ... use the face index
    a = [i for i in [0, 1, 2, 3] if i not in [b, c, d]].pop()  # ... whatever is left

    # similarly in tet1 - B gives "q".  Let "r" be the edge sharing a
    # pi with "q".  Let "s" be the vertex not meeting the given face.
    # Let "p" be the remaining vertex.

    q = embed1.vertices()[B]
    r = shares_pi_with(angle[tet1.index()], q)
    s = embed1.vertices()[3]  # ... use the face index
    p = [i for i in [0, 1, 2, 3] if i not in [q, r,s ]].pop()  # ... whatever is left
    
    # get colour of mobius strip
    pair_a = [b, c]
    pair_a.sort()
    mob_edge_a = tet0.edge( vert_pair_to_edge_num[tuple(pair_a)] )
    pair_c = [a,b]
    pair_c.sort()
    mob_edge_c = tet0.edge( vert_pair_to_edge_num[tuple(pair_c)] )
    assert mob_edge_a == mob_edge_c

    veering_colours = is_veering(tri, angle, return_type = "veering_colours")
    assert veering_colours != False  # otherwise the triangulation is not veering
    mob_colour = veering_colours[mob_edge_a.index()]

    # Now actually do the surgery
    tet0.unjoin(d)  # same as tet1.unjoin(s)
    tet_new = tri.newTetrahedron()
    if mob_colour == "red":
        tet_new.join(0, tet_new, regina.Perm4(3, 0, 1, 2))
        tet_new.join(1, tet0, regina.Perm4(c, d, a, b))
        tet_new.join(2, tet1, regina.Perm4(q, p, s, r))
    else:
        tet_new.join(1, tet_new, regina.Perm4(1, 3, 0, 2))
        tet_new.join(2, tet0, regina.Perm4(a, b, d, c))
        tet_new.join(0, tet1, regina.Perm4(s, r, p, q))

    angle.append(0)  # this is the correct taut angle for our new tetrahedron
    assert is_taut(tri, angle)
    assert is_veering(tri, angle)
    return tri, angle, tet_new.triangle(3).index()


def veering_n_mobius_dehn_surgery(tri, angle, face_num, n):
    # apply veering mobius dehn surgery n times to the same Mobius strip
    for i in range(n):
        tri, angle, face_num = veering_mobius_dehn_surgery(tri, angle, face_num)
    return tri, angle, face_num


def do_all_veering_n_surgeries(tri, angle, n = 1):
    # apply veering mobius dehn surgery n times to all Mobius strips in tri, print the results
    out = []
    print(("mob strip faces: " + str(get_mobius_strip_indices(tri))))
    for face_num in get_mobius_strip_indices(tri):
        print(("face_num", face_num))
        tri2 = regina.Triangulation3(tri)  # make a copy
        tri_surg, angle_surg, face_num_surg = veering_n_mobius_dehn_surgery(tri2, angle, face_num, n)
        sig = isosig_from_tri_angle(tri_surg, angle_surg)
        print(sig)
        out.append(sig)
    return out


def explore_mobius_surgery_graph(tri, angle, max_tetrahedra = 4):
    # apply veering mobius dehn surgery recursively
    out = []
    frontier = [(tri, angle, face_num) for face_num in get_mobius_strip_indices(tri)]
    while len(frontier) > 0:
        tri, angle, face_num = frontier.pop()
        tri_s, angle_s, face_num_s = veering_mobius_dehn_surgery(tri, angle, face_num)
        new_isosig = isosig_from_tri_angle(tri_s, angle_s)
        if new_isosig not in out:
            out.append(new_isosig)
            if tri_s.countTetrahedra() < max_tetrahedra:
                frontier.extend( [(tri_s, angle_s, face_num) for face_num in get_mobius_strip_indices(tri_s)] )
    out.sort()
    return out
