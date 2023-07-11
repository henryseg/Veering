#
# fundamental_domain.py
#

# Building a spanning tree for, and computation of loops, in the dual
# one-skeleton.

from .taut import liberal
from .transverse_taut import is_transverse_taut, edge_side_face_collections

verbose = 0

def spanning_dual_tree(triangulation, initial_tet_num = 0):
    """
    Returns three lists - tree faces (dual to edges in the spanning
    tree), non-tree faces (dual to edges not in the spanning tree),
    and the distance of each tet (in the tree) to the root.  We use
    the regina numbering to determine the tree.
    """
    explored_tetrahedra = [initial_tet_num]
    distances_to_root = [None]*triangulation.countTetrahedra()
    distances_to_root[initial_tet_num] = 0

    frontier_tet_faces = [(initial_tet_num, 0), (initial_tet_num, 1), (initial_tet_num, 2), (initial_tet_num, 3)]
    tree_faces = []
    non_tree_faces = []
    
    while len(frontier_tet_faces) > 0:
        my_tet_num, my_face_num = frontier_tet_faces.pop()
        my_tet = triangulation.tetrahedron(my_tet_num)
        my_face = my_tet.face(2, my_face_num)
        neighbour_tet = my_tet.adjacentSimplex(my_face_num)
        neighbour_face_num = my_tet.adjacentGluing(my_face_num)[my_face_num]
        if neighbour_tet.index() in explored_tetrahedra:
            non_tree_faces.append(my_face.index())
            frontier_tet_faces.remove((neighbour_tet.index(), neighbour_face_num))
        else:
            tree_faces.append(my_face.index())
            frontier_tet_faces.extend( [(neighbour_tet.index(), i) for i in range(4) if i != neighbour_face_num] )
            explored_tetrahedra.append(neighbour_tet.index())
            distances_to_root[neighbour_tet.index()] = distances_to_root[my_tet.index()] + 1
    tree_faces.sort()
    non_tree_faces.sort()
    if verbose > 0:
        print(("tree and non-tree faces", (tree_faces, non_tree_faces)))
    return (tree_faces, non_tree_faces, distances_to_root)


def walk_towards_root(tet, tree_faces, distances_to_root):
    """
    Returns the tet closer to the root, and the face we crossed.
    """
    my_dist_to_root = distances_to_root[tet.index()]
    assert my_dist_to_root > 0 # if we are at the root, blow up
    for i in range(4):
        if tet.triangle(i).index() in tree_faces:
            neighbour = tet.adjacentSimplex(i)
            if distances_to_root[neighbour.index()] < my_dist_to_root:
                return (neighbour, tet.triangle(i).index())
    assert False

    
def non_tree_face_loops(triangulation, include_tetrahedra = False):
    """
    For every non-tree face f, find the corresponding loop in T union
    f, given as a sequence of face indices, starting with f.  If
    include_tetrahedra, we also return the interstitial tetrahedra.
    In incidence order, we have f0,t0,f1,t1,...
    """
    tree_faces, non_tree_faces, distances_to_root = spanning_dual_tree(triangulation)
    loops = []
    for orig_face_ind in non_tree_faces:
        embeddings = triangulation.triangle(orig_face_ind).embeddings()
        tet0, tet1 = (embed.simplex() for embed in embeddings)
        zero_faces = []
        one_faces = []
        zero_tetrahedra = [tet0]
        one_tetrahedra = [tet1]
        while tet0.index() != tet1.index():
            dist0 = distances_to_root[tet0.index()]
            dist1 = distances_to_root[tet1.index()]
            if dist0 < dist1:
                tet1, face1_ind = walk_towards_root(tet1, tree_faces, distances_to_root)
                one_faces.append(face1_ind)
                one_tetrahedra.append(tet1)
            else:
                tet0, face0_ind = walk_towards_root(tet0, tree_faces, distances_to_root)
                zero_faces.append(face0_ind)
                zero_tetrahedra.append(tet0)
        zero_faces.reverse()
        zero_tetrahedra.reverse()
        face_loop = [orig_face_ind] + one_faces + zero_faces
        tet_loop = one_tetrahedra[:-1] + zero_tetrahedra
        if include_tetrahedra:
            loops.append((face_loop, tet_loop))
        else:
            loops.append(face_loop)
    return loops


@liberal
def non_tree_face_loops_oriented(tri, angle):
    """
    Suppose that T is the (dual) spanning tree.  For every non-tree
    face f, there is a unique oriented loop in T union e_f whose
    orientation agrees with the orientation of f (as given by the
    transverse taut structure).  We return this as a list of faces
    (starting with f) and a list of signs (starting with +1).
    """
    loops = non_tree_face_loops(tri)
    oriented_loops = []
    all_signs = []
    coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    assert coorientations != False
    for loop in loops:
        if len(loop) == 1:
            all_signs.append([1])
        elif len(loop) == 2: 
            embeddings = tri.triangle(loop[0]).embeddings()
            tet0, _ = (embed.simplex() for embed in embeddings)
            tet0_faces = [tet0.triangle(i).index() for i in range(4)]
            local_loop0 = tet0_faces.index(loop[0])
            local_loop1 = tet0_faces.index(loop[1])
            tet0_coor = coorientations[tet0.index()]
            if tet0_coor[local_loop0] == tet0_coor[local_loop1]:
                all_signs.append([1, -1])
            else:
                all_signs.append([1, 1])
        else: # at least three
            # first decide if we need to reverse the loop
            embeddings = tri.triangle(loop[0]).embeddings()
            tet0, tet1 = (embed.simplex() for embed in embeddings)
            tet0_faces = [tet0.triangle(i).index() for i in range(4)]
            local_loop0 = tet0_faces.index(loop[0])
            tet0_coor = coorientations[tet0.index()]
            if tet0_coor[local_loop0] != 1:
                loop.reverse()
                loop = loop[-1:] + loop[:-1]
                tet0, tet1 = tet1, tet0
            # now everything is sensible: that is, tet0 is below face0
            loop_signs = []
            curr_tet = tet0
            for face_ind in loop:
                curr_tet_faces = [curr_tet.triangle(i).index() for i in range(4)]
                local_face_ind = curr_tet_faces.index(face_ind)
                curr_tet_coor = coorientations[curr_tet.index()]
                loop_signs.append(curr_tet_coor[local_face_ind])
                curr_tet = curr_tet.adjacentSimplex(local_face_ind)
            assert loop_signs[0] == 1
            all_signs.append(loop_signs)

        oriented_loops.append(loop)

    return (oriented_loops, all_signs)

@liberal
def non_tree_face_cycles(tri, angle):
    """
    Returns the one-cycles associated to the non-tree faces.

    In other words: Suppose that T is the dual spanning tree.  Then,
    for every non-tree face f, we return a vector of length 2*n (num
    faces) with entries from {-1, 0, +1}, which tells us how the
    transverse orientation (angle) agrees or disagrees with the unique
    oriented loop in T union e_f.
    """
    
    cycles = []
    oriented_loops, all_signs = non_tree_face_loops_oriented(tri, angle)
    n = tri.countTetrahedra()
    for loop, signs in zip(oriented_loops, all_signs):
        cycle = [0]*2*n
        for face, sign in zip(loop, signs):
            cycle[face] = sign
        cycles.append(cycle)
    return cycles
