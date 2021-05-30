#
# fundamental_domain.py
#

# Computation of loops in pi_1(M) by building a fundamental domain realised as a spanning tree in the dual 1 skeleton.

verbose = 0

def spanning_dual_tree(triangulation):
    """
    Returns three lists - (dual) edges in the spanning tree, (dual)
    edges not in the spanning tree, and the distance of each tet (in
    the tree) to the root. We use the regina numbering to determine
    the tree.
    """
    explored_tetrahedra = [0]
    distances_to_root = [0] + [None]*(triangulation.countTetrahedra() - 1)
    frontier_tet_faces = [(0, 0), (0, 1), (0, 2), (0, 3)]
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
    returns the tet closer to the root, and the face we crossed
    """
    my_dist_to_root = distances_to_root[tet.index()]
    assert my_dist_to_root > 0 # if we are at the root, blow up
    for i in range(4):
        if tet.triangle(i).index() in tree_faces:
            neighbour = tet.adjacentSimplex(i)
            if distances_to_root[neighbour.index()] < my_dist_to_root:
                return (neighbour, tet.triangle(i).index())
    assert False

    
def non_tree_edge_loops(triangulation):
    """
    For every non-tree (dual) edge e, find the corresponding loop in T union e.
    """
    tree_faces, non_tree_faces, distances_to_root = spanning_dual_tree(triangulation)
    loops = []
    for orig_face_ind in non_tree_faces:
        embeddings = triangulation.triangle(orig_face_ind).embeddings()
        tet0, tet1 = (embed.simplex() for embed in embeddings)
        zero_faces = []
        one_faces = []
        while tet0.index() != tet1.index():
            dist0 = distances_to_root[tet0.index()]
            dist1 = distances_to_root[tet1.index()]
            if dist0 < dist1:
                tet1, face1_ind = walk_towards_root(tet1, tree_faces, distances_to_root)
                one_faces.append(face1_ind)
            else:
                tet0, face0_ind = walk_towards_root(tet0, tree_faces, distances_to_root)
                zero_faces.append(face0_ind)
        zero_faces.reverse()
        loops.append([orig_face_ind] + one_faces + zero_faces)
    return loops