
import regina
import veering 

def there_is_a_pi_here(angle_struct, embed):
    """given an embedding of an edge in a tetrahedron, tells us if there is a pi at that edge"""
    tet = embed.tetrahedron()
    vert_perm = embed.vertices() 
    vert_nums = [vert_perm[0], vert_perm[1]]
    vert_nums.sort()
    in_tet_edge_num = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)].index(tuple(vert_nums))
    if in_tet_edge_num <= 2:
        edgepair = in_tet_edge_num
    else:
        edgepair = 5 - in_tet_edge_num
    return angle_struct[tet.index()] == edgepair

def edge_equation_matrix_taut(triangulation, angle_struct):
    """for each edge, find the face numbers on either side of the pis, put +1 for one side and -1 for the other"""
    assert triangulation.isOriented()
    matrix = []
    for edge_num in range(triangulation.countEdges()):
        edge = triangulation.edge(edge_num)
        embeddings = edge.embeddings()

        triangle_num_sets = [[],[]]
        which_set_to_add_to = 0
        for embed in embeddings:
            tet = embed.tetrahedron()
            vert_perm = embed.vertices() 
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            #as we walk around the edge, leading is in front of us, trailing is behind us
            #see http://regina.sourceforge.net/engine-docs/classregina_1_1NTetrahedron.html#a54d99721b2ab2a0a0a72b6216b436440
            triangle_num_sets[which_set_to_add_to].append(tet.triangle(leading_vert_num).index())
            if there_is_a_pi_here(angle_struct, embed):
                which_set_to_add_to = (which_set_to_add_to + 1) % 2

        row = [0]*triangulation.countTriangles()
        for i in triangle_num_sets[0]:
            row[i] = row[i] + 1
        for i in triangle_num_sets[1]:
            row[i] = row[i] - 1
        matrix.append(row)
    return matrix



