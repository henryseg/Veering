from basic_math import matrix, CP1

def developed_position(A1, A2, A3, z): #use Feng's "solving Thurston's equations in a commutative ring"
    # print A1, A2, A3, z
    epsilon = 0.00000001
    a1,b1 = A1 ## eg 1, 0
    a2,b2 = A2 ## eg 0, 1
    a3,b3 = A3 ## eg 1, 1
    Xinv = matrix((a1, a2, b1, b2))
    assert abs(Xinv.det()) > epsilon, Xinv
    X = Xinv.inverse()
    a3p, b3p = X * A3
    assert abs(b3p) > epsilon and abs(a3p) > epsilon, X * A3
    out = Xinv * CP1((-z/b3p, -1/a3p))
    return out.preferred_rep() 

unknown_vert_to_known_verts_ordering = {0:(3, 2, 1), 1:(2, 3, 0), 2:(1, 0, 3), 3:(0, 1, 2)}
### if we don't know the position of vert i, this gives the order of the vertices to put in for inf, 0, 1
### Note that the complex angles we have are associated to the edges 01|23
### This is consistent with the orientation given by veering 
    
        ###   0--------inf   1--------0
        ###   |`.    ,'|     |`.    ,'|     
        ###   |  ` ,'  |     |  ` ,'  |  
        ###   |  ,' .  |     |  ,' .  |   
        ###   |,'    `.|     |,'    `.|
        ###   1--------z     2--------3   

def develop_verts_pos(vertex_posns, gluing, face_vertex, tet_shape):
    """Get vert posns for a new tetrahedron given vert posns for an existing tetrahedron and gluing info"""
    next_tet_vertex_posns = []
    current_vertex_nums = [0, 1, 2, 3]
    current_vertex_nums.remove(face_vertex)
    next_tet_vertex_posns = [None, None, None, None]
    for i in current_vertex_nums:
        next_tet_vertex_posns[gluing[i]] = vertex_posns[i]
    #get shape of next tet
    next_tet_new_vertex_num = gluing[face_vertex] #the number of the vertex we don't have the position of
    ordering = unknown_vert_to_known_verts_ordering[next_tet_new_vertex_num]
    new_position = developed_position(next_tet_vertex_posns[ordering[0]], next_tet_vertex_posns[ordering[1]], next_tet_vertex_posns[ordering[2]], tet_shape)
    next_tet_vertex_posns[next_tet_new_vertex_num] = new_position
    return next_tet_vertex_posns

