from basic_math import matrix, CP1

def developed_position(A1,A2,A3,z): #use Feng's "solving Thurston's equations in a commutative ring"
    # print A1, A2, A3, z
    a1,b1 = A1 ## eg 1,0
    a2,b2 = A2 ## eg 0,1
    a3,b3 = A3 ## eg 1,1
    assert (a1*b2-b1*a2) != complex(0,0), ",".join([str(a1), str(a2), str(b1), str(b2)]) 
    invdet = 1.0/(a1*b2-b1*a2) ## eg 1
    X = [[b2*invdet, -a2*invdet],[-b1*invdet, a1*invdet]] ## eg [[1,0],[0,1]]
    [a3p,b3p] = matrix_mult_vector(X,[a3,b3]) ## eg [1,1]
    Xinv = matrix2_inv(X)
    assert abs(b3p) > 0.00000001 and abs(a3p) > 0.00000001, str(b3p) + ' ' + str(a3p) + ' ' + str([A1,A2,A3])
    out = matrix_mult_vector(Xinv, [-z/b3p,-1/a3p]) ## eg [-z,-1] = [z,1]
    return preferred_rep(out) 

unknown_vert_to_known_verts_ordering = {0:(3,2,1), 1:(2,3,0), 2:(1,0,3), 3:(0,1,2)}
### if we don't know the position of vert i, this gives the order of the vertices to put in for inf, 0, 1
### Note that the complex angles we have are associated to the edges 01|23
### This is consistent with the orientation given by veering 
    
        ###   0--------inf   1--------0
        ###   |`.    ,'|     |`.    ,'|     
        ###   |  ` ,'  |     |  ` ,'  |  
        ###   |  ,' .  |     |  ,' .  |   
        ###   |,'    `.|     |,'    `.|
        ###   1--------z     2--------3   

def develop_verts_CP1(vertex_posns, gluing, face_vertex, tet_shape):
    """Get vert posns for a new tetrahedron given vert posns for an existing tetrahedron and gluing info"""
    next_tet_vertex_posns = []
    current_vertex_nums = [0,1,2,3]
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

