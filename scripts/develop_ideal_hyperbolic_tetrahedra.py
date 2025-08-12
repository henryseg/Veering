#
# develop_ideal_hyperbolic_tetrahedra.py
#

from veering.basic_math import matrix, KP1

def developed_position(A1, A2, A3, z, field = None): #use Feng's "solving Thurston's equations in a commutative ring"
    # print A1, A2, A3, z
    # print('field', field)
    if field == None: ### using floating point (otherwise we are using algebraic numbers and don't need to check anything!)
        epsilon = 1e-32

    a1,b1 = A1 ## eg 1, 0
    a2,b2 = A2 ## eg 0, 1
    a3,b3 = A3 ## eg 1, 1
    
    Xinv = matrix((a1, a2, b1, b2))
    
    if field == None: 
        try:
            assert abs(Xinv.det()) > epsilon
        except:
            print(Xinv) 
            print((abs(Xinv.det())))
            print(('A1', A1))
            print(('A2', A2))
            print(('A3', A3))
            raise 

    X = Xinv.inverse() 
    A3p = X * A3
    a3p, b3p = A3p

    if field == None: ### using floating point
        try:
            assert abs(b3p) > epsilon and abs(a3p) > epsilon
        except: 
            print(('A3p', X * A3))
            print(('A1', A1))
            print(('A2', A2))
            print(('A3', A3))
            raise

    # B = KP1((-z/b3p, -1/a3p)).preferred_rep_saul()
    # B = KP1((-z/b3p, -1/a3p)).preferred_rep()
    if field == None:
        num_type = z.parent()
        B = KP1((-z/b3p, -num_type(1)/a3p))   ### Danger: if a3p is an integer, the second coordinate is a float
    else:
        B = KP1((-z/b3p, -field(1)/a3p))
    
    A4 = Xinv * B
    
    if field == None: ### using floating point
        a4, b4 = A4
        try:
            # assert (A1 - A3)(A2 - A4) / (A2 - A3)(A1 - A4) == z
            # assert (a1/b1 - a3/b3)(a2/b2 - a4/b4) / (a2/b2 - a3/b3)(a1/b1 - a4/b4) == z
            # assert (a1*b2*b3*b4 - a3*b1*b2*b4)*(a2*b1*b3*b4 - a4*b1*b2*b3) / (a2*b1*b3*b4 - a3*b1*b2*b4)(a1*b2*b3*b4 - a4*b1*b2*b3) == z
            # assert b2*b4*(a1*b3 - a3*b1)*b1*b3*(a2*b4 - a4*b2) / b1*b4*(a2*b3 - a3*b2)*b2*b3*(a1*b4 - a4*b1) == z
            # assert (a1*b3 - a3*b1)*(a2*b4 - a4*b2) / (a2*b3 - a3*b2)*(a1*b4 - a4*b1) == z
            # assert (a1*b3 - a3*b1)*(a2*b4 - a4*b2) - (a2*b3 - a3*b2)*(a1*b4 - a4*b1) * z == 0
            assert abs( (a1*b3 - a3*b1)*(a2*b4 - a4*b2) - (a2*b3 - a3*b2)*(a1*b4 - a4*b1) * z ) < epsilon
        except:
            print(('cross ratio', (a1*b3 - a3*b1)*(a2*b4 - a4*b2) / ( (a2*b3 - a3*b2)*(a1*b4 - a4*b1) )))
            print(('z', z))
            print(('error', abs( (a1*b3 - a3*b1)*(a2*b4 - a4*b2) - (a2*b3 - a3*b2)*(a1*b4 - a4*b1) * z )))
            raise

        try:
            assert A4.is_infinity() or abs(A4.complex()) < 10000.0
        except:
            print(('abs(A4.complex())', abs(A4.complex()) ))
            print(Xinv) 
            print((abs(Xinv.det())))
            print(('A1', A1))
            print(('A2', A2))
            print(('A3', A3))
            raise 

    # print('developed_position', A4)
    return A4

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

def develop_verts_pos(vertex_posns, gluing, face_vertex, tet_shape, field = None):
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
    new_position = developed_position(next_tet_vertex_posns[ordering[0]], next_tet_vertex_posns[ordering[1]], next_tet_vertex_posns[ordering[2]], tet_shape, field = field)
    next_tet_vertex_posns[next_tet_new_vertex_num] = new_position
    return next_tet_vertex_posns

