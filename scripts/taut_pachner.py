#
# taut_pachner.py
#

# Pachner moves that respect taut structures

import regina
from taut import isosig_to_tri_angle, unsorted_vert_pair_to_edge_pair, is_taut, there_is_a_pi_here

def twoThreeMove(tri, angle, face_num, perform = True, return_edge = False):
    """Apply a 2-3 move to a taut triangulation, if possible. 
    If perform = False, returns if the move is possible.
    If perform = True, modifies tri, returns (tri, angle) for the performed move"""
    
    face = tri.triangle(face_num)

    embed0 = face.embedding(0)
    tet0 = embed0.simplex()
    tet_num0 = tet0.index()
    tet_0_face_num = embed0.face()
    vertices0 = embed0.vertices() # Maps vertices (0,1,2) of face to the corresponding vertex numbers of tet0

    embed1 = face.embedding(1)
    tet1 = embed1.simplex()
    tet_num1 = tet1.index()
    tet_1_face_num = embed1.face()
    vertices1 = embed1.vertices() # Maps vertices (0,1,2) of face to the corresponding vertex numbers of tet1

    if tet0 == tet1:  ### Cannot perform a 2-3 move across a self-gluing
        return False

    ### taut 2-3 move is valid if the pis are on different edges of face
    ### this never happens if we start with a veering triangulation.
    ### for veering, the two-tetrahedron ball is always a continent.

    for i in range(3):
        j = (i+1) % 3
        k = (i+2) % 3
        if angle[tet_num0] == unsorted_vert_pair_to_edge_pair[(vertices0[j], vertices0[k])]:
            pi_num_0 = i
        if angle[tet_num1] == unsorted_vert_pair_to_edge_pair[(vertices1[j], vertices1[k])]:
            pi_num_1 = i

    if pi_num_0 == pi_num_1:
        return False
    if perform == False:
        return True

    ### check we do the same as regina... 
    tri2 = regina.Triangulation3(tri)  ## make a copy
    tri2.pachner(tri2.triangle(face_num))

    ### We have to implement twoThreeMove ourselves. e.g. we do a 2-3 move to canonical fig 8 knot complement triangulation. 
    ### All of the original tetrahedra are removed. I don't see any way to carry the angle structure through without knowing
    ### exactly how Ben's implementation works.

    ## record the tetrahedra and gluings adjacent to tet0 and tet1

    tets = [tet0, tet1]
    vertices = [vertices0, vertices1]

    # print('2-3 vertices signs')
    # print([v.sign() for v in vertices])

    gluings = [] 
    for i in range(2):
        tet_gluings = []
        for j in range(3):
            tet_gluings.append( [ tets[i].adjacentTetrahedron(vertices[i][j]),  tets[i].adjacentGluing(vertices[i][j])] )
            # if tets[i].adjacentTetrahedron(vertices[i][j]) in tets:
            #     print('self gluing')
        gluings.append(tet_gluings)

    ### add new tetrahedra
    new_tets = []
    for i in range(3):
        new_tets.append(tri.newTetrahedron())

    ### glue around degree 3 edge
    for i in range(3):
        new_tets[i].join(2, new_tets[(i+1)%3], regina.Perm4(0,1,3,2))

    ### replace mapping info with corresponding info for the 3 tet. Self gluings will be annoying...

    ### write verticesi[j] as vij

    ###                  tet0                                    new_tet0
    ###                _________                                 _________
    ###              ,'\       /`.                             ,'\`.   ,'/`.
    ###            ,'   \ v03 /   `.                         ,'   \ `0' /   `. 
    ###          ,'      \   /      `.                     ,'      \ | /      `.
    ###         / \       \ /       / \                   /|\       \|/       /|\
    ###        /v02\       *       /v01\                 / | \       *       / | \
    ###       /    _\..... | ...../_    \               /  | 3\..... | ...../2 |  \ 
    ###      /_--"" /      *      \ ""--_\             /2 ,'  /      *      \  `. 3\
    ###      \`.v12/      / \      \v11,'/      `.     \`.|  /      /|\      \  |,'/ 
    ###       \ `./      /   \      \,' /     ----}     \ `./      / | \      \,' /
    ###        \ /`.    / v00 \    ,'\ /        ,'       \|/`.    /  |  \    ,'\|/
    ###         \   `. /       \ ,'   /                   \   `. /   |   \ ,'   /
    ###          \    `---------'    /                     \    * 3  |  2 *    /
    ###           \    \       /    /                       \    \   |   /    /
    ###            \    \ v10 /    /               new_tet1  \    \  |  /    /  new_tet2
    ###             \    \   /    /                           \    \ | /    /  
    ###              \    \ /    /                             \    \|/    /
    ###               \    *    /                               \    *    /
    ###         tet1   \...|.../                                 \...|.../
    ###                 \  |  /                                   \`.|.'/
    ###                  \v13/                                     \ 1 /
    ###                   \|/                                       \|/
    ###                    *                                         *

    # permutations taking the vertices for a face of the 3-tet ball to the 
    # vertices of the same face for the 2-tet ball

    # these should be even in order to preserve orientability.
    # exactly one of vertices[0] and vertices[1] is even, but it seems to depend on the face.

    # perms = [[regina.Perm4( vertices[0][3], vertices[0][0], vertices[0][1], vertices[0][2] ),   ### opposite v00
    #           regina.Perm4( vertices[0][3], vertices[0][1], vertices[0][2], vertices[0][0] ),   ### opposite v01
    #           regina.Perm4( vertices[0][3], vertices[0][2], vertices[0][0], vertices[0][1] )    ### opposite v02
    #           ],  
    #          [regina.Perm4( vertices[1][0], vertices[1][3], vertices[1][1], vertices[1][2] ),   ### opposite v10
    #           regina.Perm4( vertices[1][1], vertices[1][3], vertices[1][2], vertices[1][0] ),   ### opposite v11
    #           regina.Perm4( vertices[1][2], vertices[1][3], vertices[1][0], vertices[1][1] )    ### opposite v12
    #           ]
    #         ]

    perms = [[vertices[0] * regina.Perm4( 3,0,1,2 ),   ### opposite v00
              vertices[0] * regina.Perm4( 3,1,2,0 ),   ### opposite v01
              vertices[0] * regina.Perm4( 3,2,0,1 )    ### opposite v02
              ],  
             [vertices[1] * regina.Perm4( 0,3,1,2 ),   ### opposite v10
              vertices[1] * regina.Perm4( 1,3,2,0 ),   ### opposite v11
              vertices[1] * regina.Perm4( 2,3,0,1 )    ### opposite v12
              ]
            ]
    flip = perms[0][0].sign() == -1
    if flip:  #then all of the signs are wrong, switch 0 and 1 on input
        perms = [[p * regina.Perm4( 1,0,2,3 ) for p in a] for a in perms]

    # print('2-3 perms signs')
    # print([[p.sign() for p in a] for a in perms])

    for i in range(2):
        for j in range(3):
            gluing = gluings[i][j]
            if gluing != None:
                if gluing[0] not in tets:  ### not a self gluing
                    gluing[1] = gluing[1] * perms[i][j]
                else:
                    i_other = tets.index( gluing[0] )
                    otherfacenum = gluing[1][vertices[i][j]]
                    j_other = [vertices[i_other][k] for k in range(4)].index(otherfacenum)
                    assert gluings[i_other][j_other][0] == tets[i]
                    assert gluings[i_other][j_other][1].inverse() == gluings[i][j][1]

                    gluings[i_other][j_other] = None ### only do a self gluing from one side 
                    gluing[0] = new_tets[j_other]
                    gluing[1] = perms[i_other][j_other].inverse() * gluing[1] * perms[i][j] 

    ### unglue two tetrahedra
    tet0.isolate()
    tet1.isolate()

    ### remove the tetrahedra
    tri.removeSimplex(tet0)
    tri.removeSimplex(tet1)

    ### make the gluings on the boundary of the new ball
    for i in range(2):
        for j in range(3):
            if gluings[i][j] != None:
                if flip:
                    new_tets[j].join(i, gluings[i][j][0], gluings[i][j][1])
                else:
                    new_tets[j].join(1 - i, gluings[i][j][0], gluings[i][j][1])

    assert tri.isIsomorphicTo(tri2)
    assert tri.isOriented()

    ### update the angle structure
    tet_indices = [tet_num0, tet_num1]
    tet_indices.sort()
    angle.pop(tet_indices[1])
    angle.pop(tet_indices[0])  ## remove from the list in the correct order!

    new_angle = [None, None, None]
    new_angle[pi_num_0] = 0
    new_angle[pi_num_1] = 0 ### these two tetrahedra have their pi's on the new degree three edge

    third_index = 3 - (pi_num_0 + pi_num_1)
    if (pi_num_0 - third_index) % 3 == 1:
        new_angle[third_index] = 1
    else:
        assert (pi_num_0 - third_index) % 3 == 2
        new_angle[third_index] = 2
    if flip:
        new_angle[third_index] = 3 - new_angle[third_index]
    
    angle.extend(new_angle)

    assert is_taut(tri, angle)

    if not return_edge:
        return [ tri, angle ]
    else:
        return [ tri, angle, new_tets[0].edge(0).index() ]

def threeTwoMove(tri, angle, edge_num, perform = True, return_triangle = False):
    """Apply a 3-2 move to a taut triangulation, if possible. 
    If perform = False, returns if the move is possible.
    If perform = True, modifies tri, returns (tri, angle) for the performed move"""

    edge = tri.edge(edge_num)
    if edge.degree() != 3:
        return False

    tets = []
    tet_nums = []
    vertices = []
    non_pi_tet_num = None
    for i in range(3):
        embed = edge.embedding(i)
        tets.append(embed.simplex())
        tet_nums.append(tets[i].index())
        vertices.append(embed.vertices())
        if not there_is_a_pi_here(angle, embed):
            assert non_pi_tet_num == None
            non_pi_tet_num = embed.simplex().index()
            local_non_pi_tet_num = i

    if len(set([tet.index() for tet in tets])) != 3: 
        return False  ### tetrahedra must be distinct

    if not perform:
        return True  ### taut 3-2 move is always possible if the 3-2 move is.

    ### record the "slope" of the pis on the non_pi_tet. This is a boolean
    non_pi_tet_positive = unsorted_vert_pair_to_edge_pair[ ( vertices[local_non_pi_tet_num][0], vertices[local_non_pi_tet_num][2] ) ]
    is_positive_slope = (angle[non_pi_tet_num] == non_pi_tet_positive)
     
    ### check we do the same as regina... 
    tri2 = regina.Triangulation3(tri)  ## make a copy
    tri2.pachner(tri2.edge(edge_num))

    ## record the tetrahedra and gluings adjacent to the tets 

    gluings = [] 
    for i in range(3):
        tet_gluings = []
        for j in range(2):
            tet_gluings.append( [ tets[i].adjacentTetrahedron(vertices[i][j]),  tets[i].adjacentGluing(vertices[i][j])] )
        gluings.append(tet_gluings)

    for i in range(3):
        assert tets[i].adjacentTetrahedron(vertices[i][2]) == tets[(i+1)%3]  ### The edge embeddings should be ordered this way...

    ### add new tetrahedra
    new_tets = []
    for i in range(2):
        new_tets.append(tri.newTetrahedron())

    ### glue across face
    new_tets[0].join(3, new_tets[1], regina.Perm4(0,2,1,3))

    ### replace mapping info with corresponding info for the 2 tet. Self gluings will be annoying...

    ### write vertices[i][j] as vij

    ###                 tets[0]                                   new_tet1
    ###                _________                                 _________
    ###              ,'\`.v00,'/`.                             ,'\       /`.
    ###            ,'   \ `.' /   `.                         ,'   \  3  /   `. 
    ###          ,'   v10\ | /v20   `.                     ,'      \   /      `.
    ###         /|\       \|/       /|\                   / \       \ /       / \
    ###        / | \       *       / | \                 /   \       *       /   \
    ###    v12/  |  \..... | ...../  |  \v23            /  1 _\..... | ...../_ 2  \ 
    ###      /  ,'  /      *      \  `.  \             /_--"" /      *      \ ""--_\
    ###      \`.|  /v03   /|\   v02\  |,'/      `.     \`. 2 /      / \      \ 1 ,'/ 
    ###       \ `./      / | \      \,' /     ----}     \ `./      /   \      \,' /
    ###        \|/`.    /  |  \    ,'\|/        ,'       \ /`.    /  0  \    ,'\ /
    ###         \   `. /   |   \ ,'   /                   \   `. /       \ ,'   /
    ###          \    * v13|v22 *    /                     \    `---------'    /
    ###           \    \   |   /    /                       \    \       /    /
    ###            \    \  |  /    /                         \    \  0  /    /
    ###             \    \ | /    /                           \    \   /    /  
    ###    tets[1]   \    \|/    /   tets[2]                   \    \ /    /
    ###               \    *    /                               \    *    / new_tet0
    ###                \..v01../                                 \...|.../
    ###                 \`.|.'/                                   \  |  /
    ###               v11\ | /v21                                  \ 3 /
    ###                   \|/                                       \|/
    ###                    *                                         *

    # permutations taking the vertices for a face of the 2-tet ball to the 
    # vertices of the same face for the 3-tet ball

    # these should be even in order to preserve orientability.

    # perms = [[regina.Perm4( vertices[0][0], vertices[0][2], vertices[0][3], vertices[0][1] ),   ### opposite v00
    #           regina.Perm4( vertices[0][1], vertices[0][3], vertices[0][2], vertices[0][0] )    ### opposite v01
    #           ],
    #          [regina.Perm4( vertices[1][3], vertices[1][0], vertices[1][2], vertices[1][1] ),   ### opposite v10
    #           regina.Perm4( vertices[1][3], vertices[1][2], vertices[1][1], vertices[1][0] )    ### opposite v11
    #           ],
    #          [regina.Perm4( vertices[2][2], vertices[2][3], vertices[2][0], vertices[2][1] ),   ### opposite v20
    #           regina.Perm4( vertices[2][2], vertices[2][1], vertices[2][3], vertices[2][0] )    ### opposite v21
    #           ]
    #         ]

    perms = [[vertices[0] * regina.Perm4( 0, 2, 3, 1 ),   ### opposite v00
              vertices[0] * regina.Perm4( 1, 3, 2, 0 )    ### opposite v01
              ],
             [vertices[1] * regina.Perm4( 3, 0, 2, 1 ),   ### opposite v10
              vertices[1] * regina.Perm4( 3, 2, 1, 0 )    ### opposite v11
              ],
             [vertices[2] * regina.Perm4( 2, 3, 0, 1 ),   ### opposite v20
              vertices[2] * regina.Perm4( 2, 1, 3, 0 )    ### opposite v21
              ]
            ]

    for i in range(3):
        for j in range(2):
            gluing = gluings[i][j]
            if gluing != None:
                if gluing[0] not in tets:  ### not a self gluing
                    gluing[1] = gluing[1] * perms[i][j]
                else:
                    i_other = tets.index( gluing[0] )
                    otherfacenum = gluing[1][vertices[i][j]]
                    j_other = [vertices[i_other][k] for k in range(4)].index(otherfacenum) 
                    assert gluings[i_other][j_other][0] == tets[i]
                    assert gluings[i_other][j_other][1].inverse() == gluings[i][j][1]

                    gluings[i_other][j_other] = None ### only do a self gluing from one side 
                    gluing[0] = new_tets[j_other]  ### j refers to the vertex on the same 3 side
                    gluing[1] = perms[i_other][j_other].inverse() * gluing[1] * perms[i][j] 

    ### unglue three tetrahedra
    for tet in tets:
        tet.isolate()

    ### remove the tetrahedra
    for tet in tets:
        tri.removeSimplex(tet)

    ### make the gluings on the boundary of the new ball
    for i in range(3):
        for j in range(2):
            if gluings[i][j] != None:
                if j == 0 or i == 0:
                    assert new_tets[j].adjacentTetrahedron(i) == None ## not glued
                    assert gluings[i][j][0].adjacentTetrahedron(gluings[i][j][1][i]) == None
                    new_tets[j].join(i, gluings[i][j][0], gluings[i][j][1])
                else:
                    assert new_tets[j].adjacentTetrahedron(3 - i) == None ## not glued
                    assert gluings[i][j][0].adjacentTetrahedron(gluings[i][j][1][3 - i]) == None
                    new_tets[j].join(3 - i, gluings[i][j][0], gluings[i][j][1])  ## swap 1 and 2

    assert tri.isIsomorphicTo(tri2)
    assert tri.isOriented()

    # ### update the angle structure
    tet_nums.sort()
    angle.pop(tet_nums[2])
    angle.pop(tet_nums[1])
    angle.pop(tet_nums[0])  ## remove from the list in the correct order!

    if local_non_pi_tet_num == 0:
        if is_positive_slope:
            new_angle = [0, 0]
        else:
            new_angle = [1, 1]
    elif local_non_pi_tet_num == 1:
        if is_positive_slope:
            new_angle = [2, 1]
        else:
            new_angle = [0, 2]
    else:
        assert local_non_pi_tet_num == 2
        if is_positive_slope:
            new_angle = [1, 2]
        else:
            new_angle = [2, 0]
    
    angle.extend(new_angle)

    assert is_taut(tri, angle)

    if not return_triangle:
        return (tri, angle)    
    else:
        return (tri, angle, new_tets[0].triangle(3).index())    

    ### Extra ASCII art:

    ###                _________                                 _________
    ###              ,'\       /`.                             ,'\`.   ,'/`.
    ###            ,'   \     /   `.                         ,'   \ `.' /   `. 
    ###          ,'      \   /      `.                     ,'      \ | /      `.
    ###         / \       \ /       / \                   /|\       \|/       /|\
    ###        /   \       *       /   \                 / | \       *       / | \
    ###       /    _\..... | ...../_    \               /  |  \..... | ...../  |  \ 
    ###      /_--"" /      *      \ ""--_\             /  ,'  /      *      \  `.  \
    ###      \`.   /      / \      \   ,'/      `.     \`.|  /      /|\      \  |,'/ 
    ###       \ `./      /   \      \,' /     ----}     \ `./      / | \      \,' /
    ###        \ /`.    /     \    ,'\ /        ,'       \|/`.    /  |  \    ,'\|/
    ###         \   `. /       \ ,'   /                   \   `. /   |   \ ,'   /
    ###          \    `---------'    /                     \    *    |    *    /
    ###           \    \       /    /                       \    \   |   /    /
    ###            \    \     /    /                         \    \  |  /    /
    ###             \    \   /    /                           \    \ | /    /  
    ###              \    \ /    /                             \    \|/    /
    ###               \    *    /                               \    *    /
    ###                \...|.../                                 \...|.../
    ###                 \  |  /                                   \`.|.'/
    ###                  \ | /                                     \ | /
    ###                   \|/                                       \|/
    ###                    *                                         *


