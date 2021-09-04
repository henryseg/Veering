#
# taut_pachner.py
#

# Pachner moves that respect taut structures

import regina
from taut import unsorted_vert_pair_to_edge_pair, is_taut

def twoThreeMove(tri, angle, face_num, perform = True):
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

    if tet0 == tet1:
        return False

    ### taut 2-3 move is valid if the pis are on different edges of face

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
    tri2.twoThreeMove(tri2.triangle(face_num))

    ### We have to implement twoThreeMove ourselves. e.g. we do a 2-3 move to canonical fig 8 knot complement triangulation. 
    ### All of the original tetrahedra are removed. I don't see any way to carry the angle structure through without knowing
    ### exactly how Ben's implementation works.

    ## record the tetrahedra and gluings adjacent to tet0 and tet1

    tets = [tet0, tet1]
    vertices = [vertices0, vertices1]

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

    perms = [[regina.Perm4(vertices[0][3], vertices[0][0], vertices[0][1], vertices[0][2] ),   ### opposite v00
              regina.Perm4(vertices[0][3], vertices[0][1], vertices[0][2], vertices[0][0] ),   ### opposite v01
              regina.Perm4(vertices[0][3], vertices[0][2], vertices[0][0], vertices[0][1] )    ### opposite v02
              ],  
             [regina.Perm4(vertices[1][0], vertices[1][3], vertices[1][1], vertices[1][2] ),   ### opposite v10
              regina.Perm4(vertices[1][1], vertices[1][3], vertices[1][2], vertices[1][0] ),   ### opposite v11
              regina.Perm4(vertices[1][2], vertices[1][3], vertices[1][0], vertices[1][1] )    ### opposite v12
              ]
            ]

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
                new_tets[j].join((i + 1)%2, gluings[i][j][0], gluings[i][j][1])

    assert tri.isIsomorphicTo(tri2)

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
    
    angle.extend(new_angle)

    assert is_taut(tri, angle)

    return (tri, angle)


def main():
    t = regina.Triangulation3.fromIsoSig('jLLAvQQbcdeihhiihtsfxedxhdt')
    for i in range(t.countTriangles()):
        t2 = regina.Triangulation3(t)
        twoThreeMove(t2,[2,0,1,0,2,1,2,0,1],i) 


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


