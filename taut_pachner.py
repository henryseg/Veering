#
# taut_pachner.py
#

# Pachner moves that respect taut structures

import regina
from taut import unsorted_vert_pair_to_edge_pair

def twoThreeMove(tri, angle, face_num, perform = True):
    """Apply a 2-3 move to a taut triangulation, if possible. 
    If perform = False, returns if the move is possible.
    If perform = True, modifies tri, returns (tri, angle) for the performed move"""
    
    # tri = regina.Triangulation(tri)  ## make a copy

    face = tri.face(face_num)
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
            tet_gluing.append( tets[i].adjacentTetrahedron(vertices[i][j]).index(),  tets[i].adjacentGluing(vertices[i][j]) )
        gluings.append(tet_gluings)

    ### add new tetrahedra
    new_tets = []
    for i in range(3):
        new_tets.append(tri.newTetrahedron())

    ### glue around degree 3 edge
    for i in range(3):
        new_tets[i].join(2, new_tets[(i+1)%3], regina.Perm4(0,1,3,2))

    ### unglue two tetrahedra
    tet0.isolate()
    tet1.isolate()

    ### replace mapping info with corresponding info for the 3 tet. Self gluings will be annoying...

    ### write verticesi[j] as vij

    ###                  tet0                                     new_tet0
    ###                _________                                 _________
    ###              ,'\       /`.                             ,'\`.   ,'/`.
    ###            ,'   \ v03 /   `.                         ,'   \ `0' /   `. 
    ###          ,'      \   /      `.                     ,'      \ | /      `.
    ###         / \       \ /       / \                   /|\       \|/       /|\
    ###        /v02\       |       /v01\                 / | \       |       / | \
    ###       /    __----- | -----__    \               /  |  \----- | -----/  |  \ 
    ###      /_--"" /      |      \ ""--,\             /2 J 3 /      |      \ 2 L 3\
    ###      \`.v12/      / \      \v11,'/      `.     \`.|  /      /|\      \  |,'/ 
    ###       \ `./      /   \      \,' /     ---->     \ `./      / | \      \,' /
    ###        \ /`.    / v00 \    ,'\ /        ,'       \|/`.    /  |  \    ,'\|/
    ###         \   `. /       \ ,'   /                   \   `. /   |   \ ,'   /
    ###          \    `---------'    /                     \    `  3 | 2  '    /
    ###           \    \       /    /                       \    \   |   /    /
    ###            \    \ v10 /    /               new_tet1  \    \  |  /    /  new_tet2
    ###             \    \   /    /                           \    \ | /    /  
    ###              \    \ /    /                             \    \|/    /
    ###               \    |    /                               \    |    /
    ###         tet1   \___|___/                                 \___|___/
    ###                 \  |  /                                   \`.|.'/
    ###                  \v13/                                     \ 1 /
    ###                   \|/                                       \|/






    ### Extra ASCII art:

    ###                _________                                  _________
    ###              ,'\       /`.                             ,'\`.   ,'/`.
    ###            ,'   \     /   `.                         ,'   \ `.' /   `. 
    ###          ,'      \   /      `.                     ,'      \ | /      `.
    ###         / \       \ /       / \                   /|\       \|/       /|\
    ###        /   \       |       /   \                 / | \       |       / | \
    ###       /    __----- | -----__    \               /  |  \----- | -----/  |  \ 
    ###      /_--"" /      |      \ ""--,\             /  J   /      |      \   L  \
    ###      \`.   /      / \      \   ,'/      `.     \`.|  /      /|\      \  |,'/ 
    ###       \ `./      /   \      \,' /     ---->     \ `./      / | \      \,' /
    ###        \ /`.    /     \    ,'\ /        ,'       \|/`.    /  |  \    ,'\|/
    ###         \   `. /       \ ,'   /                   \   `. /   |   \ ,'   /
    ###          \    `---------'    /                     \    `    |    '    /
    ###           \    \       /    /                       \    \   |   /    /
    ###            \    \     /    /                         \    \  |  /    /
    ###             \    \   /    /                           \    \ | /    /  
    ###              \    \ /    /                             \    \|/    /
    ###               \    |    /                               \    |    /
    ###                \___|___/                                 \___|___/
    ###                 \  |  /                                   \`.|.'/
    ###                  \ | /                                     \ | /
    ###                   \|/                                       \|/

    ###
    ###                _________
    ###              ,'\       /`.
    ###            ,'   \     /   `.  
    ###          ,'      \   /      `.
    ###         / \       \ /       / \
    ###        /   \       |       /   \
    ###       /    __----- | -----__    \ 
    ###      /_--"" /      |      \ ""--,\
    ###      \`.   /      / \      \   ,'/
    ###       \ `./      /   \      \,' /
    ###        \ /`.    /     \    ,'\ /
    ###         \   `. /       \ ,'   / 
    ###          \    `---------'    /
    ###           \    \       /    /
    ###            \    \     /    /
    ###             \    \   /    /
    ###              \    \ /    /
    ###               \    |    /
    ###                \___|___/
    ###                 \  |  /
    ###                  \ | /
    ###                   \|/
    ###
    ###                    |
    ###                  \ | /      
    ###                   \|/
    ###                    V
    ###
    ###                _________
    ###              ,'\`.   ,'/`.
    ###            ,'   \ `.' /   `.  
    ###          ,'      \ | /      `.
    ###         /|\       \|/       /|\
    ###        / | \       |       / | \
    ###       /  |  \----- | -----/  |  \ 
    ###      /  J   /      |      \   L  \
    ###      \`.|  /      /|\      \  |,'/
    ###       \ `./      / | \      \,' /
    ###        \|/`.    /  |  \    ,'\|/
    ###         \   `. /   |   \ ,'   / 
    ###          \    `    |    '    /
    ###           \    \   |   /    /
    ###            \    \  |  /    /
    ###             \    \ | /    /
    ###              \    \|/    /
    ###               \    |    /
    ###                \___|___/
    ###                 \`.|.'/
    ###                  \ | /
    ###                   \|/




