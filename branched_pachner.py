#
# branched_pachner.py
#

# Pachner moves that respect branched surfaces (smoothings of the dual spine)

import regina
from taut import isosig_to_tri_angle, unsorted_vert_pair_to_edge_pair
from branched_surface import large_edge_of_face, determine_possible_branch_given_two_faces, determine_possible_branch_given_three_faces, is_branched
from branched_surface import all_branched_surfaces, lex_smallest_branched_surface, apply_isom_to_branched_surface

def twoThreeMove(tri, branch, face_num, perform = True, return_edge = False):
    """Apply a 2-3 move to a triangulation with a branched surface, if possible. 
    If perform = False, returns if the move is possible.
    If perform = True, modifies tri, returns (tri, possible_branches) for the performed move"""

    ### possible_branches is a list
    
    assert is_branched(tri, branch)

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

    ### are all moves valid for the branched surface?
    ### for now, lets assume yes

    ### check we do the same as regina... 
    tri2 = regina.Triangulation3(tri)  ## make a copy
    tri2.twoThreeMove(tri2.triangle(face_num))

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

    large_edges = []  ### record large edge info for the outer faces
    for i in range(2):
        this_tet_large_edges = []
        for j in range(3):
            this_tet_large_edges.append(large_edge_of_face( branch[tets[i].index()], vertices[i][j] ))
        large_edges.append(this_tet_large_edges)

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
    #     print('flip')
    # else:
    #     print('no flip')

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

    ### update the branched surface

    ### for each of the three new tetrahedra, figure out what their outer face train tracks are

    large_edges_new = []  ### record large edge info for the outer faces
    for i in range(2):
        this_tet_large_edges_new = []
        for j in range(3):
            new_large_edge = perms[i][j].inverse()[ large_edges[i][j] ]
            if flip:
                assert new_large_edge != i ### the face number cannot be the large vertex for that face
            else:
                assert new_large_edge != 1 - i 
            this_tet_large_edges_new.append( new_large_edge )
        large_edges_new.append(this_tet_large_edges_new)

    candidate_branches = []
    for j in range(3):
        if flip:   
            candidate_branches.append( determine_possible_branch_given_two_faces((0,1), (large_edges_new[0][j], large_edges_new[1][j]) ) )
        else:
            candidate_branches.append( determine_possible_branch_given_two_faces((1,0), (large_edges_new[0][j], large_edges_new[1][j]) ) )

    ### update the branch structure, many possible ways
    tet_indices = [tet_num0, tet_num1]
    tet_indices.sort()
    branch.pop(tet_indices[1])
    branch.pop(tet_indices[0])  ## remove from the list in the correct order!

    out = []
    for cand0 in candidate_branches[0]:
        for cand1 in candidate_branches[1]:
            for cand2 in candidate_branches[2]:
                candidate = branch[:] + [cand0, cand1, cand2]
                # print('candidate', candidate)
                if is_branched(tri, candidate):
                    out.append(candidate)
    assert len(out) > 0
    if not return_edge:
        return (tri, out)
    else:
        return (tri, out, new_tets[0].edge(0).index())


def threeTwoMove(tri, branch, edge_num, return_triangle = False):
    """Apply a 3-2 move to a triangulation with a branched surface, if possible. 
    If perform = False, returns if the move is possible.
    modifies tri, returns (tri, branch) for the performed move"""

    ### note if this function does not return False, then there is only one 
    ### possible branch so we just return it rather than a list

    assert is_branched(tri, branch)

    edge = tri.edge(edge_num)
    if edge.degree() != 3:
        return False

    tets = []
    tet_nums = []
    vertices = []
    for i in range(3):
        embed = edge.embedding(i)
        tets.append(embed.simplex())
        tet_nums.append(tets[i].index())
        vertices.append(embed.vertices())

    if len(set(tets)) != 3: 
        return False  ### tetrahedra must be distinct
     
    ### check we do the same as regina... 
    tri2 = regina.Triangulation3(tri)  ## make a copy
    tri2.threeTwoMove(tri2.edge(edge_num))

    ## record the tetrahedra and gluings adjacent to the tets 

    gluings = [] 
    for i in range(3):
        tet_gluings = []
        for j in range(2):
            tet_gluings.append( [ tets[i].adjacentTetrahedron(vertices[i][j]),  tets[i].adjacentGluing(vertices[i][j])] )
        gluings.append(tet_gluings)

    for i in range(3):
        assert tets[i].adjacentTetrahedron(vertices[i][2]) == tets[(i+1)%3]  ### The edge embeddings should be ordered this way...

    large_edges = []  ### record large edge info for the outer faces
    for i in range(3):
        this_tet_large_edges = []
        for j in range(2):
            this_tet_large_edges.append(large_edge_of_face( branch[tets[i].index()], vertices[i][j] ))
        large_edges.append(this_tet_large_edges)

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

    # ### update the branched surface
    tet_nums.sort()
    branch.pop(tet_nums[2])
    branch.pop(tet_nums[1])
    branch.pop(tet_nums[0])  ## remove from the list in the correct order!

    ### update the branched surface

    ### for each of the two new tetrahedra, figure out what their outer face train tracks are

    large_edges_new = []  ### record large edge info for the outer faces
    for i in range(3):
        this_tet_large_edges_new = []
        for j in range(2):
            new_large_edge = perms[i][j].inverse()[ large_edges[i][j] ]
            # assert new_large_edge != i ??
            this_tet_large_edges_new.append( new_large_edge )
        large_edges_new.append(this_tet_large_edges_new)

    large_edges_new_transposed = [list(i) for i in zip(*large_edges_new)]

    branch0 = determine_possible_branch_given_three_faces([0,1,2], large_edges_new_transposed[0])
    branch1 = determine_possible_branch_given_three_faces([0,2,1], large_edges_new_transposed[1])
    if branch0 == None or branch1 == None:
        return False
    large_edge_for_new_tet0 = large_edge_of_face( branch0, 3 )
    large_edge_for_new_tet1 = large_edge_of_face( branch1, 3 )
    if large_edge_for_new_tet0 == large_edge_for_new_tet1:
        if large_edge_for_new_tet0 != 0:
            return False
    else:
        if large_edge_for_new_tet0 + large_edge_for_new_tet1 != 3:  ### one must be 1, one must be 2
            return False

    branch.extend([branch0, branch1])

    if not is_branched(tri, branch):
        return False

    if not return_triangle:
        return (tri, branch)    
    else:
        return (tri, branch, new_tets[0].triangle(3).index())    


def main():
    # for i in range(4): 
    #     print(i)
    #     # tri, angle = isosig_to_tri_angle('cPcbbbiht_12') 
    # sigs = ['dLQacccjsnk_200', 'dLQbccchhfo_122','dLQbccchhsj_122']
    # for sig in sigs:
    #     print(sig)
    #     tri, angle = isosig_to_tri_angle(sig)   
    #     for branch in all_branched_surfaces(tri):
    #         print(lex_smallest_branched_surface(tri, branch))

    sig = 'dLQacccjsnk_200'
    for i in range(6):
        # print(i)
        tri, angle = isosig_to_tri_angle(sig)  
        tri_original = regina.Triangulation3(tri) #copy

        out = twoThreeMove(tri, [4,11,0], i, return_edge = True)
        if out != False:
            tri, possible_branches, edge_num = out
            # print('possible_branches', possible_branches)
            # print('all branches', all_branched_surfaces(tri)) 
            tri, branch = threeTwoMove(tri, possible_branches[0], edge_num)
            all_isoms = tri.findAllIsomorphisms(tri_original)
            all_branches = [apply_isom_to_branched_surface(branch, isom) for isom in all_isoms]
            assert [4,11,0] in all_branches

# 0
# False
# [(4, 11, 0), (6, 0, 11)]
# 1
# [[11, 3, 4, 2]]
# [(0, 6, 0, 9), (4, 0, 0, 0), (11, 3, 4, 2)]
# 2
# [[0, 4, 3, 0]]
# [(0, 4, 3, 0), (6, 0, 0, 0), (11, 9, 8, 6)]
# 3
# [[4, 9, 6, 1]]
# [(4, 9, 6, 1), (6, 1, 4, 3)]
# 4
# [[4, 4, 3, 1]]
# [(4, 4, 3, 1), (6, 6, 1, 9), (10, 0, 0, 0)]
# 5
# [[4, 6, 0, 9]]
# [(4, 6, 0, 9), (6, 0, 4, 3)]

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


