#
# winding_bistellar_flip.py
#

# Pachner moves that respect taut structures

import regina
import snappy

from veering import veering_census
from veering.taut import unsorted_vert_pair_to_edge_pair
from veering.transverse_taut import is_transverse_taut
from winding import is_winding, is_flat, apply_isom_to_winding, flat_reduced_windings, winding_to_preangle, analyse_angle_structs, can_deal_with_reduced_angles
from pachner_graph_path import search_Pachner_graph_for_shortest_path, build_Pachner_graph, get_path

def twoThreeMove(tri, wind, face_num, twist = 0, perform = True, test_ambiguous = False, check_twist_is_good_assuming_non_ambiguous = False, return_edge = False, verbose = 0):
    """Apply a 2-3 move to a triangulation with a winding, if possible. 
    If perform = False, returns if the move is possible.
    If perform = True, modifies tri, returns (tri, winding) for the performed move"""
    
    assert face_num < tri.countTriangles()
    M = snappy.Manifold(tri)
    assert is_winding(M, wind), (M, tri, wind)
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
    if perform == False:
        return True

    ### check we do the same as regina... 
    tri2 = regina.Triangulation3(tri)  ## make a copy
    tri2.pachner(tri2.triangle(face_num))

    ## record the tetrahedra and gluings adjacent to tet0 and tet1

    tets = [tet0, tet1]
    tet_nums = [tet_num0, tet_num1]
    vertices = [vertices0, vertices1]

    gluings = [] 
    for i in range(2):
        tet_gluings = []
        for j in range(3):
            tet_gluings.append( [ tets[i].adjacentTetrahedron(vertices[i][j]),  tets[i].adjacentGluing(vertices[i][j])] )
        gluings.append(tet_gluings)

    ### add new tetrahedra
    new_tets = []
    for i in range(3):
        new_tets.append(tri.newTetrahedron())

    ### glue around degree 3 edge
    for i in range(3):
        new_tets[i].join(2, new_tets[(i+1)%3], regina.Perm4(0,1,3,2))

    ### record winding values for tet0 and tet1
    w = []
    for t in range(2):
        w.append([None, None, None])

    # print('wind', wind)

    for i in range(2):
        for j in range(3):
            w[i][j] = wind[ 3*tet_nums[i] + unsorted_vert_pair_to_edge_pair[(vertices[i][0], vertices[i][j+1])] ]
        
    # print('w', w)

    is_ambiguous = ( (w[1][2] + w[0][2]) % 2 == (w[1][1] + w[0][1]) % 2 == (w[1][0] + w[0][0]) % 2 )
    ### ambiguous if all are even, i.e. the pi's meet at a single equatorial edge 
    if verbose > 3:
        print('is_ambiguous', is_ambiguous)
    if test_ambiguous:
        return is_ambiguous
    if is_ambiguous:
        non_ambiguous_index = None
        non_ambiguous_tet_num = None
    else:
        for i in range(3):
            if (w[1][i] + w[0][i]) % 2 == 0:  ### will hit once
                non_ambiguous_index = i  
                non_ambiguous_tet_num = 2-i

    ### replace mapping info with corresponding info for the 3 tet. Self gluings will be annoying...

    ### write verticesi[j] as vij    Paper notation is in parentheses
    ### old windings: wij, new windings: nij

    ###                  tet0                                    new_tet0
    ###                _________                                 _________
    ###    (y_3) w00 ,'\       /`. w01 (y_2)                   ,'\`.   ,'/`.  
    ###            ,'   \ v03 /   `.                         ,'   \ `0' /   `. 
    ###          ,'      \   /      `.                     ,'      \ | /      `.          n00 = w12 + w02   (a_1 = x_1 + y_1)
    ###         / \       \ /       / \                   /|\       \|/       /|\         n10 = w11 + w01   (b_1 = x_2 + y_2)
    ###        /v02\       *       /v01\                 / | \       *       / | \        n20 = w10 + w00   (c_1 = x_3 + y_3)
    ###       /    _\.(y_1)|w02.../_    \               /  | 3\..... | ...../2 |  \ 
    ###      /_--"" /      *      \ ""--_\             /2 ,'  /      *      \  `. 3\      pick a value for twist
    ###      \`.v12/      / \      \v11,'/      `.     \`.|  /      /|\      \  |,'/      n22 = twist       (c_3 = k)
    ###       \ `./      /   \      \,' /     ----}     \ `./      / | \      \,' /       n11 = w12 - n22   (b_2 = x_1 - c_3)
    ###        \ /`.    / v00 \    ,'\ /        ,'       \|/`.    /  |  \    ,'\|/        n02 = w00 - n11   (a_3 = y_3 - b_2)
    ###         \   `. /       \ ,'   /                   \   `. /   |   \ ,'   /         n21 = w11 - n02   (c_2 = x_2 - a_3)
    ###          \    `---------'    /                     \    * 3  |  2 *    /          n12 = w02 - n21   (b_3 = y_1 - c_2)
    ###           \    \       /    /                       \    \   |   /    /           n01 = w10 - n12   (a_2 = x_3 - b_3)
    ###            \    \ v10 /    /               new_tet1  \    \  |  /    /  new_tet2 (n22 = w01 - n01)
    ###   (x_3) w10 \    \   /    / w11 (x_2)                 \    \ | /    /  
    ###              \    \ /    /                             \    \|/    /
    ###               \    *    /                               \    *    /
    ###         tet1   \...|w12/ (x_1)                           \...|.../
    ###                 \  |  /                                   \`.|.'/
    ###                  \v13/                                     \ 1 /
    ###                   \|/                                       \|/
    ###                    *                                         *

    n = []
    for t in range(3):
        n.append([None, None, None])
    n[0][0] = w[1][2] + w[0][2]
    n[1][0] = w[1][1] + w[0][1]
    n[2][0] = w[1][0] + w[0][0]

    n[2][2] = twist
    n[1][1] = w[1][2] - n[2][2] 
    n[0][2] = w[0][0] - n[1][1]
    n[2][1] = w[1][1] - n[0][2]
    n[1][2] = w[0][2] - n[2][1]
    n[0][1] = w[1][0] - n[1][2]

    ### See Fig 5.14(a). In the following block of code, w[0][pi_index_on_tet0] is the edge XC
    ### and n[non_ambiguous_tet_num][(2 - non_ambiguous_tet_num - pi_index_on_tet0) % 3] is the 
    ### edge of the new tetrahedron containing BC that gets glued to XC
    if (not is_ambiguous) and check_twist_is_good_assuming_non_ambiguous:
        for i in range(3):
            if w[0][i] % 2 == 1:  ### will hit once
                pi_index_on_tet0 = i
        # if pi_index_on_tet0 == 0 and non_ambiguous_tet_num == 0:
        #     return n[0][2] % 2 == 1
        # elif pi_index_on_tet0 == 0 and non_ambiguous_tet_num == 1:
        #     return n[1][1] % 2 == 1
        # elif pi_index_on_tet0 == 1 and non_ambiguous_tet_num == 2:
        #     return n[2][2] % 2 == 1
        # elif pi_index_on_tet0 == 1 and non_ambiguous_tet_num == 0:
        #     return n[0][1] % 2 == 1
        # elif pi_index_on_tet0 == 2 and non_ambiguous_tet_num == 1:
        #     return n[1][2] % 2 == 1
        # elif pi_index_on_tet0 == 2 and non_ambiguous_tet_num == 2:
        #     return n[2][1] % 2 == 1
        if not n[non_ambiguous_tet_num][(2 - non_ambiguous_tet_num - pi_index_on_tet0) % 3] % 2 == 1:
            return False ### twist is bad

    # print('n', n)
    # permutations taking the vertices for a face of the 3-tet ball to the 
    # vertices of the same face for the 2-tet ball

    # these should be even in order to preserve orientedness.
    # exactly one of vertices[0] and vertices[1] is even, but it seems to depend on the face.

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
    # print('flip', flip)
    if flip:  #then all of the signs are wrong, switch 0 and 1 on input
        perms = [[p * regina.Perm4( 1,0,2,3 ) for p in a] for a in perms]
        for i in range(3):
            n[i][1], n[i][2] = n[i][2], n[i][1]

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

    ### update the winding
    tet_nums.sort()
    for i in range(2):
        for j in range(3):
            wind.pop( 3*tet_nums[1-i] + (2-j) )  ## remove from the list in the correct order!

    wind.extend( [n[0][0], n[0][1], n[0][2], 
                  n[1][0], n[1][1], n[1][2],
                  n[2][0], n[2][1], n[2][2]] )

    M = snappy.Manifold(tri)
    assert is_winding(M, wind)

    if not return_edge:
        return [ tri, wind ]
    else:
        return [ tri, wind, new_tets[0].edge(0).index() ]

def threeTwoMove(tri, wind, edge_num, perform = True, return_triangle = False):
    """Apply a 3-2 move to a triangulation with a winding, if possible. 
    If perform = False, returns if the move is possible.
    If perform = True, modifies tri, returns (tri, winding) for the performed move"""

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

    if len(set([tet.index() for tet in tets])) != 3: 
        return False  ### tetrahedra must be distinct

    if not perform:
        return True  ### taut 3-2 move is always possible if the 3-2 move is.

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


    ### record winding values for tets
    w = []
    for t in range(3):
        w.append([None, None, None])

    # print('wind', wind)
    for i in range(3):
        for j in range(3):
            w[i][j] = wind[ 3*tet_nums[i] + unsorted_vert_pair_to_edge_pair[(vertices[i][0], vertices[i][j+1])] ]

    ### replace mapping info with corresponding info for the 2 tet. Self gluings will be annoying...

    ### write vertices[i][j] as vij

    ###                 tets[0]                                   new_tet1
    ###                _________                                 _________
    ###              ,'\`.v00,'/`.                             ,'\       /`.
    ###            ,'   \ `.' /   `.                      n11,'   \  3  /   `.n10 
    ###          ,'   v10\ | /v20   `.                     ,'      \   /      `.
    ###         /|\       \|/       /|\                   / \       \ /       / \
    ###        / | \       *       / | \                 /   \       *       /   \
    ###    v12/  |  \..... | ...../  |  \v23            /  1 _\..... |n12.../_ 2  \ 
    ###      /  ,'  /      *      \  `.  \             /_--"" /      *      \ ""--_\
    ###      \`.|  /v03   /|\   v02\  |,'/      `.     \`. 2 /      / \      \ 1 ,'/ 
    ###       \ `./      / | \      \,' /     ----}     \ `./      /   \      \,' /
    ###        \|/`.    /  |  \    ,'\|/        ,'       \ /`.    /  0  \    ,'\ /
    ###         \   `. /   |   \ ,'   /                   \   `. /       \ ,'   /
    ###          \    * v13|v22 *    /                     \    `---------'    /
    ###           \    \   |   /    /                       \    \       /    /
    ###            \    \  |  /    /                         \    \  0  /    /
    ###             \    \ | /    /                       n00 \    \   /    / n01
    ###    tets[1]   \    \|/    /   tets[2]                   \    \ /    /
    ###               \    *    /                               \    *    / new_tet0
    ###                \..v01../                                 \...|n02/
    ###                 \`.|.'/                                   \  |  /
    ###               v11\ | /v21                                  \ 3 /
    ###                   \|/                                       \|/
    ###                    *                                         *

    n = []
    for t in range(2):
        n.append([None, None, None])
    n[0][0] = w[0][1] + w[1][2]
    n[0][1] = w[2][1] + w[0][2] 
    n[0][2] = w[1][1] + w[2][2]
    
    n[1][0] = w[2][2] + w[0][1]
    n[1][1] = w[0][2] + w[1][1] 
    n[1][2] = w[1][2] + w[2][1]    

    # permutations taking the vertices for a face of the 2-tet ball to the 
    # vertices of the same face for the 3-tet ball

    # these should be even in order to preserve orientability.

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

    # ### update the winding
    tet_nums.sort()
    for i in range(3):
        for j in range(3):
            wind.pop( 3*tet_nums[2-i] + (2-j) )  ## remove from the list in the correct order!
    
    wind.extend( [n[0][0], n[0][1], n[0][2], 
                  n[1][0], n[1][1], n[1][2]] )

    M = snappy.Manifold(tri)
    assert is_winding(M, wind)

    if not return_triangle:
        return (tri, wind)    
    else:
        return (tri, wind, new_tets[0].triangle(3).index())    

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


def carry_windings(start_winds, path_of_moves, verbose = 0):
    """Given a path of moves, carry windings as you go (eliminating non flat windings)"""
    ### each move is a pair of a sig and either a face number or edge number, written 'fn' or 'en'
    
    tri = regina.Triangulation3(path_of_moves[0][0])
    tri.orient()
    before_winds = start_winds[:]
    if verbose > 2:
        for wind in before_winds:
            print('before wind', wind)
    for move in path_of_moves[:-1]: ### except the last in the path, which only gives the final sig
        if verbose > 2:
            print('move', move)
        after_winds = []
        sig, flip = move
        assert tri.isoSig() == sig
        cell_index = int(flip[1:])
        if flip[0] == 'f':
            for wind in before_winds:
                for twist in range(2):
                    tri_temp = regina.Triangulation3(tri)
                    wind_temp = wind[:]
                    out = twoThreeMove(tri_temp, wind_temp, cell_index, twist = twist, check_twist_is_good_assuming_non_ambiguous = True, verbose = verbose)
                    if out != False:
                        new_tri, new_wind = out 
                        if is_flat(new_wind):
                            after_winds.append(new_wind)
        else:
            assert flip[0] == 'e'
            for wind in before_winds:
                tri_temp = regina.Triangulation3(tri)
                wind_temp = wind[:]
                new_tri, new_wind = threeTwoMove(tri_temp, wind_temp, cell_index)
                if is_flat(new_wind):
                    after_winds.append(new_wind)
        if verbose > 2:
            for wind in after_winds:
                print('after_wind', wind)
        if len(after_winds) == 0:
            if verbose > 1:
                print('no windings')
            return None ### no point in continuing along the path of moves.
        before_winds = after_winds
        tri = new_tri
    sig, _ = path_of_moves[-1]
    if verbose > 2:
        print('move', path_of_moves[-1])
    assert tri.isoSig() == sig, (sig, tri.isoSig())
    return after_winds, tri

def windings_diaspora(start_sig, start_winds = None, target_sigs = None, search_depth = 1000, ceiling = 5, verbose = 0):
    if start_winds == None:
        tri = regina.Triangulation3(start_sig)
        tri.orient()
        M = snappy.Manifold(tri)
        start_winds = flat_reduced_windings(M, return_lifts = True)[1]
    G = build_Pachner_graph(start_sig, search_depth = 1000, ceiling = ceiling)
    if target_sigs == None:
        target_sigs = []
        for node in G:
            if len(G[node].neighbour_moves_down_tri) == 0:
                if verbose > 2:
                    print('node', node)
                target_sigs.append(node)  ### look for local minima
    assert all(ts in G for ts in target_sigs)
    survivors = []
    junk = []
    for wind in start_winds:
        survived = True
        if verbose > 3:
            print('start winding', wind)
        for target_sig in target_sigs:
            if target_sig != start_sig:
                # print('target sig', target_sig)
                path_of_moves = get_path(target_sig, G) 
                tri = regina.Triangulation3(start_sig)
                tri.orient()
                carry_wind = wind[:] ### make a copy
                carried = carry_windings([carry_wind], path_of_moves, verbose = verbose)
                if carried == None:
                    if verbose > 0:
                        print('winding', wind, 'became', carried, 'at', target_sig)
                    survived = False
                    junk.append(wind)
                    break
                else:  
                    carried_winds, after_tri = carried
                    if verbose > 3:
                        print('sig', target_sig)
                        print('carried winds')
                        for cwind in carried_winds:
                            print(cwind)
                    for cwind in carried_winds: ### sanity check: should still have trivial internal z2 cohomology
                        cangle = winding_to_preangle(cwind)
                        assert is_transverse_taut(after_tri, cangle, check_is_taut = False)
        if survived:
            survivors.append(wind)
    if verbose > 0:
        print('survivors', survivors)   ### may be some that map to the same on some triangulation
        print('junk', junk)
    return (survivors, junk)


def search_census_for_transport_wins(census, ceiling_overhead = 1):
    for input_sig in census:
        if '_' in input_sig:    
            sig = input_sig.split('_')[0]
        else:
            sig = input_sig
        tri = regina.Triangulation3(sig)
        tri.orient()
        M = snappy.Manifold(tri)
        assert M.triangulation_isosig(decorated = False) == sig
        n = tri.countTetrahedra()
        survivors, junk = windings_diaspora(sig, ceiling = n + ceiling_overhead)
        if len(junk) > 0:  ### if no junk then we already checked these examples before
            survivors = [winding_to_preangle(wind) for wind in survivors]
            # analyse_angle_structs(M, angles = survivors, use_flipper = False, apply_symmetries = True)
            can_deal = can_deal_with_reduced_angles(M, angles = survivors)
            if can_deal[0]:
                print(input_sig)
                print(survivors)
                print(can_deal)

#### searched first 100000 of snappy census, found 
# gLLPQbefefefhhhhhha
# [[0, 2, 2, 2, 0, 2], [1, 1, 0, 0, 1, 1]]
# (True, ['veering+layered', 'veering+layered'])
# ivLAPQcedfhgghhgkjahmhtsn
# [[0, 2, 2, 1, 1, 0, 2, 0], [0, 1, 1, 2, 2, 0, 0, 1]]
# (True, ['layered', 'layered'])

#### veering census search ran out of memory somewhere in [500:2000], before then found
# fLLQcbeddeehhbghh_01110
# [[0, 2, 2, 2, 0], [1, 1, 0, 0, 1]]
# (True, ['veering+layered', 'layered'])
# gLLPQbefefefhhhhhha_011102
# [[0, 2, 2, 2, 0, 2], [1, 1, 0, 0, 1, 1]]
# (True, ['veering+layered', 'veering+layered'])
# ivvPQQcfhgghfhgfaaahhahha_10000122
# [[2, 0, 0, 0, 0, 2, 1, 1], [1, 0, 0, 0, 0, 1, 2, 2]]
# (True, ['veering+layered', 'veering+layered'])



