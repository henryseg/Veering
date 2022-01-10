#
# drill.py
#

### Given a triangulation and a loop of triangles, drill along that loop

import regina
from taut import isosig_to_tri_angle, reverse_tet_orientation, is_taut
from transverse_taut import is_transverse_taut

### anatomy of a loop of triangles:
###
### It is a list of tuples (tri_index, Perm3(v0, v1, v2), <add ordering info for multiple loop triangles carried by one triangulation triangle>)
### where tri_index is the index of the triangle the loop passes through
### v0 is the trailing vertex
### v1 is the pivot vertex
### v2 is the leading vertex


def drill(tri, loop, angle = None, branch = None, sig = None): # sig just for diagnostics
    if angle != None:
        face_coorientations = is_transverse_taut(tri, angle, return_type = "face_coorientations")
        assert face_coorientations != False
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
        assert tet_vert_coorientations != False

    original_tri = regina.Triangulation3(tri)

    original_countTetrahedra = tri.countTetrahedra()
    original_countBoundaryComponents = tri.countBoundaryComponents()
    ### add new tetrahedra
    new_0_tets = [] 
    new_1_tets = [] ## both relative to regina's two embeddings for the face
    for i in range(len(loop)):
        new_0_tets.append(tri.newTetrahedron())
        new_1_tets.append(tri.newTetrahedron())

    ### we will glue tetrahedra together with a numbering that is convenient but 
    ### unfortunately not oriented. We will orient later.

    #           1  pivot
    #          /|\
    #         / | \
    #        / ,3. \
    #       /,'   `.\
    #      0---------2 leading
    # trailing

    for i in range(len(loop)):
        new_0_tets[i].join(1, new_1_tets[i], regina.Perm4(0,1,2,3))

    ### now glue along the loop path, need to worry about regina's embeddings for neighbouring triangles

    loop_face_tets0 = []
    loop_face_vertices0 = []  
    loop_face_tets1 = []
    loop_face_vertices1 = []   ### need to store these because they cannot be recomputed once we start ungluing faces

    for i in range(len(loop)):
        # print('collect info: i', i)
        face_data = loop[i]
        face_index = face_data[0]
        vert_nums = face_data[1]
        face = tri.triangles()[face_index]
        face_embed0 = face.embedding(0) 
        face_tet = face_embed0.simplex()
        face_vertices = face_embed0.vertices()
        face_opposite_vert = face_vertices[3]
        face_other_non_edge_vert = face_vertices[vert_nums[0]] ## opposite trailing vertex

        face_data_next = loop[(i+1)%len(loop)]
        face_index_next = face_data_next[0]
        vert_nums_next = face_data_next[1]
        face_next = tri.triangles()[face_index_next]
        face_next_embed0 = face_next.embedding(0) 
        face_next_tet = face_next_embed0.simplex()
        face_next_vertices = face_next_embed0.vertices()
        face_next_opposite_vert = face_next_vertices[3]
        face_next_other_non_edge_vert = face_next_vertices[vert_nums_next[2]] ## opposite leading vertex

        ### things to store for later

        loop_face_tets0.append(face_tet) # store for later
        loop_face_vertices0.append(face_vertices) #store for later
        face_embed1 = face.embedding(1) 
        loop_face_tets1.append( face_embed1.simplex() )
        loop_face_vertices1.append( face_embed1.vertices() )

        edge = face.edge(vert_nums[0]) ## opposite trailing vertex
        # print('edge index', edge.index())
        assert edge == face_next.edge(vert_nums_next[2]) ## opposite leading vertex

        edgemapping = face.faceMapping(1,vert_nums[0])
        next_edgemapping = face_next.faceMapping(1,vert_nums_next[2])
        face_gluing_regina_numbering = next_edgemapping * (edgemapping.inverse()) ### maps vertices 0,1,2 on face to corresponding vertices on face_next
        assert face_gluing_regina_numbering[3] == 3
        face_gluing_regina_numbering = regina.Perm3(face_gluing_regina_numbering[0], face_gluing_regina_numbering[1], face_gluing_regina_numbering[2])
        assert face_gluing_regina_numbering[vert_nums[0]] == vert_nums_next[2]
        face_gluing = vert_nums_next.inverse() * face_gluing_regina_numbering * vert_nums
        assert face_gluing[0] == 2

        signs = []
        edge_embeddings = edge.embeddings()
        for embed in edge_embeddings:
            if embed.simplex() == face_tet:
                if set([embed.vertices()[2], embed.vertices()[3]]) == set([face_opposite_vert, face_other_non_edge_vert]):
                    # print('embed data face_tet', embed.simplex().index(), embed.vertices())
                    # print('embed edge vert nums', embed.vertices()[0], embed.vertices()[1])
                    signs.append( face_opposite_vert == embed.vertices()[2] )
            if embed.simplex() == face_next_tet:
                if set([embed.vertices()[2], embed.vertices()[3]]) == set([face_next_opposite_vert, face_next_other_non_edge_vert]):
                    # print('embed data face_next_tet', embed.simplex().index(), embed.vertices())
                    # print('embed edge vert nums', embed.vertices()[0], embed.vertices()[1])
                    signs.append( face_next_opposite_vert == embed.vertices()[2] )
        # print('signs', signs)
        assert len(signs) == 2
        if signs[0] == signs[1]:  ### coorientations are same around the edge (not a transverse taut coorientation!)
            new_0_tets[i].join(0, new_1_tets[(i+1)%len(loop)], regina.Perm4(2,face_gluing[1],face_gluing[2],3))
            new_1_tets[i].join(0, new_0_tets[(i+1)%len(loop)], regina.Perm4(2,face_gluing[1],face_gluing[2],3))
        else:
            new_0_tets[i].join(0, new_0_tets[(i+1)%len(loop)], regina.Perm4(2,face_gluing[1],face_gluing[2],3))
            new_1_tets[i].join(0, new_1_tets[(i+1)%len(loop)], regina.Perm4(2,face_gluing[1],face_gluing[2],3))

    ### now unglue tri along the loop and glue in the new tetrahedra

    for i in range(len(loop)):
        # print('modify triangulation: i', i)
        vert_nums = loop[i][1]

        face_tet0 = loop_face_tets0[i]
        face_vertices0 = loop_face_vertices0[i] 
        face_tet1 = loop_face_tets1[i]
        face_vertices1 = loop_face_vertices1[i] 
    
        face_opposite_vert0 = face_vertices0[3]
        face_tet0.unjoin(face_opposite_vert0)
    ### glue torus shell to the old tetrahedra

        vert_nums_Perm4 = regina.Perm4(vert_nums[0], vert_nums[1], vert_nums[2], 3)
        new_0_tets[i].join(3, face_tet0, face_vertices0 * vert_nums_Perm4)
        new_1_tets[i].join(3, face_tet1, face_vertices1 * vert_nums_Perm4)

    assert tri.isValid()
    assert tri.countBoundaryComponents() == original_countBoundaryComponents + 1

    if angle != None:
        for i in range(len(loop)):
            face_data = loop[i]
            face_index = face_data[0]
            face = tri.triangles()[face_index]
            vert_nums = face_data[1]

            # flow_agrees_with_regina_numbers = (((vert_nums[0] + 2) % 3) == vert_nums[2])
            # face_cor_agrees_with_regina_numbers = (face_coorientations[face_index] == +1)

            face_embed0 = face.embedding(0) 
            face_tet = face_embed0.simplex()
            face_vertices = face_embed0.vertices()
            face_opposite_vert = face_vertices[3]
            coor_points_out_of_tet0 = (tet_vert_coorientations[face_tet.index()][face_opposite_vert] == +1)

            # if (flow_agrees_with_regina_numbers != face_cor_agrees_with_regina_numbers) != coor_points_out_of_tet0:
            if coor_points_out_of_tet0:
                angle.extend([0,2])
            else:
                angle.extend([2,0])
        # assert is_taut(tri, angle)
        # print('taut?', is_taut(tri, angle), angle)
        
        if not is_taut(tri, angle):
            print(sig, 'our angle', angle, 'is taut', is_taut(tri, angle))
            neighbouring_tets = []
            for i in range(len(loop)):
                face_data = loop[i]
                face_index = face_data[0]
                face = original_tri.triangles()[face_index]
                neighbouring_tets.append( (face.embedding(0).simplex().index(), face.embedding(1).simplex().index()) )
            print(neighbouring_tets)
            tri.save(sig + '_' + str(loop) + '_not_taut.rga')
            # assert False


        # angle.extend([0,2] * len(loop))  ## wrong: lower and upper need to be relative to the taut coorientation, not the embed ordering



    # ### now orient
    # swaps = [regina.Perm4()] * original_countTetrahedra ### identity permutations
    # for i in range(len(loop)):
    #     if new_lower_tets[i].adjacentGluing(3).sign() == 1:
    #         swaps.append( reverse_tet_orientation(tri, new_lower_tets[i], 0) )  ### pi_location = 0 is an arbitrary choice
    #     else:
    #         swaps.append( regina.Perm4() )
    #     if new_upper_tets[i].adjacentGluing(3).sign() == 1:
    #         swaps.append( reverse_tet_orientation(tri, new_upper_tets[i], 0) ) ### pi_location = 0 is an arbitrary choice
    #     else:
    #         swaps.append( regina.Perm4() )

    


def test():
    # tri, angle = isosig_to_tri_angle('cPcbbbiht_12')
    # tri.save('cPcbbbiht_12.rga')
    # loop = [tet_to_face_data(tri, 0, 3, [2,0,1]), tet_to_face_data(tri, 0, 1, [0,2,3])]
    # print( loop )
    # drill(tri, loop)
    # tri.save('cPcbbbiht_12_drilled_short.rga')

    # loop = [tet_to_face_data(tri, 0, 3, [0,2,1]), tet_to_face_data(tri, 0, 0, [2,1,3]), tet_to_face_data(tri, 1, 3, [1,2,0]), tet_to_face_data(tri, 1, 1, [2,0,3])]
    # print( loop )  ### goes through a face twice
    # drill(tri, loop)
    # tri.save('cPcbbbiht_12_drilled_upwards.rga') 

    # tri, angle = isosig_to_tri_angle('cPcbbbdxm_10')
    # loop = [tet_to_face_data(tri, 0, 3, [0,2,1])]
    # print( loop )
    # drill(tri, loop)
    # tri.save('cPcbbbdxm_10_drilled_mobius.rga')

    tri, angle = isosig_to_tri_angle('gLLAQbecdfffhhnkqnc_120012')
    loop = [(6, regina.Perm3(0,1,2)), (2, regina.Perm3(0,1,2))]
    drill(tri, loop, angle = angle)
    tri.save('s227_drilled_upwards.rga')


