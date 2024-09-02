#
# veering_drill_midsurface_bdy.py
#

### given a veering triangulation, produce an ideal triangulation that is the result of drilling out the 
### boundary curves of the midsurfaces. The fan and toggle tetrahedra have different decompositions after drilling, 
### see the file drill_fans_and_toggles_along_bdyS.3dm

import regina #needed inside of imported files

from .transverse_taut import is_transverse_taut, get_top_and_bottom_vert_nums
from .taut import liberal, is_taut, unsorted_vert_pair_to_edge_pair
from .veering_tri import is_veering, veering_triangulation

@liberal
def drill_midsurface_bdy(tri, angle):
    vt = veering_triangulation(tri, angle)
    veering_colours = vt.veering_colours  ## "red" or "blue"
    tet_types = vt.tet_types  ### "toggle", "red" or "blue"
    tet_vert_coors = vt.coorientations

    drilledTri = regina.Triangulation3()
    subtet_indices = []
    subtet_addresses = []

    num_toggles = tet_types.count("toggle")
    drilled_num_tet = 8 * num_toggles + 4 * (tri.countTetrahedra() - num_toggles)

    meridians = [] 
    ### in each toggle, we record the meridians for the two boundary components we are drilling
    ### note that this repeats the same meridian many times - the construction here is purely local: it is much
    ### easier to figure out which are duplicates of which later

    for i in range(tri.countTetrahedra()):
        tet = tri.tetrahedron(i)
        first_subtet_index = drilledTri.countTetrahedra()
        if tet_types[i] == 'toggle':  
            num_to_add = 8
        else:
            num_to_add = 4
        for j in range(num_to_add):
            drilledTri.newTetrahedron()

        subtet_indices.append(range(first_subtet_index, first_subtet_index + num_to_add))

        ### glue sub-tetrahedra together. Note that some tetrahedra will be negatively oriented

        [(t0,t1), (b0,b1)] = get_top_and_bottom_vert_nums(tet_vert_coors, i)
        this_tet_subtet_addresses = {}
        if tet_types[i] == 'toggle':  
            ## first two subtetrahedra are top two, second two are bottom two, then the four side tetrahedra
            tet_t0, tet_t1 = drilledTri.tetrahedron(first_subtet_index), drilledTri.tetrahedron(first_subtet_index + 1)
            # tet_ti is opposite vertex bi
            tet_b0, tet_b1 = drilledTri.tetrahedron(first_subtet_index + 2), drilledTri.tetrahedron(first_subtet_index + 3)
            # tet_bi is opposite vertex ti
            tet_s00 = drilledTri.tetrahedron(first_subtet_index + 4)
            tet_s01 = drilledTri.tetrahedron(first_subtet_index + 5)
            tet_s10 = drilledTri.tetrahedron(first_subtet_index + 6)
            tet_s11 = drilledTri.tetrahedron(first_subtet_index + 7)
            # tet_sij meets vertices ti and bj
            
            ## keys are (vert not on face, vert not on edge), returns the subtet which meets the face, edge 
            this_tet_subtet_addresses[(b0, b1)] = tet_t0
            this_tet_subtet_addresses[(b1, b0)] = tet_t1
            this_tet_subtet_addresses[(t0, t1)] = tet_b0
            this_tet_subtet_addresses[(t1, t0)] = tet_b1
            this_tet_subtet_addresses[(t1, b1)] = tet_s00
            this_tet_subtet_addresses[(b1, t1)] = tet_s00
            this_tet_subtet_addresses[(t1, b0)] = tet_s01
            this_tet_subtet_addresses[(b0, t1)] = tet_s01
            this_tet_subtet_addresses[(t0, b1)] = tet_s10
            this_tet_subtet_addresses[(b1, t0)] = tet_s10
            this_tet_subtet_addresses[(t0, b0)] = tet_s11
            this_tet_subtet_addresses[(b0, t0)] = tet_s11

            tet_t0.join(b0, tet_t1, regina.Perm4(b0, b1, b1, b0, t0, t0, t1, t1))  ## b0 <-> b1, t0 and t1 fixed
            tet_b0.join(t0, tet_b1, regina.Perm4(t0, t1, t1, t0, b0, b0, b1, b1))  ## t0 <-> t1, b0 and b1 fixed

            tet_s00.join(b0, tet_t1, regina.Perm4(t0, t0, b1, b1, t1, b0, b0, t1))
            tet_s00.join(t0, tet_b1, regina.Perm4(b0, b0, t1, t1, b1, t0, t0, b1))

            tet_s01.join(b1, tet_t0, regina.Perm4(t0, t0, b0, b0, t1, b1, b1, t1))
            tet_s01.join(t0, tet_b1, regina.Perm4(b1, b1, t1, t1, b0, t0, t0, b0))

            tet_s10.join(b0, tet_t1, regina.Perm4(t1, t1, b1, b1, t0, b0, b0, t0))
            tet_s10.join(t1, tet_b0, regina.Perm4(b0, b0, t0, t0, b1, t1, t1, b1))

            tet_s11.join(b1, tet_t0, regina.Perm4(t1, t1, b0, b0, t0, b1, b1, t0))
            tet_s11.join(t1, tet_b0, regina.Perm4(b1, b1, t0, t0, b0, t1, t1, b0))

            ### meridian around upper hole
            meridian = [0]*(3*drilled_num_tet)
            meridian[ 3 * tet_t1.index() + unsorted_vert_pair_to_edge_pair[b0, t0] ] = -1
            meridian[ 3 * tet_s00.index() + unsorted_vert_pair_to_edge_pair[b1, t1] ] = 1
            meridian[ 3 * tet_b1.index() + unsorted_vert_pair_to_edge_pair[t0, t1] ] = 1
            meridian[ 3 * tet_s01.index() + unsorted_vert_pair_to_edge_pair[b0, t1] ] = 1
            meridian[ 3 * tet_t0.index() + unsorted_vert_pair_to_edge_pair[b1, t0] ] = -1
            meridians.append(meridian)

            ### meridian around lower hole - swap b with t everywhere. Note this also swaps s01 with s10
            meridian = [0]*(3*drilled_num_tet)
            meridian[ 3 * tet_b1.index() + unsorted_vert_pair_to_edge_pair[t0, b0] ] = -1
            meridian[ 3 * tet_s00.index() + unsorted_vert_pair_to_edge_pair[t1, b1] ] = 1
            meridian[ 3 * tet_t1.index() + unsorted_vert_pair_to_edge_pair[b0, b1] ] = 1
            meridian[ 3 * tet_s10.index() + unsorted_vert_pair_to_edge_pair[t0, b1] ] = 1
            meridian[ 3 * tet_b0.index() + unsorted_vert_pair_to_edge_pair[t1, b0] ] = -1
            meridians.append(meridian)

        else: ## fan
            # first two tetrahedra are the top and bottom
            tet_t, tet_b = drilledTri.tetrahedron(first_subtet_index), drilledTri.tetrahedron(first_subtet_index + 1)
            # second two tetrahedra are the sides, s0 meets vertex t0, s1 meets vertex t1
            tet_s0, tet_s1 = drilledTri.tetrahedron(first_subtet_index + 2), drilledTri.tetrahedron(first_subtet_index + 3)

            # tet_types[i] == 'blue'
            this_tet_subtet_addresses[(b0, b1)] = tet_t
            this_tet_subtet_addresses[(b1, b0)] = tet_t
            this_tet_subtet_addresses[(t0, t1)] = tet_b
            this_tet_subtet_addresses[(t1, t0)] = tet_b
            ### if this is a blue fan tet then each side tet meet a blue edge
            # tet_types[i] could be 'blue' or 'red'

            ### find which bottom vertex, when linked to t0, gives an edge of the correct colour
            if vt.get_edge_between_verts_colour(i, (t0, b0)) == tet_types[i]:
                s0b = b0  ## bottom vert of s0
                s1b = b1  ## bottom vert of s1
                this_tet_subtet_addresses[(t0, b0)] = tet_s1
                this_tet_subtet_addresses[(b0, t0)] = tet_s1
                this_tet_subtet_addresses[(t1, b1)] = tet_s0
                this_tet_subtet_addresses[(b1, t1)] = tet_s0
            else:
                s0b = b1  ## bottom vert of s0
                s1b = b0  ## bottom vert of s1
                this_tet_subtet_addresses[(t0, b1)] = tet_s1
                this_tet_subtet_addresses[(b1, t0)] = tet_s1
                this_tet_subtet_addresses[(t1, b0)] = tet_s0
                this_tet_subtet_addresses[(b0, t1)] = tet_s0  

            tet_s0.join(s0b, tet_t, regina.Perm4(t0, t0, s0b, t1, s1b, s1b, t1, s0b))
            tet_s0.join(t0,  tet_b, regina.Perm4(s0b, s0b, t0, s1b, t1, t1, s1b, t0))
            tet_s1.join(s1b, tet_t, regina.Perm4(t1, t1, s1b, t0, s0b, s0b, t0, s1b))
            tet_s1.join(t1,  tet_b, regina.Perm4(s1b, s1b, t1, s0b, t0, t0, s0b, t1))

        subtet_addresses.append(this_tet_subtet_addresses)

    ### now glue subtetrahedra from different original tetrahedra together
    ### fan tetrahedra only glue two of three subtetrahedra on each face. 
    ### toggles glue one of their subfaces all the way through to the next toggle...

    unglued_flags = []
    for f in range(tri.countTriangles()):
        for e in range(3):
            unglued_flags.append((f,e))

    while unglued_flags != []:
        (f,e) = unglued_flags.pop()
        # print f,e
        face = tri.triangle(f)
        embed0 = face.embedding(0)
        embed1 = face.embedding(1)
        tet_index_0 = embed0.simplex().index()
        tet_index_1 = embed1.simplex().index()
        face0 = embed0.face()
        face1 = embed1.face()
        vertperm0 = embed0.vertices()
        vertperm1 = embed1.vertices()
        
        edge0 = vertperm0[e]
        edge1 = vertperm1[e]

        otherverts0 = [0,1,2,3]
        otherverts0.remove(face0)
        otherverts0.remove(edge0)
        otherverts1 = [0,1,2,3]
        otherverts1.remove(face1)
        otherverts1.remove(edge1)

        if (face0, edge0) in subtet_addresses[tet_index_0]:
            subtet0 = subtet_addresses[tet_index_0][(face0, edge0)]
        else:
            subtet0 = None
        if (face1, edge1) in subtet_addresses[tet_index_1]:
            subtet1 = subtet_addresses[tet_index_1][(face1, edge1)]
        else:
            subtet1 = None

        if subtet0 == None and subtet1 == None:
            pass ### both are fans, with no subtet, skip
        elif subtet0 != None and subtet1 != None:
            ### glue
            u, v = otherverts0
            tet0perm = regina.Perm4(face0, edge0, edge0, face0, u, u, v, v)
            u, v = otherverts1
            tet1perm = regina.Perm4(face1, edge1, edge1, face1, u, u, v, v)
            gluing = embed0.simplex().adjacentGluing(face0)
            subtet0.join(edge0, subtet1, tet1perm*gluing*tet0perm) ### perms act on left

        else: ### we have to walk around to find the right place to glue this toggle subtet, and remove the other unglued flag from the list
            if subtet1 == None:
                assert tet_types[tet_index_0] == 'toggle'
                toggle_tet_index = tet_index_0
                toggle_face = face0
                toggle_edge = edge0
                subtet = subtet0
            else:
                assert tet_types[tet_index_1] == 'toggle'
                toggle_tet_index = tet_index_1
                toggle_face = face1
                toggle_edge = edge1
                subtet = subtet1
            edge_verts = [0,1,2,3]
            edge_verts.remove(toggle_edge)
            edge_verts.remove(toggle_face)
            toggle_e0, toggle_e1 = edge_verts

            e0, e1 = toggle_e0, toggle_e1
            leading_vertex = toggle_edge
            trailing_vertex = toggle_face
            tet = tri.tetrahedron(toggle_tet_index)
            while True:
                gluing = tet.adjacentGluing(trailing_vertex)
                tet = tet.adjacentTetrahedron(trailing_vertex)
                e0, e1 = gluing[e0], gluing[e1]
                leading_vertex, trailing_vertex = gluing[trailing_vertex], gluing[leading_vertex]
                if tet_types[tet.index()] == "toggle":
                    break
            other_toggle_tet_index = tet.index()
            other_toggle_face = leading_vertex
            other_toggle_edge = trailing_vertex
            other_toggle_e0 = e0
            other_toggle_e1 = e1
            other_subtet = subtet_addresses[other_toggle_tet_index][(other_toggle_face, other_toggle_edge)]
            tetperm = regina.Perm4(toggle_face, toggle_edge, toggle_edge, toggle_face, toggle_e0, toggle_e0, toggle_e1, toggle_e1)
            other_tetperm = regina.Perm4(other_toggle_face, other_toggle_edge, other_toggle_edge, other_toggle_face, other_toggle_e0, other_toggle_e0, other_toggle_e1, other_toggle_e1)
            gluing = regina.Perm4(toggle_e0, other_toggle_e0, toggle_e1, other_toggle_e1, toggle_face, other_toggle_face, toggle_edge, other_toggle_edge)
            subtet.unjoin(toggle_edge)  # 2023-05-28
            # other_subtet.unjoin(other_tetperm*gluing*tetperm[toggle_edge])  # seems not to be needed... 
            subtet.join(toggle_edge, other_subtet, other_tetperm*gluing*tetperm)

    return drilledTri, meridians




if __name__ == '__main__':
    sig = 'cPcbbbiht_12'
    # sig = 'cPcbbbdxm_10'
    # sig = 'dLQacccjsnk_200'
    # sig = 'dLQbccchhfo_122'
    # sig = 'dLQbccchhsj_122'
    # sig = 'eLMkbcddddedde_2100'
    # sig = 'eLMkbcdddhxqlm_1200'
    # sig = 'fLLQcbcdeeelonlel_02211'  ## should have three different midsurface bdy components, so four total after drilling
    # sig = 'fLLQcbddeeehrcdui_12000' ## five boundary components
    # sig = 'qvvLPQwAPLQkfhgffijlkmnopoppoaaaaaaaaaaaaadddd_1020212200111100'
    drilled_tri, meridians = drill_midsurface_bdy(sig)
    print(drilled_tri.isoSig())
    print(meridians)
    drilled_tri.save('drilled_tri_' + sig + '.rga')
    drilled_tri.saveSnapPea('drilled_tri_' + sig + '.tri')


