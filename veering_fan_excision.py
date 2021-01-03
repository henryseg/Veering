#
# veering_drill_midsurface_bdy.py
#

### given a veering triangulation, produce an ideal triangulation that is the result of drilling out the 
### boundary curves of the midsurfaces. The fan and toggle tetrahedra have different decompositions after drilling, 
### see the file drill_fans_and_toggles_along_bdyS.3dm

import regina #needed inside of imported files
from transverse_taut import get_top_and_bottom_nums
from taut import liberal, isosig_to_tri_angle
from veering import veering_triangulation

@liberal
def fan_stacks(tri, angle):
    """Find lists of tetrahedra that make up the fans, including the toggles at either end"""
    vt = veering_triangulation(tri, angle)
    veering_colours = vt.veering_colours  ## "red" or "blue"
    tet_types = vt.tet_types  ### "toggle", "red" or "blue"
    tet_vert_coors = vt.coorientations
    toggle_nums = [i for i in range(len(tet_types)) if tet_types[i] == "toggle"]
    
    out = []
    for toggle_num in toggle_nums:
        tops, bottoms = get_top_and_bottom_nums(tet_vert_coors, toggle_num)
        for i in range(2):

            tet_nums = [toggle_num]
            trailing_vertex = bottoms[i]
            other_b = bottoms[(i+1)%2]
            t0, t1 = tops
            if vt.get_edge_between_verts_colour(toggle_num, (t0, t1)) == vt.get_edge_between_verts_colour(toggle_num, (t0, other_b)):
                leading_vertex = t0
            else:
                leading_vertex = t1
            
            tet = tri.tetrahedron(toggle_num)
            while True:
                gluing = tet.adjacentGluing(trailing_vertex)
                tet = tet.adjacentTetrahedron(trailing_vertex)
                # e0, e1 = gluing[e0], gluing[e1]
                leading_vertex, trailing_vertex = gluing[trailing_vertex], gluing[leading_vertex]
                tet_nums.append(tet.index())
                if tet_types[tet.index()] == "toggle":
                    break
            out.append(tet_nums)
    return out

@liberal
def excise_fans(tri, angle, fan_nums = None):
    vt = veering_triangulation(tri, angle)
    veering_colours = vt.veering_colours  ## "red" or "blue"
    tet_types = vt.tet_types  ### "toggle", "red" or "blue"
    tet_vert_coors = vt.coorientations
    if fan_nums == None: ## do all fans
        fan_nums = [n for n in range(len(tet_types)) if tet_types[n] != "toggle"]
    
    excisedAngle = angle[:]
    for fan_num in sorted(fan_nums, reverse = True):
        del excisedAngle[fan_num]

    minority_edge_pairs = []
    for fan_num in fan_nums:
        assert tet_types[fan_num] != "toggle"
        tops, bottoms = get_top_and_bottom_nums(tet_vert_coors, fan_num)
        if 0 in tops:
            pi_pair = list(tops)
        else:
            pi_pair = list(bottoms)
        other = pi_pair[(pi_pair.index(0) + 1) % 2]
        for i in range(1,4):
            if i != other:
                if vt.get_edge_between_verts_colour(fan_num, (0,i)) != tet_types[fan_num]: ### "red" or "blue"
                    last = 6 - i - other ## the fourth vertex
                    minority_edge_pairs.append([(0,i),(other,last)])
                    break
    excisedTri = regina.Triangulation3(tri) ### copy

    fan_tets = [excisedTri.tetrahedron(fan_num) for fan_num in fan_nums]
    for k, tet in enumerate(fan_tets):
        # print 'k', k, 'tet.index', tet.index()
        minority_edge_pair = minority_edge_pairs[k]
        # print minority_edge_pair
        ### record gluings for neighbours
        neighbours = [tet.adjacentSimplex(i) for i in range(4)]
        gluings = [tet.adjacentGluing(i) for i in range(4)]

        # print [neighbour.index() for neighbour in neighbours]
        # print gluings

        tet.isolate()


        ## now glue neighbours to each other

        if tet not in neighbours:
            to_glue = [0,1,2,3]
            while to_glue != []:
                i = to_glue.pop()
                if i in minority_edge_pair[0]:
                    j = minority_edge_pair[0][(minority_edge_pair[0].index(i) + 1) % 2]
                else:
                    j = minority_edge_pair[1][(minority_edge_pair[1].index(i) + 1) % 2]
                
                teti = neighbours[i]
                tetj = neighbours[j]
                teti.join( gluings[i][i], tetj, gluings[j] * regina.Perm4(j,i) * (gluings[i].inverse()) )
                to_glue.remove(j)
                # print i,j
        else: ### tet is glued to itself, which makes it trickier to remove
            ### first, find self-gluings
            self_gluings = []
            self_gluings = [neighbours.index(tet)]
            self_gluings.append(neighbours.index(tet, self_gluings[0] + 1))  ### add second entry

            other_gluings = [0,1,2,3]
            other_gluings.remove(self_gluings[0])
            other_gluings.remove(self_gluings[1])
            i, j = other_gluings
            teti = neighbours[i]
            tetj = neighbours[j]
            if i in minority_edge_pair[0]:
                p = minority_edge_pair[0][(minority_edge_pair[0].index(i) + 1) % 2]
                q = minority_edge_pair[1][(minority_edge_pair[1].index(j) + 1) % 2]
            else:
                p = minority_edge_pair[1][(minority_edge_pair[1].index(i) + 1) % 2]
                q = minority_edge_pair[0][(minority_edge_pair[0].index(j) + 1) % 2]
            teti.join( gluings[i][i], tetj, gluings[j] * regina.Perm4(j,q) * gluings[p] * regina.Perm4(p,i) * (gluings[i].inverse()) )

        excisedTri.removeTetrahedron(tet)

    return excisedTri, excisedAngle




if __name__ == '__main__':
    # sig = 'cPcbbbiht_12'
    # sig = 'cPcbbbdxm_10'
    # sig = 'dLQacccjsnk_200'
    # sig = 'dLQbccchhfo_122'
    # sig = 'dLQbccchhsj_122'
    # sig = 'eLMkbcddddedde_2100'
    # sig = 'eLMkbcdddhxqlm_1200'
    # sig = 'eLAkaccddjsnak_2001'
    # sig = 'eLAkbbcdddhwqj_2102'
    # sig = 'fLLQcbcdeeelonlel_02211'  ## should have three different midsurface bdy components, so four total after drilling
    # sig = 'fLLQcbddeeehrcdui_12000' ## five boundary components
    # sig = 'fLAMcbbcdeedhwqqs_21020'
    # sig = 'qvvLPQwAPLQkfhgffijlkmnopoppoaaaaaaaaaaaaadddd_1020212200111100'
    # t = excise_fans(sig)
    # t.save('excise_fans_' + sig + '.rga')

    # t0, _ = isosig_to_tri_angle('cPcbbbdxm_10')
    # t1, _ = isosig_to_tri_angle('cPcbbbiht_12')

    # print t.isIsomorphicTo(t0) != None or t.isIsomorphicTo(t1) != None

    # print fan_stacks(sig)

    import veering
    import transverse_taut
    import taut
    from file_io import parse_data_file
    census = parse_data_file('Data/veering_census.txt')

    for sig in census[:200]:
        tri, angle = excise_fans(sig)
        print( (sig, tri.countTetrahedra(), angle, taut.is_taut(tri, angle)) )
        assert veering.is_veering(tri, angle)
        assert transverse_taut.is_transverse_taut(tri, angle)
        # print sig, fan_stacks(sig)



