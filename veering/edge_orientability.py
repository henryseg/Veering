#
# edge_orientability.py
#

from sage.modules.free_module_element import vector
from sage.matrix.constructor import Matrix

import regina

from .taut import liberal, vert_pair_to_edge_num
from .transverse_taut import is_transverse_taut, get_tet_top_vert_nums
from .taut_homology import faces_in_homology
from .veering_tri import is_veering, loop_twistednesses


@liberal
def is_edge_orientable(tri, angle, return_type = "boolean"):
    """
    Determines if the veering triangulation is edge-orientable.  

    Example: 

    sage: from veering.edge_orientability import is_edge_orientable
    sage: sig = "cPcbbbdxm_10"
    sage: is_edge_orientable(sig)
    False
    sage: sig = "cPcbbbiht_12"
    sage: is_edge_orientable(sig)
    True
    """
    lts = loop_twistednesses(tri, angle)
    return all([lt == 1 for lt in lts])

# 2022-09-10 I checked that the new version gives the same answer as
# the old version (below) on all manifolds in the census.

def regina_edge_orientation_agrees(tet, vert_pair):
    """
    Given tet and an ordered pair of (regina) vert nums of that tet,
    determines if this ordering agrees with regina's ordering of the
    verts of that edge of the triangulation.
    """
    edge_num = vert_pair_to_edge_num[tuple(vert_pair)]
    mapping = tet.faceMapping(1, edge_num)
    map_order = [mapping[0], mapping[1]]
    assert set(map_order) == set(vert_pair)
    return map_order == vert_pair

@liberal
def is_edge_orientable_old(tri, angle, return_type = "boolean"):
    """
    Checks to see if this veering triangulation is edge-orientable.
    If return type is "tri_angle" it returns the edge-orientable
    double cover with its angle structure.  Note that this is
    disconnected if and only if the given triangulation is
    edge-orientable.
    """
    # return type can be "boolean", "veering_tet_vert_nums", or "tri_angle"
    n = tri.countTetrahedra()
    veering_colours = is_veering(tri, angle, return_type = "veering_colours")
    assert veering_colours != False # so we are veering
    tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")

    ### assumption: the first n tetrahedra have upper edge oriented
    ### according to Regina, the last n have it against Regina

    ### build our own model vertex numbering for each tetrahedron as
    ### follows: the top edge e is oriented by regina numbering, and
    ### gets vert_nums 1 and 2 in our ordering an equatorial edge e'
    ### of the same colour as e shares a vertex v with it.  e and e'
    ### both point away from v or towards it. If away then the other
    ### end of e' is 3, else, the other end is 0

    veering_tet_vert_nums = []  ### will populate with regina's vert
    ### nums, but our order.  That is, veering_tet_vert_nums[1] and
    ### [2] will be the regina vert nums for the ends of the top edge

    for i in range(n):
        tet = tri.tetrahedron(i)
        veering_vert_nums = [None, None, None, None]
        top_vert_pair = get_tet_top_vert_nums(tet_vert_coorientations, i)
        # print(i, top_vert_pair)
        if not regina_edge_orientation_agrees(tet, top_vert_pair):
            top_vert_pair.reverse()
            assert regina_edge_orientation_agrees(tet, top_vert_pair)
        veering_vert_nums[1] = top_vert_pair[0]
        veering_vert_nums[2] = top_vert_pair[1]

        top_edge_num = vert_pair_to_edge_num[tuple(top_vert_pair)]
        top_edge_col = veering_colours[ tet.edge(top_edge_num).index() ]

        bottom_vert_pair = list(set(range(4)) - set(top_vert_pair))
        bv, bv2 = bottom_vert_pair ## choose arbitrarily, now find
                                   ## which edge from the top vertices
                                   ## has same colour as bv
        for j, tv in enumerate(top_vert_pair):
            edge_col = veering_colours[ tet.edge(vert_pair_to_edge_num[(bv, tv)]).index() ]
            if edge_col == top_edge_col:
                if j == 0: ### bv shares an edge of same colour as top with tail of top edge
                    veering_vert_nums[3] = bv
                    veering_vert_nums[0] = bv2
                else: ### bv shares an edge of same colour as top with head of top edge
                    veering_vert_nums[0] = bv
                    veering_vert_nums[3] = bv2
        veering_tet_vert_nums.append(veering_vert_nums)
    # print('veering_tet_vert_nums', veering_tet_vert_nums)
    if return_type == "veering_tet_vert_nums":
        return veering_tet_vert_nums
    
    ### Now, when we glue two tetrahedra together along a face, the
    ### first of the three vertices in the veering_vert_num order on
    ### that tet's face glues to the first of the three vertices on
    ### the other tet's face, or to the third.  depending on this, we
    ### go to the other part of the double cover, or not

    cover_tri = regina.Triangulation3()  ## starts empty
    for i in range(2*n):
        cover_tri.newTetrahedron()
    tet_faces = []
    for i in range(n):
        for j in range(4):
            tet_faces.append((i,j))
    while len(tet_faces) > 0:
        i, j = tet_faces.pop()
        tet = tri.tetrahedron(i)
        adjtet, adjgluing = tet.adjacentTetrahedron(j), tet.adjacentGluing(j)
        iN, jN = adjtet.index(), adjgluing[j]
        tet_faces.remove( (iN, jN) )

        cover_tets = [cover_tri.tetrahedron(i), cover_tri.tetrahedron(i+n)]
        cover_tetsN = [cover_tri.tetrahedron(iN), cover_tri.tetrahedron(iN+n)]
        ### find the veering indices for the verts in the gluing on
        ### this tet and on adjtet

        face_veering_nums = veering_tet_vert_nums[i][:]
        face_veering_nums.remove(j)
        a, b, c = face_veering_nums
        neighbour_face_veering_nums = veering_tet_vert_nums[iN][:]
        neighbour_face_veering_nums.remove(jN)
        aN, bN, cN = neighbour_face_veering_nums
        assert adjgluing[b] == bN ### middles should match
        if adjgluing[a] == aN:  ### veering orderings agree across the
                                ### gluing
            assert adjgluing[c] == cN
            for k in range(2):
                cover_tets[k].join(j, cover_tetsN[k], adjgluing)
        else: ### veering orderings disagree across the gluing
            assert adjgluing[a] == cN and adjgluing[c] == aN
            for k in range(2):
                cover_tets[k].join(j, cover_tetsN[(k+1)%2], adjgluing)
    assert not cover_tri.hasBoundaryFacets()
    assert is_veering(cover_tri, angle + angle)

    if return_type == "boolean":
        return not cover_tri.isConnected() ### not connected if the
                                           ### original veering
                                           ### triangulation is
                                           ### edge-orientable
    else:
        assert return_type == "tri_angle"
        return cover_tri, angle+angle

@liberal
def is_max_fab_edge_orientable(tri, angle, return_type = "boolean"):
    """
    Determines if lift of the veering triangulation, to the maximal
    free abelian cover, is edge-orientable.  

    Example: 

    sage: from veering.edge_orientability import is_fab_edge_orientable
    sage: sig = "cPcbbbdxm_10"
    sage: is_fab_edge_orientable(sig)
    True
    sage: sig = "cPcbbbiht_12"
    sage: is_fab_edge_orientable(sig)
    True
    sage: sig = "gLLAQbecdfffhhnkqnc_120012"
    sage: is_fab_edge_orientable(sig)
    False
    """
    # See the 2022-09-06 email from Anna Parlak

    # The twist of the dual loop.
    lt = loop_twistednesses(tri, angle)
    # Additive notation is more useful here, so:
    lt = [int((t - 1)/(-2)) for t in lt]
    lt = vector(lt)

    # image of the dual loop in H_1/torsion
    fh = faces_in_homology(tri, angle, [])
    FH = Matrix(fh).transpose()

    # The (right) kernel of FH is the lattice of "loops" that do not
    # unwrap in max-fab.
    K = FH.right_kernel()

    # There is an element of K which is twisted if and only if there
    # is a generator in K.gens() which is twisted.  Also, an element k
    # is untwisted if and only if the sum of its twists (lt*k), modulo
    # two, is zero. So:
    return all([((lt*k % 2) == 0) for k in K.gens()])
