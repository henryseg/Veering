#
# flow_cycles.py
#

# Given a veering triangulation, find all simple flow cycles

import regina
from veering import is_veering
from branched_surface import upper_branched_surface, branch_num_to_large_edge_and_mixed_edge_pair_num, isosig_from_tri_angle_branch
from transverse_taut import get_tet_top_and_bottom_edges, get_tet_above_edge, is_transverse_taut
from taut import isosig_to_tri_angle, edge_num_to_vert_pair
from drill import drill
from taut_polytope import is_layered

from file_io import parse_data_file

### format for loops: it is a list of tuples, 
### each tuple is (tet_index, edge_index within this tet that we exit through)

def lex_smallest_loop(loop):
    ind = loop.index(min(loop))
    loop = loop[ind:] + loop[:ind]
    return loop

def add_flow_edge(tri, branch, tet_with_this_edge_large, loop_so_far, current_tet_index, start_edge_index, available_edge_indices, found_loops):
    large_edge_num, mixed_edge_pair_num, small_edge_pair_num = branch_num_to_large_edge_and_mixed_edge_pair_num(branch[current_tet_index], return_small = True)
    tet = tri.tetrahedron(current_tet_index)
    for i in range(3): ### three choices of how to add a flow edge
        if i == 0: ### 0th small edge
            new_edge = tet.edge(small_edge_pair_num)
            if new_edge.index() in available_edge_indices:
                loop_so_far2 = loop_so_far[:]
                loop_so_far2.append((current_tet_index, small_edge_pair_num))
                available_edge_indices2 = available_edge_indices[:]
                available_edge_indices2.remove(new_edge.index())
                if new_edge.index() == start_edge_index:
                    loop_so_far2 = lex_smallest_loop(loop_so_far2)
                    if loop_so_far2 not in found_loops:
                        found_loops.append(loop_so_far2)
                else:
                    add_flow_edge(tri, branch, tet_with_this_edge_large, loop_so_far2, tet_with_this_edge_large[new_edge.index()], start_edge_index, available_edge_indices2, found_loops)

        elif i == 1: ### 1th small edge
            new_edge = tet.edge(5 - small_edge_pair_num)
            if new_edge.index() in available_edge_indices:
                loop_so_far2 = loop_so_far[:]
                loop_so_far2.append((current_tet_index, 5 - small_edge_pair_num))
                available_edge_indices2 = available_edge_indices[:]
                available_edge_indices2.remove(new_edge.index())
                if new_edge.index() == start_edge_index:
                    loop_so_far2 = lex_smallest_loop(loop_so_far2)
                    if loop_so_far2 not in found_loops:
                        found_loops.append(loop_so_far2)
                else:
                    add_flow_edge(tri, branch, tet_with_this_edge_large, loop_so_far2, tet_with_this_edge_large[new_edge.index()], start_edge_index, available_edge_indices2, found_loops)
        elif i == 2: ### tiny edge
            new_edge = tet.edge(5 - large_edge_num) ### top edge
            if new_edge.index() in available_edge_indices:
                loop_so_far2 = loop_so_far[:]
                loop_so_far2.append((current_tet_index, 5 - large_edge_num))
                available_edge_indices2 = available_edge_indices[:]
                available_edge_indices2.remove(new_edge.index())
                if new_edge.index() == start_edge_index:
                    loop_so_far2 = lex_smallest_loop(loop_so_far2)
                    if loop_so_far2 not in found_loops:
                        found_loops.append(loop_so_far2)
                else:
                    add_flow_edge(tri, branch, tet_with_this_edge_large, loop_so_far2, tet_with_this_edge_large[new_edge.index()], start_edge_index, available_edge_indices2, found_loops)

def make_list_of_tet_with_this_edge_large(tri, branch):
    out = [None] * tri.countEdges()
    for i in range(tri.countTetrahedra()):
        large_edge, _ = branch_num_to_large_edge_and_mixed_edge_pair_num(branch[i])
        large_edge_index = tri.tetrahedron(i).edge(large_edge).index()
        out[large_edge_index] = i
    return out

def find_flow_cycles(tri, branch):
    """Given a branched surface dual to a triangulation, find all simple flow cycles in the upper flow graph"""
    tet_with_this_edge_large = make_list_of_tet_with_this_edge_large(tri, branch)
    found_loops = []
    for start_tet_index in range(tri.countTetrahedra()):
        large_edge, _ = branch_num_to_large_edge_and_mixed_edge_pair_num(branch[start_tet_index])
        large_edge_index = tri.tetrahedron(start_tet_index).edge(large_edge).index()
        add_flow_edge(tri, branch, tet_with_this_edge_large, [], start_tet_index, large_edge_index, list(range(tri.countEdges())), found_loops)
    return found_loops

def tet_to_face_data(tri, tet_num, face_num, vertices): ### vertices is list in order (trailing, pivot, leading) in numbering for the tet
    """given a tetrahedron, face_num and vertices in the labelling for the tetrahedron,
       convert to face_index and vertices in the labelling for the face.
       This is the correct format to feed into drill."""

    assert face_num not in vertices
    # print('tet_num, face_num, vertices (trailing, pivot, leading)', tet_num, face_num, vertices)
    tet = tri.tetrahedron(tet_num)
    facemapping = tet.faceMapping(2,face_num)
    face = tet.triangle(face_num)
    new_vertices = facemapping.inverse() * regina.Perm4(vertices[0], vertices[1], vertices[2], face_num)
    return (face.index(), regina.Perm3(new_vertices[0], new_vertices[1], new_vertices[2]))

def neighbouring_edges_to_triangle_loop_element(start_edge_num, end_edge_num, tri, tet_index):
    # print('start_edge_num, end_edge_num, tet_index', start_edge_num, end_edge_num, tet_index)
    start_vert_pair = edge_num_to_vert_pair[start_edge_num]
    end_vert_pair = edge_num_to_vert_pair[end_edge_num]
    start_vert = ( set(start_vert_pair).difference(set(end_vert_pair)) ).pop()
    pivot_vert = ( set(start_vert_pair).intersection(set(end_vert_pair)) ).pop()
    end_vert = ( set(end_vert_pair).difference(set(start_vert_pair)) ).pop()
    face_num = {0,1,2,3}.difference( {start_vert, pivot_vert, end_vert} ).pop()
    # print('neighbouring_edges_to_triangle_loop_element', tet_index, face_num, [start_vert, pivot_vert, end_vert])
    return tet_to_face_data(tri, tet_index, face_num, [start_vert, pivot_vert, end_vert]) 

def flow_cycle_to_triangle_loop(tri, branch, cycle):
    available_edges = list(range(tri.countEdges()))
    out = []
    for tet_index, edge_num in cycle:
        # print('tet_index, edge_num', tet_index, edge_num)
        large_edge_num, mixed_edge_pair_num, small_edge_pair_num = branch_num_to_large_edge_and_mixed_edge_pair_num(branch[tet_index], return_small = True)
        # print('large_edge_num, mixed_edge_pair_num, small_edge_pair_num', large_edge_num, mixed_edge_pair_num, small_edge_pair_num)
        assert edge_num != large_edge_num ### should be leaving large edge and going to the edge_num
        assert edge_num != mixed_edge_pair_num
        assert 5 - edge_num != mixed_edge_pair_num
        
        edge_index = tri.tetrahedron(tet_index).edge(edge_num).index()
        # print('edge_index', edge_index)
        if not edge_index in available_edges:
            return False ### give up, our path meets an edge more than once

        if edge_num == small_edge_pair_num or edge_num == 5 - small_edge_pair_num:
            out.append( neighbouring_edges_to_triangle_loop_element(large_edge_num, edge_num, tri, tet_index) )
            
        elif edge_num == 5 - large_edge_num:
            ### try going through both mixed edges, if we find one that works, go with it. Possibly this choice 
            ### messes up some later edge, so be it. We're just trying to generate examples, don't need to get everything possible
            mixed_edge_1_num = mixed_edge_pair_num
            mixed_edge_2_num = 5 - mixed_edge_pair_num
            mixed_edge_1_index = tri.tetrahedron(tet_index).edge(mixed_edge_1_num).index()
            mixed_edge_2_index = tri.tetrahedron(tet_index).edge(mixed_edge_2_num).index()
            if mixed_edge_1_index in available_edges and mixed_edge_1_index != edge_index:
                used_mixed_edge_num = mixed_edge_1_num
                used_mixed_edge_index = mixed_edge_1_index
            elif mixed_edge_2_index in available_edges and mixed_edge_2_index != edge_index:
                used_mixed_edge_num = mixed_edge_2_num
                used_mixed_edge_index = mixed_edge_2_index
            else:
                return False ### give up
            out.append( neighbouring_edges_to_triangle_loop_element(large_edge_num, used_mixed_edge_num, tri, tet_index) )
            out.append( neighbouring_edges_to_triangle_loop_element(used_mixed_edge_num, edge_num, tri, tet_index) )
            available_edges.remove(used_mixed_edge_index)
        else:
            assert False
        available_edges.remove(edge_index)
    return out

def tri_loop_is_boundary_parallel(tri_loop, tri):
    for i in range(len(tri_loop)):
        face_data = tri_loop[i]
        face_index = face_data[0]
        vert_nums = face_data[1]
        face = tri.triangles()[face_index]

        face_data_next = tri_loop[(i+1)%len(tri_loop)]
        face_index_next = face_data_next[0]
        vert_nums_next = face_data_next[1]
        face_next = tri.triangles()[face_index_next]

        edge = face.edge(vert_nums[0]) ## opposite trailing vertex
        assert edge == face_next.edge(vert_nums_next[2]) ## opposite leading vertex

        edgemapping = face.faceMapping(1,vert_nums[0])
        next_edgemapping = face_next.faceMapping(1,vert_nums_next[2])
        face_gluing_regina_numbering = next_edgemapping * (edgemapping.inverse()) ### maps vertices 0,1,2 on face to corresponding vertices on face_next
        # print('face_gluing_regina_numbering', face_gluing_regina_numbering)
        if face_gluing_regina_numbering[vert_nums[1]] != vert_nums_next[1]: ## pivot is not sent to pivot
            return False
    return True

def test():
    # sig = 'cPcbbbiht_12'
    # sig = 'dLQacccjsnk_200'
    # sig = 'dLQbccchhsj_122'
    # sig = 'eLAkaccddjsnak_2001'
    # sig = 'eLAkbccddhhsqs_1220'
    # sig = 'eLMkbcddddedde_2100'
    # sig = 'eLMkbcdddhhhdu_1221'
    # sig = 'eLMkbcdddhhhml_1221'
    # sig = 'eLMkbcdddhhqqa_1220'
# eLMkbcdddhhqxh_1220
# eLMkbcdddhxqdu_1200
# eLMkbcdddhxqlm_1200
# eLPkaccddjnkaj_2002
# eLPkbcdddhrrcv_1200
    # sig = 'gLLAQbecdfffhhnkqnc_120012'

    sigs = parse_data_file('Data/veering_census.txt')

    for j, sig in enumerate(sigs[:5]):
        if j%100 == 0:
            print(j)
        tri, angle = isosig_to_tri_angle(sig)
        # tri.save(sig + '.rga')
        branch = upper_branched_surface(tri, angle) ### also checks for veering and transverse taut
        found_loops = find_flow_cycles(tri, branch)
        # print(len(found_loops))
        # for loop in found_loops:
        #   print(loop)

        # print('found_loops', found_loops)
        # print(sig)
        for loop in found_loops:
            tri, angle = isosig_to_tri_angle(sig)
            branch = upper_branched_surface(tri, angle) 
            tri_loop = flow_cycle_to_triangle_loop(tri, branch, loop)
            if tri_loop != False:
                if not tri_loop_is_boundary_parallel(tri_loop, tri):
                    print('sig', isosig_from_tri_angle_branch(tri, angle, branch), 'loop', loop, 'tri_loop', tri_loop) 
                    drill(tri, tri_loop, angle = angle, branch = branch, sig = sig)
                    print('new angle, branch', angle, branch)
                    print(isosig_from_tri_angle_branch(tri, angle, branch))
                    # tri.save('drilled_' + sig + '_' + str(tri_loop) + '.rga')
                    # print(tri.countTetrahedra())

        
def test_drillings(sig):
    tri, angle = isosig_to_tri_angle(sig)
    branch = upper_branched_surface(tri, angle)
    loops = find_flow_cycles(tri, branch)
    tri_loops = [flow_cycle_to_triangle_loop(tri, branch, loop) for loop in loops]
    
    for tri_loop in tri_loops:
        if tri_loop != False: # False means that tri_loop goes more than once  along the same triangle - not currently implemented
            tri, angle = isosig_to_tri_angle(sig)
            if tri_loop_is_boundary_parallel(tri_loop, tri) == False: # if a loop is boundary parallel then we don't drill
                tri, angle = isosig_to_tri_angle(sig)
                branch = upper_branched_surface(tri, angle)
                print ("drilling", sig, "along", tri_loop)
                drill(tri, tri_loop, angle, branch)
                print("drilled:", tri.isoSig(), angle, branch)
                print("is layered:", is_layered(tri, angle))


