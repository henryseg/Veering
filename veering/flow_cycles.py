#
# flow_cycles.py
#

# Given a veering triangulation, find all simple flow cycles

import regina
from functools import reduce

from .branched_surface import upper_branched_surface, branch_num_to_large_edge_and_mixed_edge_pair_num, isosig_from_tri_angle_branch, has_non_sing_semiflow
from .transverse_taut import get_tet_top_and_bottom_edges, get_tet_above_edge, is_transverse_taut, edge_side_face_collections
from .taut import isosig_to_tri_angle, edge_num_to_vert_pair
from .veering_tri import is_veering, is_AB_turn
from .drill import drill

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

def make_list_of_tet_with_this_edge_large(tri, branch, return_reverse_map = False):
    out = [None] * tri.countEdges()
    out2 = []
    for i in range(tri.countTetrahedra()):
        large_edge, _ = branch_num_to_large_edge_and_mixed_edge_pair_num(branch[i])
        large_edge_index = tri.tetrahedron(i).edge(large_edge).index()
        out2.append(large_edge_index)
        out[large_edge_index] = i
    if return_reverse_map:
        return (out, out2)
    else:
        return out

def generate_simple_flow_cycles(tri, branch):
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

def find_tri_loops(sig):
    tri, angle = isosig_to_tri_angle(sig)
    branch = upper_branched_surface(tri, angle)
    loops = generate_simple_flow_cycles(tri, branch)
    return [flow_cycle_to_triangle_loop(tri, branch, loop) for loop in loops]

def is_twisted(vt, flow_cycle):
    """Is this flow cycle in the veering triangulation twisted?"""

    triangles = edge_side_face_collections(vt.tri, vt.angle, tet_vert_coorientations = vt.coorientations)

    ### Go through the list of edges in the flow_cycle and write down entering/exiting faces around that edge.
    ### These are the top face of the entering tet and the bottom face of the exiting tet

    ### Then calculate AB_turning for each tet

    entering_face_nums = []
    exiting_face_nums = []
    for tet_num, edge_index in flow_cycle:

        ### find a face on top of the tet that is next to the edge we are exiting through
        neighbouring_faces = [i for i in [0,1,2,3] if i not in edge_num_to_vert_pair[edge_index]]
        tet = vt.tri.tetrahedron(tet_num)
        e = tet.edge(edge_index)
        left_triangles_for_e, right_triangles_for_e = triangles[e.index()]

        for k, i in enumerate(neighbouring_faces):
            if vt.coorientations[tet_num][i] == +1: ### must be at least one, just pick the first one
                j = neighbouring_faces[1 - k] ### the other neighbouring face num
                f = tet.triangle(i)
                f_index = f.index()
                entering_face_nums.append(f_index)
                ### now find the exiting face num
                fmapping = tet.faceMapping(2, i)
                index_of_opposite_vert_in_f = fmapping.inverse()[j]
                triangle_in_list = (f_index, index_of_opposite_vert_in_f)

                if triangle_in_list in left_triangles_for_e:
                    correct_triangle_list = left_triangles_for_e
                else:
                    assert triangle_in_list in right_triangles_for_e
                    correct_triangle_list = right_triangles_for_e
                exiting_face_nums.append(correct_triangle_list[-1][0])      
                break

    # print('entering_face_nums', entering_face_nums)
    # print('exiting_face_nums', exiting_face_nums)
    num_AB_turns = len(flow_cycle)
    for i in range(len(entering_face_nums)):
        if is_AB_turn(vt, exiting_face_nums[i], entering_face_nums[(i+1)%len(entering_face_nums)], +1, +1): ## always go up so +1, +1
            num_AB_turns += 1
    return num_AB_turns % 2 == 1

def ternary(n, num_digits):  ### modified from https://stackoverflow.com/questions/34559663/convert-decimal-to-ternarybase3-in-python
    nums = []
    for i in range(num_digits):
        n, r = divmod(n, 3)
        nums.append(r)
    return nums

def rotate_to_lex_smallest(inputs):
    rotations = []
    for i in range(len(inputs)):
        rotations.append(inputs[i:] + inputs[:i])
    return min(rotations)

def factors(n): ### from https://stackoverflow.com/questions/6800193/what-is-the-most-efficient-way-of-finding-all-the-factors-of-a-number-in-python
    facs = list(set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0))))
    facs.sort()
    return facs

def is_power(flow_cycle):
    n = len(flow_cycle)
    for f in factors(n)[:-1]: ### don't check full length for repeats
        g, r = divmod(n, f)
        assert r == 0 
        sub_cycles = [flow_cycle[i*f:(i+1)*f] for i in range(g)]
        if sub_cycles.count(sub_cycles[0]) == g: ### all the same as the first
            return True
    return False

### new version generates non-simple cycles as well. Is not efficient, but we will not be able to drill long cycles quickly anyway
def generate_flow_cycles(sig, max_length = 5, include_num_steps_up = False, monochromatic_only = False, max_length_only = False):
    tri, angle = isosig_to_tri_angle(sig)
    edge_colours = is_veering(tri, angle, return_type = "veering_colours")
    branch = upper_branched_surface(tri, angle)
    tet_with_this_edge_large, edge_at_bottom_of_tet = make_list_of_tet_with_this_edge_large(tri, branch, return_reverse_map = True)
    # print('tet_with_this_edge_large', tet_with_this_edge_large)
    exit_edges_for_each_tet = []
    for i in range(tri.countTetrahedra()):
        large_edge_num, mixed_edge_pair_num, small_edge_pair_num = branch_num_to_large_edge_and_mixed_edge_pair_num(branch[i], return_small = True)
        exit_edges_for_each_tet.append([small_edge_pair_num, 5 - small_edge_pair_num, 5 - large_edge_num])

    # print('exit_edges_for_each_tet', exit_edges_for_each_tet)

    flow_cycles = set([])
    for i in range(tri.countTetrahedra()):
        original_tet_index = i
        # print('original_tet_index', original_tet_index)
        for n in range(3**max_length):
            current_tet_index = original_tet_index
            first_edge_colour = edge_colours[edge_at_bottom_of_tet[i]]
            digits = ternary(n, max_length)
            # print(digits)
            path = []
            num_steps_straight_up = 0
            for j in digits:    
                # print('path', path, 'j', j, 'current tet index', current_tet_index)
                exit_edge_num = exit_edges_for_each_tet[current_tet_index][j]
                path.append((current_tet_index, exit_edge_num))
                if j == 2:
                    num_steps_straight_up += 1

                current_tet = tri.tetrahedron(current_tet_index)
                next_edge_index = current_tet.edge(exit_edge_num).index()
                if monochromatic_only:
                    if edge_colours[next_edge_index] != first_edge_colour:
                        break
                # print('next edge index', next_edge_index)
                next_tet_index = tet_with_this_edge_large[next_edge_index]
                # print('next_tet_index', next_tet_index)
                if next_tet_index == original_tet_index:
                    if not is_power(path):
                        # print('add path', path)
                        if include_num_steps_up:
                            flow_cycles.add( (tuple(rotate_to_lex_smallest(path)), num_steps_straight_up))
                        else:
                            flow_cycles.add(tuple(rotate_to_lex_smallest(path)))
                current_tet_index = next_tet_index
    flow_cycles = list(flow_cycles)
    flow_cycles.sort()
    if include_num_steps_up:
        flow_cycles.sort(key = lambda x: len(x[0]))
        if max_length_only:
            flow_cycles = [c for c in flow_cycles if len(c[0]) == max_length]
    else:
        flow_cycles.sort(key = lambda x: len(x))
        if max_length_only:
            flow_cycles = [c for c in flow_cycles if len(c) == max_length]
    return flow_cycles








