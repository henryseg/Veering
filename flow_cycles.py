#
# flow_cycles.py
#

# Given a veering triangulation, find all simple flow cycles

from veering import is_veering
from branched_surface import upper_branched_surface, branch_num_to_large_edge_and_mixed_edge_pair_num
from transverse_taut import get_tet_top_and_bottom_edges, get_tet_above_edge, is_transverse_taut
from taut import isosig_to_tri_angle


### format for loops: it is a list of tuples, each tuple is (tet_index, edge_index within this tet that we exit through)

def lex_smallest_loop(loop):
	ind = loop.index(min(loop))
	loop = loop[ind:] + loop[:ind]
	return loop

def add_flow_edge(tri, branch, tet_with_this_edge_large, loop_so_far, current_tet_index, start_edge_index, available_edge_indices, found_loops):
	large_edge_num, mixed_edge_pair_num, small_edge_pair_num = branch_num_to_large_edge_and_mixed_edge_pair_num(branch[current_tet_index], return_small = True)
	tet = tri.tetrahedron(current_tet_index)
	for i in range(3): ### four choices of how to add a flow edge
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

### previous idea: have vertical flow edges detour through a mixed edge. Let's do this later when we drill
		# else:
		# 	top_edge = tet.edge(5 - large_edge_num)
		# 	if top_edge.index() in available_edge_indices:
		# 		if i == 2: ### top edge through 0th mixed edge
		# 			new_edge = tet.edge(mixed_edge_pair_num)
		# 			if new_edge.index() in available_edge_indices and new_edge.index() != start_edge_index:
		# 				loop_so_far2 = loop_so_far[:]
		# 				loop_so_far2.append((current_tet_index, mixed_edge_pair_num, 5 - large_edge_num))
		# 				available_edge_indices2 = available_edge_indices[:]
		# 				# available_edge_indices2.remove(new_edge.index())
		# 				available_edge_indices2.remove(top_edge.index())
		# 				if top_edge.index() == start_edge_index:
		# 					found_loops.append(loop_so_far2)
		# 				else:
		# 					add_flow_edge(tri, branch, loop_so_far2, get_tet_above_edge(tri, angle, top_edge).index(), start_edge_index, available_edge_indices2, found_loops)

		# 		elif i == 3: ### top edge through 1th mixed edge
		# 			new_edge = tet.edge(5 - mixed_edge_pair_num)
		# 			if new_edge.index() in available_edge_indices and new_edge.index() != start_edge_index:
		# 				loop_so_far2 = loop_so_far[:]
		# 				loop_so_far2.append((current_tet_index, 5 - mixed_edge_pair_num, 5 - large_edge_num))
		# 				available_edge_indices2 = available_edge_indices[:]
		# 				# available_edge_indices2.remove(new_edge.index())
		# 				available_edge_indices2.remove(top_edge.index())
		# 				if top_edge.index() == start_edge_index:
		# 					found_loops.append(loop_so_far2)
		# 				else:
		# 					add_flow_edge(tri, branch, loop_so_far2, get_tet_above_edge(tri, angle, top_edge).index(), start_edge_index, available_edge_indices2, found_loops)

def make_list_of_tet_with_this_edge_large(tri, branch):
	out = [None] * tri.countEdges()
	for i in range(tri.countTetrahedra()):
		large_edge, _ = branch_num_to_large_edge_and_mixed_edge_pair_num(branch[i])
		large_edge_index = tri.tetrahedron(i).edge(large_edge).index()
		out[large_edge_index] = i
	return out

def flow_cycles(tri, branch):
	"""Given a branched surface dual to a triangulation, find all simple flow cycles in the upper flow graph"""
	tet_with_this_edge_large = make_list_of_tet_with_this_edge_large(tri, branch)
	found_loops = []
	for start_tet_index in range(tri.countTetrahedra()):
		large_edge, _ = branch_num_to_large_edge_and_mixed_edge_pair_num(branch[start_tet_index])
		large_edge_index = tri.tetrahedron(start_tet_index).edge(large_edge).index()
		add_flow_edge(tri, branch, tet_with_this_edge_large, [], start_tet_index, large_edge_index, list(range(tri.countEdges())), found_loops)
	return found_loops

def test():
	# tri, angle = isosig_to_tri_angle("cPcbbbiht_12")
	# tri, angle = isosig_to_tri_angle('gLLAQbecdfffhhnkqnc_120012')
	tri, angle = isosig_to_tri_angle('dLQacccjsnk_200')
	branch = upper_branched_surface(tri, angle) ### also checks for veering and transverse taut
	found_loops = flow_cycles(tri, branch)
	print(len(found_loops))
	for loop in found_loops:
		print(loop)
		

