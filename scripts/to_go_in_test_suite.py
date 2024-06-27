

import random
import time

from veering.file_io import veering_census, parse_data_file, read_from_pickle

def test_drilling_methods_agree(veering_isosigs, num_to_check, smaller_num_to_check):
    import regina
    import snappy
    from veering.taut import isosig_to_tri_angle
    from veering.branched_surface import upper_branched_surface    
    from veering.flow_cycles import find_flow_cycles, flow_cycle_to_triangle_loop, tri_loop_is_boundary_parallel, generate_flow_cycles 
    from veering.drill import drill
    from drill_flow_cycle import drill_flow_cycle
    from pachner_graph_path import search_Pachner_graph_for_shortest_path

    for sig in random.sample(veering_isosigs[100:200], num_to_check): ### only try small examples
        tri, angle = isosig_to_tri_angle(sig)
        branch = upper_branched_surface(tri, angle)   
        simple_flow_cycles = find_flow_cycles(tri, branch) ### only simple flow cycles
        for fc in simple_flow_cycles:
        	tri_copy = regina.Triangulation3(tri)
        	angle_copy = angle[:]
        	branch_copy = branch[:]
        	tl = flow_cycle_to_triangle_loop(tri_copy, branch, fc)
        	if tl != False:
        		if not tri_loop_is_boundary_parallel(tl, tri):
        			drill(tri_copy, tl, angle = angle_copy, branch = branch_copy)
        			tri_2, _, found_parallel = drill_flow_cycle(sig, fc, return_tri_angle = True, return_found_parallel = True) 
        			if not found_parallel:
        				print('testing drilling', sig, 'flow cycle', fc)
	        			### now tri_copy and tri_2 should be the same manifold.
	        			sig1 = tri_copy.isoSig()
	        			sig2 = tri_2.isoSig()
	        			# print('drilling flow cycle', fc, 'sigs', sig1, sig2)
	        			
	        			#### search pachner graph is too slow
	        			# ceiling = 2 + max(tri_copy.countTetrahedra(), tri_2.countTetrahedra())
	        			# search_Pachner_graph_for_shortest_path(sig1, sig2, search_depth = 10, ceiling = ceiling)

	        			### but snappy often fails to find an isometry, or canonize the tri_copy triangulation.
	        			M1 = snappy.Manifold(tri_copy.snapPea())
	        			M2 = snappy.Manifold(tri_2.snapPea())
	        			M1.simplify()
	        			M2.simplify()
	        			# print('volumes', M1.volume(), M2.volume())
	        			assert M1.is_isometric_to(M2), sig

#### is_isometric failed on hvLPQkcdegffggfsshqqho_0221100 flow cycle [(1, 2)]??

# def run_tests(num_to_check = 20, smaller_num_to_check = 10):
#     veering_isosigs = veering_census()

def run_tests(num_to_check = 20, smaller_num_to_check = 10):
	veering_isosigs = veering_census()
	test_drilling_methods_agree(veering_isosigs, 5, 3)   
