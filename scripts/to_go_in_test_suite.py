

import random
import time

from veering.file_io import veering_census, parse_data_file, read_from_pickle

def test_drilling_methods_agree(veering_isosigs, num_to_check, smaller_num_to_check):
    import regina
    import snappy
    from veering.taut import isosig_to_tri_angle
    from veering.branched_surface import upper_branched_surface    
    from veering.flow_cycles import generate_simple_flow_cycles, flow_cycle_to_triangle_loop, tri_loop_is_boundary_parallel, generate_flow_cycles, is_twisted 
    from veering.drill import drill
    from veering.veering_tri import veering_triangulation
    from drilling_flow_cycle import drill_flow_cycle
    from pachner_graph_path import search_Pachner_graph_for_shortest_path

    for sig in random.sample(veering_isosigs, num_to_check): ### only try small examples
        tri, angle = isosig_to_tri_angle(sig)
        vt = veering_triangulation(tri, angle)
        branch = upper_branched_surface(tri, angle)   
        simple_flow_cycles = generate_simple_flow_cycles(tri, branch) ### only simple flow cycles
        for fc in random.sample(simple_flow_cycles, min(num_to_check, len(simple_flow_cycles))):
            tri_copy = regina.Triangulation3(tri)
            angle_copy = angle[:]
            branch_copy = branch[:]
            fc_is_twisted = is_twisted(vt, fc)
            tl = flow_cycle_to_triangle_loop(tri_copy, branch, fc)
            if tl != False:
                if not tri_loop_is_boundary_parallel(tl, tri):
                    drill(tri_copy, tl, angle = angle_copy, branch = branch_copy)
                    veering_sig_2, tri_2, angle_2, found_parallel = drill_flow_cycle(sig, fc, return_tri_angle = True, return_found_parallel = True)    
                    if not found_parallel:
                        print('testing drilling', sig, 'flow cycle', fc, 'is_twisted', fc_is_twisted)
                        ### now tri_copy and tri_2 should be the same manifold.
                        sig1 = tri_copy.isoSig()
                        sig2 = tri_2.isoSig()
                        
                        if not fc_is_twisted:
                            veering_sig_3, tri_3, angle_3, found_parallel = drill_flow_cycle(sig, fc, return_tri_angle = True, return_found_parallel = True, use_untwisted_speed_up = False)    
                            assert veering_sig_2 == veering_sig_3, 'untwisted speedup gave a different answer ' + sig + ' ' + fc

                        # print('drilling flow cycle', fc, 'sigs', sig1, sig2)
                        
                        #### search pachner graph is too slow
                        # ceiling = 2 + max(tri_copy.countTetrahedra(), tri_2.countTetrahedra())
                        # search_Pachner_graph_for_shortest_path(sig1, sig2, search_depth = 10, ceiling = ceiling)

                        M1 = snappy.Manifold(tri_copy.snapPea())
                        M2 = snappy.Manifold(tri_2.snapPea())
                        M1.simplify()
                        M2.simplify()
                        # print('volumes', M1.volume(), M2.volume())
                        found_isometry = False
                        for i in range(10):
                            if M1.is_isometric_to(M2):  ### from the docstring for is_isometric_to:
                            ### The answer True is rigorous, but the answer False may
                            ### not be as there could be numerical errors resulting in finding
                            ### an incorrect canonical triangulation.
                                found_isometry = True
                                break
                        if not found_isometry:
                            print('checking with verified isometry_signature', sig, fc)
                            isomsig1 = M1.isometry_signature(verified = True)
                            isomsig2 = M2.isometry_signature(verified = True)
                            assert not isomsig1 == None, 'isom signature failed ' + sig + ' ' + fc
                            assert not isomsig2 == None, 'isom signature failed ' + sig + ' ' + fc
                            assert isomsig1 == isomsig2, sig + '_' + fc


def test_drilling_properties(veering_isosigs, num_to_check, smaller_num_to_check, flow_cycle_max_length = 5):
    import regina
    import snappy
    from veering.taut import isosig_to_tri_angle   
    from veering.flow_cycles import generate_flow_cycles, is_twisted 
    from veering.veering_tri import veering_triangulation
    from veering.taut_polytope import is_layered
    from veering.edge_orientability import is_edge_orientable
    from drilling_flow_cycle import drill_flow_cycle
    from boundary_triangulation import generate_boundary_triangulation

    for sig in random.sample(veering_isosigs, num_to_check): ### only try small examples
        tri, angle = isosig_to_tri_angle(sig)
        original_is_layered = is_layered(tri, angle)
        original_is_edge_orientable = is_edge_orientable(tri, angle)
        vt = veering_triangulation(tri, angle)
        flow_cycles = generate_flow_cycles(sig, max_length = flow_cycle_max_length)
        for fc in random.sample(flow_cycles, min(num_to_check, len(flow_cycles))):
            out = drill_flow_cycle(sig, fc, return_tri_angle = True, return_found_parallel = True, return_cusp_mapping = True)  
            if out != None:
                drilled_sig, drilled_tri, drilled_angle, found_parallel, isosig_to_original_cusp_mapping = out
                fc_is_twisted = is_twisted(vt, fc)
                print('testing drilling', sig, 'flow cycle', fc, 'is_twisted', fc_is_twisted)

                if not fc_is_twisted:
                    drilled_sig_2, _, _, _ = drill_flow_cycle(sig, fc, return_tri_angle = True, return_found_parallel = True, use_untwisted_speed_up = False)   
                    assert drilled_sig == drilled_sig_2, 'untwisted speedup gave a different answer ' + sig + ' ' + fc

                # print(sig, '->', drilled_sig, 'cusp_mapping', isosig_to_original_cusp_mapping)
                if original_is_layered:
                    # print('testing layered')
                    assert is_layered(drilled_tri, drilled_angle), 'drilling a layered tri produced a non layered tri'

                assert original_is_edge_orientable == is_edge_orientable(drilled_tri, drilled_angle), 'drilling changed edge edge orientability'

                original_b = generate_boundary_triangulation(sig, draw = False)
                original_ladder_counts = original_b.ladder_counts()
                drilled_b = generate_boundary_triangulation(drilled_sig, draw = False)
                drilled_ladder_counts = drilled_b.ladder_counts()

                # print('original_ladder_counts', original_ladder_counts)
                # print('drilled_ladder_counts', drilled_ladder_counts)

                for i in range(drilled_tri.countVertices()):
                    corresponding_original_cusp_index = isosig_to_original_cusp_mapping[i]
                    if corresponding_original_cusp_index == tri.countVertices(): ### this is the new cusp
                        if not found_parallel: 
                            if fc_is_twisted:
                                assert drilled_ladder_counts[i] == 2, 'drilled ladder count wrong, twisted'
                            else:
                                assert drilled_ladder_counts[i] == 4, 'drilled ladder count wrong, untwisted'
                        else: ### we drill along the core of a Mobius strip. Example: dLQacccjsnk_200 flow cycle ((0, 3),)
                            assert drilled_ladder_counts[i] == 2, 'drilled ladder count wrong, Mobius strip core'
                    else:  ### means that fc is the boundary of a Mobius strip, and we actually drill the core of the strip
                        # print('i', i, 'corresponding_original_cusp_index', corresponding_original_cusp_index)
                        assert drilled_ladder_counts[i] == original_ladder_counts[corresponding_original_cusp_index], 'cusp ' + str(corresponding_original_cusp_index) + ' changed number of ladders when we drilled'
                    
# def run_tests(num_to_check = 20, smaller_num_to_check = 10):
#     veering_isosigs = veering_census()

def run_tests(num_to_check = 20, smaller_num_to_check = 10):
    veering_isosigs = veering_census()
    # test_drilling_methods_agree(veering_isosigs, 5, 3)   
    test_drilling_properties(veering_isosigs, 5, 3)
