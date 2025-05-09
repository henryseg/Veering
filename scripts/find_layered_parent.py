from veering.file_io import veering_census, parse_data_file, read_from_pickle
from veering.taut import isosig_to_tri_angle  
from veering.flow_cycles import generate_flow_cycles
from veering.taut_polytope import is_layered, taut_rays
from veering.carried_surface import genus_punctures
from drilling_flow_cycle import drill_flow_cycle

def try_to_find_layered_parent(sig, flow_cycle_max_length = 5, flow_cycle_min_length = None, verbose = 0):
    tri, angle = isosig_to_tri_angle(sig)
    if is_layered(tri, angle):
        print('already layered!')
        return None

    flow_cycles = generate_flow_cycles(sig, max_length = flow_cycle_max_length, min_length = flow_cycle_min_length)
    for fc in flow_cycles:
        if verbose > 2:
            print(fc)
        out = drill_flow_cycle(sig, fc, return_tri_angle = True, return_found_parallel = True, return_cusp_mapping = True)  
        if out != None:
            drilled_sig, drilled_tri, drilled_angle, found_parallel, isosig_to_original_cusp_mapping = out
            if is_layered(drilled_tri, drilled_angle):
                if verbose > 1:
                    print(sig, fc, drilled_sig)
                return(fc, drilled_sig)
    return False ### did not find a layered parent

### difficult examples:

# find_layered_parent.try_to_find_layered_parent('kLLALLQkccedhgjijjilnhchvhhrwr_1200111222', flow_cycle_max_length = 7, flow_cycle_min_length = 7)                                                                                                                               
# (((0, 0), (3, 0), (3, 0), (3, 0), (3, 0), (3, 1), (1, 4)), '8LLPvvLzvAPPQMQLPvzQvvMwvMzQQwwPQQvQvQvQQcdcfhokrmursrpqoqwtvyzxDECBKLJOMRQQTTPPQUXYYXWXVY0Z1342675677hsraqqfaoahhhxxsoqfjaahaahoaahaaaapppiiiaaaxxxxaxaaxaaxaaxnnn_122201021122110000011122220011111112220000000011111000111000')

# find_layered_parent.try_to_find_layered_parent('kLLAwAAkccedfhghijjlnxkxtigqqj_1200112202', flow_cycle_max_length = 6, flow_cycle_min_length = 6)                                                                                                                               
# (((0, 4), (2, 1), (6, 0), (4, 5), (0, 5), (0, 5)), 'KLLwQLwALwPwLwvvQLMLQQPMQbeefgehijkmnpqotuBCyBDzEDGFDCHIEIJJIJhhaaaaaaaaagaaaaaaajgagqaxhaghajgqhbh_011100011001112001111222222001000012')

def append_to_file(output_filename, string):
    # lock.acquire()
    output_file = open(output_filename, 'a')  #append mode
    output_file.write(string+'\n')
    output_file.close()
    # lock.release()

def search_for_layered_parents(veering_isosigs, flow_cycle_max_length = 5, flow_cycle_min_length = 0, verbose = 0):
    win_filename = 'data/layered_parents_' + str(flow_cycle_max_length) + '_' + str(flow_cycle_min_length) + '.txt'
    win_file = open(win_filename, 'w')  #write mode, clear any existing file
    win_file.close()
    fail_filename = 'data/could_not_find_layered_parents_' + str(flow_cycle_max_length) + '_' + str(flow_cycle_min_length) + '.txt'
    fail_file = open(fail_filename, 'w')  #write mode, clear any existing file
    fail_file.close()

    for sig in veering_isosigs:
        tri, angle = isosig_to_tri_angle(sig)
        if not is_layered(tri, angle):
            if verbose > 1:
                print(sig)
            out = try_to_find_layered_parent(sig, flow_cycle_max_length = flow_cycle_max_length, flow_cycle_min_length = flow_cycle_min_length, verbose = verbose)
            if out == False:
                append_to_file(fail_filename, sig)
            else:
                fc, drilled_sig = out
                append_to_file(win_filename, sig + ' ' + str(fc) + ' ' + drilled_sig)

def run():
    # veering_isosigs = veering_census()
    # search_for_layered_parents(veering_isosigs, flow_cycle_max_length = 5, flow_cycle_min_length = 0)
    # veering_isosigs = parse_data_file('data/could_not_find_layered_parents_5_0.txt')
    # search_for_layered_parents(veering_isosigs, flow_cycle_max_length = 8, flow_cycle_min_length = 6, verbose = 3)
    # veering_isosigs = parse_data_file('data/could_not_find_layered_parents_8_6.txt')
    # search_for_layered_parents(veering_isosigs, flow_cycle_max_length = 9, flow_cycle_min_length = 9, verbose = 3)
    # veering_isosigs = parse_data_file('data/could_not_find_layered_parents_9_9.txt')
    # search_for_layered_parents(veering_isosigs, flow_cycle_max_length = 10, flow_cycle_min_length = 10, verbose = 3)
    # veering_isosigs = parse_data_file('data/could_not_find_layered_parents_10_10.txt')
    # search_for_layered_parents(veering_isosigs, flow_cycle_max_length = 11, flow_cycle_min_length = 11, verbose = 3)
    veering_isosigs = parse_data_file('data/could_not_find_layered_parents_11_11.txt')
    search_for_layered_parents(veering_isosigs, flow_cycle_max_length = 12, flow_cycle_min_length = 12, verbose = 3)

def find_genus_bound(sig, num_rays_to_try = 10):
    """Given a sig for a fibered veering triangulation, get upper bound on the minimal genus of the fiber"""
    tri, angle = isosig_to_tri_angle(sig)
    rays = taut_rays(tri, angle)
    rays_with_num_tris = [(ray, sum(ray)) for ray in rays]
    rays_with_num_tris.sort(key = lambda x: x[1])
    rays_to_try = rays_with_num_tris[:num_rays_to_try]
    rays_to_try = [ray for (ray, num) in rays_to_try]
    # print(rays_to_try)
    outputs = []
    for ray in rays_to_try:
        g, p = genus_punctures(tri, angle, ray)
        outputs.append((g,p))
    if len(outputs) == 0:
        return None
    else:
        outputs.sort(key = lambda x: x[0])
        return outputs[0]

def add_genus_bounds_to_layered_parents_file():
    data = parse_data_file('data/layered_parents.txt')
    out_filename = 'data/layered_parents_g_p.txt'
    out_file = open(out_filename, 'w')  #write mode, clear any existing file
    out_file.close()

    for line in data[-10:]:
        parent_sig = line.split(' ')[-1]
        # print(parent_sig)
        g_p = find_genus_bound(parent_sig)
        new_line = line + " " + str(g_p)
        print(new_line)
        append_to_file(out_filename, new_line)
    print('done')






