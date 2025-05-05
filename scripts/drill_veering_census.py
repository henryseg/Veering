from veering.file_io import parse_data_file, write_data_file
from veering.taut import isosig_to_tri_angle
from veering.veering_tri import veering_triangulation
from veering.flow_cycles import generate_flow_cycles, is_twisted 
from veering.edge_orientability import is_edge_orientable

from drilling_flow_cycle import drill_flow_cycle
from boundary_triangulation import generate_boundary_triangulation

from multiprocessing import Process, Lock, Queue

def append_to_file(lock, output_filename, string):
    lock.acquire()
    output_file = open(output_filename, 'a')  #append mode
    output_file.write(string+'\n')
    output_file.close()
    lock.release()

def printlock(lock, string):
    lock.acquire()
    print(string)
    lock.release()

# def init(l):
#     global lock
#     lock = l

def make_many_drillings(veering_isosig, max_length, min_length):
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle)
    cycles = generate_flow_cycles(veering_isosig, max_length = max_length, min_length = min_length)
    drillings = {}
    for cycle in cycles:
        ### We could check that edge orientablility is preserved
        # if edge_orientable:
        #     assert not is_twisted(vt, cycle) ### edge orientable should be preserved
        ### We could check that the ladder counts for the different cusps are consistent
        drilled = drill_flow_cycle(veering_isosig, cycle, return_found_parallel = True)
        if drilled != None:
            drilled_sig, found_parallel = drilled
            if not drilled_sig in drillings.keys():
                drillings[drilled_sig] = cycle
            else:
                if cycle < drillings[drilled_sig] and len(cycle) <= len(drillings[drilled_sig]):
                    ### update cycle for this drilling if it is lex earlier and not longer
                    drillings[drilled_sig] = cycle
    drillings_list = [(sig, drillings[sig]) for sig in drillings.keys()]
    drillings_list.sort()
    return drillings_list

def make_many_drillings_with_lock(lock, output_filename, instructions_packet):
    veering_isosig, max_length, min_length = instructions_packet
    drillings_list = make_many_drillings(veering_isosig, max_length, min_length)
    result_string = veering_isosig + '_' + str(drillings_list)
    # printlock(lock, result_string)
    append_to_file(lock, output_filename, result_string)
    return (veering_isosig, drillings_list)

# def worker_func(lock, request_queue):
#     while True:
#         instructions_packet = request_queue.get()
#         if instructions_packet != None:
#             make_many_drillings_with_lock(lock, ..., instructions_packet)
#             # result = do_the_big_math(file_name_to_process)
#             # result_queue.put(result)
#         else: # None means that there are no more files to process
#             # result_queue.put(None)
#             break

def get_census_subset(desired_ladder_count = 4, edge_orientable = True, census_cap = None):
    census_data = parse_data_file('veering_census_with_data.txt') 
    if census_cap != None:
        census = census_data[:census_cap]
    else:
        census = census_data
    subset = []
    for i, line in enumerate(census):
        # if i % 50 == 0:
        #     print(float(i)/float(len(census)))
        line = line.split(' ')
        veering_isosig = line[0]
        if line[5] == 'E':
            is_edge_ori = True
        else:
            assert line[5] == 'N'
            is_edge_ori = False
        ladder_counts = line[7]
        ladder_counts = ladder_counts[1:-1]
        ladder_counts = ladder_counts.split(',')
        ladder_counts = [int(c) for c in ladder_counts]
        # tri, angle = isosig_to_tri_angle(veering_isosig)
        # bdry_triang = generate_boundary_triangulation(veering_isosig, draw = False)
        # ladder_counts = bdry_triang.ladder_counts()
        # if is_edge_orientable(tri, angle) == edge_orientable:
        if is_edge_ori == edge_orientable:
            if all(((count == desired_ladder_count) or (count == 4)) for count in ladder_counts):
                subset.append(veering_isosig)
    return subset


def generate_drillings_census(desired_ladder_count = 4, edge_orientable = True, max_length = 2, min_length = 1, census_cap = 100, start_after = None, output_filename_base = 'data/drillings_census', parallel = 0):
    output_filename = output_filename_base + 'drillings_census_ladders_' + str(desired_ladder_count) + '_edge_orientable_' + str(edge_orientable) + '_max_len_' + str(max_length) + '_min_len_' + str(min_length) + '.txt'
    # census_data = parse_data_file('veering_census.txt') 
    census_data = parse_data_file('veering_census_with_data.txt') 
    # subset = get_census_subset(desired_ladder_count = desired_ladder_count, edge_orientable = edge_orientable, census_cap = census_cap)

    ### hack: get the remaining ones done
    subset1 = get_census_subset(desired_ladder_count = 4, edge_orientable = True, census_cap = None)
    subset2 = get_census_subset(desired_ladder_count = 2, edge_orientable = False, census_cap = None)
    print(len(subset1), len(subset2))
    census = parse_data_file('veering_census.txt') 
    subset = set(census)
    subset.difference_update(set(subset1))
    subset.difference_update(set(subset2))
    subset = list(subset)
    subset.sort()

    instructions_packets = []
    for veering_isosig in subset:
        instructions_packet = (veering_isosig, max_length, min_length)
        instructions_packets.append(instructions_packet)
    print('full instructions_packets length', len(instructions_packets))
    
    if start_after != None:
        ind = instructions_packets.index((start_after, max_length, min_length))
        instructions_packets = instructions_packets[ind+1:]
        print('start after instructions_packets length', len(instructions_packets))



    lock = Lock()
    if parallel == 0:
        for instructions_packet in instructions_packets:            
            sig, drillings_list = make_many_drillings_with_lock(lock, output_filename, instructions_packet)
            print(sig, drillings_list)
    else:
        pass ### multiprocessing seems to be broken... ModuleNotFoundError.
        # request_queue = Queue()
        # result_queue = Queue()

        # for instructions_packet in instructions_packets:
        #     request_queue.put(instructions_packet)
        # for i in range(parallel):
        #     request_queue.put(None) ## tell them to stop

        # for i in range(parallel):
        #     handle = Process(target = worker_func, args = (lock, request_queue))
        #     handle.start()

        # pool = Pool(initializer=init, initargs=(l,), processes = parallel)
        # pool.map_async(make_many_drillings, instructions_packets, chunksize = 1)
        # pool.close()
        # pool.join()  

def tri_size_from_sig(sig):
    return len(sig.split('_')[1])

def min_size_tri(collection):
    return min(collection, key = tri_size_from_sig)

def recursive_drilling_search(starting_sig_list, target_sig_class, drilling_pattern = (2,2)):
    starting_sig_list = list(starting_sig_list)
    frontier_sigs = starting_sig_list[:]
    visited_sigs_and_paths = {}
    for sig in starting_sig_list:
        # print(sig)
        # assert type(sig) == str, type(sig)
        visited_sigs_and_paths[sig] = [sig] ### path starts with initial sig. Will add list of drilling flow cycles
    for i, length in enumerate(drilling_pattern):
        print('pattern step', i, 'of', len(drilling_pattern), 'length', length, 'size of frontier', len(frontier_sigs))
        new_frontier = []
        for sig in frontier_sigs:
            print(sig, visited_sigs_and_paths[sig])
            drillings_list = make_many_drillings(sig, length, 1)
            for new_sig, cycle in drillings_list:
                if new_sig not in visited_sigs_and_paths.keys():
                    path = visited_sigs_and_paths[sig] + [cycle]
                    visited_sigs_and_paths[new_sig] = path
                    if new_sig in target_sig_class:
                        print('found connection', (new_sig, path))
                        return (new_sig, path)
                    new_frontier.append(new_sig)
        frontier_sigs = new_frontier
        # print(visited_sigs_and_paths)
    print('no connection found')




# ModuleNotFoundError: No module named 'drill_veering_census'
# Traceback (most recent call last):
#   File "<string>", line 1, in <module>
#   File "/private/var/tmp/sage-10.4-current/local/var/lib/sage/venv-python3.12.4/lib/python3.12/multiprocessing/spawn.py", line 122, in spawn_main
#     exitcode = _main(fd, parent_sentinel)
#                ^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/private/var/tmp/sage-10.4-current/local/var/lib/sage/venv-python3.12.4/lib/python3.12/multiprocessing/spawn.py", line 132, in _main
#     self = reduction.pickle.load(from_parent)
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

                    