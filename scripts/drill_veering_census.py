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

def make_many_drillings(lock, output_filename, instructions_packet):
    veering_isosig, max_length = instructions_packet
    tri, angle = isosig_to_tri_angle(veering_isosig)
    vt = veering_triangulation(tri, angle)
    cycles = generate_flow_cycles(veering_isosig, max_length = max_length)
    drillings = {}
    for cycle in cycles:
        assert not is_twisted(vt, cycle)
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
    result_string = veering_isosig + '_' + str(drillings_list)
    printlock(lock, result_string)
    append_to_file(lock, output_filename, result_string)
    return (veering_isosig, drillings_list)

def worker_func(lock, request_queue):
    while True:
        instructions_packet = request_queue.get()
        if instructions_packet != None:
            make_many_drillings(lock, instructions_packet)
            # result = do_the_big_math(file_name_to_process)
            # result_queue.put(result)
        else: # None means that there are no more files to process
            # result_queue.put(None)
            break

def generate_drillings_census(max_length = 2, parallel = 0):
    output_filename = 'data/drillings_census.txt'
    census_data = parse_data_file('veering_census.txt') 
    instructions_packets = []
    census = census_data[:50]
    output_file = open(output_filename, 'w') 
    for i, veering_isosig in enumerate(census):
        if i % 50 == 0:
            print(float(i)/float(len(census)))
        tri, angle = isosig_to_tri_angle(veering_isosig)
        bdry_triang = generate_boundary_triangulation(veering_isosig, draw = False)
        ladder_counts = bdry_triang.ladder_counts()
        if is_edge_orientable(tri, angle):
            if all(count == 4 for count in ladder_counts):
                instructions_packet = (veering_isosig, max_length)
                instructions_packets.append(instructions_packet)
    print('instructions_packets done')
    
    lock = Lock()
    if parallel == 0:
        for instructions_packet in instructions_packets:            
            sig, drillings_list = make_many_drillings(lock, output_filename, instructions_packet)
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



# ModuleNotFoundError: No module named 'drill_veering_census'
# Traceback (most recent call last):
#   File "<string>", line 1, in <module>
#   File "/private/var/tmp/sage-10.4-current/local/var/lib/sage/venv-python3.12.4/lib/python3.12/multiprocessing/spawn.py", line 122, in spawn_main
#     exitcode = _main(fd, parent_sentinel)
#                ^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/private/var/tmp/sage-10.4-current/local/var/lib/sage/venv-python3.12.4/lib/python3.12/multiprocessing/spawn.py", line 132, in _main
#     self = reduction.pickle.load(from_parent)
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

                    