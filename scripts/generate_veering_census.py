import regina
import os # for sound
import time
import datetime
import cProfile
from multiprocessing import Pool, Lock

from veering.file_io import parse_data_file, write_data_file
from veering.taut import isosig_from_tri_angle
from veering.veering_tri import is_veering

### Usage: make sure you have a "/data" folder inside your veering/scripts directory.
### Generate a census with the following syntax:
### main_function(n = 4, verbose = 2, parallel = 5)

### n is the number of tetrahedra to search for.
### verbose is how much you want printed to the screen and log file
### parallel is how many processes you want to use. Set it to 0 for no parallel processing.

def append_to_file(output_filename, string):
    lock.acquire()
    output_file = open(output_filename, 'a')  #append mode
    output_file.write(string+'\n')
    output_file.close()
    lock.release()

def printlock(string):
    lock.acquire()
    print(string)
    lock.release()

def read_veering_data_file(output_filename):
    data_file = open(output_filename, 'r')  #read mode
    veering_isosigs = []
    regina_isosigs = []
    output_data = []
    for line in data_file:
        data_line = line[:-1] # remove '\n'
        entries = data_line.split("_")
        regina_isosig = entries[0]
        angles_string = entries[1]
        veering_isosig = regina_isosig + "_" + angles_string
        veering_isosigs.append(veering_isosig)
        regina_isosigs.append(regina_isosig)
        output_data.append(data_line)
    data_file.close()
    return veering_isosigs, regina_isosigs, output_data

def torus_gluings_to_triangulation(R,B,faces,signs):
    r = sum(R)
    triang = regina.Triangulation3()
    for i in range(sum(R)+sum(B)):
        triang.newTetrahedron()

    # faces and signs correspond to existing_gluings and existing_signs in triangulation_project_foo

    # join command: https://regina-normal.github.io/engine-docs/classregina_1_1Face_3_013_00_013_01_4.html#ad1f0a3046ec46f13f270bbc82193c5cc
    # join(facet_of_simplex (ranges from 0 to 3), simplex_we_glue_to, gluing_permutation (permutation according to adjacent_gluing))
    # adjacent_gluing info: https://regina-normal.github.io/engine-docs/classregina_1_1detail_1_1SimplexBase.html#a94011a912f21f381fb96f8f4e5be1ddc

    ### faces refers to squares, facets refers to tetrahedra

    for i in range(len(faces)):
        #print('i is ', i, 'faces[i] is', faces[i], 'signs[i] is ', signs[i])

        ## red top
        if (i < r) and (faces[i] < r):
             # glue red to red
             if signs[i] == 1:
                 # without rotation
                 triang.tetrahedron(i).join(1,triang.tetrahedron(faces[i]),regina.Perm4(0,2,1,3))
                 offset_face = offset_in_solid_torus(faces[i], R, B, -1)
                 triang.tetrahedron(i).join(3,triang.tetrahedron( offset_face) ,regina.Perm4(3,1,2,0)) # offset ok
             else: # signs[i] == -1:
                 # rotate by pi
                 offset_face = offset_in_solid_torus(faces[i], R, B, -1)
                 triang.tetrahedron(i).join(1,triang.tetrahedron( offset_face ),regina.Perm4(2,0,3,1))
                 triang.tetrahedron(i).join(3,triang.tetrahedron(faces[i]),regina.Perm4(1,3,0,2)) # offset ok
        if (i < r) and (faces[i] >=r):
            # glue red to blue
            if signs[i] == 1:
                # rotate cw by pi/2
                triang.tetrahedron(i).join(1,triang.tetrahedron( faces[i] ),regina.Perm4(0,1,3,2))
                offset_face = offset_in_solid_torus(faces[i], R,B, -1)
                triang.tetrahedron(i).join(3,triang.tetrahedron( offset_face ),regina.Perm4(1,0,2,3)) # offset ok
            else: # signs[i] == -1:
                # rotate ccw by pi/2
                offset_face = offset_in_solid_torus(faces[i], R, B, -1)
                triang.tetrahedron(i).join(1,triang.tetrahedron( offset_face ),regina.Perm4(2,3,1,0))
                triang.tetrahedron(i).join(3,triang.tetrahedron( faces[i] ),regina.Perm4(3,2,0,1)) # offset ok

        ## blue top
        if (i >=r) and (faces[i] < r):
            # glue blue to red
            if signs[i] == 1:
                # rotate cw by pi/2   ### note, switched these two from versions previous to build_veering_census5
                triang.tetrahedron(i).join(0, triang.tetrahedron( faces[i] ), regina.Perm4(2,0,3,1))
                offset_face = offset_in_solid_torus(faces[i], R, B, -1)
                triang.tetrahedron(i).join(2, triang.tetrahedron( offset_face ), regina.Perm4(1,3,0,2)) # offset ok
            else: # signs[i] == -1:
                # rotate ccw by pi/2  ### this makes the code for edge glomming simpler 
                offset_face = offset_in_solid_torus(faces[i], R,B, -1)
                triang.tetrahedron(i).join(0, triang.tetrahedron( offset_face ), regina.Perm4(0,2,1,3))
                triang.tetrahedron(i).join(2, triang.tetrahedron(faces[i]), regina.Perm4(3,1,2,0)) # offset ok 

        if (i >= r) and (faces[i] >= r):
             # glue blue to blue
            if signs[i] == 1:
                 # without rotation
                 triang.tetrahedron(i).join(0,triang.tetrahedron( faces[i] ),regina.Perm4(1,0,2,3))
                 offset_face = offset_in_solid_torus(faces[i], R, B, -1) 
                 triang.tetrahedron(i).join(2,triang.tetrahedron( offset_face ),regina.Perm4(0,1,3,2)) # offset ok
            else: # signs[i] == -1:
                 # rotate by pi
                 offset_face = offset_in_solid_torus(faces[i], R, B, -1)
                 triang.tetrahedron(i).join(0,triang.tetrahedron( offset_face ),regina.Perm4(3,2,0,1))
                 triang.tetrahedron(i).join(2,triang.tetrahedron( faces[i] ),regina.Perm4(2,3,1,0)) # offset ok
    return triang

def list_to_isosig_verify(R,B, existing_gluings, existing_signs, found_data, instructions_packet, verbose = 0):
    # the combination R,B,existing_gluings,existing_sings determine a manifold and a triangulation
    # existing_xxx is the gluing data, which relies on the edge color data from R,B
    pair, output_filename, log_filename, verbose = instructions_packet
    tri = torus_gluings_to_triangulation(R,B,existing_gluings,existing_signs)
    angle = [1]*tri.countTetrahedra()
    veering_isosig = isosig_from_tri_angle(tri, angle)
    if veering_isosig not in found_data[0]: 
        check = is_veering(tri, angle)
        assert check
        found_data[0].append(veering_isosig)
        output_data = veering_isosig + '_' + str(R) + '_' + str(B) + '_' + str(existing_gluings) + '_' + str(existing_signs) ##+ '_' + str(os.getpid())
        found_data[2].append(output_data)
        if verbose > 2:
            printlock(output_data)
        append_to_file(output_filename, output_data)
        # os.system('afplay /System/Library/Sounds/Submarine.aiff')
    else:
        if verbose > 5:
            printlock('duplicate isosig: ' + str(veering_isosig))

class edgeClass:
    def __init__(self, model_edge): 
        self.edges = []
        self.angle = 0
        # self.free_faces = []
        if model_edge != None:
            self.edges.append(model_edge)
            self.angle = model_edge % 2

    def eat(self, other_edge_class):
        self.edges.extend(other_edge_class.edges)
        self.angle += other_edge_class.angle

    def uneat(self, other_edge_class):
        del self.edges[-len(other_edge_class.edges):] ### delete the stuff we extended by
        self.angle -= other_edge_class.angle

    def __repr__(self):
        return str(self.edges) + '_' +str(id(self))

class connectedComponent:
    def __init__(self, first_square_index, length):  # initialises as a single solid torus
        if first_square_index == None:
            self.unglued_upper_faces = set([])
            self.unglued_lower_faces = set([])
        else:
            self.unglued_upper_faces = set(range(first_square_index, first_square_index + length))
            self.unglued_lower_faces = set(self.unglued_upper_faces)

    def self_gluing(self, upper, lower):
        self.unglued_upper_faces.remove(upper)
        self.unglued_lower_faces.remove(lower)

    def other_gluing(self, this_upper, other_lower, other_connected_component):
        self.unglued_upper_faces = self.unglued_upper_faces.union(other_connected_component.unglued_upper_faces)
        self.unglued_lower_faces = self.unglued_lower_faces.union(other_connected_component.unglued_lower_faces)
        self.self_gluing(this_upper, other_lower)

    def is_closed(self):
        return (len(self.unglued_upper_faces) == 0) and (len(self.unglued_lower_faces) == 0)

    def copy(self):
        newCC = connectedComponent(None, None)
        newCC.unglued_upper_faces = self.unglued_upper_faces.copy()
        newCC.unglued_lower_faces = self.unglued_lower_faces.copy()
        return newCC

## checked it
def offset_in_solid_torus(l, R, B, offset):
    '''Given global index of a square of a solid torus, find the global index of
    the offset of this square in its solid torus'''
    colour, torus_index_in_colour, place = find_torus(l, R, B)
    r = sum(R)
    # if l < r:  # we are red
    if colour == 'R':
        new_place = (place + offset) % R[torus_index_in_colour]
        global_place = 0
        for i in range(torus_index_in_colour):
            global_place = global_place + R[i]
        global_place = global_place + new_place
        return global_place
    else: # we are blue
        new_place = (place + offset) % B[torus_index_in_colour]
        global_place = r
        for i in range(torus_index_in_colour):
            global_place = global_place + B[i]
        global_place = global_place + new_place
        return global_place

def glom(edge_classes, i, j):
    # find the edge classes that contain model edges i and j, combine them into one edge class
    ECI = edge_classes[i]
    ECJ = edge_classes[j]
    if ECI != ECJ:
        ECI.eat(ECJ)
        for k in ECJ.edges:
            edge_classes[k] = ECI
        return (ECI.angle <= 2, ECJ) # boolean: are we are still ok?, and data to let us unglom
    else: # we have closed up around the edge
        return (ECI.angle == 2, None) # boolean, are we still ok?

def glom_connected_components(connected_components, upper_face, lower_face):
    for conn_comp in connected_components:
        if upper_face in conn_comp.unglued_upper_faces:
            upper_conn_comp = conn_comp
        if lower_face in conn_comp.unglued_lower_faces:
            lower_conn_comp = conn_comp
    if upper_conn_comp == lower_conn_comp:
        upper_conn_comp.self_gluing(upper_face, lower_face)
    else:
        upper_conn_comp.other_gluing(upper_face, lower_face, lower_conn_comp)
        connected_components.remove(lower_conn_comp)
    return upper_conn_comp.is_closed() ### this component may have become closed as a result of the glomming

## is called in offset_in_solid_torus
def find_torus(l, R, B):
    # finds which torus and where in that torus the global square index l is.
    torus_index_in_colour = 0
    r = sum(R)
    if l >= r:
        colour = 'B'
        l = l - r
        # we check which blue torus the new l lies in
        for i in range(len(B)):
            if l - B[i] >= 0:
                torus_index_in_colour = torus_index_in_colour + 1
                l = l - B[i]
            else:
                break
    else:
        colour = 'R'
        for i in range(len(R)):
            if l - R[i] >= 0:
                torus_index_in_colour = torus_index_in_colour + 1
                l = l - R[i]
            else:
                break
    place = l
    # colour is either 'R' or 'B'
    # torus_index_in_colour tells us which torus in the appropriate one of R or B l lies in
    # place is the local index of the square l in the torus R[torus_index_in_colour] or B[torus_index_in_colour]
    return colour, torus_index_in_colour, place 

def first_torus_index_in_colour_list_of_my_length(my_length, colour, R, B):
    if colour == 'R':
        torus_length_list = R
    else:
        torus_length_list = B
    for i, torus_length in enumerate(torus_length_list):
        if torus_length == my_length:
            return i

def get_first_square_index(torus_colour, torus_index_in_colour, R, B):
    if torus_colour == 'R':
        return sum(R[:torus_index_in_colour])
    else: 
        assert torus_colour == 'B'
        return sum(R) + sum(B[:torus_index_in_colour])

def edge_in_torus_to_edge_class_init(j,R,B,index,torus_color):
    # 0 => red, 1 => blue
    global_place = 0
    r = sum(R)
    if torus_color == 0: # in red
        for i in range(index):   # index-1?     # range of a negative number is empty, this way we avoid two cases for index >,= 0
            global_place = global_place + R[i]*4
            #print('edge_in_torus', global_place)
        global_place = global_place + j
    if torus_color == 1: # in blue 
        global_place = sum(R)*4
        for i in range(index):
            global_place = global_place + B[i]*4
        global_place = global_place + j
    return global_place

def make_edge_classes_init(hf):
    edge_classes_init = {}
    for i in range(4 * hf):
        edge_classes_init[i] = edgeClass(i)
    return edge_classes_init

def make_connected_components_init(R, B):
    connected_components = set([])
    first_square_index = 0
    for solid_torus_length in (R + B):
        connected_components.add(connectedComponent(first_square_index, solid_torus_length))
        first_square_index += solid_torus_length
    return connected_components

def try_to_build_veering(instructions_packet):
    pair, output_filename, log_filename, verbose = instructions_packet
    # print pair, os.getpid()
    R, B = pair
    if verbose > 1:
        log_string = str(datetime.datetime.now()) + ' | Process ' + str(os.getpid()) + ' starting search: R=' + str(R) + ' B=' + str(B)
        printlock(log_string)
        append_to_file(log_filename, log_string)
    start_time = time.time()
    r = sum(R)
    hf = r + sum(B)
    edge_classes_init = make_edge_classes_init(hf)
    unglued_red_torus_indices = range(1,len(R)) ## first torus is always assumed to be glued
    unglued_blue_torus_indices = range(len(B))
    connected_components = make_connected_components_init(R, B)
    
    found_data = [[],[],[]]  
    ### In parallel processing each process gets a solid torus decomposition to check. 
    ### This is an invariant of the veering triangulation so there is no way that different processes could find the same veering triangulation.
    ### Therefore each process has its own found_data to check for duplicates.
    recursive_try_to_build_veering([], range(hf), [], R, B, hf, r, edge_classes_init, unglued_red_torus_indices, unglued_blue_torus_indices, connected_components, found_data, instructions_packet, verbose = verbose)
    if verbose > 1:
        log_string = str(datetime.datetime.now()) + ' | Process ' + str(os.getpid()) + ' finished search: R=' + str(R) + ' B=' + str(B) + ' | found veering isosigs:' + str(len(found_data[0])) + ' | done in time:' + str(time.time() - start_time) 
        printlock(log_string)
        append_to_file(log_filename, log_string)
    return found_data

def recursive_try_to_build_veering(existing_gluings, unused_gluings, existing_signs, R, B, hf, r, existing_edge_classes, unglued_red_torus_indices, unglued_blue_torus_indices, connected_components, found_data, instructions_packet, verbose = 0):
    if len(existing_gluings) >= hf:  ### we have won
        ## now we make the regina triangulation from existing_gluings, existing_signs and print out the isoSig
        list_to_isosig_verify(R, B, existing_gluings, existing_signs, found_data, instructions_packet, verbose = verbose)
        return True 
    for next_gluing in unused_gluings:
        torus_colour, torus_index_in_colour, place = find_torus(next_gluing, R, B)
        skip = False
        require_sign = False
        if len(existing_gluings) < R[0]: 
            ### we are gluing off the first torus
            if len(existing_gluings) == 0:  # this is the first gluing
                if next_gluing < R[0]: 
                    ### we may assume that the first gluing off of the first torus is not a self gluing
                    skip = True
            elif next_gluing >= R[0]: # this is another non-self gluing off the first torus
            ### we may assume that the first gluing goes to the smallest index non-self gluing
                if next_gluing < existing_gluings[0]:
                    skip = True

        if (torus_colour == 'R' and torus_index_in_colour in unglued_red_torus_indices) or (torus_colour == 'B' and torus_index_in_colour in unglued_blue_torus_indices):
            ## then we are about to glue to a brand new torus - so we should only glue to the first square with sign 1, say
            ### skip if we are not the smallest index unglued torus of this colour and length
            if torus_colour == 'R':
                torus_length_list = R
                unglued_list = unglued_red_torus_indices
            else:
                torus_length_list = B
                unglued_list = unglued_blue_torus_indices
            my_length = torus_length_list[torus_index_in_colour]
            first_ind = first_torus_index_in_colour_list_of_my_length(my_length, torus_colour, R, B)
            for i in range(first_ind, torus_index_in_colour):
                if i in unglued_list:
                    skip = True

            if not skip:  # if we aren't already skipping for the above reason, we might also skip for the following reason
                require_sign = True
                first_square_index = get_first_square_index(torus_colour, torus_index_in_colour, R, B)
                if next_gluing != first_square_index:  ### skip
                    skip = True 

        if not skip:
            for next_sign in (-1,1):   # would like to change range to [1,-1], to be consistent with handwritten notes
                # print 'next_sign', next_sign
                if not require_sign or next_sign == 1: ### dont skip  
                    connected_components_copy = set([])
                    for conn_comp in connected_components:
                        connected_components_copy.add(conn_comp.copy())
                    upper_face = len(existing_gluings)
                    lower_face = next_gluing
                    
                    made_closed_component = glom_connected_components(connected_components_copy, upper_face, lower_face)
                    if not made_closed_component or upper_face == hf-1:  # don't continue if we made a closed component, unless this is last gluing
                        unglued_red_torus_indices_copy = list(unglued_red_torus_indices)
                        unglued_blue_torus_indices_copy = list(unglued_blue_torus_indices)
                        #### update the unglued torus indices
                        if torus_colour == 'R' and torus_index_in_colour in unglued_red_torus_indices_copy:
                            unglued_red_torus_indices_copy.remove(torus_index_in_colour)
                        if torus_colour == 'B' and torus_index_in_colour in unglued_blue_torus_indices_copy:
                            unglued_blue_torus_indices_copy.remove(torus_index_in_colour)
                        ### if the torus we are gluing from (on the bottom of the gluing) has no previous gluing 
                        ### then that also takes it off the unglued list
                        lower_torus_colour, lower_torus_index_in_colour, lowerplace = find_torus(len(existing_gluings), R, B)
                        if lower_torus_colour == 'R' and lower_torus_index_in_colour in unglued_red_torus_indices_copy:
                            unglued_red_torus_indices_copy.remove(lower_torus_index_in_colour)
                        if lower_torus_colour == 'B' and lower_torus_index_in_colour in unglued_blue_torus_indices_copy:
                            unglued_blue_torus_indices_copy.remove(lower_torus_index_in_colour)
                        # continue the recursion
                        unglue_data = edge_glom_one_gluing(list(existing_gluings), list(unused_gluings), list(existing_signs), next_gluing, next_sign, R, B, hf, r, existing_edge_classes, unglued_red_torus_indices_copy, unglued_blue_torus_indices_copy, connected_components_copy, found_data, instructions_packet, verbose = verbose)
                        # now undo the changes made by edge_glom_one_gluing
                        unglue_data.reverse()
                        for x, ECJ in unglue_data:
                            existing_edge_classes[x].uneat(ECJ)
                            for k in ECJ.edges:
                                existing_edge_classes[k] = ECJ
    return False

## checked
def edge_glom_one_gluing(existing_gluings, unused_gluings, existing_signs, new_gluing, new_sign, R, B, hf, r, existing_edge_classes, unglued_red_torus_indices, unglued_blue_torus_indices, connected_components, found_data, instructions_packet, verbose = 0):
    # print 'edge_glom_one_gluing'
    colour1, a, b = find_torus(len(existing_gluings), R, B)       # returns the index j in the list R or B, and the place in the torus R[j] or B[j], depending on whether j<r
    colour2, c, d = find_torus(new_gluing, R, B)

    if (len(existing_gluings) < r):
        torus_top_color = 0 # 0 => red, 1 => blue
        torus_top_length = R[a]
    else:
        torus_top_color = 1 # 0 => red, 1 => blue
        torus_top_length = B[a]

    x1 = edge_in_torus_to_edge_class_init((4*b-1)%(4*torus_top_length),R,B,a,torus_top_color)
    x2 = edge_in_torus_to_edge_class_init((4*b+0)%(4*torus_top_length),R,B,a,torus_top_color)
    x3 = edge_in_torus_to_edge_class_init((4*b+2)%(4*torus_top_length),R,B,a,torus_top_color)
    x4 = edge_in_torus_to_edge_class_init((4*b+3)%(4*torus_top_length),R,B,a,torus_top_color)

    if new_gluing < r:
        torus_bottom_color = 0
        torus_bottom_length = R[c]
    else:
        torus_bottom_color = 1
        torus_bottom_length = B[c]

    if (torus_top_color + torus_bottom_color) % 2 == 0:  # same colour gluing
        if new_sign == 1:
            # no rotation
            y1 = edge_in_torus_to_edge_class_init((4*d-3)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y2 = edge_in_torus_to_edge_class_init((4*d+0)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y3 = edge_in_torus_to_edge_class_init((4*d-2)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y4 = edge_in_torus_to_edge_class_init((4*d+1)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
        else: #if new_sign == -1:
            # rotate by pi
            y1 = edge_in_torus_to_edge_class_init((4*d+1)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y2 = edge_in_torus_to_edge_class_init((4*d-2)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y3 = edge_in_torus_to_edge_class_init((4*d+0)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y4 = edge_in_torus_to_edge_class_init((4*d-3)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
    else: # opposite colour gluing
        if new_sign == 1:
            # rotate cw by pi/2 for gluing red to blue, ccw for gluing blue to red
            y1 = edge_in_torus_to_edge_class_init((4*d-2)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y2 = edge_in_torus_to_edge_class_init((4*d+1)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y3 = edge_in_torus_to_edge_class_init((4*d-3)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y4 = edge_in_torus_to_edge_class_init((4*d+0)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
        else: #if new_sign == -1:
            # rotate ccw by pi/2 for gluing red to blue, cw for gluing blue to red
            y1 = edge_in_torus_to_edge_class_init((4*d+0)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y2 = edge_in_torus_to_edge_class_init((4*d-3)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y3 = edge_in_torus_to_edge_class_init((4*d+1)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
            y4 = edge_in_torus_to_edge_class_init((4*d-2)%(4*torus_bottom_length),R,B,c,torus_bottom_color)
    
    xy_pairs = [(x1,y1),(x2,y2),(x3,y3),(x4,y4)]
    unglue_data = []
    for x,y in xy_pairs:
        still_ok, ECJ = glom(existing_edge_classes,x,y)
        if ECJ != None:
            unglue_data.append([x, ECJ])
        if not still_ok:
            return unglue_data  ## pairs of [model edge, the thing it ate]

    existing_gluings = existing_gluings + [new_gluing]  # we already made copies of these before entering the function, so ok to mangle them
    existing_signs = existing_signs + [new_sign]        ## do not use existing_signs.append(new_sign)
    unused_gluings.remove(new_gluing)   ## syntax?    ## remove() removes the first instance of the value new_gluing
 
    recursive_try_to_build_veering(existing_gluings, unused_gluings, existing_signs, R,B,hf,r, existing_edge_classes, unglued_red_torus_indices, unglued_blue_torus_indices, connected_components, found_data, instructions_packet, verbose = verbose)
    return unglue_data

   # adding up to an integer: en.wikipedia.org/wiki/Partition_(number_theory)
   # cite this: https://stackoverflow.com/questions/10035752/elegant-python-code-for-integer-partitioning
   # we used skovorodkin's answer
def partition(number):
    answer = set([])
    answer.add((number, ))
    for x in range(1, number):
        for y in partition(number - x):
            answer.add(tuple(sorted((x, ) + y)))
    answer = list(answer)
    answer.sort()
    return answer

def generate_census(hf, output_filename, log_filename, verbose = 0, parallel = 0, restart_at = None):
    if hf%2 == 0:
        top = int(hf/2.0)+1
    else:
        top = int((hf-1)/2.0)+1
    partition_pairs = []
    for i in range(1,top): # +int(math.ceil(float(hf)/2.0))):
        xa = partition(i)     # list of possible red solid torus lengths
        xb = partition(hf-i)  # list of possible blue solid torus lengths
        temp_pairs = []
        for j1 in xa:           # will create duplicates when len(xa) == len(xb)
            for j2 in xb:
                temp_pairs.append([j1, j2])
        if i == hf - i:  # remove duplicates
            for pair in temp_pairs:
                pair.sort()
            newlist = [pair for n,pair in enumerate(temp_pairs) if pair not in temp_pairs[:n]]
            temp_pairs = newlist
        partition_pairs.extend(temp_pairs)
    partition_pairs.reverse() # better use of parallel processing this way?: fast pairs are now at the end, so processes that finish have something else to do.
    if restart_at != None:
        partition_pairs = partition_pairs[partition_pairs.index(restart_at):]
    instructions_packets = []
    for pair in partition_pairs:
        instructions_packets.append( (pair, output_filename, log_filename, verbose) )
    l = Lock()
    if parallel == 0: 
        init(l)
        for instructions_packet in instructions_packets:
            found_data = try_to_build_veering(instructions_packet)
            # all_data[0].extend(found_data[0]) # veering_isosigs
            # all_data[1].extend(found_data[1]) # regina_isosigs
            # all_data[2].extend(found_data[2]) # full output string
    else:
        pool = Pool(initializer=init, initargs=(l,), processes = parallel)
        # we are recording to the output file as we go anyway - main_function might as well use output_file also.
        pool.map_async(try_to_build_veering, instructions_packets, chunksize = 1)
        pool.close()
        pool.join()  
        # for found_data in returned_data:
        #     all_data[0].extend(found_data[0]) # veering_isosigs
        #     all_data[1].extend(found_data[1]) # regina_isosigs
        #     all_data[2].extend(found_data[2]) # full output string
    # return all_data

def init(l):
    global lock
    lock = l

def main_function(n = 3, verbose = 0, parallel = 0, restart_at = None):
    # finds isosigs given a specific number of tetrahedra. restart_at is the last R,B partition we got to (in case we want to stop and restart for some reason)
    version_num = str(10)
    main_start_time = time.time()
    if parallel == 0:
        output_filename = 'Data/veering_isoSigs_' + str(n) + '_tetrahedra_census_' + version_num + '.txt'
    else:
        output_filename = 'Data/veering_isoSigs_' + str(n) + '_tetrahedra_census_p' + version_num + '.txt'
    # if os.path.exists(output_filename): # delete any existing file:
    #     os.remove(output_filename)
    if restart_at == None:
        output_file = open(output_filename, 'w')  #write mode, clear any existing file
        output_file.close()

    log_filename = 'Data/generate_veering_census_log_file_detailed_' + version_num + '.txt'

    generate_census(n, output_filename, log_filename, verbose = verbose, parallel = parallel, restart_at = restart_at)

    l = Lock()
    init(l)

    veering_isosigs, regina_isosigs, output_data = read_veering_data_file(output_filename)
    # regina_isosigs may have duplicates
    unduplicated_regina_isosigs = set(regina_isosigs)
    output_data.sort()  ## now its sorted, let's write to file again
    output_filename_sorted = output_filename[:-4] + '_sorted.txt'
    output_file_sorted = open(output_filename_sorted, 'w')  #write mode
    for output_datum in output_data:
        output_file_sorted.write(output_datum+'\n')
    output_file_sorted.close()

    log_string = str(datetime.datetime.now()) + ' | search n=' + str(n) + ' parallel=' + str(parallel) + ' complete in time: ' + str(time.time() - main_start_time) + ' | num veering isosigs:' + str(len(veering_isosigs)) + ' | num duplicate regina isosigs:' + str(len(veering_isosigs)-len(unduplicated_regina_isosigs))

    if verbose > 0:
        print(log_string)
        log_file = open('Data/generate_veering_census_log_file_' + version_num + '.txt', 'a')
        log_file.write(log_string + '\n')
        log_file.close()

        log_file_detail = open(log_filename, 'a')
        log_file_detail.write(log_string + '\n')
        log_file_detail.close()

    # os.system('afplay /System/Library/Sounds/Ping.aiff')

if __name__ == '__main__':
    # for n in range(12,15):
    #     main(n = n)
    # cProfile.run('main(n=6)')

    # for n in range(14,17):
    #     # main_function(n=n, verbose = 0.5)
    #     main_function(n=n, verbose = 2, parallel = 5)

    # cProfile.run('main_function(n=9, verbose = 0.5)')

    # main_function(n=15, verbose = 2, parallel = 5, restart_at = [(4,), (1, 1, 1, 1, 1, 1, 1, 1, 3)])

    main_function(n=4, verbose = 2, parallel = 5)
