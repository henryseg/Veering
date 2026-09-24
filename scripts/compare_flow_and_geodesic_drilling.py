from veering.taut import isosig_from_tri_angle, isosig_to_tri_angle
from veering.flow_cycles import generate_flow_cycles, flow_cycle_to_dual_edge_loop
from veering.drill_flow_cycle import drill_flow_cycles
from veering.file_io import veering_census, parse_data_file
from snappy_drill_homotopic import drill_tet_and_face_indices, tet_and_face_indices_to_word
from snappy.drilling.exceptions import GeodesicSystemNotSimpleError
from snappy.geometric_structure.geodesic.exceptions import WordAppearsToBeParabolic
from snappy.verify.exceptions import ShapePositiveImaginaryPartNumericalVerifyError
import snappy
import sys
sys.setrecursionlimit(1000000)


def append_to_file(output_filename, string):
    output_file = open(output_filename, 'a')  #append mode
    output_file.write(string+'\n')
    output_file.close()

def orbit_drilled_different_from_geodesic_drilled(sig, fc, drilled_sig, output_filename = None, verbose = 0):
    if verbose > 0:
        print('drilled sig', drilled_sig)
    tri, angle = isosig_to_tri_angle(sig) 
    drilled_tri, _ = isosig_to_tri_angle(drilled_sig)
    drilled_M = snappy.Manifold(drilled_tri) 

    orig_M = snappy.Manifold(tri) 
    dual_loop = flow_cycle_to_dual_edge_loop(tri, angle, fc) 
    if verbose > 1: 
        print('dual_loop', dual_loop, 'word', tet_and_face_indices_to_word(orig_M, dual_loop)) 
    try:
        ### May need to sys.setrecursionlimit(1000000) to make this work
        snappy_drilled_M = drill_tet_and_face_indices(orig_M, dual_loop, verified = True, bits_prec = 2000) 
    except GeodesicSystemNotSimpleError as e:
        if verbose > 1:
            print('not simple error |', e)
        return False
    except ShapePositiveImaginaryPartNumericalVerifyError as e:
        if verbose > 0:
            print('shape error |', e)              
        return False
    except WordAppearsToBeParabolic as e:
        if verbose > 0:
            print('parabolic error |', e)    
        assert False ### we should never get here
    except Exception as e:
        print(type(e), e)  
        return False
    snappy_drilled_M.simplify()
    drilled_M.simplify()

    ### try to find an isometry
    for i in range(10):
        if drilled_M.is_isometric_to(snappy_drilled_M):  ### from the docstring for is_isometric_to:
        ### The answer True is rigorous, but the answer False may
        ### not be as there could be numerical errors resulting in finding
        ### an incorrect canonical triangulation.
            return False

    # print('checking with verified isometry_signature', sig, fc)
    isomsig1 = drilled_M.isometry_signature(verified = True)
    isomsig2 = snappy_drilled_M.isometry_signature(verified = True)
    assert not isomsig1 == None, 'isom signature failed ' + sig + ' ' + fc
    assert not isomsig2 == None, 'isom signature failed ' + sig + ' ' + fc
    # assert isomsig1 == isomsig2, sig + '_' + fc
    if isomsig1 != isomsig2:
        print('drilling', sig, 'along', fc, 'gives different results', isomsig1, drilled_M.identify(), isomsig2, snappy_drilled_M.identify())
        out_line = sig + '|' + str(fc) + '|' + str(isomsig1) + '|' + str(isomsig2)
        if output_filename != None:
            append_to_file(output_filename, out_line)
        return True
    else:
        return False

def compare_flow_and_geodesic_drilling_script_search(sig, output_filename = None, max_length = 5, min_length = 1, quit_after_finding_one = True, verbose = 0):  
    cycles = generate_flow_cycles(sig, max_length = max_length, min_length = min_length)
    for fc in cycles:  
        if verbose > 0:
            print('flow cycle', fc)
        out = drill_flow_cycles(sig, [fc], return_isosig_tri_angle = True) 
        drilled_sig, drilled_tri, drilled_angle = out 
        if drilled_sig != sig:  ### This happens if you try to drill a peripheral flow cycle 
            different = orbit_drilled_different_from_geodesic_drilled(sig, fc, drilled_sig, output_filename = output_filename, verbose = verbose)
            if different and quit_after_finding_one:
                return True
    if quit_after_finding_one:
        print('no differences found for', sig)
        return False

def census_compare_flow_and_geodesic(max_length = 5, min_length = 1, filename_suffix = "", census_start = 0, census_end = -1, verbose = 0):
    # output_filename = "../../../Dropbox/Data/compare_flow_and_geodesic" + filename_suffix + ".txt"
    # fail_filename = "../../../Dropbox/Data/compare_flow_and_geodesic_fail" + filename_suffix + ".txt"
    output_filename = "data/compare_flow_and_geodesic" + filename_suffix + ".txt"
    fail_filename = "data/compare_flow_and_geodesic_fail" + filename_suffix + ".txt"

    output_file = open(output_filename, 'w')  #write mode, clear any existing file
    output_file.close()
    fail_file = open(fail_filename, 'w')  #write mode, clear any existing file
    fail_file.close()

    if census_end != -1:
        census = veering_census()[census_start:census_end]
    else:
        census = veering_census()[census_start:]
    print(len(census))
    win = []
    lose = []
    for sig in census:
        print(sig)
        if compare_flow_and_geodesic_drilling_script_search(sig, output_filename = output_filename, max_length = 5, min_length = 1, verbose = verbose):
            win.append(sig)
        else:
            lose.append(sig)
            append_to_file(fail_filename, sig)
    return (win, lose)

def get_sig(line):
    parts = line.split('_')
    return parts[0] + '_' + parts[1]

def get_drillings(line):
    ind = line.index('[')
    return eval(line[ind:])

def census_from_data_compare_flow_and_geodesic(census_data_filename = 'data/drillings_census_eo_4_ladders_cycle_len_up_to_5.txt', 
        max_length = 5, min_length = 1, filename_suffix = "", census_start = 0, census_end = -1, quit_after_finding_one = True, verbose = 0):
    # output_filename = "../../../Dropbox/Data/compare_flow_and_geodesic" + filename_suffix + ".txt"
    # fail_filename = "../../../Dropbox/Data/compare_flow_and_geodesic_fail" + filename_suffix + ".txt"
    output_filename = "data/compare_flow_and_geodesic" + filename_suffix + ".txt"
    fail_filename = "data/compare_flow_and_geodesic_fail" + filename_suffix + ".txt"

    output_file = open(output_filename, 'w')  #write mode, clear any existing file
    output_file.close()
    fail_file = open(fail_filename, 'w')  #write mode, clear any existing file
    fail_file.close()

    census_data = parse_data_file(census_data_filename)

    if census_end != -1:
        census_data = census_data[census_start:census_end]
    else:
        census_data = census_data[census_start:]

    for line in census_data:
        sig = get_sig(line)
        print(sig)
        drillings = get_drillings(line)  ### these are sorted by drilled_sig, not by the flow cycles. We fix that to match the direct search
        drillings.sort(key = lambda x: (len(x[1]), x[1]))
        found_drilling = False
        for drilling in drillings:
            drilled_sig, fc = drilling
            if verbose > 0:
                print('flow cycle', fc)
            
            assert drilled_sig != sig  ### This happens if you try to drill a peripheral flow cycle 
            different = orbit_drilled_different_from_geodesic_drilled(sig, fc, drilled_sig, output_filename = output_filename, verbose = verbose)
            if verbose > 2:
                print('different', different)
            if different and quit_after_finding_one:
                found_drilling = True
                break
        if quit_after_finding_one and not found_drilling:
            print('no differences found for', sig)
            append_to_file(fail_filename, sig)


def compare_flow_and_geodesic_drilling_script_specific():  

    # sig = 'cPcbbbdxm_10'
    sig = 'cPcbbbiht_12'
    # sig = 'dLQacccjsnk_200' 
    # sig = 'dLQbccchhfo_122'
    # sig = 'dLQbccchhsj_122'

    fc = [(0, 0), (0, 0), (0, 5), (0, 0), (0, 5)]
    print(fc)
    out = drill_flow_cycles(sig, [fc], return_isosig_tri_angle = True) 
    drilled_sig, drilled_tri, drilled_angle = out 
    if drilled_sig != sig:
        different = orbit_drilled_different_from_geodesic_drilled(sig, fc, drilled_sig, output_filename = None, verbose = 0)
              