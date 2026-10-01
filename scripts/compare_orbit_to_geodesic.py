import snappy
from snappy_drill_homotopic import tet_and_face_indices_to_word, drill_tet_and_face_indices
from snappy.drilling.exceptions import GeodesicSystemNotSimpleError
from snappy.geometric_structure.geodesic import geodesic_start_point_info
from veering.taut import isosig_to_tri_angle
from veering.flow_cycles import flow_cycle_to_dual_edge_loop
from veering.drill_flow_cycle import drill_flow_cycles

import sys
sys.setrecursionlimit(1000000)



def compare_orbit_to_geodesic(sig, flow_cycle, verbose=True):
    tri, angle = isosig_to_tri_angle(sig)

    orbit_sig, orbit_tri, orbit_angle, cusp_mapping = drill_flow_cycles(
        sig,
        [flow_cycle],
        return_isosig_tri_angle=True,
        return_cusp_mapping=True)

    if verbose:
        print(orbit_sig, cusp_mapping)

    orbit_mfd = snappy.Manifold(orbit_tri)

    if verbose:
        print(orbit_mfd.identify())

    tet_and_face_indices = flow_cycle_to_dual_edge_loop(
        tri, angle, flow_cycle)

    mfd = snappy.Manifold(tri)
    word = tet_and_face_indices_to_word(
        mfd, tet_and_face_indices)

    if verbose:
        print(sig, flow_cycle, word)
    else:
        print(flow_cycle)

    try:
        drilled = drill_tet_and_face_indices(
            mfd, tet_and_face_indices, verified = True, bits_prec = 1000
        )
    except Exception as e:
        if verbose:
            print("SKIPPED: geodesic drilling failed")
            print("sig:", sig)
            print("flow_cycle:", flow_cycle)
            print("word:", word)
            print("error:", type(e).__name__)
        return None

    if verbose:
        print(drilled.identify())

    try:
        different = not orbit_mfd.is_isometric_to(drilled)
    except RuntimeError as e: ### when snappy cannot determine if they are isometric
        if verbose:
            print("SKIPPED: isometry test failed")
            print("sig:", sig)
            print("flow_cycle:", flow_cycle)
            print("error:", e)
        return None

	### If orbit_mfd.is_isometric_to(drilled) returns True, that certifies that the manifolds are isometric
	### If it returns False, that does not necessarily certify that they are non-isometric - it might just failed to find an isometry
	### so if different == True, use verified = True

    if different: 
    	isomsig1 = orbit_mfd.isometry_signature(verified=True)
    	isomsig2 = drilled.isometry_signature(verified=True)
    	assert isomsig1 is not None, 'isom signature failed ' + sig + ' ' + flow_cycle
    	assert isomsig2 is not None, 'isom signature failed ' + sig + ' ' + flow_cycle
    	different = (isomsig1 != isomsig2)
    
    if different:
        if verbose:
            print("DIFFERENT:", sig, flow_cycle)
        else:
            print("orbit volume:", orbit_mfd.volume())
            print("geodesic volume:", drilled.volume())

    return different

def main():
    ### examples where drilling the flow cycle and drilling the geodesic give different answers:
    ### drilling cPcbbbiht_12 along ((0, 0), (0, 0), (0, 5), (0, 0), (0, 5)) gives different results jLLwQLQbeefgehiiixxxaaxxxcv [o9_40888(0,0)(0,0)] nLvALzAAQkbeffhhikjlkmmmhaihggfhujcvcf [L14n33639(0,0)(0,0)]
    ### snappy word: 'bbCabCa'

    ### drilling dLQacccjsnk_200 along ((0, 4), (2, 2), (2, 5), (1, 1)) gives different results pLLPwvAPPAQccdfejhmjklnmnooqffaakvachckcvhw [o9_43267(0,0)(0,0)] oLLzMLLzQQcaceefiljkmnlnmnjkxccnabqqarggr [o9_41941(0,0)(0,0)]



    sig = 'cPcbbbiht_12'
    flow_cycle = ((0, 0), (0, 0), (0, 5), (0, 0), (0, 5))
    tri, angle = isosig_to_tri_angle(sig)


    orbit_sig, orbit_tri, orbit_angle, cusp_mapping = drill_flow_cycles(sig, [flow_cycle], return_isosig_tri_angle = True, return_cusp_mapping = True)
    print(orbit_sig, cusp_mapping)
    orbit_mfd = snappy.Manifold(orbit_tri)
    print(orbit_mfd.identify())

    tet_and_face_indices = flow_cycle_to_dual_edge_loop(tri, angle, flow_cycle)
    mfd = snappy.Manifold(tri)
    word = tet_and_face_indices_to_word(mfd, tet_and_face_indices)
    print(sig, flow_cycle, word)
    drilled = drill_tet_and_face_indices(mfd, tet_and_face_indices)
    print(drilled.identify())

    # output: cPcbbbiht_12 ((0, 0), (0, 0), (0, 5), (0, 0), (0, 5)) bbCabCa
    #         [L14n33639(0,0)(0,0)]

    return orbit_tri
