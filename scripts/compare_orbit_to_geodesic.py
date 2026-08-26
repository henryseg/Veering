
import snappy
from snappy_drill_homotopic import tet_and_face_indices_to_word, drill_tet_and_face_indices
from veering.taut import isosig_to_tri_angle
from veering.flow_cycles import flow_cycle_to_dual_edge_loop
from veering.drill_flow_cycle import drill_flow_cycles


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
