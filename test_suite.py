#
# test_suite.py
#

from file_io import parse_data_file
from random import random

import taut
import transverse_taut
import veering
import veering_dehn_surgery

import snappy


# def run_tests(num_to_check = 10):
def run_tests(num_to_check=1000):

    veering_isosigs = parse_data_file("Data/veering_census.txt")

    print("testing is_taut")
    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert taut.is_taut(tri, angle)

    print("testing isosig round trip")
    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        recovered_sig = taut.isosig_from_tri_angle(tri, angle)
        assert sig == recovered_sig
        # we only test this round trip - the other round trip does not
        # make sense because tri->isosig is many to one.

    print("testing is_transverse_taut")
    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert transverse_taut.is_transverse_taut(tri, angle)

    print("testing not is_transverse_taut")
    non_transverse_taut_isosigs = parse_data_file("Data/veering_non_transverse_taut_examples.txt")
    for sig in non_transverse_taut_isosigs:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert not transverse_taut.is_transverse_taut(tri, angle)

    print("testing is_veering")
    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert veering.is_veering(tri, angle)

    # tri, angle = taut.isosig_to_tri_angle("cPcbbbdxm_10")
    # explore_mobius_surgery_graph(tri, angle, max_tetrahedra = 12)
    # # tests to see that it makes only veering triangulations as it goes

    print("testing veering_dehn_surgery")
    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        for face_num in veering_dehn_surgery.get_mobius_strip_indices(tri):
            (tri_s, angle_s, face_num_s) = veering_dehn_surgery.veering_mobius_dehn_surgery(tri, angle, face_num)
            assert veering.is_veering(tri_s, angle_s)

    print("all tests not depending on sage passed")

    veering_polys = {
        "cPcbbbiht_12": "a^3 - 4*a^2 + 4*a - 1",
        "eLMkbcddddedde_2100": "a^6*b - a^6 - 2*a^5*b - a^4*b^2 + a^5 + 2*a^4*b + a^3*b^2 - 2*a^3*b + a^3 + 2*a^2*b + a*b^2 - a^2 - 2*a*b - b^2 + b",
        "gLLAQbecdfffhhnkqnc_120012": "a^7 + a^6 + a^5 + a^4 - a^3 - a^2 - a - 1",
        "gLLPQcdfefefuoaaauo_022110": "a^12*b^3 - a^11*b^2 - a^10*b^3 - a^10*b^2 - a^7*b^3 - a^7*b^2 - a^6*b^3 + a^7*b + a^5*b^2 - a^6 - a^5*b - a^5 - a^2*b - a^2 - a*b + 1",
    }

    taut_polys = {
        "cPcbbbiht_12": "a^2 - 3*a + 1",
        "eLMkbcddddedde_2100": "a^2*b - a^2 - a*b - b^2 + b",
        "iLLAwQcccedfghhhlnhcqeesr_12001122": "0",
    }

    torus_bundles = [
        "cPcbbbiht_12",
        "eLMkbcdddhhqqa_1220",
        "gLMzQbcdefffhhqqqdl_122002",
    ]

    measured = [
        "gLLAQbecdfffhhnkqnc_120012",
        "iLLALQcccedhgghhlnxkxrkaa_12001112",
        "iLLAwQcccedfghhhlnhcqeesr_12001122",
    ]

    empties = [
        "fLAMcaccdeejsnaxk_20010",
        "gLALQbcbeeffhhwsras_211220",
        "hLALAkbcbeefgghhwsraqj_2112202",
    ]
    
    try:
        from sage.rings.integer_ring import ZZ
        import taut_polytope
        import taut_polynomial
        import veering_polynomial

        sage_working = True
    except:
        print("failed to import from sage?")
        sage_working = False

    if sage_working:
        for sig in veering_isosigs[:17]:
            assert taut_polytope.is_layered(sig)
        for sig in veering_isosigs[17:21]:
            assert not taut_polytope.is_layered(sig)

    if sage_working:
        for sig in veering_polys:
            print("testing veering", sig)
            p = veering_polynomial.veering_polynomial(sig)
            assert p.__repr__() == veering_polys[sig]
        for sig in taut_polys:
            print("testing taut", sig)
            p = taut_polynomial.taut_polynomial_via_tree(sig)
            assert p.__repr__() == taut_polys[sig]
        for i in range(3):
            j = int(5000 * random())
            sig = veering_isosigs[j]
            print("testing divide", sig)
            p = veering_polynomial.veering_polynomial(sig)
            q = taut_polynomial.taut_polynomial_via_tree(sig)
            if q == 0:
                assert p == 0
            else:
                assert q.divides(p)

    if sage_working:
        for i in range(3):
            j = int(5000 * random())
            sig = veering_isosigs[j]
            print("testing alex", sig)
            snap_sig = sig.split("_")[0]
            M = snappy.Manifold(snap_sig)
            if M.homology().betti_number() == 1:
                assert taut_polynomial.taut_polynomial_via_tree(sig, mode="alexander") == M.alexander_polynomial()

    if sage_working:
        for sig in torus_bundles:
            print("testing torus bundle", sig)
            assert taut_polytope.is_torus_bundle(sig)

    if sage_working:
        for sig in torus_bundles:
            print("testing layered", sig)
            assert taut_polytope.is_layered(sig)
        for sig in measured:
            print("testing measured", sig)
            assert taut_polytope.LMN_tri_angle(sig) == "M"
        for sig in empties:
            print("testing empty", sig)
            assert taut_polytope.LMN_tri_angle(sig) == "N"

    if sage_working:  # warning - this takes random amounts of time!
        for i in range(3):
            j = int(10000 * random())
            sig = veering_isosigs[j]
            print(("testing hom dim", sig))
            # dimension = zero if and only if nothing is carried.
            assert (taut_polytope.taut_cone_homological_dim(sig) == 0) == (taut_polytope.LMN_tri_angle(sig) == "N")


    import taut_carried     

    boundary_cycles = {
        ('eLMkbcddddedde_2100',(2,5,5,1,3,4,7,1)): "((-7, -7, 0, 0, 4, -3, 7, 0), (7, 7, 0, 0, -4, 3, -7, 0))",
        ('iLLLQPcbeegefhhhhhhahahha_01110221',(0,1,0,0,0,1,0,0,0,0,0,0,1,0,1,0)): "((0, 0, -1, 1, 1, 0, 1, 1, -1, 0, 0, 0, 0, 1, 0, 1), (0, 0, 1, -1, -1, 0, -1, -1, 1, 0, 0, 0, 0, -1, 0, -1))",
        ('ivvPQQcfhghgfghfaaaaaaaaa_01122000',(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)): "((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2), (1, 1, 0, 2, -1, 0, -3, 1, 2, -1, -2, 0, 3, -2, -1, 0), (-2, 0, -3, 1, 2, -1, 0, 2, -1, 0, 3, 1, -2, 1, 0, -1), (0, -2, 1, -3, 0, -1, 2, 0, -1, 2, -1, 1, 0, 1, -2, 3))",
    }

    taut_polys_with_cycles = {
        ('eLMkbcddddedde_2100',((7, 7, 0, 0, -4, 3, -7, 0),)): "a^14 - a^8 - a^7 - a^6 + 1",
        ('iLLLQPcbeegefhhhhhhahahha_01110221', ((0, 0, 1, -1, -1, 0, -1, -1, 1, 0, 0, 0, 0, -1, 0, -1),)): "a^2 + 2*a + 1",
        ('ivvPQQcfhghgfghfaaaaaaaaa_01122000', ((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2), (1, 1, 0, 2, -1, 0, -3, 1, 2, -1, -2, 0, 3, -2, -1, 0))): "a*b^2 - a^2 - 4*a*b - b^2 + a",
    }

    alex_polys_with_cycles = {
        ('eLMkbcddddedde_2100',((7, 7, 0, 0, -4, 3, -7, 0),)): "a^15 - a^14 + a^9 - 2*a^8 + 2*a^7 - a^6 + a - 1",
        ('iLLLQPcbeegefhhhhhhahahha_01110221', ((0, 0, 1, -1, -1, 0, -1, -1, 1, 0, 0, 0, 0, -1, 0, -1),)): "3*a^3 - a^2 + a - 3",
        ('ivvPQQcfhghgfghfaaaaaaaaa_01122000', ((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2), (1, 1, 0, 2, -1, 0, -3, 1, 2, -1, -2, 0, 3, -2, -1, 0))): "a*b^2 - a^2 - b^2 + a",
    }

    if sage_working:
        for sig, surface in boundary_cycles:
            print("testing boundary cycles", sig, surface)
            surface_list = list(surface)
            cycles = taut_carried.boundary_cycles_from_surface(sig,surface_list)
            cycles = tuple(tuple(cycle) for cycle in cycles)
            assert cycles.__repr__() == boundary_cycles[(sig, surface)]

    if sage_working:
        for sig, cycles in taut_polys_with_cycles:
            print("testing taut with cycles", sig, cycles)
            cycles_in = [list(cycle) for cycle in cycles]
            p = taut_polynomial.taut_polynomial_via_tree(sig, cycles_in)
            assert p.__repr__() == taut_polys_with_cycles[(sig, cycles)]

    if sage_working:
        for sig, cycles in alex_polys_with_cycles:
            print("testing Alex with cycles", sig, cycles)
            cycles_in = [list(cycle) for cycle in cycles]
            p = taut_polynomial.taut_polynomial_via_tree(sig, cycles_in, mode='alexander')
            assert p.__repr__() == alex_polys_with_cycles[(sig, cycles)]

    if sage_working:
        for i in range(3):
            j = int(5000 * random())
            sig = veering_isosigs[j]
            print(("testing euler and edge orientability", sig))
            # Theorem: If (tri, angle) is edge orientable then e = 0.
            assert not ( edge_orientability.is_edge_orientable(sig) and
                         (taut_euler_class.order_of_euler_class_wrapper(sig) == 2) )

    if sage_working:
        # Theorem: If (tri, angle) is edge orientable then taut poly = alex poly.
        # taut_polynomial.taut_polynomial_via_tree(sig, mode = "alexander") ==
        #      taut_polynomial.taut_polynomial_via_tree(sig, mode = "taut")
        pass
            
    if sage_working:
        print("all tests depending on sage passed")


if __name__ == "__main__":
    run_tests()
