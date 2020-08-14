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

    print "testing is_taut"
    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert taut.is_taut(tri, angle)

    print "testing isosig round trip"
    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        recovered_sig = taut.isosig_from_tri_angle(tri, angle)
        assert sig == recovered_sig
        # we only test this round trip - the other round trip does not
        # make sense because tri->isosig is many to one.

    print "testing is_transverse_taut"
    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert transverse_taut.is_transverse_taut(tri, angle)

    print "testing not is_transverse_taut"
    non_transverse_taut_isosigs = parse_data_file("Data/veering_non_transverse_taut_examples.txt")
    for sig in non_transverse_taut_isosigs:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert not transverse_taut.is_transverse_taut(tri, angle)

    print "testing is_veering"
    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert veering.is_veering(tri, angle)

    # tri, angle = taut.isosig_to_tri_angle("cPcbbbdxm_10")
    # explore_mobius_surgery_graph(tri, angle, max_tetrahedra = 12)
    # # tests to see that it makes only veering triangulations as it goes

    print "testing veering_dehn_surgery"
    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        for face_num in veering_dehn_surgery.get_mobius_strip_indices(tri):
            (tri_s, angle_s, face_num_s) = veering_dehn_surgery.veering_mobius_dehn_surgery(tri, angle, face_num)
            assert veering.is_veering(tri_s, angle_s)

    print "all tests not depending on sage passed"

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

    try:
        from sage.rings.integer_ring import ZZ
        import taut_polytope
        import taut_polynomial
        import veering_polynomial

        sage_working = True
    except:
        print "failed to import from sage?"
        sage_working = False

    if sage_working:
        for sig in veering_isosigs[:17]:
            assert taut_polytope.is_layered(sig)
        for sig in veering_isosigs[17:21]:
            assert not taut_polytope.is_layered(sig)

    if sage_working:
        for sig in veering_polys:
            print "testing veering", sig
            p = veering_polynomial.veering_polynomial(sig)
            assert p.__repr__() == veering_polys[sig]
        for sig in taut_polys:
            print "testing taut", sig
            p = taut_polynomial.taut_polynomial_via_tree(sig)
            assert p.__repr__() == taut_polys[sig]
        for i in range(3):
            j = int(5000 * random())
            sig = veering_isosigs[j]
            print "testing divide", sig
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
            print "testing alex", sig
            snap_sig = sig.split("_")[0]
            M = snappy.Manifold(snap_sig)
            if M.homology().betti_number() == 1:
                assert taut_polynomial.taut_polynomial_via_tree(sig, mode="alexander") == M.alexander_polynomial()

    if sage_working:
        for sig in torus_bundles:
            print "testing torus bundle", sig
            assert taut_polytope.is_torus_bundle(sig)

    if sage_working:
        for i in range(3):
            j = int(10000 * random())
            sig = veering_isosigs[j]
            print "testing hom dim", sig
            assert (taut_polytope.taut_cone_homological_dim(sig) == 0 and taut_polytope.LMN_tri_angle(sig) == "N") or (
                taut_polytope.taut_cone_homological_dim(sig) != 0 and taut_polytope.LMN_tri_angle(sig) != "N"
            )

    if sage_working:
        print "all tests depending on sage passed"


if __name__ == "__main__":
    run_tests()
