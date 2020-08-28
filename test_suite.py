#
# test_suite.py
#

from random import random
from hashlib import md5
from os import remove

import snappy

from file_io import parse_data_file, read_from_pickle

import taut
import transverse_taut
import veering
import veering_dehn_surgery

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

    print("all tests depending on regina/snappy passed")

    try:
        import pyx
        from boundary_triangulation import draw_triangulation_boundary_from_veering_isosig
        pyx_working = True
    except:
        print("failed to import from pyx?")
        pyx_working = False

    ladders_style_sigs = {
        "cPcbbbiht_12": "f34c1fdf65db9d02994752814803ae01",
        "gLLAQbecdfffhhnkqnc_120012": "091c85b4f4877276bfd8a955b769b496",
        "kLALPPzkcbbegfhgijjhhrwaaxnxxn_1221100101": "a0f15a8454f715f492c74ce1073a13a4",
    }

    geometric_style_sigs = {
        "cPcbbbiht_12": "1e74d0b68160c4922e85a5adb20a0f1d",
        "gLLAQbecdfffhhnkqnc_120012": "856a1fce74eb64f519bcda083303bd8f",
        "kLALPPzkcbbegfhgijjhhrwaaxnxxn_1221100101": "33bd23b34c5d977a103fa50ffe63120a",
    }

    args = {
        "draw_boundary_triangulation":True,
        "draw_triangles_near_poles": False,
        "ct_depth":-1,
        "ct_epsilon":0.03,
        "global_drawing_scale": 4,
        "delta": 0.2,
        "ladder_width": 10.0,
        "ladder_height": 20.0,
        "draw_labels": True,
    }

    shapes_data = read_from_pickle("Data/veering_shapes_up_to_ten_tetrahedra.pkl")

    if pyx_working:
        for sig in ladders_style_sigs:
            print("testing boundary triangulation pictures, ladder style", sig)
            args["tet_shapes"] = shapes_data[sig]
            args["style"] = "ladders"
            file_name = draw_triangulation_boundary_from_veering_isosig(sig, args = args) 
            f = open(file_name, "rb")
            file_hash = md5(f.read())
            assert file_hash.hexdigest() == ladders_style_sigs[sig]
            f.close()
            remove(file_name)
        
    if pyx_working:
        for sig in geometric_style_sigs:
            print("testing boundary triangulation pictures, ladder style", sig)
            args["tet_shapes"] = shapes_data[sig]
            args["style'="] = "geometric"
            file_name = draw_triangulation_boundary_from_veering_isosig(sig, args = args) 
            f = open(file_name, "rb")
            file_hash = md5(f.read())
            assert file_hash.hexdigest() == geometric_style_sigs[sig]
            f.close()
            remove(file_name)

    if pyx_working: 
        print("all tests depending on pyx passed")

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
        print("testing is_layered")
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
            print("testing hom dim", sig)
            assert (taut_polytope.taut_cone_homological_dim(sig) == 0) == (taut_polytope.LMN_tri_angle(sig) == "N")  # that is, iff

    taut_polys_with_cycles = {
        ('eLMkbcddddedde_2100',((7, 7, 0, 0, -4, 3, -7, 0),)): "a^14 - a^8 - a^7 - a^6 + 1",
    }

    if sage_working:
        for sig, cycles in taut_polys_with_cycles:
            print("testing taut with cycles", sig, cycles)
            cycles_in = [list(cycle) for cycle in cycles]
            p = taut_polynomial.taut_polynomial_via_tree(sig, cycles_in)
            assert p.__repr__() == taut_polys_with_cycles[(sig, cycles)]
            
    if sage_working:
        print("all tests depending on sage passed")


if __name__ == "__main__":
    run_tests()
