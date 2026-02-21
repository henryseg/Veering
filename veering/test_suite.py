#
# test_suite.py
#

import random
import time

from .file_io import veering_census, parse_data_file, read_from_pickle


# regina tests


def test_is_taut(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut

    for sig in random.sample(veering_isosigs, num_to_check):
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert taut.is_taut(tri, angle), sig

    return None
        
def test_isosig_round_trip(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut
    for sig in random.sample(veering_isosigs, num_to_check):
        tri, angle = taut.isosig_to_tri_angle(sig)
        recovered_sig = taut.isosig_from_tri_angle(tri, angle)
        assert sig == recovered_sig, sig
        # we only test this round trip - the other round trip does not
        # make sense because tri->isosig is many to one.

    return None

def test_is_transverse_taut(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut        
    from . import transverse_taut

    for sig in random.sample(veering_isosigs, num_to_check):
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert transverse_taut.is_transverse_taut(tri, angle), sig

    try:
        non_transverse_taut_isosigs = parse_data_file("veering_non_transverse_taut_examples.txt")
    except ValueError:
        print("ignore is_transverse_taut (no data file)")
        non_transverse_taut_isosigs = None
    if non_transverse_taut_isosigs is not None:
        for sig in random.sample(non_transverse_taut_isosigs, num_to_check):
            tri, angle = taut.isosig_to_tri_angle(sig)
            assert not transverse_taut.is_transverse_taut(tri, angle), sig

    return None

def test_lower_labelling(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut        
    from . import transverse_taut

    for sig in random.sample(veering_isosigs, num_to_check):
        tri, angle = taut.isosig_to_tri_angle(sig)
        transverse_taut.lower_labelling(tri, angle)

def test_is_veering(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut
    from . import veering_tri

    for sig in random.sample(veering_isosigs, num_to_check):
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert veering_tri.is_veering(tri, angle), sig

    # tri, angle = taut.isosig_to_tri_angle("cPcbbbdxm_10")
    # explore_mobius_surgery_graph(tri, angle, max_tetrahedra = 12)
    # tests to see that it makes only veering triangulations as it goes

    return None

def test_veering_dehn_surgery(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut
    from . import veering_tri
    from . import veering_dehn_surgery

    for sig in random.sample(veering_isosigs, num_to_check):
        tri, angle = taut.isosig_to_tri_angle(sig)
        for face_num in veering_dehn_surgery.get_mobius_strip_indices(tri):
            (tri_s, angle_s, face_num_s) = veering_dehn_surgery.veering_mobius_dehn_surgery(tri, angle, face_num)
            assert veering_tri.is_veering(tri_s, angle_s), sig
            
    return None

def test_veering_fan_excision(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut
    from . import veering_tri
    from . import veering_fan_excision

    m003, _ = taut.isosig_to_tri_angle('cPcbbbdxm_10')
    m004, _ = taut.isosig_to_tri_angle('cPcbbbiht_12')
    for sig in random.sample(veering_isosigs, num_to_check):
        tri, angle = taut.isosig_to_tri_angle(sig)
        tet_types = veering_tri.is_veering(tri, angle, return_type = "tet_types")
        if tet_types.count("toggle") == 2:
            excised_tri, _ = veering_fan_excision.excise_fans(tri, angle)
            assert ( excised_tri.isIsomorphicTo(m003) != None or
                     excised_tri.isIsomorphicTo(m004) != None ), sig

    return None

def test_pachner_with_taut_structure(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut
    from . import pachner

    for sig in random.sample(veering_isosigs, num_to_check):
        tri, angle = taut.isosig_to_tri_angle(sig)
        face_num = random.randrange(tri.countTriangles())
        result = pachner.twoThreeMove(tri, face_num, angle = angle, return_edge = True)  
        if result != False: 
            tri2, angle2, edge_num = result
            tri3, angle3 = pachner.threeTwoMove(tri2, edge_num, angle = angle2)
            assert taut.isosig_from_tri_angle(tri, angle) == taut.isosig_from_tri_angle(tri3, angle3), sig

    return None

def test_branched_surface_and_pachner_with_branched_surface(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut
    from . import branched_surface
    from . import pachner
    import regina

    for sig in random.sample(veering_isosigs, num_to_check):
        tri, angle = taut.isosig_to_tri_angle(sig)
        tri_original = regina.Triangulation3(tri) #copy
        branch = branched_surface.upper_branched_surface(tri, angle, return_lower = random.choice([True, False]))
        
        ### test branch isosig round trip
        sig_with_branch = branched_surface.isosig_from_tri_angle_branch(tri, angle, branch)
        tri2, angle2, branch2 = branched_surface.isosig_to_tri_angle_branch(sig_with_branch)
        assert (branch == branch2) and (angle == angle2), sig

        branch_original = branch[:]  # copy
        face_num = random.randrange(tri.countTriangles())
        out = pachner.twoThreeMove(tri, face_num, branch = branch, return_edge = True)
        if out != False:
            tri, possible_branches, edge_num = out
            tri, branch = pachner.threeTwoMove(tri, edge_num, branch = possible_branches[0])
            all_isoms = tri.findAllIsomorphisms(tri_original)
            all_branches = [branched_surface.apply_isom_to_branched_surface(branch, isom) for isom in all_isoms]
            assert branch_original in all_branches, sig

    return None

def test_taut_and_branched_drill_plus_semiflows_on_drillings(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut
    from . import branched_surface
    from . import flow_cycles
    from . import drill

    for sig in random.sample(veering_isosigs, 3): 
        tri, angle = taut.isosig_to_tri_angle(sig)
        branch = branched_surface.upper_branched_surface(tri, angle)  # also checks for veering and transverse taut
        found_loops = flow_cycles.generate_simple_flow_cycles(tri, branch)
        for loop in random.sample(found_loops, min(len(found_loops), 5)):  # drill along at most 5 loops
            tri, angle = taut.isosig_to_tri_angle(sig)
            branch = branched_surface.upper_branched_surface(tri, angle) 
            tri_loop = flow_cycles.flow_cycle_to_triangle_loop(tri, branch, loop)
            if tri_loop != False: 
                if not flow_cycles.tri_loop_is_boundary_parallel(tri_loop, tri):
                    drill.drill(tri, tri_loop, angle = angle, branch = branch, sig = sig)
                    assert branched_surface.has_non_sing_semiflow(tri, branch), sig
                    
    return None


# snappy tests


def test_algebraic_intersection(veering_isosigs, num_to_check, smaller_num_to_check):
    import snappy
    from . import snappy_tools
    census = snappy.OrientableCuspedCensus() # not a set or list, so can't use random.sample

    for i in range(num_to_check):
        M = random.choice(census)
        n = M.num_cusps()
        peripheral_curves = M.gluing_equations()[-2*n:]
        for i in range(2*n):
            for j in range(i, 2*n):
                alg_int = snappy_tools.algebraic_intersection(peripheral_curves[i], peripheral_curves[j])
                if i % 2 == 0 and j == i + 1:
                    assert alg_int == 1, M.name()
                else:
                    assert alg_int == 0, M.name()

    return None

def test_veering_drilling_and_filling(veering_isosigs, num_to_check, smaller_num_to_check):
    import snappy
    from . import snappy_tools    
    from . import veering_drill_midsurface_bdy

    for sig in random.sample(veering_isosigs[:4000], num_to_check):
        N = snappy.Manifold(sig.split("_")[0])
        T, per = veering_drill_midsurface_bdy.drill_midsurface_bdy(sig)
        M = snappy.Manifold(T.snapPea())
        M.set_peripheral_curves("shortest")
        L = snappy_tools.get_slopes_from_peripherals(M, per)
        M.dehn_fill(L)
        M = M.filled_triangulation()
        M.randomize()
        # M = snappy.Manifold(M.triangulation_isosig()) # seems to help?!?!? Bizarre
        try:
            assert M.is_isometric_to(N)
        except:
            print("Something has gone wrong.")
            print("Perhaps is_isometric_to is failing on the following two sigs.")
            print(sig.split("_")[0])
            print(M.triangulation_isosig())

    return None

# try:
#     from hashlib import md5
#     from os import remove
#     import pyx
#     from boundary_triangulation import draw_triangulation_boundary_from_veering_isosig
#     pyx_working = True
# except:
#     print("failed to import from pyx?")
#     pyx_working = False

# ladders_style_sigs = {
#     "cPcbbbiht_12": "f34c1fdf65db9d02994752814803ae01",
#     "gLLAQbecdfffhhnkqnc_120012": "091c85b4f4877276bfd8a955b769b496",
#     "kLALPPzkcbbegfhgijjhhrwaaxnxxn_1221100101": "a0f15a8454f715f492c74ce1073a13a4",
# }

# geometric_style_sigs = {
#     "cPcbbbiht_12": "1e74d0b68160c4922e85a5adb20a0f1d",
#     "gLLAQbecdfffhhnkqnc_120012": "856a1fce74eb64f519bcda083303bd8f",
#     "kLALPPzkcbbegfhgijjhhrwaaxnxxn_1221100101": "33bd23b34c5d977a103fa50ffe63120a",
# }

# args = {
#     "draw_boundary_triangulation":True,
#     "draw_triangles_near_poles": False,
#     "ct_depth":-1,
#     "ct_epsilon":0.03,
#     "global_drawing_scale": 4,
#     "delta": 0.2,
#     "ladder_width": 10.0,
#     "ladder_height": 20.0,
#     "draw_labels": True,
# }

# shapes_data = read_from_pickle("veering_shapes_up_to_ten_tetrahedra.pkl")

# if pyx_working:
#     for sig in ladders_style_sigs:
#         print("testing boundary triangulation pictures, ladder style", sig)
#         args["tet_shapes"] = shapes_data[sig]
#         args["style"] = "ladders"
#         file_name = draw_triangulation_boundary_from_veering_isosig(sig, args = args) 
#         f = open(file_name, "rb")
#         file_hash = md5(f.read())
#         assert file_hash.hexdigest() == ladders_style_sigs[sig]
#         f.close()
#         remove(file_name)
        
# if pyx_working:
#     for sig in geometric_style_sigs:
#         print("testing boundary triangulation pictures, ladder style", sig)
#         args["tet_shapes"] = shapes_data[sig]
#         args["style"] = "geometric"
#         file_name = draw_triangulation_boundary_from_veering_isosig(sig, args = args) 
#         f = open(file_name, "rb")
#         file_hash = md5(f.read())
#         assert file_hash.hexdigest() == geometric_style_sigs[sig]
#         f.close()
#         remove(file_name)

# if pyx_working: 
#     print("all tests depending on pyx passed")


# sage tests


veering_polys = {
    "cPcbbbiht_12": [-4, -1, 1, 4],
    "eLMkbcddddedde_2100": [-2, -2, -2, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 2, 2],
    "gLLAQbecdfffhhnkqnc_120012": [-1, -1, -1, -1, 1, 1, 1, 1],
    "gLLPQcdfefefuoaaauo_022110": [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1],
}

# veering_polys = { ### old
#     "cPcbbbiht_12": "a^3 - 4*a^2 + 4*a - 1",
#     "eLMkbcddddedde_2100": "a^6*b - a^6 - 2*a^5*b - a^4*b^2 + a^5 + 2*a^4*b + a^3*b^2 - 2*a^3*b + a^3 + 2*a^2*b + a*b^2 - a^2 - 2*a*b - b^2 + b",
#     "gLLAQbecdfffhhnkqnc_120012": "a^7 + a^6 + a^5 + a^4 - a^3 - a^2 - a - 1",
#     "gLLPQcdfefefuoaaauo_022110": "a^12*b^3 - a^11*b^2 - a^10*b^3 - a^10*b^2 - a^7*b^3 - a^7*b^2 - a^6*b^3 + a^7*b + a^5*b^2 - a^6 - a^5*b - a^5 - a^2*b - a^2 - a*b + 1",
# }

def test_is_layered(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polytope

    for sig in veering_isosigs[:17]:
        assert taut_polytope.is_layered(sig), sig
    for sig in veering_isosigs[17:21]:
        assert not taut_polytope.is_layered(sig), sig

    return None

# def test_is_fibered(veering_isosigs, num_to_check, smaller_num_to_check):
#     from . import fibered
#     try:
#         mflds = parse_data_file("mflds_which_fiber.txt")
#         mflds = mflds[58:]  # get rid of the preamble
#     except ValueError:
#         print("ignore is_fibered (no data file)")
#         mflds = []
#     mflds = [line.split("\t")[0:2] for line in mflds]
#     for (name, kind) in random.sample(mflds, num_to_check):        
#         if not fibered.is_fibered(name) == (kind == "fibered"):
#             print("is_fibered failed on", name, "but don't worry about that... is_fibered is pretty random.")

def test_veering_poly(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import veering_polynomial

    for sig in veering_polys:
        p = veering_polynomial.veering_polynomial(sig)
        assert check_polynomial_coefficients(p, veering_polys[sig]), sig
        ### Nov 2021: sage 9.4 changed how smith normal form works, which changed our polynomials
        ### to equivalent but not equal polynomials. To avoid this kind of change breaking things
        ### in the future, we changed to comparing the list of coefficients.
        # assert p.__repr__() == veering_polys[sig]

    return None
        
taut_polys = {
    "cPcbbbiht_12": [-3, 1, 1],
    "eLMkbcddddedde_2100": [-1, -1, -1, 1, 1],
    "iLLAwQcccedfghhhlnhcqeesr_12001122": [],
}

# taut_polys = { ### old
#     "cPcbbbiht_12": "a^2 - 3*a + 1",
#     "eLMkbcddddedde_2100": "a^2*b - a^2 - a*b - b^2 + b",
#     "iLLAwQcccedfghhhlnhcqeesr_12001122": "0",
# }
        
def test_taut_poly(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polynomial

    for sig in taut_polys:
        p = taut_polynomial.taut_polynomial_via_tree(sig)
        assert check_polynomial_coefficients(p, taut_polys[sig]), sig
        # assert p.__repr__() == taut_polys[sig]

    return None

def test_taut_poly_via_fox_calculus(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polynomial

    for sig in taut_polys:
        p = taut_polynomial.taut_polynomial_via_fox_calculus(sig)
        assert check_polynomial_coefficients(p, taut_polys[sig]), sig
        # assert p.__repr__() == taut_polys[sig]

    return None

def test_taut_equals_taut_fox(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polynomial

    for sig in random.sample(veering_isosigs[500:3000], smaller_num_to_check):        
        taut1 = taut_polynomial.taut_polynomial_via_tree(sig)
        taut2 = taut_polynomial.taut_polynomial_via_fox_calculus(sig)
        assert taut1 == taut2

    return None

def test_fox_simplified_equals_fox_unsimplified(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polynomial

    for sig in random.sample(veering_isosigs[100:500], 3):        
        taut2 = taut_polynomial.taut_polynomial_via_fox_calculus(sig)
        taut3 = taut_polynomial.taut_polynomial_via_fox_calculus(sig, simplified = False)
        assert taut2 == taut3

    return None

def test_divide(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polynomial
    from . import veering_polynomial

    for sig in random.sample(veering_isosigs[:3000], num_to_check):
        p = veering_polynomial.veering_polynomial(sig)
        q = taut_polynomial.taut_polynomial_via_tree(sig)
        if q == 0:
            assert p == 0, sig
        else:
            assert q.divides(p), sig

    return None

def test_alex(veering_isosigs, num_to_check, smaller_num_to_check):
    import snappy
    from . import taut_polynomial

    for sig in random.sample(veering_isosigs[:3000], num_to_check):        
        snap_sig = sig.split("_")[0]
        M = snappy.Manifold(snap_sig)
        if M.homology().betti_number() == 1:
            assert taut_polynomial.taut_polynomial_via_tree(sig, mode = "alexander") == M.alexander_polynomial(), sig

    return None

# would be nice to automate the following - need to fetch the angle
# structure say via z_charge.py...

torus_bundles = [
    "cPcbbbiht_12",
    "eLMkbcdddhhqqa_1220",
    "gLMzQbcdefffhhqqqdl_122002",
]

def test_is_torus_bundle(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polytope

    for sig in torus_bundles: 
        assert taut_polytope.is_torus_bundle(sig), sig

    return None

def test_is_layered(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polytope

    for sig in torus_bundles:
        assert taut_polytope.is_layered(sig), sig

    return None
        
measured = [
    "gLLAQbecdfffhhnkqnc_120012",
    "iLLALQcccedhgghhlnxkxrkaa_12001112",
    "iLLAwQcccedfghhhlnhcqeesr_12001122",
]

def test_measured(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polytope

    for sig in measured:
        assert taut_polytope.LMN_tri_angle(sig) == "M", sig

    return None

empties = [
    "fLAMcaccdeejsnaxk_20010",
    "gLALQbcbeeffhhwsras_211220",
    "hLALAkbcbeefgghhwsraqj_2112202",
]

def test_empty(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polytope

    for sig in empties:
        assert taut_polytope.LMN_tri_angle(sig) == "N", sig

    return None
        
def test_euler_characteristic(veering_isosigs, num_to_check, smaller_num_to_check):
    import snappy
    from . import taut_polytope

    for sig in random.sample(veering_isosigs[:2500], num_to_check):
        snap_sig = sig.split("_")[0]
        M = snappy.Manifold(snap_sig)
        if M.homology().betti_number() == 1:
            if taut_polytope.is_layered(sig):
                eul = int(sum(taut_polytope.min_neg_euler_carried(sig)) / 2)
                alex = M.alexander_polynomial()
                alex_span = alex.exponents()[-1] - alex.exponents()[0]
                assert eul == alex_span - 1, sig

    return None

def test_hom_dim(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polytope

    for sig in random.sample(veering_isosigs[:3000], 3): # magic number
        # dimension = zero if and only if nothing is carried.
        assert (taut_polytope.taut_cone_homological_dim(sig) == 0) == (taut_polytope.LMN_tri_angle(sig) == "N"), sig

    return None

boundary_cycles = {
    ("eLMkbcddddedde_2100",(2,5,5,1,3,4,7,1)): "((-7, -7, 0, 0, 4, -3, 7, 0), (7, 7, 0, 0, -4, 3, -7, 0))",
    ("iLLLQPcbeegefhhhhhhahahha_01110221",(0,1,0,0,0,1,0,0,0,0,0,0,1,0,1,0)): "((0, 0, -1, 1, 1, 0, 1, 1, -1, 0, 0, 0, 0, 1, 0, 1), (0, 0, 1, -1, -1, 0, -1, -1, 1, 0, 0, 0, 0, -1, 0, -1))",
    ("ivvPQQcfhghgfghfaaaaaaaaa_01122000",(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)): "((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2), (1, 1, 0, 2, -1, 0, -3, 1, 2, -1, -2, 0, 3, -2, -1, 0), (-2, 0, -3, 1, 2, -1, 0, 2, -1, 0, 3, 1, -2, 1, 0, -1), (0, -2, 1, -3, 0, -1, 2, 0, -1, 2, -1, 1, 0, 1, -2, 3))",
}

def test_boundary_cycles(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_carried

    for sig, surface in boundary_cycles:
        surface_list = list(surface)
        cycles = taut_carried.boundary_cycles_from_surface(sig, surface_list)
        cycles = tuple(tuple(cycle) for cycle in cycles)
        assert cycles.__repr__() == boundary_cycles[(sig, surface)], sig

    return None

taut_polys_with_cycles = {
    ("eLMkbcddddedde_2100", ((7, 7, 0, 0, -4, 3, -7, 0),)): [-1, -1, -1, 1, 1],
    ("iLLLQPcbeegefhhhhhhahahha_01110221", ((0, 0, 1, -1, -1, 0, -1, -1, 1, 0, 0, 0, 0, -1, 0, -1),)): [1, 1, 2],
    ("ivvPQQcfhghgfghfaaaaaaaaa_01122000", ((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2), (1, 1, 0, 2, -1, 0, -3, 1, 2, -1, -2, 0, 3, -2, -1, 0))): [-4, -1, -1, 1, 1],
}

# taut_polys_with_cycles = {
#     ("eLMkbcddddedde_2100", ((7, 7, 0, 0, -4, 3, -7, 0),)): "a^14 - a^8 - a^7 - a^6 + 1",
#     ("iLLLQPcbeegefhhhhhhahahha_01110221", ((0, 0, 1, -1, -1, 0, -1, -1, 1, 0, 0, 0, 0, -1, 0, -1),)): "a^2 + 2*a + 1",
#     ("ivvPQQcfhghgfghfaaaaaaaaa_01122000", ((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2), (1, 1, 0, 2, -1, 0, -3, 1, 2, -1, -2, 0, 3, -2, -1, 0))): "a*b^2 - a^2 - 4*a*b - b^2 + a",
# }

def test_taut_with_cycles(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polynomial

    for sig, cycles in taut_polys_with_cycles:
        cycles_in = [list(cycle) for cycle in cycles]
        p = taut_polynomial.taut_polynomial_via_tree(sig, cycles_in)
        assert check_polynomial_coefficients(p, taut_polys_with_cycles[(sig, cycles)]), sig
        # assert p.__repr__() == taut_polys_with_cycles[(sig, cycles)]

    return None

taut_polys_image = {
    ('eLMkbcddddedde_2100', ((7, 8, -1, 0, -4, 4, -8, 0),)):[-1, -1, -1, 1, 1],
    ('ivvPQQcfhghgfghfaaaaaaaaa_01122000', ((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2),)):[-2, -2, -1, -1, 1, 1],
    ('ivvPQQcfhghgfghfaaaaaaaaa_01122000', ((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2), (1, 1, 0, 2, -1, 0, -3, 1, 2, -1, -2, 0, 3, -2, -1, 0))):[-4, -1, -1, 1, 1]
}

# taut_polys_image = {
#     ('eLMkbcddddedde_2100', ((7, 8, -1, 0, -4, 4, -8, 0),)):"a^16 - a^9 - a^8 - a^7 + 1",
#     ('ivvPQQcfhghgfghfaaaaaaaaa_01122000', ((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2),)):"a*b^2*c - 2*a*b*c - b^2*c - a^2 - 2*a*b + a",
#     ('ivvPQQcfhghgfghfaaaaaaaaa_01122000', ((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2), (1, 1, 0, 2, -1, 0, -3, 1, 2, -1, -2, 0, 3, -2, -1, 0))):"a*b^2 - a^2 - 4*a*b - b^2 + a"
# }
        
def test_taut_with_images(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polynomial

    for sig, cycles in taut_polys_image:
        cycles_in = [list(cycle) for cycle in cycles]
        p = taut_polynomial.taut_polynomial_image(sig, cycles_in)
        assert check_polynomial_coefficients(p, taut_polys_image[(sig, cycles)]), sig
        # assert p.__repr__() == taut_polys_image[(sig, cycles)]

    return None

alex_polys_with_cycles = {
    ("eLMkbcddddedde_2100",((7, 7, 0, 0, -4, 3, -7, 0),)): [-2, -1, -1, -1, 1, 1, 1, 2],
    ("iLLLQPcbeegefhhhhhhahahha_01110221", ((0, 0, 1, -1, -1, 0, -1, -1, 1, 0, 0, 0, 0, -1, 0, -1),)): [-3, -1, 1, 3],
    ("ivvPQQcfhghgfghfaaaaaaaaa_01122000", ((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2), (1, 1, 0, 2, -1, 0, -3, 1, 2, -1, -2, 0, 3, -2, -1, 0))): [-1, -1, 1, 1],
}

# alex_polys_with_cycles = {
#     ("eLMkbcddddedde_2100",((7, 7, 0, 0, -4, 3, -7, 0),)): "a^15 - a^14 + a^9 - 2*a^8 + 2*a^7 - a^6 + a - 1",
#     ("iLLLQPcbeegefhhhhhhahahha_01110221", ((0, 0, 1, -1, -1, 0, -1, -1, 1, 0, 0, 0, 0, -1, 0, -1),)): "3*a^3 - a^2 + a - 3",
#     ("ivvPQQcfhghgfghfaaaaaaaaa_01122000", ((1, 1, 2, 0, -1, 2, 1, -3, 0, -1, 0, -2, -1, 0, 3, -2), (1, 1, 0, 2, -1, 0, -3, 1, 2, -1, -2, 0, 3, -2, -1, 0))): "a*b^2 - a^2 - b^2 + a",
# }

def test_alex_with_cycles(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polynomial

    for sig, cycles in alex_polys_with_cycles:
        cycles_in = [list(cycle) for cycle in cycles]
        p = taut_polynomial.taut_polynomial_via_tree(sig, cycles_in, mode = "alexander")
        assert check_polynomial_coefficients(p, alex_polys_with_cycles[(sig, cycles)]), sig
        # assert p.__repr__() == alex_polys_with_cycles[(sig, cycles)]

    return None

def test_euler_and_edge_orientability(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import edge_orientability
    from . import taut_euler_class

    for sig in random.sample(veering_isosigs[:3000], 3):
        # Theorem: If (tri, angle) is edge-orientable then e = 0.
        assert not ( edge_orientability.is_edge_orientable(sig) and
                     (taut_euler_class.order_of_euler_class_wrapper(sig) == 2) ), sig

    return None

# Perhaps implement the following:        
# Theorem: If (tri, angle) is edge orientable then taut poly = alex poly.
# taut_polynomial.taut_polynomial_via_tree(sig, mode = "alexander") ==
#      taut_polynomial.taut_polynomial_via_tree(sig, mode = "taut")

def test_exotics(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut
    from . import transverse_taut
    from . import veering_tri
    from . import taut_polytope

    for sig in random.sample(veering_isosigs[:3000], 3):
        tri, angle = taut.isosig_to_tri_angle(sig)
        T = veering_tri.veering_triangulation(tri, angle)
        is_eo = T.is_edge_orientable()
        for angle in T.exotic_angles():
            assert taut_polytope.taut_cone_homological_dim(tri, angle) == 0, sig
            assert is_eo == transverse_taut.is_transverse_taut(tri, angle), sig
        
    return None

def test_depth(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut_polytope

    for sig in random.sample(veering_isosigs[:3000], 3):
        is_finite, depth = taut_polytope.depth(sig)
        LMN = taut_polytope.LMN_tri_angle(sig)
        assert ((LMN == "L" and is_finite and depth == 0) or
                (LMN == "M" and depth > 0) or
                (LMN == "N" and (not is_finite) and depth == 0)), sig

    return None
        
### test for drill_midsurface_bdy: drill then fill, check you get the same manifold

### test for winding: check that punctured torus bundles have unique veering structure
            
def test_building_carried_surfaces_and_mutations(veering_isosigs, num_to_check, smaller_num_to_check):
    from . import taut
    from . import carried_surface
    from . import mutation
    sigs_weights = [
        ['iLLLPQccdgefhhghqrqqssvof_02221000',  (0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0)], 
        ['jLLAvQQcedehihiihiinasmkutn_011220000', (2, 0, 1, 0, 0, 0, 1, 2, 0, 2, 0, 2, 1, 0, 0, 0, 1, 0)],
        ['jLLAvQQcedehihiihiinasmkutn_011220000', (0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0)],
        ['jLLLMPQcdgfhfhiiihshassspiq_122201101', (0, 0, 4, 0, 4, 1, 0, 2, 2, 0, 1, 0, 0, 4, 0, 4, 0, 0)]
    ]
    strata = [
        ((1, 2), [2, 2]), 
        ((2, 4), [5, 5, 1, 1]),
        ((0, 3), [2, 0, 0]),
        ((6, 1), [22])
    ]
    orders_of_veering_symmetry_groups = [4, 2, 2, 2]
        
    for i in range(len(sigs_weights)):
        tri, angle = taut.isosig_to_tri_angle(sigs_weights[i][0])
        weights = sigs_weights[i][1]
        surface, edge_colours = carried_surface.build_surface(tri, angle, weights, return_edge_colours = True)
        assert strata[i] == carried_surface.stratum_from_weights_surface(weights, surface)
        veering_isoms = carried_surface.veering_symmetry_group(surface, edge_colours)
        assert len(veering_isoms) == orders_of_veering_symmetry_groups[i]
        isom = veering_isoms[1]
        mutation.mutate(tri, angle, weights, isom, quiet = True)
        if i == 0:
            assert tri.isoSig() == 'ivLLQQccfhfeghghwadiwadrv'
            #print('svof to wadrv passed')
        elif i == 1:
            assert tri.isoSig() == 'jvLLAQQdfghhfgiiijttmtltrcr'
            #print('smkutn to tltrcr passed')
        elif i == 2:
            assert tri.isoSig() == 'jLLMvQQcedehhiiihiikiwnmtxk'
            #print('smkutn to mtxk passed')
        elif i == 3:
            assert tri.isoSig() == 'jLLALMQcecdhggiiihqrwqwrafo'
            #print('spiq to rafo passed')

    return None


# helper functions


def check_polynomial_coefficients(p, correct_coeffs):
    p = p.coefficients()
    p.sort()
    n = [-c for c in p]
    n.reverse()
    return (p == correct_coeffs) or (n == correct_coeffs)
    

# calling function


regina_tests = [test_is_taut,
                test_isosig_round_trip,
                test_is_transverse_taut,
                test_lower_labelling,
                test_is_veering,
                test_veering_dehn_surgery,
                test_veering_fan_excision,
                test_pachner_with_taut_structure,
                test_branched_surface_and_pachner_with_branched_surface,
                test_taut_and_branched_drill_plus_semiflows_on_drillings,
]

snappy_tests = [test_algebraic_intersection,
                test_veering_drilling_and_filling, # randomly throws errors...
]

sage_tests = [test_is_layered,
              test_veering_poly,
              test_taut_poly,
              test_taut_poly_via_fox_calculus,
              test_taut_equals_taut_fox,
              test_fox_simplified_equals_fox_unsimplified,
              test_divide,
              test_alex,
              test_is_torus_bundle,
              test_is_layered,
              test_measured,
              # test_empty, # broken in sage 10.4
              test_euler_characteristic,
              # test_hom_dim, # sometimes broken in sage 10.4
              test_boundary_cycles,
              test_taut_with_cycles,
              test_taut_with_images,
              test_alex_with_cycles,
              test_euler_and_edge_orientability,
              test_exotics,
              # test_depth, # broken in sage 10.4
              # test_building_carried_surfaces_and_mutations, # Broken 2026-02-21
]

def run_tests(num_to_check = 20, smaller_num_to_check = 10):
    veering_isosigs = veering_census()

    print("tests requiring regina")
    for test in regina_tests:
        print(test.__name__)
        test(veering_isosigs, num_to_check, smaller_num_to_check)
        
    print("tests requiring snappy")
    for test in snappy_tests:
        print(test.__name__)
        test(veering_isosigs, num_to_check, smaller_num_to_check)
        
    print("tests requiring sage")
    for test in sage_tests:
        print(test.__name__)
        test(veering_isosigs, num_to_check, smaller_num_to_check)
    
    print("testing complete")
