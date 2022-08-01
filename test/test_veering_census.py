import pytest
import random

@pytest.fixture
def veering_census():
    from veering.file_io import parse_data_file
    return parse_data_file("veering_census.txt")

def test_is_taut(veering_census):
    from veering import taut
    for sig in random.sample(veering_census, 10):
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert taut.is_taut(tri, angle), sig

def test_is_transverse_taut(veering_census):
    from veering import taut, transverse_taut
    for sig in random.sample(veering_census, 10):
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert transverse_taut.is_transverse_taut(tri, angle), sig

def test_is_veering(veering_census):
    from veering import taut, veering
    for sig in random.sample(veering_census, 10):
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert veering.is_veering(tri, angle), sig

def test_dehn_surgery(veering_census):
    from veering import taut, veering, veering_dehn_surgery
    for sig in random.sample(veering_census, 10):
        tri, angle = taut.isosig_to_tri_angle(sig)
        for face_num in veering_dehn_surgery.get_mobius_strip_indices(tri):
            (tri_s, angle_s, face_num_s) = veering_dehn_surgery.veering_mobius_dehn_surgery(tri, angle, face_num)
            assert veering.is_veering(tri_s, angle_s), sig
 
def test_fan_excision(veering_census):
    from veering import taut, veering, veering_fan_excision
    m003, _ = taut.isosig_to_tri_angle('cPcbbbdxm_10')
    m004, _ = taut.isosig_to_tri_angle('cPcbbbiht_12')
    for sig in random.sample(veering_census, 10):
        tri, angle = taut.isosig_to_tri_angle(sig)
        tet_types = veering.is_veering(tri, angle, return_type = "tet_types")
        if tet_types.count("toggle") == 2:
            excised_tri, _ = veering_fan_excision.excise_fans(tri, angle)
            assert ( excised_tri.isIsomorphicTo(m003) != None or
                     excised_tri.isIsomorphicTo(m004) != None ), sig




