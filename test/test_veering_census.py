import pytest
import random

@pytest.fixture
def veering_sample():
    from veering import file_io, taut
    veering_census = file_io.parse_data_file("veering_census.txt")
    result = []
    for sig in random.sample(veering_census, 10):
        result.append(taut.isosig_to_tri_angle(sig))
    return result

def test_is_taut(veering_sample):
    from veering import taut
    for tri, angle in veering_sample:
        assert taut.is_taut(tri, angle)

def test_is_transverse_taut(veering_sample):
    from veering import taut, transverse_taut
    for tri, angle in veering_sample:
        assert transverse_taut.is_transverse_taut(tri, angle)

def test_is_veering(veering_sample):
    from veering import veering_tro
    for tri, angle in veering_sample:
        assert veering_tri.is_veering(tri, angle)

def test_dehn_surgery(veering_sample):
    from veering import taut, veering_tri, veering_dehn_surgery
    for tri, angle in veering_sample:
        for face_num in veering_dehn_surgery.get_mobius_strip_indices(tri):
            (tri_s, angle_s, face_num_s) = veering_dehn_surgery.veering_mobius_dehn_surgery(tri, angle, face_num)
            assert veering_tri.is_veering(tri_s, angle_s)
 
def test_fan_excision(veering_sample):
    from veering import taut, veering_tri, veering_fan_excision
    m003, _ = taut.isosig_to_tri_angle('cPcbbbdxm_10')
    m004, _ = taut.isosig_to_tri_angle('cPcbbbiht_12')
    for tri, angle in veering_sample:
        tet_types = veering_tri.is_veering(tri, angle, return_type = "tet_types")
        if tet_types.count("toggle") == 2:
            excised_tri, _ = veering_fan_excision.excise_fans(tri, angle)
            assert ( excised_tri.isIsomorphicTo(m003) != None or
                     excised_tri.isIsomorphicTo(m004) != None )




