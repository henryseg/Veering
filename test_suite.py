#
# test_suite.py
#

from file_io import parse_data_file
from random import random

import taut
import transverse_taut
import veering
import veering_dehn_surgery

def run_tests(num_to_check = 1000):

    veering_isosigs = parse_data_file('Data/veering_census.txt')

    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert taut.is_taut(tri, angle)

    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        recovered_sig = taut.isosig_from_tri_angle(tri, angle)
        assert sig == recovered_sig
        # we only test this round trip - the other round trip does not
        # make sense because tri->isosig is many to one.

    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert transverse_taut.is_transverse_taut(tri, angle)

    non_transverse_taut_isosigs = parse_data_file('Data/veering_non_transverse_taut_examples.txt')
    for sig in non_transverse_taut_isosigs:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert not transverse_taut.is_transverse_taut(tri, angle)

    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        assert veering.is_veering(tri, angle)

    # tri, angle = taut.isosig_to_tri_angle('cPcbbbdxm_10')
    # explore_mobius_surgery_graph(tri, angle, max_tetrahedra = 12)
    # ### tests to see that it makes only veering triangulations as it goes

    for sig in veering_isosigs[:num_to_check]:
        tri, angle = taut.isosig_to_tri_angle(sig)
        for face_num in veering_dehn_surgery.get_mobius_strip_indices(tri):
            tri_s, angle_s, face_num_s = veering_dehn_surgery.veering_mobius_dehn_surgery(tri, angle, face_num)
            assert veering.is_veering(tri_s, angle_s)

    big_polys = {'cPcbbbiht_12':'a^3 - 4*a^2 + 4*a - 1',
                 'eLMkbcddddedde_2100':'a^6*b - a^6 - 2*a^5*b - a^4*b^2 + a^5 + 2*a^4*b + a^3*b^2 - 2*a^3*b + a^3 + 2*a^2*b + a*b^2 - a^2 - 2*a*b - b^2 + b',
                 'gLLAQbecdfffhhnkqnc_120012':'a^7 + a^6 + a^5 + a^4 - a^3 - a^2 - a - 1',
                 'gLLPQcdfefefuoaaauo_022110':'a^12*b^3 - a^11*b^2 - a^10*b^3 - a^10*b^2 - a^7*b^3 - a^7*b^2 - a^6*b^3 + a^7*b + a^5*b^2 - a^6 - a^5*b - a^5 - a^2*b - a^2 - a*b + 1'}

    small_polys = {'cPcbbbiht_12':'a^2 - 3*a + 1',
                   'eLMkbcddddedde_2100':'a^2*b - a^2 - a*b - b^2 + b',
                   'iLLAwQcccedfghhhlnhcqeesr_12001122':'0'}
    
    try:
        import veering_polynomial
        poly_working = True
    except:
        print 'failed to import veering_polynomial - probably not running in sage'
        poly_working = False

    if poly_working:
        for sig in big_polys:
            print 'testing big', sig
            p = veering_polynomial.big_polynomial(sig)
            assert p.__repr__() == big_polys[sig]
        for sig in small_polys:
            print 'testing small', sig
            p = veering_polynomial.small_polynomial(sig)
            assert p.__repr__() == small_polys[sig]
        for i in range(5):
            j = int(200*random())
            sig = veering_isosigs[j]
            print 'testing divide', sig
            p = veering_polynomial.big_polynomial(sig)
            q = veering_polynomial.small_polynomial(sig)
            if q == 0:
                assert p == 0
            else:
                assert q.divides(p)

    if poly_working:
        print 'all tests passed'
    else:
        print 'all tests not depending on veering_polynomial passed'

                
if __name__ == '__main__':
    run_tests()
