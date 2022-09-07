#
# dilatation.py
#

# Goal - compute (normalised) dilatations of (triangulation, fibre) pairs

from veering.file_io import parse_data_file, write_data_file

from veering.taut import liberal, isosig_to_tri_angle
from veering.taut_polynomial import taut_polynomial_via_fox_calculus
from veering.taut_polytope import min_carried_neg_euler

# Let's start with the case of b_1 = 1

@liberal
def dilatation_betti_one(tri, angle, normalised = "True"):
    assert tri.homology().rank() == 1  # rank of H_1 / torsion
    p = taut_polynomial_via_fox_calculus(tri, angle)
    print(p.span())
    # p = taut_polynomial_via_tree(tri, angle)
    # p = taut_polynomial_via_tree_and_smith(tri, angle)
    R = p.parent()
    a = R('a')
    q = p.polynomial(a)
    dil = max(q.real_roots())
    # print(dil)
    euler = min_carried_neg_euler(tri, angle)
    print(euler)
    # return dil**euler
