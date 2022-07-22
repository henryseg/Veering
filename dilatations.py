#
# dilatation.py
#

# Goal - compute (normalised) dilatations of (triangulation, fibre) pairs

from file_io import parse_data_file, write_data_file

from taut import liberal, isosig_to_tri_angle
from taut_polynomial import taut_polynomial_via_tree, taut_polynomial_via_tree_and_smith
from taut_polytope import min_carried_neg_euler

# Let's start with the case of b_1 = 1

@liberal
def dilatation_betti_one(tri, angle, normalised = "True"):
    assert tri.homology().rank() == 1  # rank of H_1 / torsion
    p = taut_polynomial_via_tree(tri, angle)
    # p = taut_polynomial_via_tree_and_smith(tri, angle)
    R = p.parent()
    a = R('a')
    q = p.polynomial(a)
    dil = max(q.real_roots())
    euler = min_carried_neg_euler(tri, angle)
    return dil**euler

def main(report = 50):
    data = parse_data_file("Data/veering_census_with_data.txt")
    out_filename = "Data/betti_one_dilatations_tree.txt"
    # out_filename = "Data/betti_one_dilatations_tree_and_smith.txt"
    out = [] 
    for i, line in enumerate(data): 
        line = line.split(" ") 
        sig = line[0] 
        tri, angle = isosig_to_tri_angle(sig) 
        if tri.homology().rank() == 1 and line[1] == "F0": 
            out.append( [sig, str(dilatation_betti_one(tri, angle))] ) 
        if i % report == 0: 
            print(i, len(out)) 
            write_data_file(out, out_filename)
    write_data_file(out, out_filename)
