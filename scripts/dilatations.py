#
# dilatations.py
#

# Goal - compute (normalised) dilatations of (triangulation, fibre) pairs

from sage.functions.log import log

from snappy import Manifold

from veering.file_io import parse_data_file, write_data_file
from veering.taut import liberal, isosig_to_tri_angle
from veering.taut_polynomial import taut_polynomial_via_fox_calculus
from veering.taut_polytope import min_neg_euler_carried


@liberal
def dilatation_betti_one_fibred(tri, angle):
    # Old way, with linear programming
    # if V_fab.is_edge_orientable():
    #     euler = span - 1  # See Parlak 2021
    # else:
    #     euler = int(sum(min_neg_euler_carried(tri, angle, solver))/2)

    M = Manifold(tri)
    q = M.alexander_polynomial()
    q_span = q.exponents()[-1] - q.exponents()[0]
    euler = q_span - 1

    p = taut_polynomial_via_fox_calculus(tri, angle)
    p_span = p.exponents()[0][0] - p.exponents()[-1][0]  # NB
    R = p.parent()
    a = R('a')
    p = p.polynomial(a)
    dil = max(p.real_roots())  # rigour?
    if euler != p_span - 1 and euler != p_span:
        print(sig, p_span, euler)
    return dil, euler


def dilatation_script(report = 100, start = 0, end = 87047):  #, solver = "GPLK"):
    data = parse_data_file("veering_census_with_data.txt")
    data = data[start:end]
    out_filename = "betti_one_dilatations_fox_alex.txt"
    # out_filename = "betti_one_dilatations_tree_fox_cplex.txt"
    # out_filename = "betti_one_dilatations_tree_and_smith.txt"
    out = [] 
    for i, line in enumerate(data): 
        line = line.split(" ") 
        sig = line[0]
        if i % report == 0:
            print(i, sig, len(out))
        if line[1] == "F0":  # fibered
            tri, angle = isosig_to_tri_angle(sig) 
            if tri.homology().rank() == 1:  # b_1 = 1
                dil, euler = dilatation_betti_one_fibred(tri, angle)
                out.append( [sig, str(dil), str(euler), str(log(dil)*euler)] ) 
        if i % (10*report) == 0 and len(out) > 0: 
            write_data_file(out, out_filename)
    write_data_file(out, out_filename)
    
