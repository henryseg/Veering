#
# veering_knots.py
#

# Find all triangulations in the veering census that are
# triangulations of knot complements (in the three-sphere).

import snappy
from veering import veering_census

cen = veering_census()

# filter to get the one-cusped manifolds

def get_one_cusped(cen):
    cen_one = []
    for sig in cen:
        sig_sp = sig.split("_")[0]
        M = snappy.Manifold(sig_sp)
        if M.num_cusps() == 1:
            cen_one.append(sig)
    return cen_one
            
# this is less than a minute
# 59114 manifolds with one cusp.

def get_positives(cen_one):
    pos_uns = []
    for sig in cen_one:
        sig_sp = sig.split("_")[0]
        M = snappy.Manifold(sig_sp)
        i = 0
        while M.solution_type() != "all tetrahedra positively oriented":
            M.randomize()
            i = i + 1
            if (i + 1) % 2000:
                print(sig, "is having a hard time...")
        pos_uns.append((sig, M.triangulation_isosig()))
    return pos_uns

# Also less than a minute!
# 59114 still 

def get_slopes():
    pos_plus_slopes = []
    bad_uns = []
    i = 0
    for sig, pos_sig in pos_uns:
        i = i + 1
        if i % 500 == 0:
            print(i, len(pos_plus_slopes), len(bad_uns))
        M = snappy.Manifold(pos_sig)
        try:
            slopes = M.short_slopes()[0]
            pos_plus_slopes.append((sig, pos_sig, slopes))
        except:
            N = M.high_precision()
            try:
                slopes = N.short_slopes()[0]
                pos_plus_slopes.append((sig, pos_sig, slopes))
            except:
                bad_uns.append((sig, pos_sig))
    assert bad_uns == 0
    return pos_plus_slopes

# this is a bit slow, but it works.

def get_hom_trivial_slope(pos_plus_slopes):
    pos_plus_special_slope = []
    i = 0
    for sig, pos_sig, slopes in pos_plus_slopes:
        i = i + 1
        if (i + 1) % 500 == 0:
            print(sig, len(pos_plus_special_slope))
        for slope in slopes:
            M = snappy.Manifold(pos_sig)
            M.dehn_fill(slope)
            if str(M.homology()) == "0":
                pos_plus_special_slope.append((sig, pos_sig, slope))
    return pos_plus_special_slope
    
# 1504 of these

def first_sort(pos_plus_special_slope):
    prob_hyp = []
    prob_not_hyp = []
    for sig, pos_sig, slope in pos_plus_special_slope:
        M = snappy.Manifold(pos_sig)
        M.dehn_fill(slope, 0)
        if M.solution_type() == "all tetrahedra positively oriented":
            prob_hyp.append((sig, pos_sig, slope))
        else:
            prob_not_hyp.append((sig, pos_sig, slope))
    return prob_hyp, prob_not_hyp

# 289 and 1215

def deal_with_hyp(prob_hyp):
    def_hyp = []
    for sig, pos_sig, slope in prob_hyp:
        M = snappy.Manifold(pos_sig)
        M.dehn_fill(slope, 0)
        N = M.high_precision()
        if N.verify_hyperbolicity()[0]:
            def_hyp.append((sig, pos_sig, slope))
        else:
            print(sig, pos_sig, slope, "hyp check failed")
    return def_hyp

# This catches all 289

def actual_knots(prob_not_hyp):
    def_knot = []
    prob_not_knot = []
    for sig, pos_sig, slope in prob_not_hyp:
        win = False
        M = snappy.Manifold(pos_sig)
        M.dehn_fill(slope, 0)
        T = M.filled_triangulation()
        for i in range(100):
            if T.num_tetrahedra() == 1:
                G = T.fundamental_group()
                if G.num_generators() == 0:
                    win = True
                    break
                else:
                    print(sig, pos_sig, slope, "is super weird!")
            T.randomize()
        if win:
            def_knot.append((sig, pos_sig, slope))
        else:
            prob_not_knot.append((sig, pos_sig, slope))
    return def_knot, prob_not_knot

# Should have 609 knots and 606 almost certainly not knots.

def get_rid_of_non_knots(prob_not_knot):
    hard_cases = []
    def_not_knot = []
    for sig, pos_sig, slope in prob_not_knot:
        win = False
        M = snappy.Manifold(pos_sig)
        M.dehn_fill(slope, 0)
        for i in range(2, 8):
            if len(M.covers(i)) > 0:
                win = True
                break
        if win:
            def_not_knot.append((sig, pos_sig, slope))
        else:
            hard_cases.append((sig, pos_sig, slope))
    return def_not_knot, hard_cases

def get_rid_of_more_non_knots(hard_cases):
    hard_but_def_hyp = []
    hard_cases_two = []
    for sig, pos_sig, slope in hard_cases:
        M = snappy.Manifold(pos_sig)
        M.dehn_fill(slope, 0)
        win = False
        T = M.filled_triangulation()
        for i in range(200):
            T.randomize()
            N = T.with_hyperbolic_structure()
            if N.solution_type() == "all tetrahedra positively oriented":
                if N.verify_hyperbolicity()[0]:
                    win = True
                    break
        if win:
            hard_but_def_hyp.append((sig, pos_sig, slope))
        else:
            hard_cases_two.append((sig, pos_sig, slope))
    return hard_but_def_hyp, hard_cases_two

# 59 and 55 - need to deal with the latter.
