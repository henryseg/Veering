#
# sage_tools.py
#

# The code below is copied and modified (with permission) from
# https://github.com/3-manifolds/SnapPy/blob/master/python/snap/nsagetools.py

from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import Matrix


def join_lists(list_of_lists):
    out = []
    for l in list_of_lists:
        out = out + l
    return out


def uniform_exponents(poly):
    return [list(e) if hasattr(e, "__getitem__") else (e,) for e in poly.exponents()]


def monomial_multiplier(elts, ZH):
    if all(elt == 0 for elt in elts):
        # Zero (unlike other constants) has valuation -\infty.  This
        # can show up as an issue when computing the big polynomial.
        # For an example, compute ET for
        # "kLLLMPPkcdgfehijjijhshassqhdqr_1222011022"
        return ZH(1)
    elts = [ZH(elt) for elt in elts]
    A =  Matrix(ZZ, join_lists([ uniform_exponents(p) for p in elts]))
    min_exp = tuple( [min(row) for row in A.transpose()] )
    return ZH( {min_exp:1} )


def laurent_to_poly(elt, P):
    if type(elt) is int:
        return P(elt)
    return P( elt.dict() )


def matrix_laurent_to_poly(M, ZH, P):
    # convert to polynomials after shifting rows
    muls = [ monomial_multiplier(row, ZH) for row in M ]
    return Matrix( [ [ laurent_to_poly(p / mul, P) for p in row ] for row, mul in zip(M, muls) ] )


def normalise_poly(poly, ZH, P):
    if poly == 0:
        return poly
    mul = monomial_multiplier([poly], ZH)
    poly = laurent_to_poly(poly / mul, P)
    # if poly.coefficients()[-1] < 0:  # I don't trust this.
    # if str(poly)[0] == "-":  # This is a hack, of course.
    if poly.lt().coefficients()[0] < 0:  # lt = leading term
        poly = -poly
    return poly


### end of copied/modified code
