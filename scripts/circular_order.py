#
# circular_order.py
#

# def enforce_fraction(a):   ### other implementations of rational numbers we have tried include QQ and Fraction. They seem slow.
#     if type(a) != tuple:
#         return (a, 1)
#     else:
#         return a

def are_anticlockwise(a, b, c): ### given indices in this order, are they anticlockwise on the circle?
    return (b - a)*(c - b)*(a - c) < 0  ### either one or two are negative
    ### If you give it two of the same index, the answer is False

def are_anticlockwise_pairs(a, b, c): ### given indices in this order, are they anticlockwise on the circle?
    # an, ad = enforce_fraction(a)
    # bn, bd = enforce_fraction(b)
    # cn, cd = enforce_fraction(c)
    an, ad = a
    bn, bd = b
    cn, cd = c
    return (bn*ad - an*bd) * (cn*bd - bn*cd) * (an*cd - cn*ad) < 0  ### either one or two are negative
    ### If you give it two of the same index, the answer is False

def are_linking(a1, a2, b1, b2): ### given indices on a circle, do the a's link the b's?
    return (a1 - b1)*(a1 - b2)*(a2 - b1)*(a2 - b2) < 0 ### same as the sign of the cross-ratio

def are_linking_pairs(a1, a2, b1, b2): ### given indices on a circle, do the a's link the b's?
    # a1n, a1d = enforce_fraction(a1) 
    # b1n, b1d = enforce_fraction(b1)
    # a2n, a2d = enforce_fraction(a2)
    # b2n, b2d = enforce_fraction(b2)
    a1n, a1d = a1
    b1n, b1d = b1
    a2n, a2d = a2
    b2n, b2d = b2
    return (a1n*b1d - b1n*a1d)*(a1n*b2d - b2n*a1d)*(a2n*b1d - b1n*a2d)*(a2n*b2d - b2n*a2d) < 0 ### same as the sign of the cross-ratio
