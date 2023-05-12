#
# circular_order.py
#

def are_anticlockwise(a, b, c): ### given indices in this order, are they anticlockwise on the circle?
    return (b - a)*(c - b)*(a - c) < 0  ### either one or two are negative

def are_linking(a1, a2, b1, b2): ### given indices on a circle, do the a's link the b's?
    return (a1 - b1)*(a1 - b2)*(a2 - b1)*(a2 - b2) < 0 ### same as the sign of the cross-ratio
