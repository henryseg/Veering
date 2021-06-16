#
# snappy_util.py
# 

# Functions that call upon snappy.

import snappy

from file_io import output_to_pickle
from taut import isosig_to_tri_angle


# Shapes

def tet_norm(z):
    if abs(z) < 1 and abs(1-z) < 1:   
        return z        
    elif z.real() > 0.5:
        return (z - 1)/z
    else:   
        return 1/(1 - z)
    
# From a given collection of isosigs, build the snappy shapes and put
# them in a pickled dictionary.

def shapes_to_pickle(isosigs, filename, progress = 100):
    shapes = {}
    for i, sig in enumerate(isosigs):
        if i % progress == 0: print((i, sig))

        tri, angle = isosig_to_tri_angle(sig)

        N = snappy.Manifold(tri)
        N_shapes = [complex(shape['rect']) for shape in N.tetrahedra_shapes()]
        # N = snappy.ManifoldHP(tri) 
        # N_shapes = [shape['rect'] for shape in N.tetrahedra_shapes()]
        
        shapes[sig] = N_shapes

    output_to_pickle(shapes, filename)
    return None

# Peripheral intersection numbers

def triangle_sum(x, y):
    x0, x1, x2 = x
    y0, y1, y2 = y
    
    return ( x2 * y1 + x0 * y2 + x1 * y0 ) - \
           ( x1 * y2 + x2 * y0 + x0 * y1 )

def algebraic_intersection(a, b):
    assert len(a) == len(b)
    assert len(a) % 3 == 0

    num_tets = len(a) // 3
    
    # See Neumann-Zagier equation 26. There they use powers of z and 1 - z.
    # We instead use powers of z_0, z_1, and z_2.

    return sum([triangle_sum(a[3*i:3*i + 3], b[3*i:3*i + 3]) for i in range(num_tets)]) // 2

def cusp_slope(m, l, a):
    p, q = (algebraic_intersection(a, l), algebraic_intersection(m, a))
    if p > 0 or (p == 0 and q > 0):
        return (p, q)
    return (-p, -q)

def get_slopes_from_peripherals(M, peripherals):
    # given a snappy manifold and list of peripheral curves, return a
    # list "slopes", where slopes[i] is the slope in cusp i.  If a cusp
    # is not visited then slopes[i] == (0, 0)

    slopes = {}
    n = M.num_cusps()
    frames = M.gluing_equations()[-2*n:]
    merids = frames[::2]
    longs  = frames[1::2]
    for curve in peripherals:
        curr_slopes = [cusp_slope(merids[i], longs[i], curve) for i in range(n)]
        curr_booleans = [curr_slope == (0, 0) for curr_slope in curr_slopes]
        assert curr_booleans.count(False) == 1
        ind = curr_booleans.index(False)
        slope = curr_slopes[ind]
        if ind in slopes:
            assert slopes[ind] == slope
        else:
            slopes[ind] = slope
    
    for i in range(n):
        if i not in slopes:
            slopes[i] = (0, 0) 

    return [slopes[i] for i in range(n)]
