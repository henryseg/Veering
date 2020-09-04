#
# shapes.py
# 

# From a given collection of isosigs, build the snappy shapes and put
# them in a pickled dictionary.

import snappy

from file_io import output_to_pickle
from taut import isosig_to_tri_angle

def shapes_to_pickle(isosigs, filename, progress = 100):
    shapes = {}
    for i, sig in enumerate(isosigs):
        if i % progress == 0: print(i, sig)

        tri, angle = isosig_to_tri_angle(sig)

        N = snappy.Manifold(tri.snapPea()) 
        N_shapes = [complex(shape['rect']) for shape in N.tetrahedra_shapes()]
        
        shapes[sig] = N_shapes

    output_to_pickle(shapes, filename)
    return None
