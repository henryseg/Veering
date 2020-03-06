#
#  veering_triangulation.py
# 

from transverse_taut import is_transverse_taut
from veering import is_veering

class veering_triangulation():
    """Container class for a triangulation with transverse veering structure, possibly with hyperbolic shapes"""
    def __init__(self, tri, angle, tet_shapes = None):
        self.tri = tri
        self.angle = angle
        self.coorientations = is_transverse_taut(tri, angle, return_type = 'tet_vert_coorientations')
        assert self.coorientations != False
        self.veering_colours = is_veering(tri, angle, return_type = 'veering_colours')
        assert self.veering_colours != False
        self.tet_shapes = tet_shapes