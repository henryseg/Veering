#
# fibered.py
#
# Goal: Best effort to detect fibered manifolds

import regina
import snappy

from .taut import taut_regina_angle_struct_to_taut_struct
from .taut_polytope import is_layered
from .transverse_taut import is_transverse_taut


def get_many_sigs(snappy_name, tries):
    """
    Given a snappy name (isosig, etc), tries to find many sigs 
    for the manifold.
    """
    M = snappy.Manifold(snappy_name)
    sigs = []
    for i in range(tries):
        sig = M.triangulation_isosig(decorated = False)
        if sig not in sigs: sigs.append(sig)
        M.randomize()
    return sigs
        

def is_fibered(snappy_name, tries = 1000, with_data = False):
    """
    Given an snappy_name, tries to find a layered transverse taut
    triangulation for the manifold via random search.  Thus an answer
    of "False" is not rigorous.
    """
    sigs = get_many_sigs(snappy_name, tries)
    for sig in sigs:
        tri = regina.Triangulation3(sig)
        tri.orient()
        angle_strs = regina.AngleStructures(tri, True)
        angle_strs = [taut_regina_angle_struct_to_taut_struct(angle_str) for angle_str in angle_strs]
        for angle in angle_strs:
            # we are not interested in semi-fibered manifolds so:
            if is_transverse_taut(tri, angle) and is_layered(tri, angle):
                if with_data:
                    return (True, tri.isoSig(), angle)
                return True
    return False
