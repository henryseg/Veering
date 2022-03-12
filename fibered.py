#
# fibered.py
#

import regina
import snappy

from taut import taut_regina_angle_struct_to_taut_struct
from taut_polytope import is_layered
from transverse_taut import is_transverse_taut

# Best effort attempt to detect fibered manifolds

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
    triangulation for the manifold via random search.
    """
    sigs = get_many_sigs(snappy_name, tries)
    for sig in sigs:
        T = regina.Triangulation3(sig)
        T.orient()
        angle_strs = list(regina.AngleStructures.enumerate(T, True))
        angle_strs = [taut_regina_angle_struct_to_taut_struct(angle_str) for angle_str in angle_strs]
        for angle in angle_strs:
            # we are not interested in semi-fibered manifolds so:
            if is_transverse_taut(T, angle) and is_layered(T, angle):
                if with_data:
                    return (True, T.isoSig(), angle)
                return True

    return False
