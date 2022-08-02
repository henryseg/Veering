#
# normal_surfaces.py
#

import regina

from veering.file_io import parse_data_file

from veering.taut import isosig_to_tri_angle

def count_quads(surf):
    count = 0
    for i in range(surf.triangulation().countTetrahedra()):
        for j in range(3):
            count += surf.quads(i, j)
    return count

def count_quad_types(surf):
    count = 0
    for i in range(surf.triangulation().countTetrahedra()):
        for j in range(3):
            if surf.quads(i, j) > 0:
                count += 1
    return count

def analyze_sig(sig):
    # print(sig)
    tri, angle = isosig_to_tri_angle(sig)
    surfs = regina.NormalSurfaces.enumerate(tri, regina.NS_QUAD_CLOSED, regina.NS_FUNDAMENTAL)
    if surfs != None:
        two_quad_type_surfs = []
        for i in range(surfs.size()):
            surf = surfs.surface(i)
            if count_quad_types(surf) <= 2:
                two_quad_type_surfs.append(surf)
                # print(count_quads(surf), sig, surf)
        if len(two_quad_type_surfs) > 2:
            print(sig)
            for surf in two_quad_type_surfs:
                print(surf)
