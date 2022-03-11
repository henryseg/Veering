#
# normal_surfaces.py
#

from file_io import parse_data_file
from taut import isosig_to_tri_angle
import regina


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
                

def main(num_to_check = 10):
    lines = parse_data_file('Data/veering_census.txt')
    for line in lines[:num_to_check]:
        sig = line.strip()
        analyze_sig(sig)
