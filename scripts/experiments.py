#
# experiments.py
#

import snappy

from veering.snappy_tools import tet_norm, get_slopes_from_peripherals
from veering.taut import isosig_to_tri_angle
from veering.veering_tri import veering_triangulation
from veering.veering_drill_midsurface_bdy import drill_midsurface_bdy

def get_sub_shapes(census):

    sub_tog_tb = []
    sub_tog_sides = []
    sub_fan_tb = []
    sub_fan_sides = []

    for sig in census:
        tri, angle = isosig_to_tri_angle(sig)
        vt = veering_triangulation(tri, angle)
        tet_types = vt.tet_types
        # num_tog = tet_types.count("toggle") 

        subtet_start_index = 0
        subtet_indices = []
        for i in range(tri.countTetrahedra()):
            tet = tri.tetrahedron(i)
            if tet_types[i] == 'toggle':
                num_to_add = 8
            else:
                num_to_add = 4
            subtet_indices.append(range(subtet_start_index, subtet_start_index + num_to_add))
            subtet_start_index = subtet_start_index + num_to_add

        T, _ = drill_midsurface_bdy(tri, angle)
        M = snappy.Manifold(T)
        new_tet_shapes = M.tetrahedra_shapes("rect")
        new_tet_shapes = [tet_norm(z) for z in new_tet_shapes]
        for z in new_tet_shapes:
            if z.imag() > 1:
                print z, sig
        for i, tet_type in enumerate(tet_types):
            if tet_type == "toggle":
               for k, j in enumerate(subtet_indices[i]):
                   if k < 4:
                       sub_tog_tb.append(new_tet_shapes[j])
                   else:
                       sub_tog_sides.append(new_tet_shapes[j])
            else:
               for k, j in enumerate(subtet_indices[i]):
                   if k < 2:
                       sub_fan_tb.append(new_tet_shapes[j])
                   else:
                       sub_fan_sides.append(new_tet_shapes[j])

    A = list_plot(sub_tog_tb, aspect_ratio = 1, color = "red")
    B = list_plot(sub_tog_sides, aspect_ratio = 1, color = "blue")
    C = list_plot(sub_fan_tb, aspect_ratio = 1, color = "green")
    D = list_plot(sub_fan_sides, aspect_ratio = 1, color = "purple")

    return A, B, C, D

def get_merid_lengths(census, report = 100):

    old_cusp_diffs = []
    new_cusp_shapes = []
    holonomies = []

    for k, sig in enumerate(census):
        if (k % report) == 0:
            print k, sig
        tri, angle = isosig_to_tri_angle(sig)
        T, merids = drill_midsurface_bdy(tri, angle)
        M = snappy.Manifold(T)
        M.set_peripheral_curves("shortest")
        num_cusps = M.num_cusps()
        curr_slopes = get_slopes_from_peripherals(M, merids)
        assert len(curr_slopes) == num_cusps

        temp_shapes = []
        for i, slope in enumerate(curr_slopes):
            if slope == (0, 0):
                temp_shapes.append(M.cusp_info(i).shape)
            else:
                new_cusp_shapes.append(M.cusp_info(i).shape)

        for i, slope in enumerate(curr_slopes):
            if slope != (0, 0):
                C = M.cusp_neighborhood()
                # maximise the ith cusp
                try:
                    C.set_displacement(C.stopping_displacement(i), i)
                except:
                    continue
                # measure the length of merids[i]
                p, q = slope
                m, l = C.translations(i)
                holonomies.append(p*m + q*l)
                
        # fill and compute diffs
        temps_two = []
        M.dehn_fill(curr_slopes)
        for i, slope in enumerate(curr_slopes):
            if slope == (0, 0):
                diff = temp_shapes.pop(0)/M.cusp_info(i).shape
                temps_two.append(diff)

        old_cusp_diffs.extend(temps_two)

    A = list_plot(old_cusp_diffs, aspect_ratio = 1, color = "red")
    B = list_plot(new_cusp_shapes, aspect_ratio = 1, color = "blue")
    C = list_plot(holonomies, aspect_ratio = 1, color = "green")

    return A, B, C
