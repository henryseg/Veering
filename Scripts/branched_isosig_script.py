from taut import isosig_to_tri_angle, isosig_from_tri_angle, apply_isom_to_angle_struct_list, is_taut
from branched_surface import upper_branched_surface, isosig_from_tri_angle_branch, isosig_to_tri_angle_branch, apply_isom_to_branched_surface, is_branched
from flow_cycles import flow_cycle_to_triangle_loop
from drill import drill 
import regina

def main():
    tri, angle = isosig_to_tri_angle('cPcbbbdxm_10')
    branch = upper_branched_surface(tri, angle)
    print(branch)
    tl = flow_cycle_to_triangle_loop(tri, branch, [(0,2)])
    drill(tri, tl, angle = angle, branch = branch)
    print(branch)

    # sig_taut, isom1 = isosig_from_tri_angle(tri, angle, return_isom = True)
    # tri2, angle2, isom2 = isosig_to_tri_angle(sig_taut, return_isom = True)

    # combined_isom = isom2 * isom1
    # tri3 = combined_isom.apply(tri)

    # print(isom1)
    # print(isom2)
    # print(combined_isom)

    # assert tri.isIsomorphicTo(tri3)

    # sig = isosig_from_tri_angle_branch(tri, angle, branch)
    # # print(sig, angle, branch)

    # # tri3, angle3, branch3 = isosig_to_tri_angle_branch(sig)


    # tri, angle = isosig_to_tri_angle('cPcbbbdxm_10')
    # branch = upper_branched_surface(tri, angle)

    # print('original tri', tri, tri.countTetrahedra())
    # print('original angle, branch', angle, branch)
    # assert is_taut(tri, angle)
    # assert is_branched(tri, branch)

    # all_isoms = tri.findAllIsomorphisms(tri)
    # for isom in all_isoms:
    #     print(isom)
    #     new_angle = apply_isom_to_angle_struct_list(angle, isom)
    #     new_branch = apply_isom_to_branched_surface(branch, isom)
    #     new_tri = isom.apply(tri) # not in place
    #     print('new_angle, new_branch', new_angle, new_branch)
    #     assert is_taut(tri, new_angle)
    #     assert is_taut(new_tri, new_angle)
    #     assert is_branched(tri, new_branch)
    #     assert is_branched(new_tri, new_branch)

    isosig, isosig_isom = tri.isoSigDetail()  # isom is the mapping from the original triangulation to the isosig triangulation
    isosig_tri = regina.Triangulation3.fromIsoSig(isosig)
    isosig_branch = apply_isom_to_branched_surface(branch, isosig_isom)
    assert is_branched(isosig_tri, isosig_branch)