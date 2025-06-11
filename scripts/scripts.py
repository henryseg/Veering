#
# scripts.py
#

# Examples of usage for some of the modules - most need cleaning before they will work. 

from sage.functions.log import log

from veering.file_io import parse_data_file, write_data_file
from veering.taut import liberal, isosig_to_tri_angle
from veering.snappy_tools import shapes
from veering.taut_polynomial import taut_polynomial_via_fox_calculus
from veering.taut_polytope import min_neg_euler_carried

def boundary_triangulation_script(veering_isosig):
    from boundary_triangulation import draw_triangulation_boundary_from_veering_isosig
    # Set 'ct_depth': <some non-negative integer> to do cannon-thurston
    args = {'draw_boundary_triangulation':True, 'draw_triangles_near_poles': False, 'ct_depth':-1, 'ct_epsilon':0.03, 'global_drawing_scale': 12, 'delta': 0.2, 'ladder_width': 10.0, 'ladder_height': 20.0, 'draw_labels': True}
    ### for standard ladder picture, set 'draw_triangles_near_poles' = False. Set True for CT pictures

    # num_to_draw = 87 ## up to 6 tet
    # # num_to_draw = 5699 ## up to 12 tet  
    # args['style'] = 'ladders'
    args['style'] = 'geometric'

    # args = {'draw_boundary_triangulation':False, 'draw_triangles_near_poles': False, 'ct_depth':-1, 'ct_epsilon':0.03, 'global_drawing_scale': 4, 'delta': 0.2, 'ladder_width': 10.0, 'ladder_height': 20.0, 'draw_labels': True}
    # args['style'] = 'geometric'
    # num_to_draw = 1000
    # draw_triangulations_from_veering_isosigs_file('veering_census.txt', 'Images/Boundary_triangulations/Geometric', args = args, num_to_draw = num_to_draw)

    name = veering_isosig

    # name = 'cPcbbbdxm_10'
    # name = 'cPcbbbiht_12'
    # name = 'eLMkbcddddedde_2100'
    # name = 'eLAkbccddhhsqs_1220'
    # name = 'fLAMcaccdeejsnaxk_20010'
    # # # # name = 'fLLQcbecdeepuwsua_20102'
    # # # # name = 'fLLQcbeddeehhbghh_01110'
    # # # # name = 'jLAwwAQbcbdfghihihhwhnaaxrn_211211021'
    # name = 'eLAkaccddjsnak_2001'
    # name = 'eLAkbbcdddhwqj_2102'
    # name = 'dLQacccjsnk_200' 
    # name = 'dLQbccchhsj_122'
    # name = 'iLLLAQccdffgfhhhqgdatgqdm_21012210'
    # name = 'gLvQQadfedefjaaajkk_200211'

    # name = 'gLPLQbcdfeffhrraarw_120011'
    # name = 'qvvLPQMvQLQkfgfhhgfknlmoppopohahhaaahaqqaqqaaa_1222211100222200'

    # shapes_data = read_from_pickle('veering_shapes_up_to_ten_tetrahedra.pkl')
    # args['style'] = 'geometric'
    # args['tet_shapes'] = shapes_data[name]
    # # args['tet_shapes'] = None

    draw_triangulation_boundary_from_veering_isosig(name, output_filename = name + '_' + args['style'] + '.pdf', args = args) 

    # B = generate_boundary_triangulation(name)
    # print(B.ladder_counts())

def branched_pachner_script():
    from branched_pachner import twoThreeMove, threeTwoMove
    # for i in range(4): 
    #     print(i)
    #     # tri, angle = isosig_to_tri_angle('cPcbbbiht_12') 
    # sigs = ['dLQacccjsnk_200', 'dLQbccchhfo_122','dLQbccchhsj_122']
    # for sig in sigs:
    #     print(sig)
    #     tri, angle = isosig_to_tri_angle(sig)   
    #     for branch in all_branched_surfaces(tri):
    #         print(lex_smallest_branched_surface(tri, branch))

    sig = 'dLQacccjsnk_200'
    for i in range(6):
        # print(i)
        tri, angle = isosig_to_tri_angle(sig)  
        tri_original = regina.Triangulation3(tri) #copy

        out = twoThreeMove(tri, [4,11,0], i, return_edge = True)
        if out != False:
            tri, possible_branches, edge_num = out
            # print('possible_branches', possible_branches)
            # print('all branches', all_branched_surfaces(tri)) 
            tri, branch = threeTwoMove(tri, possible_branches[0], edge_num)
            all_isoms = tri.findAllIsomorphisms(tri_original)
            all_branches = [apply_isom_to_branched_surface(branch, isom) for isom in all_isoms]
            assert [4,11,0] in all_branches


def draw_continent_script():
    from draw_continent import draw_continent

    # draw_args = {'draw_boundary_triangulation':True, 'draw_labels': False, 'only_draw_ladderpoles': True, 'ct_lw': 0.002, 'global_drawing_scale': 4, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps
    # draw_args = {'draw_boundary_triangulation':True, 'draw_labels': True, 'only_draw_ladderpoles': False, 'ct_lw': 0.02, 'global_drawing_scale': 4, 'style': 'geometric', 'draw_triangles_near_poles': False, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps
    draw_args = {'draw_boundary_triangulation':False, 'draw_labels': False, 'only_draw_ladderpoles': True, 'ct_lw': 0.02, 'global_drawing_scale': 4, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps

    # max_num_tetrahedra = 40
    # max_num_tetrahedra = 5000
    # max_num_tetrahedra = 50000
    # max_num_tetrahedra = 100000
    max_num_tetrahedra = 400000
    # max_num_tetrahedra = 2000000
    # max_num_tetrahedra = 50000000
    # max_length = 0.4
    # max_length = 0.3
    # max_length = 0.2
    # max_length = 0.15
    # max_length = 0.14
    max_length = 0.1
    # max_length = 0.09
    # max_length = 0.07
    # max_length = 0.06
    # max_length = 0.05
    # max_length = 0.02
    # max_length = 0.015
    # max_length = 0.01
    # max_length = 0.005

    # draw_args['ct_lw'] = 0.2 * max_length 

    # build_type = 'build_naive'
    # build_type = 'build_on_coast'
    # build_type = 'build_make_long_descendant_edges_internal'
    # build_type = 'build_explore_prongs'
    build_type = 'build_long_and_mid'

    veering_isosig = 'cPcbbbiht_12'
    # # # # veering_isosig = 'cPcbbbdxm_10'
    # veering_isosig = 'dLQacccjsnk_200'
    # veering_isosig = 'eLMkbcddddedde_2100'
    # veering_isosig = 'eLAkaccddjsnak_2001'
    # veering_isosig = 'gLAMPbbcdeffdhwqqqj_210202'
    # veering_isosig = 'fLLQcbecdeehhnkei_12001'
    # veering_isosig = 'gLLAQbecdfffhhnkqnc_120012'
    # veering_isosig = 'iLLLAQccdffgfhhhqgdatgqdm_21012210' ## no symmetry - helps us spot errors
    # veering_isosig = 'iLLPwQcccdfehghhhggaahhbg_20102211'
    # # veering_isosig = 'jLAwwAQbcbdfghihihhwhnaaxrn_211211021' ## first non geometric
    # # veering_isosig = 'nLLwMLPMMkbeefeihjkjlmlmhhaaaektxnaqrs_0111000011220'  ### quite big negative shape
    # veering_isosig = 'qLvPvvMQQLQkccgkgjkmlknpooppoqjaajqqhhqqaqxhhh_0222110112222211'
    # veering_isosig = 'fLLQcbeddeehhnkhh_21112'
    # veering_isosig = 'eLAkbbcdddhwqj_2102'

    # # # veering_isosig = 'mLvLLLQQQbegikhjiilkllhiardrnnkxeif_120000112222'
    # veering_isosig = 'mLvLLMMQQcehfhjlklkjlktilbbjumhtfai_011220220111'
    # veering_isosig = 'mLLvLQLQQbeffjglhlkjklxxxjsfqjhhoqo_102210101022'
    
    # veering_isosig ='kLLLAPPkcdgfehhjijjhfhaqiphffj_2010222001'
    # veering_isosig = 'jLLLLQQbegfhihihihhqakkkkoo_120011211'
    # veering_isosig = 'dLQbccchhsj_122'

    # veering_isosig = 'mvLALLMQQecfgjkkjiljllccaxvvwkfekix_100001122112'
    # veering_isosig = 'mvLALPMPQecfggjgikllklccaxxvcfaqdmo_100001122100'
    # veering_isosig = 'nLLLLzAPQkbefgjkjimlllmmxxqhxqubxtivwb_1022101110220'

    
    ### for jigsaws: 
    # draw_args['draw_CT_curve'] = False
    # draw_args['draw_lightning_curve'] = True
    # draw_args['draw_jordan_curve'] = False
    # draw_args['draw_dividers'] = True
    # draw_args['draw_landscapes'] = False
    # draw_args['draw_box_for_cohom_frac'] = False
    # draw_args['draw_alignment_dots'] = False
    # draw_args['draw_desired_vertices'] = False
    # draw_args['expand_fund_dom'] = True


    draw_args['draw_CT_curve'] = True
    draw_args['draw_lightning_curve'] = False
    draw_args['draw_jordan_curve'] = False
    draw_args['draw_dividers'] = False
    draw_args['draw_landscapes'] = False
    draw_args['draw_box_for_cohom_frac'] = False
    draw_args['draw_alignment_dots'] = False
    draw_args['draw_desired_vertices'] = False
    draw_args['expand_fund_dom'] = True  ### needed for jordan curve?

    draw_type = ''
    if draw_args['draw_CT_curve']:
        draw_type = 'CT_spectrum'
    elif draw_args['draw_jordan_curve']:
        draw_type = 'CT_jordan'

    filename = 'Images/Cannon-Thurston/' + veering_isosig + '_' + str(max_length) + '_' + build_type + '_' + draw_type + '.pdf'
    # filename = 'Images/Jigsaw/' + veering_isosig + '_' + str(max_num_tetrahedra) + '_' + str(max_length) + '_' + build_type + '.pdf'


    draw_continent( veering_isosig, max_num_tetrahedra = max_num_tetrahedra, max_length = max_length, output_filename = filename, draw_args = draw_args, build_type = build_type )
    


    ## draw many:
    # start_num = 500
    # end_num = 502
    # draw_cannon_thurston_from_veering_isosigs_file('veering_census.txt', 'Images/Jordan_curve', max_num_tetrahedra = max_num_tetrahedra, max_length = max_length, interval_to_draw = (start_num, end_num), draw_args = draw_args, build_type = build_type)
    
    ### jigsaws

    # draw_jigsaw_from_veering_isosigs_file('sigs_for_jigs_no_symm.txt', 'Images/Jigsaw', jigsaw_data_out_filename = "jigsaw_data_no_symm.pkl", max_num_tetrahedra = max_num_tetrahedra, max_length = max_length, draw_args = draw_args, interval_to_draw = (start_num, end_num))  # all up through n's is 876. The 281th has trouble developing

    # veering_isosigs_list = ['kLLLAPPkcdgfehhjijjhfhaqiphffj_2010222001', 'mvLALLMQQecfgjkkjiljllccaxvvwkfekix_100001122112', 'mvLALPMPQecfggjgikllklccaxxvcfaqdmo_100001122100']
    # draw_jigsaw_from_veering_isosigs_list(veering_isosigs_list, 'Images/Jigsaw', jigsaw_data_out_filename = "jigsaw_data.pkl", max_num_tetrahedra = 2000000, max_length = max_length, draw_args = draw_args)


def draw_continent_from_isosig(veering_isosig, max_length = 0.1, max_num_tetrahedra = 500000, use_algebraic_numbers = True):
    from draw_continent import draw_continent
    draw_args = {'draw_boundary_triangulation':False, 'draw_labels': False, 'only_draw_ladderpoles': True, 'ct_lw': 0.02, 'global_drawing_scale': 12, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps
    # draw_args = {'draw_boundary_triangulation':False, 'draw_labels': False, 'only_draw_ladderpoles': False, 'ct_lw': 0.02, 'global_drawing_scale': 4, 'style': 'geometric', 'draw_triangles_near_poles': True, 'ct_depth': -1} #ct_depth is the old way to try to build ct maps
    # build_type = 'build_make_long_descendant_edges_internal'
    build_type = 'build_long_and_mid'
    # build_type = 'build_naive'

    # max_num_tetrahedra = 40
    # max_num_tetrahedra = 5000
    # max_num_tetrahedra = 50000
    # max_num_tetrahedra = 100000
    # max_num_tetrahedra = 400000
    # max_num_tetrahedra = 2000000
    # max_num_tetrahedra = 50000000
    # max_length = 0.4
    # max_length = 0.3
    # max_length = 0.2
    # max_length = 0.15
    # max_length = 0.14
    # max_length = 0.1
    # max_length = 0.09
    # max_length = 0.07
    # max_length = 0.06
    # max_length = 0.05
    # max_length = 0.02
    # max_length = 0.015
    # max_length = 0.01
    # max_length = 0.005

    # draw_args['ct_lw'] = 0.2 * max_length 
    draw_args['draw_CT_curve'] = False
    draw_args['draw_lightning_curve'] = False
    draw_args['draw_jordan_curve'] = True
    draw_args['draw_dividers'] = False
    draw_args['draw_landscapes'] = False
    draw_args['draw_box_for_cohom_frac'] = False
    draw_args['draw_alignment_dots'] = False
    draw_args['draw_desired_vertices'] = False
    draw_args['expand_fund_dom'] = True  ### needed for jordan curve?

    draw_type = ''
    if draw_args['draw_CT_curve']:
        draw_type = 'CT_spectrum'
    elif draw_args['draw_jordan_curve']:
        draw_type = 'CT_jordan'

    string_max_length = str(max_length)
    string_max_length = string_max_length.strip('0')
    if use_algebraic_numbers:
        alg_status = 'algeb'
    else:
        alg_status = 'float'
    filename_base = veering_isosig + '_' + build_type + '_' + draw_type + '_' + string_max_length + '_' + alg_status 
    filename = 'Images/Cannon-Thurston/' + filename_base + '.pdf'
    load_continents_filename = None  ## don't load anything
    # load_continents_filename = "data/gLLAQbecdfffhhnkqnc_120012_build_long_and_mid_CT_spectrum_.5_float.pkl"
    # save_continents_filename = 'data/' + filename_base + '.pkl'
    save_continents_filename = None  ## don't save anything (pickling doesn't seem to work for large continents anyway)

    draw_continent( veering_isosig, max_num_tetrahedra = max_num_tetrahedra, max_length = max_length, use_algebraic_numbers = use_algebraic_numbers, output_filename = filename, draw_args = draw_args, build_type = build_type, load_continents_filename = load_continents_filename, expand_continents = True, save_continents_filename = save_continents_filename )
    
def draw_census_continents(census_cap = 10):
    census = parse_data_file('../veering/data/veering_census.txt')
    for sig in census[:census_cap]:
        draw_continent_from_isosig(sig)

def draw_fund_dom_continent_circle_script():
    from draw_continent_circle import draw_continent_circle, complete_tetrahedron_rectangles
    from build_continent import get_fund_domain_tetrahedra, make_continent_naive, make_continent_fund_dom
    # veering_isosig = 'cPcbbbdxm_10' 

    # veering_isosig = 'cPcbbbiht_12'

    # veering_isosig = 'dLQacccjsnk_200'

    # veering_isosig = 'eLAkaccddjsnak_2001'

    # veering_isosig = 'dLQacccjsnk_200'
   
    # veering_isosig = 'eLAkaccddjsnak_2001'

    # veering_isosig = 'eLMkbcddddedde_2100'
    # veering_isosig = 'fLLQcbecdeepuwsua_20102'
    # veering_isosig = 'gLLAQbecdfffhhnkqnc_120012'
    # veering_isosig = 'ivLLQQccdgffhghhvvaaaaavv_01111220'
    veering_isosig = 'ivvPQQcfhghgfghfaaaaaaaaa_01122000'

    con, continent_fund_dom_tets = make_continent_fund_dom(veering_isosig, maintain_coast = True, maintain_coastal_edges = True)

    print('con size', len(con.tetrahedra), con.num_tetrahedra)


    # edge_rectangles_to_draw = interval.edges

    # for e in edge_rectangles_to_draw:
        # e.ensure_continent_contains_rectangle()

    tetrahedra_to_draw = []
    tetrahedron_rectangles_to_shade = []
    # tetrahedra_to_draw.extend(continent_fund_dom_tets[:1])  ## only 0
    # tetrahedra_to_draw.extend(continent_fund_dom_tets)  ## only 0
    # complete_tetrahedron_rectangles(con, tetrahedra_to_draw)  ### broken?

    name = veering_isosig + '_fund_dom'

    edge_rectangles_to_draw = []
    draw_continent_circle(con, name = name, draw_labels = True,
        # draw_upper_landscape = True, draw_lower_landscape = True, 
        # draw_upper_landscape = False, draw_lower_landscape = False, 
        draw_coastal_edges = False, draw_all_edges = True,
        draw_cusp_leaves = False,
        shade_triangles = False, 
        draw_fund_domain = False, fund_dom_tets = continent_fund_dom_tets,
        draw_fund_domain_edges = False, 
        #draw_tetrahedron_rectangles = tetrahedra_to_draw,
        edge_rectangles_to_draw = edge_rectangles_to_draw,
        tetrahedron_rectangles_to_draw = tetrahedra_to_draw,
        tetrahedron_rectangles_to_shade = tetrahedron_rectangles_to_shade,
        triangles_to_draw = [],
        tetrahedra_to_draw = tetrahedra_to_draw,
        draw_edges_for_edge_rectangles = True,
        leaves_to_draw = [],
        intervals_to_draw = [],
        edge_thickness = 0.02,
        leaf_thickness = 0.02,
        transparency = 0.8,
        text_size = -2)


def draw_drill_continent_circle_script():
    from draw_continent_circle import draw_continent_circle, complete_tetrahedron_rectangles
    from build_continent import make_continent_drill_flow_cycle, get_fund_domain_tetrahedra, make_continent_naive, make_continent_fund_dom
    from draw_drilled_tetrahedra import draw_drilled_tetrahedra
    # veering_isosig = 'cPcbbbdxm_10' 
    # flow_cycle = [(0, 2)]

    # veering_isosig = 'cPcbbbiht_12'
    # # flow_cycle = [(0, 0)]  ### boundary parallel
    # flow_cycle = [(0, 4), (1,2)]
    # flow_cycle = [(0, 0), (0, 5)] 
    # flow_cycle = [(1, 0), (1, 5)]  ### blue edge rectangle has two punctures, red has none.

    # veering_isosig = 'dLQacccjsnk_200'
    # flow_cycle = [(1, 4)]
    # flow_cycle = [(1, 0), (2, 5)]
    # flow_cycle = [(0, 1), (2, 5), (1, 1)]  ## sideways and up

    # veering_isosig = 'eLAkaccddjsnak_2001'

    # veering_isosig = 'dLQacccjsnk_200'
   
    # veering_isosig = 'eLAkaccddjsnak_2001'
    # flow_cycle = [(0, 1), (2, 2), (3, 2)]

    veering_isosig = 'gLLAQbecdfffhhnkqnc_120012'
    # flow_cycle = [(0, 0)]
    # flow_cycle = [(0, 4), (4, 0)]
    flow_cycle = [(0, 4), (4, 5), (2, 4), (1, 2), (5, 1)]  
    # flow_cycle = [(0, 4), (4, 5), (2, 4), (1, 2), (5, 1), (0, 4), (4, 0)]  ### non simple cycle example
    # flow_cycle = [(0, 4), (4, 5), (2, 4), (1, 2), (5, 1), (0, 4), (4, 0), (0, 4), (4, 5), (2, 4), (1, 2), (5, 1)]
    
    # flow cycle which is not boundary parallel but completely in a blue region 
    # veering_isosig = 'eLAkbbcdddhwqj_2102'
    # flow_cycle = [(3, 2)]   ### this is isotopic to an LMT "AB-cycle", 
    ### so there will be distinct flow lines that fellow-travel, and we cannot separate them
    # flow_cycle = [(2, 5)] ### has same problem of being isotopic to a pi_1 translate
    # flow_cycle = [(0, 3), (1, 1)]  ### goes straight up but does not stay in one colour so we can separate

    # short cycle in a large triangulation
    # veering_isosig = 'qLLALvQLwQMkcecdhkjilmmlnoppphqrwqqaqxhhxofhxq_2010222010122221'
    # flow_cycle = [(1, 0)]
    # # flow_cycle = [(3, 5)]
    # # flow_cycle = [(4, 2)]
    # flow_cycle = [(7, 5)]
    # flow_cycle = [(9, 0)]

    # from veering.branched_surface import upper_branched_surface
    # from veering.flow_cycles import find_flow_cycles
    # tri, angle = isosig_to_tri_angle(veering_isosig)
    # branch = upper_branched_surface(tri, angle)
    # loops = find_flow_cycles(tri, branch)
    # print(loops)

    
    # flow_cycle = [(1, 4)]

    # # # for num_steps in range(10):
    
    print(veering_isosig, flow_cycle)

    out = make_continent_drill_flow_cycle(veering_isosig, flow_cycle, verbose = 0, maintain_coast = True, maintain_coastal_edges = True)
    if out == None:
        print('flow cycle is boundary parallel')
        return None
    con, tetrahedra_cusp_orders, tetrahedra_chunks, intervals_inside_tet_rectangles, all_intervals, continent_fund_dom_tets, leaves_to_draw, triangles_to_draw, found_parallel = out

    # print('interval tet indices', [t.index for t in interval.tetrahedra])
    print('con size', len(con.tetrahedra), con.num_tetrahedra)


    # t1 = con.tetrahedra[-1]
    # t2 = con.tetrahedra[-11]
    # path = t1.face_num_path_to_other_tet(t2)
    # print(path)
    # print(t1.follow_face_num_path(path) == t2)

    # edge_rectangles_to_draw = interval.edges

    # for e in edge_rectangles_to_draw:
        # e.ensure_continent_contains_rectangle()

    tetrahedra_to_draw = []
    tetrahedron_rectangles_to_shade = []
    # tetrahedra_to_draw.extend(continent_fund_dom_tets[:1])  ## only 0
    tetrahedra_to_draw.extend(continent_fund_dom_tets)  ## only 0
    # complete_tetrahedron_rectangles(con, tetrahedra_to_draw) ### broken?
    # for interval in all_intervals:
    #     tetrahedra_to_draw.extend([interval.tetrahedra[0], interval.tetrahedra[-1]])
    #     tetrahedron_rectangles_to_shade.extend([interval.tetrahedra[0], interval.tetrahedra[-1]])
    intervals_to_draw = all_intervals
    # intervals_to_draw = intervals_inside_tet_rectangles[0]
    for interval in intervals_to_draw:
        g, p = interval.crossing_leaves() ### this function can grow the continent - do so before we start drawing
        ### we should have flags for "do not change the continent, assert False if you need to" for drawing routines.


    # for interval in intervals_to_draw:  ### important to do this before we start drawing!
    #     for e in interval.tetrahedra[0].equatorial_edges:
    #         e.ensure_continent_contains_rectangle()
    #     for e in interval.tetrahedra[-1].equatorial_edges:
    #         e.ensure_continent_contains_rectangle()
    # for t in tetrahedra_to_draw:
    #     for e in t.equatorial_edges:
    #         e.ensure_continent_contains_rectangle()

    name = veering_isosig + '' + str(flow_cycle) + '_enough_interval'

    edge_rectangles_to_draw = []
    draw_continent_circle(con, name = name, draw_labels = True,
        # draw_upper_landscape = True, draw_lower_landscape = True, 
        draw_upper_landscape = False, draw_lower_landscape = False, 
        draw_coastal_edges = False, draw_all_edges = False,
        draw_cusp_leaves = False,
        shade_triangles = False, 
        draw_fund_domain = False, fund_dom_tets = continent_fund_dom_tets,
        draw_fund_domain_edges = False, 
        #draw_tetrahedron_rectangles = tetrahedra_to_draw,
        edge_rectangles_to_draw = edge_rectangles_to_draw,
        tetrahedron_rectangles_to_draw = tetrahedra_to_draw,
        tetrahedron_rectangles_to_shade = tetrahedron_rectangles_to_shade,
        triangles_to_draw = triangles_to_draw,
        tetrahedra_to_draw = tetrahedra_to_draw,
        draw_edges_for_edge_rectangles = True,
        leaves_to_draw = leaves_to_draw,
        intervals_to_draw = intervals_to_draw,
        edge_thickness = 0.02,
        leaf_thickness = 0.02,
        transparency = 0.8)

    # draw_drilled_tetrahedra(con, name = name, draw_vertex_numbers = False, tetrahedra_cusp_orders = tetrahedra_cusp_orders, intervals_inside_tet_rectangles = intervals_inside_tet_rectangles, tetrahedra_chunks = tetrahedra_chunks)


def drill_flow_cycle_script():
    from draw_continent_circle import draw_continent_circle, complete_tetrahedron_rectangles
    from build_continent import make_continent_drill_flow_cycle, get_fund_domain_tetrahedra, make_continent_naive, make_continent_fund_dom
    from draw_drilled_tetrahedra import draw_drilled_tetrahedra
    from ordered_rectangles import build_tetrahedron_rectangle_orderings, sanity_check, build_drilled_triangulation_data
    from drilling_flow_cycle import triangulation_data_to_tri_angle, drill_flow_cycle
    from veering.taut import isosig_from_tri_angle
    from boundary_triangulation import generate_boundary_triangulation

    # veering_isosig = 'cPcbbbdxm_10' 
    # flow_cycle = [(0, 2)]

    # veering_isosig = 'cPcbbbiht_12'
    # flow_cycle = [(0, 0)]  ### boundary parallel
    # flow_cycle = [(0, 4), (1, 2)]
    # flow_cycle = [(0, 0), (0, 5)]  
    # flow_cycle = [(1, 0), (1, 5)]  ### blue edge rectangle has two punctures, red has none.

    # veering_isosig = 'dLQacccjsnk_200'
    # flow_cycle = [(1, 4)]
    # flow_cycle = [(1, 0), (2, 5)]
    # flow_cycle = [(0, 1), (2, 5), (1, 1)]  ## sideways and up

    # veering_isosig = 'eLAkaccddjsnak_2001'
    # flow_cycle = [(0, 1), (2, 2), (3, 2)]

    veering_isosig = 'gLLAQbecdfffhhnkqnc_120012'
    # flow_cycle = [(0, 0)]
    # flow_cycle = [(0, 4), (4, 0)]
    # flow_cycle = [(0, 4), (4, 5), (2, 4), (1, 2), (5, 1)]  
    flow_cycle = [(0, 4), (4, 5), (2, 4), (1, 2), (5, 1), (0, 4), (4, 0)]  ### non simple cycle example
    # flow_cycle = [(0, 4), (4, 5), (2, 4), (1, 2), (5, 1), (0, 4), (4, 0), (0, 4), (4, 5), (2, 4), (1, 2), (5, 1)]  
    
    # flow cycle which is not boundary parallel but completely in a blue region 
    # veering_isosig = 'eLAkbbcdddhwqj_2102'
    # flow_cycle = [(3, 2)]   ### this is isotopic to an LMT "AB-cycle", 
    ### so there will be distinct flow lines that fellow-travel
    # flow_cycle = [(2, 5)] ### has same issue of being isotopic to a pi_1 translate
    # flow_cycle = [(0, 3), (1, 1)]  ### goes straight up but does not stay in one colour so we can separate

    # short cycle in a large triangulation
    # veering_isosig = 'qLLALvQLwQMkcecdhkjilmmlnoppphqrwqqaqxhhxofhxq_2010222010122221'
    # flow_cycle = [(1, 0)]
    # # flow_cycle = [(3, 5)]
    # # flow_cycle = [(4, 2)]
    # flow_cycle = [(7, 5)]
    # flow_cycle = [(9, 0)]
    
    print(veering_isosig, flow_cycle)

    out = drill_flow_cycle(veering_isosig, flow_cycle, return_tri_angle = True, draw_rectangles = True, return_found_parallel = True, return_cusp_mapping = True, use_untwisted_speed_up = True, verbose = 0)
    drilled_sig, drilled_tri, drilled_angle, found_parallel, cusp_mapping_isosig = out


    print(veering_isosig, '->', drilled_sig, 'cusp_mapping', cusp_mapping_isosig)

    original_b = generate_boundary_triangulation(veering_isosig, draw = False)
    print(original_b.ladder_counts())
    drilled_b = generate_boundary_triangulation(drilled_sig, draw = False)
    print(drilled_b.ladder_counts())

def compare_flow_and_geodesic_drilling_script():  
    from veering.taut import isosig_from_tri_angle
    from veering.flow_cycles import generate_flow_cycles, flow_cycle_to_dual_edge_loop
    from drilling_flow_cycle import drill_flow_cycle
    from snappy_drill_homotopic import drill_tet_and_face_indices
    from snappy.drilling.exceptions import GeodesicSystemNotSimpleError
    import snappy
# 
    # sig = 'cPcbbbdxm_10'
    sig = 'cPcbbbiht_12'
    # sig = 'dLQacccjsnk_200' 
    # sig = 'dLQbccchhfo_122'
    # sig = 'dLQbccchhsj_122'

    cycles = generate_flow_cycles(sig, max_length = 5)

    for fc in cycles:  
        print(fc)
        out = drill_flow_cycle(sig, fc, return_tri_angle = True) 
        if out != None: 
            tri, angle = isosig_to_tri_angle(sig) 
            
            drilled_sig, drilled_tri, drilled_angle = out     
            drilled_M = snappy.Manifold(drilled_tri) 

            orig_M = snappy.Manifold(tri) 
            dual_loop = flow_cycle_to_dual_edge_loop(tri, angle, fc) 
            # print(fc, dual_loop) 
            try:
                ### May need to sys.setrecursionlimit(1000000) to make this work
                snappy_drilled_M = drill_tet_and_face_indices(orig_M, dual_loop, verified = True, bits_prec = 200) 
            except GeodesicSystemNotSimpleError as e:
                print(e)
                continue
            snappy_drilled_M.simplify()
            drilled_M.simplify()
            found_isometry = False
            for i in range(10):
                if drilled_M.is_isometric_to(snappy_drilled_M):  ### from the docstring for is_isometric_to:
                ### The answer True is rigorous, but the answer False may
                ### not be as there could be numerical errors resulting in finding
                ### an incorrect canonical triangulation.
                    found_isometry = True
                    break
            if not found_isometry:
                # print('checking with verified isometry_signature', sig, fc)
                isomsig1 = drilled_M.isometry_signature(verified = True)
                isomsig2 = snappy_drilled_M.isometry_signature(verified = True)
                assert not isomsig1 == None, 'isom signature failed ' + sig + ' ' + fc
                assert not isomsig2 == None, 'isom signature failed ' + sig + ' ' + fc
                # assert isomsig1 == isomsig2, sig + '_' + fc
                if isomsig1 != isomsig2:
                    print('drilling', sig, 'along', fc, 'gives different results', isomsig1, drilled_M.identify(), isomsig2, snappy_drilled_M.identify())



def draw_veering_triangulation_and_mid_annuli_script():
    from draw_veering_triangulation_and_mid_annuli import draw_triangulation_from_veering_isosig
    # draw_triangulations_from_veering_isosigs_file('Veering_census/veering_census_up_to_12.txt', '../Veering_mid-annuli', '../Veering_tetrahedra')

    # calculate_well_framed_from_veering_isosigs_file('Veering_census/veering_census.txt', 'veering_census_with_well_framed.txt')
    # calculate_well_framed_from_veering_isosigs_file('Veering_census/veering_euler.txt', 'veering_euler_with_well_framed.txt')
    # calculate_well_framed_from_veering_isosigs_file('Veering_census/veering_covers_veering_isosigs.txt', 'veering_covers_with_well_framed.txt')


    # draw_triangulation_from_veering_isosig('dLQacccjsnk_200') 
    # draw_triangulation_from_veering_isosig('iLMzMPcbcdefghhhhhqqqqhxq_12200221')
    # draw_triangulation_from_veering_isosig('lLLAvPAQcbcdejhihkjjktsfxhxhqsjfw_20102120121')
    # draw_triangulation_from_veering_isosig('lLLvLQAQcbefgigijkikkxxxgbrgrqraq_01110220202')
    # draw_triangulation_from_veering_isosig('lLLvLQPQcbeffhgjikjkkxxxfjssfxfhq_10221010102')
    # draw_triangulation_from_veering_isosig('mvvLPQwQQfjfihgfkllklkhahaaaaqqxxaa_210012222211')
   
    # draw_triangulation_from_veering_isosig('gLLAQcdecfffhsermws_122201') 
   
    # draw_triangulation_from_veering_isosig('iLLLQPcbeefeghhhhhhqhqxxa_21112001')
    # draw_triangulation_from_veering_isosig('iLLLQPcbeefeghhhhhhqhqxxa_01110221')

    # draw_triangulation_from_veering_isosig('qLAMzMzMzMLkbcbdefghijklmnppphhwhhhhhhhhhhhnkw_2112121212121211')

    # draw_triangulations_from_veering_isosigs_file('Gadgets/double_gadget_veering_isosigs.txt', 'Gadgets/Veering_mid-annuli', 'Gadgets/Veering_tetrahedra' )

    # draw_triangulations_from_veering_isosigs_file('iLLLQPcbeefeghhhhhhqhqxxa_veer_isosigs.txt', 'iLLLQPcbeefeghhhhhhqhqxxa_mid-annuli', 'iLLLQPcbeefeghhhhhhqhqxxa_veering_tetrahedra' )

    # draw_triangulations_from_veering_isosigs_file('m016_surgeries.txt', 'm016_surgeries_mid-annuli', 'm016_surgeries_veering_tetrahedra' )



    # sigs = ['ivvPQQcfhghgfghfaaaaaaaaa_01122000', 
    #         'mvvLPQwQQfghffjikllklkaaaaaaaaaaaaa_102021111100', 
    #         'mvvLPQwQQhjhihgflkkllkaaaaaaaaaaaaa_022220001122',
    #         'oLvLAvwQQQccfghkmnkjklnmmnhialaatiqqqffff_20102220210100',
    #         'oLvLAwzPQQccfghhilkjnmmnmnhialgoiqqffffff_20102220200011']

    # flat_toggles = ['qLLvAvAMQLQkbeehklmnjnnppopooxxxahahxxxaxqxxxq_2111200221111100',
    #                 'qLLvAvPMQLQkbeehlkmnjnnpoopopxxxaaxxxxxaxaxxax_2111200221111100',
    #                 'qvvLPQMvQLQkfgfhhgfknlmoppopohahhaaahaqqaqqaaa_1222211100222200']
    # for sig in flat_toggles:
    #     draw_triangulation_from_veering_isosig(sig)


    # to_draw_for_paper = ['cPcbbbdxm_10',
    #                      'cPcbbbiht_12', 
    #                      'eLMkbcddddedde_2100', 
    #                      'gLLAQbecdfffhhnkqnc_120012', 
    #                      'fLLQccecddehqrwjj_20102']
    # for sig in to_draw_for_paper:
    #     draw_triangulation_from_veering_isosig(sig)
    
    # draw_triangulation_from_veering_isosig('qvvLPQMvQLQkfgfhhgfknlmoppopohahhaaahaqqaqqaaa_1222211100222200')

    # draw_triangulation_from_veering_isosig('cPcbbbiht_12')

    # from file_io import parse_data_file
    
    # census = parse_data_file('veering_census.txt')
    # for sig in census[:5]:
    #     draw_triangulation_from_veering_isosig(sig)

    census = parse_data_file('veering_census.txt')
    for sig in census:
        # if ord(sig[0]) <= 109: ### up to 'm'
        if ord(sig[0]) == 109: ### 'm'
            draw_triangulation_from_veering_isosig(sig, tetrahedra_filename = 'Images/Triangulation/Census/' + sig + '_tetrahedra.pdf', midannuli_filename = 'Images/Mid-annuli/Census/' + sig +'_mid-annuli.pdf')

    # only_toggles = parse_data_file('only_toggles.txt')
    # for sig in only_toggles:
    #     draw_triangulation_from_veering_isosig(sig, midannuli_filename = 'Images/Mid-annuli/only_toggles/' + sig +'_mid-annuli.pdf')

    # data = parse_data_file('low_volume_stacks.txt')
    # for line in data:
    #     split_data = line.split(" ")
    #     sig, vol = split_data[:2]
    #     stack = "".join(split_data[2:])
    #     draw_triangulation_from_veering_isosig(sig, midannuli_filename = 'Images/Mid-annuli/Low_volume_stacks/' + sig + stack +'_mid-annuli.pdf')

def normal_surfaces_script(num_to_check = 10):
    from normal_surfaces import analyze_sig
    lines = parse_data_file('veering_census.txt')
    for line in lines[:num_to_check]:
        sig = line.strip()
        analyze_sig(sig)

def pachner_graph_script():
    # from pachner_graph import *
    depth = 100
    # ceiling = 10
    extra_ceiling = 5  ### above how many tetrahedra we start with
    # print('depth', depth)
    # print('ceiling', ceiling)

    # # sig = 'cPcbbbdxm_10'
    # # sig = 'dLQacccjsnk_200'
    # # sig = 'dLQbccchhfo_122'
    # sig = 'eLAkbccddhhsqs_1220'
    # tri, angle = isosig_to_tri_angle(sig) 
    # branch = upper_branched_surface(tri, angle)
    # # tl = flow_cycle_to_triangle_loop(tri, branch, [(0, 2)]) 
    # # drilled_cusp_index = drill(tri, tl, angle = angle, branch = branch) 
    # # print('angle', angle, 'branch', branch, 'drilled_cusp_index', drilled_cusp_index)
    # # # assert has_non_sing_semiflow(tri, branch)

    # # # start veering: cPcbbbdxm_10_dl loop [(0, 2)] tri_loop [(2, 102)]
    # # # drill: eLMkbbddddhapu_2100_fjek

    # # start_isoSig = isosig_from_tri_angle_branch(tri, angle, branch)
    # # assert start_isoSig == 'eLMkbbddddhapu_2100_fjek'  
    # # tri, angle, branch = isosig_to_tri_angle_branch(start_isoSig)     ### fails?? 
    

    # # graph = search_Pachner_graph_for_shortest_path(start_isoSig, tri, angle, branch,  name=None, search_depth = depth, ceiling = ceiling, drilled_cusp_index = drilled_cusp_index, check_property = False, property = None, save_dir = None)


    # loops = find_flow_cycles(tri, branch)
    # tri_loops = [flow_cycle_to_triangle_loop(tri, branch, loop) for loop in loops]
    
    # for i, tri_loop in enumerate(tri_loops):
    #     if tri_loop != False: # False means that tri_loop goes more than once  along the same triangle - not currently implemented
    #         tri, angle = isosig_to_tri_angle(sig)
    #         if tri_loop_is_boundary_parallel(tri_loop, tri) == False: # if a loop is boundary parallel then we don't drill
    #             branch = upper_branched_surface(tri, angle)
    #             drilled_cusp_index = drill(tri, tri_loop, angle = angle, branch = branch)
    #             print('loop', loops[i], 'tri_loop', tri_loop, 'has_essential_torus', has_essential_torus(tri))
    #             # start_isoSig = isosig_from_tri_angle_branch(tri, angle, branch)
    #             # print(start_isoSig, 'angle', angle, 'branch', branch, 'drilled_cusp_index', drilled_cusp_index)
    #             # graph = search_Pachner_graph_for_shortest_path(start_isoSig, tri, angle, branch,  name=None, search_depth = depth, ceiling = ceiling, drilled_cusp_index = drilled_cusp_index, check_property = False, property = None, save_dir = None)


    sigs = parse_data_file('veering_census.txt')

    for j, sig in enumerate(sigs[:5]):
        if j%100 == 0:
            print(j)
        tri, angle = isosig_to_tri_angle(sig)
        branch = upper_branched_surface(tri, angle)
        loops = find_flow_cycles(tri, branch)
        tri_loops = [flow_cycle_to_triangle_loop(tri, branch, loop) for loop in loops]
            
        for i, tri_loop in enumerate(tri_loops):
            if tri_loop != False: # False means that tri_loop goes more than once  along the same triangle - not currently implemented
                tri, angle = isosig_to_tri_angle(sig)
                if tri_loop_is_boundary_parallel(tri_loop, tri) == False: # if a loop is boundary parallel then we don't drill
                    branch = upper_branched_surface(tri, angle)
                    drilled_cusp_index = drill(tri, tri_loop, angle = angle, branch = branch)
                    if has_essential_torus(tri):
                        print('sig', sig, 'has_essential_torus', has_essential_torus(tri), 'loop', loops[i], 'tri_loop', tri_loop)
                    else:
                        print('sig', sig, 'has_essential_torus', has_essential_torus(tri), 'loop', loops[i], 'tri_loop', tri_loop)
                        ceiling = tri.countTetrahedra() + extra_ceiling
                        print('ceiling', ceiling)
                        start_isoSig = isosig_from_tri_angle_branch(tri, angle, branch)
                        print(start_isoSig, 'angle', angle, 'branch', branch, 'drilled_cusp_index', drilled_cusp_index)
                        graph = search_Pachner_graph_for_shortest_path(start_isoSig, tri, angle, branch,  name=None, search_depth = depth, ceiling = ceiling, drilled_cusp_index = drilled_cusp_index, check_property = False, property = None, save_dir = None)

                    # start_isoSig = isosig_from_tri_angle_branch(tri, angle, branch)
                    # print(start_isoSig, 'angle', angle, 'branch', branch, 'drilled_cusp_index', drilled_cusp_index)
                    # graph = search_Pachner_graph_for_shortest_path(start_isoSig, tri, angle, branch,  name=None, search_depth = depth, ceiling = ceiling, drilled_cusp_index = drilled_cusp_index, check_property = False, property = None, save_dir = None)

    ## drilling fLLQcbeddeehhbghh_01110 along loop [(0, 2), (3, 1), (2, 3)] has essential torus - where is the mobius strip?
    ## we got through the first five manifolds and found veering for all simple drillings without essential tori.

def taut_branched_pachner_graph_script():
    # from taut_branched_pachner_graph import *
    depth = 100
    ceiling = 8

    print('depth', depth)
    print('ceiling', ceiling)

    target_isoSig = 'gLLPQceeffefiiaellu_012110'  ### drilled
    start_isoSig = 'gLLPQccdfeffhggaagb_201022'  ### veering
    tri, angle = isosig_to_tri_angle(start_isoSig)
    branch = upper_branched_surface(tri, angle)
    start_isoSig = isosig_from_tri_angle_branch(tri, angle, branch)
    graph = search_Pachner_graph_for_shortest_path(start_isoSig, tri, angle, branch, target_isoSig, name=None, search_depth = depth, ceiling = ceiling, check_property = False, property = None, save_dir = None)

def taut_pachner_script():
    # from taut_pachner import *
    tri, angle = isosig_to_tri_angle('jLLAvQQbcdeihhiihtsfxedxhdt_201021201')
    # tri, angle = isosig_to_tri_angle('cPcbbbiht_12')
    for i in range(tri.countTriangles()):
        print('triangle_num', i)
        tri2 = regina.Triangulation3(tri)
        angle2 = angle[:]
        tri3, angle3, edge_num = twoThreeMove(tri2, angle2, i, return_edge = True) 

        print('angle3', angle3, 'edge_num', edge_num)
        threeTwoMove(tri3, angle3, edge_num)

        # for e_ind in range(tri3.countEdges()):
        #     if threeTwoMove(tri3, angle3, e_ind, perform = False):
        #         print('could do 3-2 on e_ind', e_ind)
        #         tri4, angle4 = threeTwoMove(tri3, angle3, e_ind)
        #         break

def taut_pachner_graph_script():
    # from taut_pachner_graph import *
    depth = 100
    ceiling = 8

    print('depth', depth)
    print('ceiling', ceiling)
    
    target_isoSig = 'gLLPQceeffefiiaellu_012110'  ### drilled
    start_isoSig = 'gLLPQccdfeffhggaagb_201022'  ### veering
    
    graph = search_Pachner_graph_for_shortest_path(start_isoSig, target_isoSig, name=None, search_depth = depth, ceiling = ceiling, check_property = False, property = None, save_dir = None)

### Pulled out of flow_cycles.py

def test_branched_isosig_script():
    # from test_branched_isosig import *
    # tri, angle = isosig_to_tri_angle("cPcbbbdxm_10")
    tri, angle = isosig_to_tri_angle("gLAMPbbcdeffdhwqqqj_210202")

    out = explore_mobius_surgery_graph(tri, angle, max_tetrahedra = 7)
    for sig in out:
        print(sig)

def test():
    # sig = 'cPcbbbiht_12'
    # sig = 'dLQacccjsnk_200'
    # sig = 'dLQbccchhsj_122'
    # sig = 'eLAkaccddjsnak_2001'
    # sig = 'eLAkbccddhhsqs_1220'
    # sig = 'eLMkbcddddedde_2100'
    # sig = 'eLMkbcdddhhhdu_1221'
    # sig = 'eLMkbcdddhhhml_1221'
    # sig = 'eLMkbcdddhhqqa_1220'
# eLMkbcdddhhqxh_1220
# eLMkbcdddhxqdu_1200
# eLMkbcdddhxqlm_1200
# eLPkaccddjnkaj_2002
# eLPkbcdddhrrcv_1200
    # sig = 'gLLAQbecdfffhhnkqnc_120012'

    sigs = parse_data_file('veering_census.txt')

    for j, sig in enumerate(sigs[:5]):
        if j%100 == 0:
            print(j)
        tri, angle = isosig_to_tri_angle(sig)
        # tri.save(sig + '.rga')
        branch = upper_branched_surface(tri, angle) ### also checks for veering and transverse taut
        found_loops = find_flow_cycles(tri, branch)
        # print(len(found_loops))
        # for loop in found_loops:
        #   print(loop)

        # print('found_loops', found_loops)
        # print(sig)
        for loop in found_loops:
            tri, angle = isosig_to_tri_angle(sig)
            branch = upper_branched_surface(tri, angle) 
            tri_loop = flow_cycle_to_triangle_loop(tri, branch, loop)
            if tri_loop != False:
                if not tri_loop_is_boundary_parallel(tri_loop, tri):
                    print('sig', isosig_from_tri_angle_branch(tri, angle, branch), 'loop', loop, 'tri_loop', tri_loop) 
                    drill(tri, tri_loop, angle = angle, branch = branch, sig = sig)
                    print('new angle, branch', angle, branch)
                    print(isosig_from_tri_angle_branch(tri, angle, branch))
                    # tri.save('drilled_' + sig + '_' + str(tri_loop) + '.rga')
                    # print(tri.countTetrahedra())
                    
def test_layered_parents(sig, quiet = False):    
    tri, angle = isosig_to_tri_angle(sig)
    branch = upper_branched_surface(tri, angle)
    loops = find_flow_cycles(tri, branch)
    tri_loops = [flow_cycle_to_triangle_loop(tri, branch, loop) for loop in loops]
    
    no_of_layered_parents = 0 # drillings through some simple cycles is not implemented so might get 0 even if there is a simple cycle which gives a layered parent
    for tri_loop in tri_loops:
        if tri_loop != False: # False means that tri_loop goes more than once  along the same triangle - not currently implemented
            tri, angle = isosig_to_tri_angle(sig)
            if tri_loop_is_boundary_parallel(tri_loop, tri) == False: # if a loop is boundary parallel then we don't drill
                tri, angle = isosig_to_tri_angle(sig)
                branch = upper_branched_surface(tri, angle)
                if quiet == False:
                    print ("drilling", sig, "along", tri_loop)
                drill(tri, tri_loop, angle, branch)
                if quiet == False:
                    print("drilled:", tri.isoSig(), angle, branch)
                    print("is layered:", is_layered(tri, angle))
                if is_layered(tri, angle):
                    no_of_layered_parents = no_of_layered_parents + 1
    if no_of_layered_parents == 0:
        print (sig, "does not have a layered parent")

        
def test_semiflow_on_drillings(sig):
    
    tri, angle = isosig_to_tri_angle(sig)
    branch = upper_branched_surface(tri, angle)
    loops = find_flow_cycles(tri, branch)
    tri_loops = [flow_cycle_to_triangle_loop(tri, branch, loop) for loop in loops]
    
    for tri_loop in tri_loops:
        if tri_loop != False: # False means that tri_loop goes more than once  along the same triangle - not currently implemented
            tri, angle = isosig_to_tri_angle(sig)
            if tri_loop_is_boundary_parallel(tri_loop, tri) == False: # if a loop is boundary parallel then we don't drill
                tri, angle = isosig_to_tri_angle(sig)
                branch = upper_branched_surface(tri, angle)
                drill(tri, tri_loop, angle, branch)
                assert has_non_sing_semiflow(tri, branch)
                # print (tri.isoSig(), branch, "has nonsing semiflow")
