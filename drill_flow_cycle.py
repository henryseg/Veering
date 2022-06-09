from draw_continent_circle import make_continent_drill_flow_cycle, get_fund_domain_tetrahedra, complete_tetrahedron_rectangles



def main():
    veering_isosig = 'cPcbbbdxm_10' 
    flow_cycle = [(0, 2)]

    # veering_isosig = 'eLAkaccddjsnak_2001'
    # flow_cycle = [(1, 0), (2, 5)]

    # for num_steps in range(10):
    num_steps = 10
    con, flow_tetrahedra, flow_edges = make_continent_drill_flow_cycle(veering_isosig, flow_cycle, num_steps)
    fund_dom_tets = get_fund_domain_tetrahedra(con)
    # tets_to_draw = [flow_tetrahedra[0], flow_tetrahedra[-1]]
    # tets_to_draw = fund_dom_tets
    tets_to_draw = fund_dom_tets + [flow_tetrahedra[0], flow_tetrahedra[-1]]
    complete_tetrahedron_rectangles(con, tets_to_draw)
    print(len(flow_tetrahedra))

    name = veering_isosig + '_' + str(flow_cycle) + '_' + str(num_steps) + '_cusp_leaves'

### 1. build continent for the flow cycle
### 2. identify fund domain D (perhaps we dont care if we havent hit all tetrahedra - some may not be affected by the drilling)
###    The main puncture p is the intersection of the image of the flow cycle in the continent
### 3. We have the main puncture p. For each edge rectangle e of D, determine whether or not p is in e (possibly extend the flow cycle)
### 4. For each punctured edge rectangle in D, find the flow cycles through its translates in D 




### 4. Let S be a support for the flow cycle upstairs (the tetrahedra the flow cycle goes through)