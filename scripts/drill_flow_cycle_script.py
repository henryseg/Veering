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

### 0. Given a veering triangulation and a flow cycle (all in the downstairs manifold)
### 1. following a flow line (an elevation of the flow cycle), build continent (upstairs). We get a \emph{flow arc} in the continent.
### 2. identify fund domain D (perhaps we dont care if we havent hit all tetrahedra - some may not be affected by the drilling).
###    Let X be the image of D in the link space. 
###    The main puncture p is the intersection of the tetrahedron rectangles along the flow line
### 3. We have the main puncture p. For each edge rectangle E of X, determine whether or not p is in E (possibly extend the flow arc)
### 4. For each punctured edge rectangle E in X, find the flow lines through the translates of E (in X).







### ??. Let S be the support for a flow line upstairs (the tetrahedra the flow line goes through)


### function: given edge rectangle and flow line (specified by some flow arc), is the point specified by the flow line in the edge rectangle? 
### (may need to extend continent to find out)
