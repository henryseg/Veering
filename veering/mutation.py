#
# mutation.py
#

# mutating veering triangulations along carried surfaces

import regina

from .taut import isosig_to_tri_angle, isosig_from_tri_angle, is_taut
from .transverse_taut import is_transverse_taut, top_bottom_embeddings_of_faces
from .veering_tri import is_veering
from .taut_polynomial import tet_lower_upper_edges
from .taut_polytope import is_layered
from .edge_orientability import is_edge_orientable
from .carried_surface import build_surface, triangle_is_red, veering_symmetry_group, stratum

def vertex_correspondence_between_embeds(embed0, embed1):
    vertex00 = embed0.vertices()[0]
    vertex01 = embed0.vertices()[1]
    vertex02 = embed0.vertices()[2]
    vertex10 = embed1.vertices()[0]
    vertex11 = embed1.vertices()[1]
    vertex12 = embed1.vertices()[2]
    
    return [vertex00, vertex10], [vertex01, vertex11], [vertex02, vertex12]
    
def verts_of_large_small_same_coloured_edges_in_embeds(tri, angle, face_index, tet_vert_coorientations = None):
    """
    returns [vertices of large edge, vertices of small edge of the same colour] for top embedding and bottom embedding of a face
    """
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
        
    top_embed, bottom_embed = top_bottom_embeddings_of_faces(tri, angle, tet_vert_coorientations)[0][face_index], top_bottom_embeddings_of_faces(tri, angle, tet_vert_coorientations)[1][face_index]
    
    tet_below = top_embed.simplex()
    tet_above = bottom_embed.simplex()
    tet_below_coor = tet_vert_coorientations[tet_below.index()]
    tet_above_coor = tet_vert_coorientations[tet_above.index()]
    verts_of_top_diagonal_of_tet_below = [i for i in range(4) if tet_below_coor[i] == -1] # this edge is small, of the same colour as large
    verts_of_bottom_diagonal_of_tet_above = [i for i in range(4) if tet_above_coor[i] == +1] # this edge is large
    #print('small edge in top embed', verts_of_top_diagonal_of_tet_below)
    #print('large edge in bottom embed', verts_of_bottom_diagonal_of_tet_above)
    
    vertex_correspondence = vertex_correspondence_between_embeds(bottom_embed, top_embed) 
    #print('bottom top vertex correspondence', vertex_correspondence)
    verts_of_large_edge_of_top_embed = []
    for x in vertex_correspondence:
        if x[0] == verts_of_bottom_diagonal_of_tet_above[0] or x[0] == verts_of_bottom_diagonal_of_tet_above[1]:   
            verts_of_large_edge_of_top_embed.append(x[1])
    verts_of_other_edge_of_same_colour_in_bottom_embed = []
    for x in vertex_correspondence:
        if x[1] == verts_of_top_diagonal_of_tet_below[0] or x[1] == verts_of_top_diagonal_of_tet_below[1]:   
            verts_of_other_edge_of_same_colour_in_bottom_embed.append(x[0])
        
    #print('large edge in top embed', verts_of_large_edge_of_top_embed)
    #print('small edge in bottom embed', verts_of_other_edge_of_same_colour_in_bottom_embed)
    
            
    return [verts_of_large_edge_of_top_embed, verts_of_top_diagonal_of_tet_below], [verts_of_bottom_diagonal_of_tet_above, verts_of_other_edge_of_same_colour_in_bottom_embed]
       
def top_triangles_to_surface_triangles(tri, angle, weights, tet_vert_coorientations = None):
    """
    Gives a map from the top triangles of the cut 3-manifold to the surface triangles (works only for surface built using build_surface!).
    If weight >  1 then the map sends a triangle f of the 3D-triangulation to the triangle of the surface corresponding to the "lowermost" (first) copy of f in the surface.
    The output format is a list of length 2*num_tet, whose entries are either None (if the face is not unglued in the cut manifold) or a list [tet_num, [vertices of tet], triangle_num, [vertices of tri]]
    """
    
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    
    n = tri.countTetrahedra()
    
    top_embeddings, _ = top_bottom_embeddings_of_faces(tri, angle, tet_vert_coorientations)
    #print('top embeddings', top_embeddings)
    
    for i in range(2*n):
        if weights[i] > 0: # the ith face of the 3D-triangulation is unglued
            which_surface_triangle = sum(weights[:i]) # the index of the surface triangle corresponding to the "lowermost" copy of face i
            top_embed = top_embeddings[i]
            [verts_of_large_edge, verts_of_small_edge_of_same_colour] = verts_of_large_small_same_coloured_edges_in_embeds(tri, angle, i, tet_vert_coorientations)[0]
            common_vertex = [x for x in verts_of_large_edge if x in verts_of_small_edge_of_same_colour][0]
            #print('verts of large', verts_of_large_edge)
            #print('verts_of_small_edge_of_same_colour', verts_of_small_edge_of_same_colour)
            #print('common vertex', common_vertex)
            #verts_of_large_edge_copy = verts_of_large_edge.copy()
            verts_of_large_edge.remove(common_vertex)
            other_vertex_of_large = verts_of_large_edge[0]
            #print('other_vertex_of_large', other_vertex_of_large)
            verts_of_small_edge_of_same_colour.remove(common_vertex)
            other_vertex_of_small = verts_of_small_edge_of_same_colour[0]
            if triangle_is_red(tri, angle, i, tet_vert_coorientations):
                # mapping common_vertex --> 0, other_vertex_of_large --> 1, other_vertex_of_small --> 2
                vertex_mapping = [top_embed.simplex().index(), [common_vertex, other_vertex_of_large, other_vertex_of_small], which_surface_triangle,  [0,1,2]]
            else:
                # mapping common_vertex --> 1, other_vertex_of_large --> 0, other_vertex_of_small --> 2
                vertex_mapping = [top_embed.simplex().index(), [other_vertex_of_large, common_vertex, other_vertex_of_small], which_surface_triangle, [0,1,2]]
            #top_embeddings[i] = [top_embeddings[i]]
            #top_embeddings[i].append(vertex_mapping)
            top_embeddings[i] = vertex_mapping
        else:
            top_embeddings[i] =  None #so that it is clear that this face is not unglued in the cut manifold
            
    return top_embeddings
                        
def surface_triangles_to_bottom_triangles(tri, angle, weights, tet_vert_coorientations = None):
    """
    Gives a map from (some) surface triangles to bottom triangles of the cut 3-manifold.
    If weight >  1 then only the uppermost (last) copy of the 3-manifold triangle in the surface has a 3-manifold triangle assigned.
    The output format is a list of length 2*num_tet, whose entries are either None (if the face is not unglued in the cut manifold) or a list [triangle_num, [vertices of tri], tet_num, [vertices of tet]]
    """
    
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    
    n = tri.countTetrahedra()
    
    _, bottom_embeddings = top_bottom_embeddings_of_faces(tri, angle, tet_vert_coorientations)

    for i in range(2*n):
        if weights[i] > 0: # the ith face of the 3D-triangulation is unglued
            which_surface_triangle = sum(weights[:i+1]) - 1 # the "uppermost" copy of face i which appears in the surface
            bottom_embed = bottom_embeddings[i]
            [verts_of_large_edge, verts_of_small_edge_of_same_colour] = verts_of_large_small_same_coloured_edges_in_embeds(tri, angle, i, tet_vert_coorientations)[1]
            common_vertex = [x for x in verts_of_large_edge if x in verts_of_small_edge_of_same_colour][0]
            #print('verts of large', verts_of_large_edge)
            #print('verts_of_small_edge_of_same_colour', verts_of_small_edge_of_same_colour)
            #print('common vertex', common_vertex)
            verts_of_large_edge.remove(common_vertex)
            other_vertex_of_large = verts_of_large_edge[0]
            #print('other_vertex_of_large', other_vertex_of_large)
            verts_of_small_edge_of_same_colour.remove(common_vertex)
            other_vertex_of_small = verts_of_small_edge_of_same_colour[0]
            if triangle_is_red(tri, angle, i, tet_vert_coorientations):
                #mapping common_vertex --> 0, other_vertex_of_large --> 1, other_vertex_of_small --> 2
                vertex_mapping = [which_surface_triangle,  [0,1,2], bottom_embed.simplex().index(), [common_vertex, other_vertex_of_large, other_vertex_of_small]]
            else:
                #mapping common_vertex --> 1, other_vertex_of_large --> 0, other_vertex_of_small --> 2
                vertex_mapping = [which_surface_triangle, [0,1,2], bottom_embed.simplex().index(), [other_vertex_of_large, common_vertex, other_vertex_of_small]]
            #bottom_embeddings[i] = [bottom_embeddings[i]]
            #bottom_embeddings[i].append(vertex_mapping)
            bottom_embeddings[i] = vertex_mapping
        else:
            bottom_embeddings[i] =  None #so that it is clear that this face is not unglued in the cut manifold
            
    return bottom_embeddings

def face_num_in_tet(verts_of_face):
    all_verts = [0,1,2,3]
    for i in range(3):
        all_verts.remove(verts_of_face[i])
    return all_verts[0]    

def vertex_correspondence_to_perm4(x, y): 
    """
    given a correspondence between a triangle of a 3D-triangulation and a triangle of a 2D-triangulation, returns the associated permutation on 0,1,2,3
    """
    vertex_to_vertex = [[x[i], y[i]] for i in range(3)]
    vertex_to_vertex.sort(key = lambda x:x[0])
    missing_vertex = face_num_in_tet(x)
    
    if missing_vertex == 0:
        perm = regina.Perm4(3, vertex_to_vertex[0][1], vertex_to_vertex[1][1], vertex_to_vertex[2][1])
    elif missing_vertex == 1:
        perm = regina.Perm4(vertex_to_vertex[0][1], 3, vertex_to_vertex[1][1], vertex_to_vertex[2][1])
    elif missing_vertex == 2:
        perm = regina.Perm4(vertex_to_vertex[0][1], vertex_to_vertex[1][1], 3, vertex_to_vertex[2][1])
    else:
        perm = regina.Perm4(vertex_to_vertex[0][1], vertex_to_vertex[1][1], vertex_to_vertex[2][1], 3)
    
    return perm
        
def surface_isom_to_regluing_pattern(tri, angle, weights, isom, tet_vert_coorientations = None): 
    """
    Translation between a combinatorial isomorphism of a carried surface and the associated mutation pattern.
    """
    #isom must be a combinatorial isomorphism of a surface build using build_surface
    
    # 1. Get the top embeds to surface triangles map (if weights > 1 then the lowermost copy - the one which has a tet below, and not a 'prism').
    # 2. Find images of triangles and permutations on vertices under the isomorphism.
    # 3. For every simpImage check if it is in the list obtained via surface_triangles_to_bottom triangles. If so, then it has a tet above (and not a prism), so the regluing pattern is 
    # top embed --> surface triangle (lowermost over top embed) --isom--> surface triangle (uppermost over some tet) --> bottom face of the tet above
    # Otherwise, we look at the copy of the 3D-triangle in the surface which is above this surface triangle, apply isom again etc.
    
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    
    top_tri_to_surface_tri = top_triangles_to_surface_triangles(tri, angle, weights, tet_vert_coorientations)
    surface_tri_to_bottom_tri = surface_triangles_to_bottom_triangles(tri, angle, weights, tet_vert_coorientations)
    
    k = len(surface_tri_to_bottom_tri)
    auxiliary = [None]*k 
    for i in range(k):
        if surface_tri_to_bottom_tri[i] != None:
            auxiliary[i] = surface_tri_to_bottom_tri[i][0]
    # auxiliary[k] is either None or an index of a surface triangle which has a tet above -- I want to keep Nones so that later I can find the index of the bottom triangle in surface_tri_to_bottom_tri by value
            
    regluing = []
    
    for i in range(len(top_tri_to_surface_tri)):
        if top_tri_to_surface_tri[i] != None: # this triangle is unglued in the cut manifold
            #print('top triangle to surface triangle', top_tri_to_surface_tri[i])
            surface_tri_index = top_tri_to_surface_tri[i][2] 
            tri_image = isom.simpImage(surface_tri_index)
            isom_perm = isom.facetPerm(surface_tri_index)
            #print('first tri image and perm', tri_image, isom_perm)
            if tri_image in auxiliary:
                index_in_surface_to_bottom_tri = auxiliary.index(tri_image)
                #print('surface triangle', tri_image, 'has a tet above')
            while tri_image not in auxiliary: #look at further images
                #print('tri_image does not have a tet above')
                triangle_above = tri_image + 1
                tri_image = isom.simpImage(triangle_above) #look at the image of the upper copy
                isom_perm = isom.facetPerm(triangle_above)*isom_perm
                #print('next tri image and perm', tri_image, isom_perm)
                if tri_image in auxiliary:
                    index_in_surface_to_bottom_tri = auxiliary.index(tri_image)
                    #print('surface triangle', tri_image, 'has a tet above')
            top_to_surface_perm = vertex_correspondence_to_perm4(top_tri_to_surface_tri[i][1], top_tri_to_surface_tri[i][3])
            #print('top tri to surface tri perm', top_to_surface_perm)
            isom_perm4 = regina.Perm4(isom_perm[0], isom_perm[1], isom_perm[2], 3)
            #print('4perm associated to isom', isom_perm4)
            bottom_to_surface_perm = vertex_correspondence_to_perm4(surface_tri_to_bottom_tri[index_in_surface_to_bottom_tri][3], surface_tri_to_bottom_tri[index_in_surface_to_bottom_tri][1])
            #print('bottom tri to surface tri perm', bottom_to_surface_perm)
            surface_to_bottom_perm = bottom_to_surface_perm.inverse()
            #print('surface to bottom perm', surface_to_bottom_perm)
            perm = surface_to_bottom_perm*isom_perm4*top_to_surface_perm 
            #print('gluing perm', perm)
            
            tet_below_top_embed = top_tri_to_surface_tri[i][0]
            which_face = face_num_in_tet(top_tri_to_surface_tri[i][1])
            tet_above_image = surface_tri_to_bottom_tri[index_in_surface_to_bottom_tri][2]
                                         
            regluing.append([tet_below_top_embed, which_face, tet_above_image, perm])
              
    return regluing
    
def mutate(tri, angle, weights, isom, tet_vert_coorientations = None, quiet = False): 
    
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    
    regluing = surface_isom_to_regluing_pattern(tri, angle, weights, isom, tet_vert_coorientations)
    #print('regluing data:', regluing)
    
    # first unglue all faces that need to be unglued, then reglue them via isom
    # can't do both at the same step, because regina goes crazy when you try to glue a face to a face which is already glued somewhere else
    for x in regluing:
        tet_below_top_triangle = tri.tetrahedron(x[0])
        which_face = x[1]
        tet_below_top_triangle.unjoin(which_face)
        #print('unglued face', which_face, 'of tet', tet_below_top_triangle.index())
    
    for x in regluing:
        tet_below_top_triangle = tri.tetrahedron(x[0])
        which_face = x[1] 
        tet_above_image_triangle = tri.tetrahedron(x[2])
        perm = x[3]
        tet_below_top_triangle.join(which_face, tet_above_image_triangle, perm)
        #print('glued face', which_face, 'of tet', tet_below_top_triangle.index(), 'to tet', tet_above_image_triangle.index(), 'via', perm)
    
    assert tri.isValid()
    assert tri.countBoundaryFacets() == 0
    
    if quiet == False:
        print('triangulation isosig:', tri.isoSig())
        print('taut:', is_taut(tri, angle))
        if is_taut(tri, angle):
            print('taut isosig:', isosig_from_tri_angle(tri, angle))
            print('transverse taut:', is_transverse_taut(tri, angle))
            print('veering:', is_veering(tri, angle))
            print('layered:', is_layered(tri, angle))
            print('edge-orientable:', is_edge_orientable(tri, angle))
        else:
            edge_num = tri.countEdges()
            print('Got', edge_num, 'edges out of', tri.countTetrahedra())
            totals = is_taut(tri, angle, return_totals = True)
            print('Angles:', totals)
            

def perform_all_mutations(tri, angle, weights, tet_vert_coorientations = None, print_stratum = True):
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
        
    surface, edge_colours = build_surface(tri, angle, weights, return_edge_colours = True)
    
    if print_stratum == True:
        this_stratum = stratum(tri, angle, weights)
        print('stratum:', this_stratum)
        
    
    veering_isoms = veering_symmetry_group(surface, edge_colours)
    
    if len(veering_isoms) > 1:
        sig = isosig_from_tri_angle(tri, angle)
        for i in range(1, len(veering_isoms)):
            isom = veering_isoms[i]
            tri, angle = isosig_to_tri_angle(sig)
            mutate(tri, angle, weights, isom)
    else:
        print('surface has no veering symmetries')
    

        
    

        
        
            
   

    
    
