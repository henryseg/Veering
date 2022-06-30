#
# carried_surface.py
#

# build regina 2D-triangulation of a surface carried by a veering triangulation

# Conventions: We build a carried surface out of positively oriented triangles. 
# For every triangle the large edge corresponds to the edge (01)


import regina 
from transverse_taut import is_transverse_taut, edge_side_face_collections
from taut import liberal
from veering import is_veering

def where_faces_are_uppermost(tri, angle, tet_vert_coorientations = None):
    """
    returns a list whose ith elements is the index of the edge where the ith face is uppermost
    """
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
        
    n = tri.countTetrahedra()
        
    where_uppermost = [None]*2*n
    sides = edge_side_face_collections(tri, angle, tet_vert_coorientations)
    for edge in tri.edges():
        uppermost0_index = sides[edge.index()][0][-1][0]
        where_uppermost[uppermost0_index] = edge.index()
        uppermost1_index = sides[edge.index()][1][-1][0]
        where_uppermost[uppermost1_index] = edge.index()
        
    for x in where_uppermost:
        assert x != None
        
    return where_uppermost

def pairing_of_faces_across_an_edge(tri, angle, weights, edge_index, tet_vert_coorientations = None): 
    """
    given a solution to the system of branch equations gives a pairing of faces across an edge (together with the index of the edge in that face)
    """

    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    
    edge_sides = edge_side_face_collections(tri, angle, tet_vert_coorientations = tet_vert_coorientations)[edge_index]
    
    side0 = edge_sides[0]
    for i in range(len(side0)):
        face_index = side0[i][0]
        weight = weights[face_index]
        side0[i] =  [side0[i], weight] #first replace side_collections into weighted_side_collections
    
    sum_of_weights0 = sum(weight for [_, weight] in side0)
    
    side0_nonzero_weights = []
    for x in side0:
        for k in range(x[1]):
            side0_nonzero_weights.append(x[0]) 
              
    side1 = edge_sides[1]
    for i in range(len(side1)):
        face_index = side1[i][0]
        weight = weights[face_index]
        side1[i] =  [side1[i], weight]
    
    sum_of_weights1 = sum(weight for [_, weight] in side1)
    assert sum_of_weights0 == sum_of_weights1
    
    side1_nonzero_weights = []
    for x in side1:
        for k in range(x[1]):
            side1_nonzero_weights.append(x[0])
    
    return([side0_nonzero_weights, side1_nonzero_weights])

def is_uppermost(tri, angle, face_index, edge_index_in_face, edge_index, tet_vert_coorientations = None):
    """
    tells if (face_index, edge_index_in_face) is uppermost for that edge
    """
    
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    
    edge_sides = edge_side_face_collections(tri, angle, tet_vert_coorientations = tet_vert_coorientations)[edge_index]
    side0 = edge_sides[0]
    side1 = edge_sides[1]
    
    return (side0[-1] == (face_index, edge_index_in_face) or side1[-1] == (face_index, edge_index_in_face))

def triangle_is_red(tri, angle, face_index, tet_vert_coorientations = None):
    
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
        
    edge_index_where_uppermost = where_faces_are_uppermost(tri, angle, tet_vert_coorientations)[face_index] #face has the same colour as the edge where it is uppermost
    
    colour = is_veering(tri, angle, return_type = 'veering_colours')[edge_index_where_uppermost]
    
    return (colour == 'red')
    

def appears_twice_on_the_same_side(tri, angle, face_index, edge_index, tet_vert_coorientations = None):
    """
    tells if a face is attached twice on the same side of an edge (if so, then lowermost and uppermost)
    """
    
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    
    edge_sides = edge_side_face_collections(tri, angle, tet_vert_coorientations = tet_vert_coorientations)[edge_index]
    side0 = edge_sides[0]
    side1 = edge_sides[1]
        
    return side0[0][0] == side0[-1][0] and side0[-1][0] == face_index or side1[0][0] == side1[-1][0] and side1[-1][0] == face_index
    
@liberal
def build_surface(tri, angle, weights, tet_vert_coorientations = None, return_edge_colours = False): 
    
    if tet_vert_coorientations == None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    
    surface = regina.Triangulation2()
    n = tri.countTetrahedra()
      
    for i in range(2*n):
        for j in range(weights[i]):
            surface.newTriangle()
    
    colours = is_veering(tri, angle, return_type = 'veering_colours')
    
    edge_colours_pairs = [] # red, blue or black (dual to large branches)
        
    for edge in tri.edges():
        pairings = pairing_of_faces_across_an_edge(tri, angle, weights, edge.index(), tet_vert_coorientations)
        height = len(pairings[0])
        #print('height', height)
        
        #now we get (face_index, index of e in the face) so have to pick the first one
        face_indices0 = [pairings[0][i][0] for i in range(height)]
        face_indices1 = [pairings[1][i][0] for i in range(height)]
        
        for i in range(height):
            #face0_index = pairings[0][i][0] #now we get (face_index, index of e in the face) so have to pick the first one
            face0_index = face_indices0[i]
            edge_index_in_face0 = pairings[0][i][1]
            #face1_index = pairings[1][i][0]
            face1_index = face_indices1[i]
            edge_index_in_face1 = pairings[1][i][1]
            
            which_copy0 = face_indices0[:i].count(face0_index)
            # check if not overcounting (if a face is both uppermost and lowermost on the same side)
            if is_uppermost(tri, angle, face0_index, edge_index_in_face0, edge.index(), tet_vert_coorientations) and appears_twice_on_the_same_side(tri, angle, face0_index, edge.index(), tet_vert_coorientations):
                which_copy0 = which_copy0 - weights[face0_index]
            
            which_copy1 = face_indices1[:i].count(face1_index)
            if is_uppermost(tri, angle, face1_index, edge_index_in_face1, edge.index(), tet_vert_coorientations) and appears_twice_on_the_same_side(tri, angle, face1_index, edge.index(), tet_vert_coorientations):
                which_copy1 = which_copy1 - weights[face1_index]

            #which_copy0 = face_indices0[:i].count(face0_index)
            #which_copy1 = face_indices1[:i].count(face1_index)
    
            which_triangle0 = sum(weights[:face0_index]) + which_copy0
            which_triangle1 = sum(weights[:face1_index]) + which_copy1
            
            #print('which triangle 0', which_triangle0)
            #print('which triangle 1', which_triangle1)
            
            if colours[edge.index()] == 'red': # gluing triangles across a red edge
                if is_uppermost(tri, angle, face0_index, edge_index_in_face0, edge.index(), tet_vert_coorientations):
                    if is_uppermost(tri, angle, face1_index, edge_index_in_face1, edge.index(), tet_vert_coorientations):
                        surface.triangle(which_triangle0).join(2, surface.triangle(which_triangle1), regina.Perm3(1,0,2))
                        #first we append edge_colours with tuples [face_index, edge_index_in_face, colour]; after finishing the triangulation we will change it into a tuple of colours ordered by indices
                        edge_colours_pairs.append([which_triangle0, 2, 'black'])
                    else: # middle and lowermost have the same permutation
                        surface.triangle(which_triangle0).join(2, surface.triangle(which_triangle1), regina.Perm3(0,2,1)) 
                        edge_colours_pairs.append([which_triangle0, 2, 'red'])
                else: # middle or lowermost
                    if is_uppermost(tri, angle, face1_index, edge_index_in_face1, edge.index(), tet_vert_coorientations):
                        surface.triangle(which_triangle0).join(1, surface.triangle(which_triangle1), regina.Perm3(0,2,1))
                        edge_colours_pairs.append([which_triangle0, 1, 'red'])
                    else:
                        surface.triangle(which_triangle0).join(1, surface.triangle(which_triangle1), regina.Perm3(2,1,0)) 
                        edge_colours_pairs.append([which_triangle0, 1, 'red'])
            else: # gluing triangles across a blue edge
                if is_uppermost(tri, angle, face0_index, edge_index_in_face0, edge.index(), tet_vert_coorientations):
                    if is_uppermost(tri, angle, face1_index, edge_index_in_face1, edge.index(), tet_vert_coorientations):
                        surface.triangle(which_triangle0).join(2, surface.triangle(which_triangle1), regina.Perm3(1,0,2))
                        edge_colours_pairs.append([which_triangle0, 2, 'black'])
                    else: # middle and lowermost have the same permutation
                        surface.triangle(which_triangle0).join(2, surface.triangle(which_triangle1), regina.Perm3(2,1,0)) 
                        edge_colours_pairs.append([which_triangle0, 2, 'blue'])
                else: # middle or lowermost
                    if is_uppermost(tri, angle, face1_index, edge_index_in_face1, edge.index(), tet_vert_coorientations):
                        surface.triangle(which_triangle0).join(0, surface.triangle(which_triangle1), regina.Perm3(2,1,0))
                        edge_colours_pairs.append([which_triangle0, 0, 'blue'])
                    else:
                        surface.triangle(which_triangle0).join(0, surface.triangle(which_triangle1), regina.Perm3(0,2,1)) 
                        edge_colours_pairs.append([which_triangle0, 0, 'blue'])
    
    assert surface.isClosed()
    
    #print (edge_colours_pairs)
    
    if return_edge_colours == True:
        k = surface.countEdges()
        edge_colours = [None]*k
        for x in edge_colours_pairs:
            edge_index = surface.triangle(x[0]).edge(x[1]).index()
            edge_colours[edge_index] = x[2]
        return surface, edge_colours
       
    return surface
                
def genus_punctures(tri, angle, weights, return_surface = False):
        
    real_Euler_char = -1/2*sum(weights)
    surface = build_surface(tri, angle, weights)
    regina_Euler_char = surface.eulerChar()
    # regina fills in ideal vertices
    genus = int(1/2*(2 - regina_Euler_char))
    boundaries = int(2 - 2*genus - real_Euler_char)
    
    if return_surface == True:
        return genus, boundaries, surface
        
    return genus, boundaries

def genus_punctures_from_weights_surface(weights, surface): #if the surface is already built
    real_Euler_char = -1/2*sum(weights)
    regina_Euler_char = surface.eulerChar()
    genus = int(1/2*(2 - regina_Euler_char))
    boundaries = int(2 - 2*genus - real_Euler_char)
    
    return genus, boundaries
    

def count_prongs(surface):
    # for a surface built using build_surface it is enough to count how many times a vertex appears as the vertex 2 of faces in which it is embedded
    out = []
    for v in surface.vertices():
        prongs = 0
        embeds = v.embeddings()
        for embed in embeds:
            if embed.vertex() == 2:
                prongs = prongs + 1
        out.append(prongs)
    return out

def stratum(tri, angle, weights, return_surface = False):
    
    genus, punctures, surface = genus_punctures(tri, angle, weights, return_surface =  True)
    prongs = count_prongs(surface)
    
    if return_surface == True:
        return (genus, punctures), prongs, surface

    return (genus, punctures), prongs

def stratum_from_weights_surface(weights, surface): # when the surface is already built
    
    genus, punctures = genus_punctures_from_weights_surface(weights, surface)
    prongs = count_prongs(surface)
    
    return (genus, punctures), prongs
    
# below: finding combinatorial isomorphisms of the surface which are orientation preserving and preserve the stable track    

def edge_image(isom, surface, edge_index):
    """
    returns the index of the edge image under a combinatorial isomorphism --- can't find regina doing it?
    """
    edge =  surface.edge(edge_index)
    index_of_face_where_edge_is_embedded = edge.embedding(0).simplex().index()
    #print(edge_index, 'is embedded in face', index_of_face_where_edge_is_embedded)
    index_of_edge_in_that_face = edge.embedding(0).vertices()[2]
    #print (edge.embedding(0), 'is edge', index_of_edge_in_that_face, 'of', index_of_face_where_edge_is_embedded)
    
    index_of_face_image = isom.simpImage(index_of_face_where_edge_is_embedded)
    vertex_perm = isom.facetPerm(index_of_face_where_edge_is_embedded)
    image_edge_index_in_image = vertex_perm[index_of_edge_in_that_face]
    edge_image = surface.triangle(index_of_face_image).edge(vertex_perm[index_of_edge_in_that_face])
    
    return edge_image.index()

def isom_preserves_colours(isom, surface, edge_colours): # note that this is not the same as preserving the stable track!
    for edge in surface.edges():
        colour = edge_colours[edge.index()]
        edge_image_index = edge_image(isom, surface, edge.index())
        if edge_colours[edge_image_index] != colour:
            return False
            print('Edge', edge.index(), 'is', colour, 'and is mapped to edge', edge_image_index, 'which is', edge_colours[edge_image_index])
    return True
        
def isom_preserves_stable_track(isom, surface, edge_colours, quiet = True):
    # to preserve the stable track the colours need to be preserved, and for every face the image of its large edge must be large in the face image
    if not isom_preserves_colours(isom, surface, edge_colours):
        return False 
    for triangle in surface.triangles(): 
        perm = isom.facetPerm(triangle.index())
        if perm[0] == 2 or perm[1] == 2:
            if quiet != True:
                print('The large edge of face', triangle.index(), 'is mapped to a small edge', perm[2], 'of face', isom.simpImage(triangle.index()), '\n(by our conventions, the large edge of a face always has index 2 in that face)')
            return False
    return True

def veering_symmetry_group(surface, edge_colours, return_isoms = True):
    """
    returns the group of orientation-preserving symmetries of the carried surface which preserve the stable train track
    """
    veering_isoms = []
    isoms = surface.findAllIsomorphisms(surface)
    for isom in isoms:
        if isom_preserves_stable_track(isom, surface, edge_colours, quiet = True):
            veering_isoms.append(isom)
    if return_isoms == True:
        return veering_isoms
    return len(veering_isoms)

        