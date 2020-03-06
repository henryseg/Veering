###run this using regina-python to produce a .pkl file that rhino can read in and draw edges from

import cPickle
from taut import isosig_to_tri_angle
from veering_triangulation import veering_triangulation

import pyx ### vector graphics 

from develop_ideal_hyperbolic_tetrahedra import develop_vert_posns, convert_to_complex

def read_from_cpickle(filename):                                                                                                     
    f = open(filename, 'r')                                                                                                          
    data = cPickle.load(f)                                                                                                           
    f.close()                                                                                                                        
    return data

edge_vert_index_map = {(0,1):0, (0,2):1, (0,3):2, (1,2): 3, (1,3):4, (2,3):5}

class ct_edge():  ### cooriented edge is determined by the tetrahedron above it, and properties thereof
  def __init__(self, tet_num, forward_face_vertex, vert_posns, depth, vt):
    self.tet_num = tet_num
    self.forward_face_vertex = forward_face_vertex
    self.vert_posns = vert_posns
    self.depth = depth
    self.vt = vt

    self.top, self.bottom = get_top_bottom(self.vt, self.tet_num)

    self.end_complex = convert_to_complex(vert_posns[self.bottom[self.top.index(forward_face_vertex)]])
    self.start_complex = convert_to_complex(vert_posns[self.bottom[(self.top.index(forward_face_vertex)+1)%2]])

    self.length = abs(self.start_complex - self.end_complex)
    self.index = self.vt.tri.tetrahedron(self.tet_num).edge( edge_vert_index_map[ tuple(sorted(self.bottom)) ] ).index() 
    self.colour = self.vt.veering_colours[self.index]

  def __str__(self):
    return str(self.tet_num) + ',' + str(self.forward_face_vertex) + '|' + str(self.start_complex) + ',' + str(self.end_complex) + '|' + str(self.depth) + ',' + self.colour
  def __repr__(self):
    return str(self.tet_num) + ',' + str(self.forward_face_vertex) + '|' + str(self.start_complex) + ',' + str(self.end_complex) + '|' + str(self.depth) + ',' + self.colour

  def too_close_to_infty(self, max_size = 4.0):
    return abs(self.start_complex) > max_size and abs(self.end_complex) > max_size

  def develop_edge_outwards(self, depth_increment = 1, verbose = 0.0):
    """Given edge_data, find list of edges (i.e. edge_data) on far side of this edge"""
    if verbose >= 1.0: 
      print 'develop through', self
    ### if bottom edge is right veering then as we go anticlockwise around the circle, the 
    ### the new edges we get will be in the correct order as we work down from the top of the 
    ### edge. if bottom edge is left veering then we have the opposite. 
    vt = self.vt
    tet_num = self.tet_num
    forward_face_vertex = self.forward_face_vertex
    vert_posns = self.vert_posns
    depth = self.depth
    colour = self.colour

    top, bottom = self.top, self.bottom

  ###      *-----top[0] 
  ###       \`-. / | \
  ###        \  /-.|  \
  ###         \/   |`-.\
  ### bottom[0]--- | ---bottom[1]  <- red bottom edge, e
  ###          \   |   /
  ###           \  |  /
  ###           top[1] = forward_face_vertex

    ### assuming red lower edge...
    back_edge_vertex = bottom[ top.index(forward_face_vertex) ]  ### the vertex of edge e that is behind us as we go around the edge
    front_edge_vertex = bottom[ (top.index(forward_face_vertex)+1)%2 ] ### the vertex of edge e that is in front of us as we go around the edge
    if colour != 'R':
      back_edge_vertex, front_edge_vertex = front_edge_vertex, back_edge_vertex
    previous_face_vertex = top[ (top.index(forward_face_vertex)+1)%2 ]

  ###      *-----prev 
  ###       \`-. / | \
  ###        \  /-.|  \
  ###         \/   |`-.\
  ###    front --- | --- back  <- red bottom edge, e
  ###          \   |   /
  ###           \  |  /
  ###            next

    current_tet = vt.tri.tetrahedron(tet_num)
    next_face_vertex = forward_face_vertex # forward in terms of the edge coorientation is the same as going around the edge here, but not for any other tetrahedra
    
    # print 'tet_num going down', tet_num
    # print 'f,b,p,n', front_edge_vertex, back_edge_vertex, previous_face_vertex, next_face_vertex

    out = []
    ### new edge connects prev-back for the top tet, prev-next for others, and next-front for the bottom tet
    out.append( get_ct_edge_above( current_tet, vert_posns, next_face_vertex, front_edge_vertex, self, depth_increment = depth_increment, verbose = verbose ) )
    while True:
      next_tet = current_tet.adjacentTetrahedron(next_face_vertex)
      gluing = current_tet.adjacentGluing(next_face_vertex)
      # print 'next tet index', next_tet.index()
      vert_posns = develop_vert_posns(vert_posns, gluing, next_face_vertex, tet_shapes[next_tet.index()])
      back_edge_vertex = gluing[back_edge_vertex]
      front_edge_vertex = gluing[front_edge_vertex]
      next_face_vertex, previous_face_vertex = gluing[previous_face_vertex], gluing[next_face_vertex]

      current_tet = next_tet

      # print 'tet_num going down', current_tet.index()
      # print 'f,b,p,n', front_edge_vertex, back_edge_vertex, previous_face_vertex, next_face_vertex
      # print 'coor', coorientations[next_tet.index()][next_face_vertex]

      ###            next ------*
      ###            / | \ _,-'/
      ###           /  _,-'   /
      ###          /,-'|   \ /
      ###    front --------- back  <- red bottom edge, e
      ###          \   |   /
      ###           \  |  /
      ###            prev

      assert vt.coorientations[current_tet.index()][previous_face_vertex] == +1
      if vt.coorientations[current_tet.index()][next_face_vertex] != +1: # keep going around the edge
        out.append( get_ct_edge_above( current_tet, vert_posns, front_edge_vertex, back_edge_vertex, self, depth_increment = depth_increment, verbose = verbose ) )
      else: ### we are now at the bottom of the edge
        out.append( get_ct_edge_above( current_tet, vert_posns, back_edge_vertex, previous_face_vertex, self, depth_increment = depth_increment, verbose = verbose ) )
        if colour != 'R': # we meet the edges in the opposite order when going anticlockwise
          out.reverse()  ### we want these edges to be in reverse order for popping off the to_do list
        return out

def get_ct_edge_above(current_tet, vert_posns, edge_vertex, face_vertex, old_ct_edge, depth_increment = 1, verbose = 0.0):
  """move up from current_tet around the edge, developing as we go, find the edge data for this edge (meaning data for tet above the edge)."""
  ### face vertex gives a face of the tet, edge vertex gives the edge of that face that we rotate around, starting to move away from face_vertex
  coorientations = old_ct_edge.vt.coorientations
  tet_shapes = old_ct_edge.vt.tet_shapes
  depth = old_ct_edge.depth + depth_increment
  assert coorientations[current_tet.index()][face_vertex] == +1
  if verbose > 3.0:
    print 'tet_num, edge_vertex, face_vertex', current_tet.index(), edge_vertex, face_vertex
    print 'vert_posns', vert_posns
  while True:
    # print 'tet_num going up', current_tet.index()
    # print '1,2,3,f', first_edge_vertex, second_edge_vertex, third_vertex, face_vertex
    next_tet = current_tet.adjacentTetrahedron(face_vertex)
    gluing = current_tet.adjacentGluing(face_vertex)
    vert_posns = develop_vert_posns(vert_posns, gluing, face_vertex, tet_shapes[next_tet.index()])

    # first_edge_vertex = gluing[first_edge_vertex]
    # second_edge_vertex = gluing[second_edge_vertex]
    edge_vertex, face_vertex = gluing[face_vertex], gluing[edge_vertex]

    current_tet = next_tet
    if verbose > 3.0:
      print 'tet_num, edge_vertex, face_vertex', current_tet.index(), edge_vertex, face_vertex
      print 'vert_posns', vert_posns
    assert coorientations[current_tet.index()][edge_vertex] == -1
    if coorientations[current_tet.index()][face_vertex] == -1: ### we are now at the top of the edge
      # print 'tet_num going up', current_tet.index()
      # print '1,2,3,f', first_edge_vertex, second_edge_vertex, edge_vertex, face_vertex
      # print 'done going up'
      return ct_edge(current_tet.index(), face_vertex, vert_posns, depth, old_ct_edge.vt)

def get_top_bottom(vt, tet_num): ## modified from draw_veering_mid-annuli2
  ### this says which vertices are 
  ### top and bottom, and nails down the orientation of the tetrahedron in terms of the colouring
  coors = vt.coorientations[tet_num]
  tri = vt.tri
  veering_colours = vt.veering_colours
  top_vertices = []
  bottom_vertices = []
  #coor for a tet looks like [1, -1, 1, -1], is 1 for pointing out of tet, -1 for in
  for i in range(4):
    if coors[i] == 1:
      bottom_vertices.append(i) 
      # Assume coorientation is upwards. If == 1, face coorientation is outwards, 
      # means this face is on top, means corresponding vertex is on bottom
    else:
      top_vertices.append(i) 
  edge_pair = [top_vertices[0], bottom_vertices[0]]
  edge_pair.sort()
  edge_num = tri.tetrahedron(tet_num).edge(edge_vert_index_map[tuple(edge_pair)]).index()
  col = veering_colours[edge_num]
  if col == 'L':
    # print 'flip orientation'
    bottom_vertices = [bottom_vertices[1], bottom_vertices[0]]
  return (top_vertices, bottom_vertices)

###           top[0]
###          /   |   \
### bottom[0]--- | ---bottom[1]
###          \   |   /
###           top[1]

###              fwd_face
###              ,' | `.
### out_vert_num -- | --*
###              `. | ,'
###                 *

### We want to build an approximation to the Cannon-Thurston map for a veering triangulation, 
### i.e. mapping from S^1 to \bdy H^3
### Our map will be polygonal, each edge corresponding to an edge of the universal cover of the triangulation
### Each edge will be given by the tetrahedron above it (relative to the transverse taut structure)
### together with one of the two other vertices, which gives a horizontal coorientation to the edge

### We develop out from an edge by replacing it with all zero angle edges of the tetrahedra
### on the side of the edge pointed to by the coorientation

def draw_path(edges, dots, name = 'foobar', lw = 0.005, verbose = 0.0):
  canv = pyx.canvas.canvas() 
  for dot in dots:
    canv.fill(pyx.path.circle(dot.real, dot.imag, 0.02))
  comp_coords = [edges[0].start_complex] + [edge.end_complex for edge in edges]
  if verbose > 0.1:
    print 'len(comp_coords)', len(comp_coords)
  p = pyx.path.path( pyx.path.moveto(comp_coords[0].real, comp_coords[0].imag) )
  for coord in comp_coords[1:]: 
    p.append( pyx.path.lineto(coord.real, coord.imag) )
  canv.stroke(p, [pyx.style.linewidth(lw)])
  name = name + '.pdf'
  canv.writePDFfile(name)

def draw_edges(edges, dots, name = 'foobar', lw = 0.005, verbose = 0.0):
  canv = pyx.canvas.canvas() 
  for dot in dots:
    canv.fill(pyx.path.circle(dot.real, dot.imag, 0.02))
  for edge in edges:
    s, e = edge.start_complex, edge.end_complex
    canv.stroke(pyx.path.line(s.real, s.imag, e.real, e.imag), [pyx.style.linewidth(lw)])
  name = name + '.pdf'
  canv.writePDFfile(name)

def initial_path_single(tri, angle, tet_shapes, verbose = 0.0):
  vt = veering_triangulation(tri, angle, tet_shapes = tet_shapes)
  lower_faces = [i for i in range(4) if vt.coorientations[0][i] == -1]
  vert_posns = [[1,0],[0,1],[1,1],[tet_shapes[0],1]]
  return ( [ct_edge(0,lower_faces[1],vert_posns,0, vt)], [] )

def initial_path_sideways(tri, angle, tet_shapes, verbose = 0.0):
  edges, _ = initial_path_single(tri, angle, tet_shapes, verbose = 0.0)
  current_edge = edges[0]
  rho = []
  oriented_edges = []
  infinity_edge_endpoints = []
  while True:
    assert current_edge.end_complex == infty
    infinity_edge_endpoints.append(current_edge.start_complex)
    oriented_e = (current_edge.tet_num, current_edge.forward_face_vertex)
    if oriented_e in oriented_edges:
      end_ori_e = oriented_e
      break
    oriented_edges.append( oriented_e )
    if verbose > 2.0: print 'oriented edges', oriented_edges
    
    children = current_edge.develop_edge_outwards(depth_increment = 0, verbose = verbose)
    children.reverse() ### we want the opposite order in the main develop function
    if verbose > 2.0:
      print 'children'
      for child in children:
        print child.tet_num, child.forward_face_vertex, child.start_complex, child.end_complex
    current_edge = children.pop() 
    rho.append(children)
  if verbose > 2.0: print 'end_ori_e, index in oriented edges', end_ori_e, oriented_edges.index(end_ori_e)
  rho = rho[oriented_edges.index(end_ori_e):] ### throw away the initial segment of the rho, leaving the loop
  if verbose > 2.0: 
    print 'rho'
    for c in rho: print c
  # out = [edge for edge in children for children in rho] ## doesn't concat the right way
  out = []
  for c in rho:
    out.extend(c)
  print out
  return ( out, infinity_edge_endpoints )

def initial_path_one_ladderpole(tri, angle, tet_shapes, verbose = 0.0):
  edges, _ = initial_path_single(tri, angle, tet_shapes, verbose = 0.0)
  current_edge = edges[0]
  children = current_edge.develop_edge_outwards(depth_increment = 0, verbose = verbose)
  e = children[-1]
  return ( [e], [e.start_complex, e.end_complex] )

def initial_path_up_ladderpole(tri, angle, tet_shapes, verbose = 0.0):
  edges, _ = initial_path_single(tri, angle, tet_shapes, verbose = 0.0)
  current_edge = edges[0]
  if verbose > 2.0: 
    print 'current_edge', current_edge
    print 'top, bottom', current_edge.top, current_edge.bottom
  oriented_edges = []
  infinity_edge_endpoints = []
  ladderpole_edges = []
  while True:
    # assert current_edge.end_complex == infty
    # infinity_edge_endpoints.append(current_edge.start_complex)
    infinity_edge_endpoints.append(current_edge.end_complex)
    oriented_e = (current_edge.tet_num, current_edge.forward_face_vertex)
    if oriented_e in oriented_edges:
      end_ori_e = oriented_e
      break
    oriented_edges.append( oriented_e )
    if verbose > 2.0: print 'oriented edges', oriented_edges
    vt = current_edge.vt
    current_tet = vt.tri.tetrahedron(current_edge.tet_num)
###              
###                ,*.
###           B r,' | `.7 C
###            ,'   |   `.
###  new face *---^-|-A---*  new_face_vert for edge B
###  vert for  `.   |   ,'
###  edge C      `. | ,'
###                `*'
###               face_vert = new_edge_vert
    edge_vertex = current_edge.forward_face_vertex
    if verbose > 2.0: print 'edge vertex', edge_vertex
    for vertex in current_edge.bottom:
      if verbose > 2.0: print 'bottom vertex', vertex
      ct_edge = get_ct_edge_above(current_tet, current_edge.vert_posns, edge_vertex, vertex, current_edge, depth_increment = 0, verbose = verbose) 
      if verbose > 2.0: print 'ct edge', ct_edge
      if ct_edge.colour == current_edge.colour:
        next_edge = ct_edge
      else:
        ladderpole_edges.append(ct_edge)
    current_edge = next_edge
  return (ladderpole_edges, infinity_edge_endpoints)

def make_cannon_thurston(tri, angle, tet_shapes, init_function = initial_path_single, name = 'foobar', max_depth = 1, epsilon = 0.02, lw = 0.005, verbose = 0.0):
  to_do_edges, infinity_edge_endpoints = init_function(tri, angle, tet_shapes, verbose = verbose)
  if verbose > 0.0: print 'len initial to_do_edges', len(to_do_edges)
  to_do_edges.reverse() 
  if verbose > 4.5: 
    print 'initial to_do_edges'
    for e in to_do_edges:
      print e 
  done_edges = []
  while len(to_do_edges) > 0:
    if verbose > 2.0:
      print 'to do edges', len(to_do_edges)
      print 'done edges', len(done_edges)
    edge = to_do_edges.pop()
    if verbose >= 2.0:
      print 'depth', edge.depth
    if edge.depth >= max_depth or edge.length < epsilon or edge.too_close_to_infty():
      done_edges.append(edge)
    else:
      to_do_edges.extend( edge.develop_edge_outwards(verbose = verbose) )
  if verbose > 0.0: print 'len final done_edges', len(done_edges)
  # draw_path(done_edges, infinity_edge_endpoints, name = name, lw = lw, verbose = verbose)
  draw_edges(done_edges, infinity_edge_endpoints, name = name, lw = lw, verbose = verbose)


if __name__ == '__main__':
  # tet_shapes = [complex(0.5,math.sqrt(3)*0.5), complex(0.5,math.sqrt(3)*0.5)]
  # tri, angle = isosig_to_tri_angle('cPcbbbiht_12')

  # tet_shapes = [(0.9231422157742049+0.9853733543905342j), (-0.29947574053100134+0.44832432157891566j), (0.07867777692993604+1.0087070002132608j), (0.07867777692993473+1.0087070002132614j), (0.7378949835259692+0.13340505745111064j), (-0.29947574053100123+0.4483243215789155j)]
  # tri, angle = isosig_to_tri_angle('gLLAQbecdfffhhnkqnc_120012')
  # make_cannon_thurston(tri, angle, tet_shapes, max_depth = 30, epsilon = 0.02)

  data = read_from_cpickle('Data/shapes_for_cusped_374.pkl')
  names = data.keys()
  names.sort()
  # for name in names:
  name = 'gLLAQbecdfffhhnkqnc_120012'
  # name = names[1]
  print name
  tri, angle = isosig_to_tri_angle(name)
  tet_shapes = data[name]
  # print 'tet_shapes', tet_shapes
  # make_cannon_thurston(tri, angle, tet_shapes, name = 'Cannon-Thurston_images/' + name, max_depth = 25, epsilon = 0.05, verbose = 0.0)

  # sideways = initial_path_sideways(tri, angle, tet_shapes)
  # draw_path(sideways, name = name + '_sideways', verbose = 0.0)

  # make_cannon_thurston(tri, angle, tet_shapes, init_function = initial_path_sideways, name = name + '_sideways', max_depth = 0, epsilon = 0.05, verbose = 5.0)
  # make_cannon_thurston(tri, angle, tet_shapes, name = name, max_depth = 5, epsilon = 0.05, verbose = 5.0)
  # make_cannon_thurston(tri, angle, tet_shapes, init_function = initial_path_one_ladderpole, name = name + '_one_ladderpole', max_depth = 25, epsilon = 0.005, lw = 0.002, verbose = 0.0)
  make_cannon_thurston(tri, angle, tet_shapes, init_function = initial_path_up_ladderpole, name = name + '_up_ladderpole', max_depth = 8, epsilon = 0.005, lw = 0.002, verbose = 0.0)

  # vt = veering_triangulation(tri, angle, tet_shapes = tet_shapes)
  # for i in range(tri.countTetrahedra()):
  #   print get_top_bottom(vt, i)


