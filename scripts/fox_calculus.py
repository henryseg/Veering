#
#  fox_calculus.py
#

from veering.taut import liberal, vert_pair_to_edge_num
from veering.veering_tri import veering_triangulation
from veering.transverse_taut import edge_side_face_collections, top_bottom_embeddings_of_faces, get_tet_top_vert_nums
from sage.groups.free_group import FreeGroup
from veering.fundamental_domain import spanning_dual_tree, non_tree_edge_loops_oriented


def is_AB_turn(vt, top_bottom_embeddings, face0, face1, face0_dir, face1_dir):
	"""We go through face0 in direction face0_dir (+1 if with the coorientation) 
	   into a tet t. We then leave through face1 in direction face1_dir. 
	   Return True if this turn is an AB turn. That is, the triangles are
	   adjacent along an equatorial edge of t of the same colour as the top
	   diagonal of the edge. See Lemma 5.6 of https://arxiv.org/pdf/2008.04836.pdf"""
	top, bottom = top_bottom_embeddings
	if face0_dir == 1:
		embed0 = bottom[face0]
	else:
		embed0 = top[face0]
	if face1_dir == -1:
		embed1 = bottom[face1]
	else:
		embed1 = top[face1]
	t0 = embed0.tetrahedron()
	t1 = embed1.tetrahedron()
	# print(t0.index(), t1.index())
	assert(t0 == t1)
	t = t0
	if face0_dir != face1_dir:
		return False
	f0 = embed0.face()
	f1 = embed1.face()
	equatorial_nums = [0,1,2,3]
	equatorial_nums.remove(f0)
	equatorial_nums.remove(f1)
	equatorial_colour = vt.get_edge_between_verts_colour(t.index(), equatorial_nums)

	top_vert_nums = get_tet_top_vert_nums(vt.coorientations, t.index())
	top_colour = vt.get_edge_between_verts_colour(t.index(), top_vert_nums)
	return top_colour == equatorial_colour

@liberal
def taut_polynomial_via_fox_calculus(tri, angle):
	vt = veering_triangulation(tri, angle)
	top_bottom_embeddings = top_bottom_embeddings_of_faces(tri, angle)
	tree = spanning_dual_tree(tri)
	tree_edges = tree[0]
	esfc = edge_side_face_collections(tri, angle)
	F = FreeGroup(tri.countTriangles())
	gens = F.gens()
	# new_gens = []
	# for i in range(range(len(gens))):
	# 	if i in tree_edges: 
	# 		new_gens.append(F.one())
	# 	else:
	# 		new_gens.append(gens[i])
	print(tree_edges)
	# gens = [F.one() if i in tree_edges else gens[i] for i in range(len(gens))]
	relators = []
	for edge in esfc:
		rel = F.one()
		left, right = edge
		left.reverse()
		for term in left:
			rel = rel * gens[term[0]].inverse()
		for term in right:
			rel = rel * gens[term[0]]  ## are these the correct way round?
		relators.append(rel)

	gens = F.gens()
	for f in tree_edges:
		relators.append(gens[f])
	G = F/relators
	print(G)

	loop_twistednesses = {}
	oriented_loops, all_signs = non_tree_edge_loops_oriented(tri, angle)
	for i in range(len(oriented_loops)):
		loop = oriented_loops[i]
		signs = all_signs[i]
		count = 0
		# print(loop, signs)
		for j in range(len(loop)):
			f0, f1 = loop[j], loop[(j+1)%len(loop)]
			f0d, f1d = signs[j], signs[(j+1)%len(loop)]
			if is_AB_turn(vt, top_bottom_embeddings, f0, f1, f0d, f1d):
				# print('isAB', f0, f1)
				count += 1
		loop_twistednesses[loop[0]] = (-1)**(count % 2) ### first in loop is the non tree edge
	# print(loop_twistednesses)
	for i in range(len(gens)):
		if i not in loop_twistednesses:
			loop_twistednesses[i] = 1
	
	loop_twistednesses = [loop_twistednesses[i] for i in range(len(gens))]
	print(loop_twistednesses)





