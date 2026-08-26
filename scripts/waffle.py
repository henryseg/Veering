from regina import Triangulation3
from veering.taut import isosig_to_tri_angle, liberal
from veering.transverse_taut import is_transverse_taut, get_top_and_bottom_vert_nums
from veering.veering_tri import veering_triangulation

def global_tet_index(tet_num, j, k, n = 2):
	return j + k * n + tet_num * n * n 

	###   t0------b1   04  14  24  34  44    n = 5, 25 tetrahedra
	###    |`.  ,'|    03  13  23  33  43 
	###    |  `.  |    02  12  22  32  42
	###    |,'  `.|    01  11  21  31  41
	###   b0------t1   00  10  20  30  40   

	### Barycentric coords: entries in vector sum to n-1. 

	###   t0------b1   040                  barycentric coordinates for SW tri   (order b0, t0, t1)
	###    |`.  ,'|    130 031          
	###    |  `.  |    220 121 022
	###    |,'  `.|    310 211 112 013
	###   b0------t1   400 301 202 103 004  (j, k) -> (n-1-j-k, k, j) 

def SW(j, k, n = 2):
	return (n-1-j-k, k, j) 

    ###   t0------b1   040 130 220 310 400  barycentric coordinates for NE tri   (order b1, t0, t1)
	###    |`.  ,'|        031 121 211 301    
	###    |  `.  |            022 112 202
	###    |,'  `.|                013 103
	###   b0------t1                   004  (i, j) -> (-n+1+j+k, n-1-j, n-1-k)   

def NE(i, j, n = 2):
	return (-n+1+j+k, n-1-j, n-1-k)

    ###   t0------b1   004 013 022 031 040  barycentric coordinates for NW tri   (order b0, b1, t0)
	###    |`.  ,'|    103 112 121 130
	###    |  `.  |    202 211 220    
	###    |,'  `.|    301 310        
	###   b0------t1   400                  (i, j) -> (n-1-k, j, k-i). inverse (u, v, w) -> (v, v+w)   

def NW_inv(triplet):
	u, v, w = triplet
	return (v, v+w)

	###   t0------b1                   040  barycentric coordinates for SE tri   (order b0, b1, t1)
	###    |`.  ,'|                130 031
	###    |  `.  |            220 121 022
	###    |,'  `.|        310 211 112 013
	###   b0------t1   400 301 202 103 004  (i, j) -> (n-1-j, k, j-k). inverse (u, v, w) -> (v+w, v) 

def SE_inv(triplet):
	u, v, w = triplet
	return (v+w, v)  

def get_new_coords(images, barycentric_coords, t0a, t1a):
	pairs = zip(images, barycentric_coords)
	pairs.sort(key = lambda x : x[0])  ###  sort on 0th coord
	new_barycentric_coords = [pair[1] for pair in pairs]
	if not t0a in images:
		return SE_inv(new_barycentric_coords)
	else:
		assert not t1a in images:
		return NW_inv(new_barycentric_coords)

@liberal
def waffle_press(tri, angle, n = 2):
	vt = veering_triangulation(tri, angle)

	new_tri = Triangulation3()
	new_angle = []

	for i in range(tri.countTetrahedra()):
		for j in range(n*n):
			new_tri.newSimplex()
			new_angle.append(angle[i]) 
			### each new small tetrahedra inherits vert numbers (and so angle) from the original

	top_bot_verts = []
	for i in range(tri.countTetrahedra()):
		(t0,t1), (b0,b1) = get_top_and_bottom_vert_nums(vt.coorientations, i) 
		### The top vertices come in arbitrary order, as do the bottom vertices.
		### We need to normalise to fix orientation. So let's ensure that the
		### edge from b0 to t0 is red.
		if vt.get_edge_between_verts_colour(i, (b0, t0)) != "red":
			(b0, b1) = (b1, b0)
		assert vt.get_edge_between_verts_colour(i, (b0, t0)) == "red"
		top_bot_verts.append(((t0, t1), (b0, b1)))

	for i in range(tri.countTetrahedra()):
		tet = tri.tetrahedron(i)
		b0_gluing = tet.adjacentGluing(b0)
        b0_adj_tet_num = tet.adjacentTetrahedron(b0).index()
        b1_gluing = tet.adjacentGluing(b1)
        b1_adj_tet_num = tet.adjacentTetrahedron(b1).index()

		###   t0------b1   04  14  24  34  44    n = 5, 25 tetrahedra
		###    |`.  ,'|    03  13  23  33  43 
		###    |  `.  |    02  12  22  32  42
		###    |,'  `.|    01  11  21  31  41
		###   b0------t1   00  10  20  30  40   

		for j in range(n):
			for k in range(n):
				this_tet_index = global_tet_index(i, j, k)
				if j + k <= n - 1:  ### glue SW upper face using the b1 face of original tet
					images = b1_gluing[b0], b1_gluing[t0], b1_gluing[t1]
					(t0a, t1a), (b0a, b1a) = top_bot_verts[b1_adj_tet_num]
					barycentric_coords = SW(i, j, n=n)
					new_j, new_k = get_new_coords(images, barycentric_coords, t0a, t1a)
					other_tet_index = global_tet_index(b1_adj_tet_num, new_j, new_k)
					new_tri.tetrahedron(this_tet_index).join(b1, new_tri.tetrahedron(other_tet_index), b1_gluing)

				else:               ### glue SW upper face using the b0 face of original tet
					images = b0_gluing[b1], b0_gluing[t0], b0_gluing[t1]
					(t0a, t1a), (b0a, b1a) = top_bot_verts[b0_adj_tet_num]
					barycentric_coords = NE(i, j, n=n)
					new_j, new_k = get_new_coords(images, barycentric_coords, t0a, t1a)
					other_tet_index = global_tet_index(b0_adj_tet_num, new_j, new_k)
					new_tri.tetrahedron(this_tet_index).join(b0, new_tri.tetrahedron(other_tet_index), b0_gluing)

				if j + k < n - 1:   ### glue NE upper face using the b1 face of original tet

				else:				### glue NE upper face using the b0 face of original tet






