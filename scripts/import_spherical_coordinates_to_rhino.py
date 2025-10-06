import rhinoscriptsyntax as rs
import json
from math import acos

def read_from_json(filename):                           
    with open(filename) as f:
  		data = json.load(f)
    return data

def add_arc_two_points(v, w):
	u = rs.VectorUnitize(v+w)
	rs.AddArc3Pt(v, w, u)

def transform_saul_to_henry_convention(v):
	refl = rs.XformMirror([0,0,0],[1,0,0])
	rot = rs.XformRotation2(60,[0,0,1],[0,0,0])
	v = rs.VectorTransform(v, refl)
	return rs.VectorTransform(v, rot)

def spherical_distance(v, w):
	return acos( sum([a*b for (a,b) in zip(v, w)]) )

def draw_path(filename):
	path = read_from_json(filename)
	path = [rs.CreateVector(v) for v in path]

	dists = []
	n = len(path)
	for i in range(n):
		v = path[i]
		w = path[(i+1)%n]
		dists.append(spherical_distance(v, w))
		add_arc_two_points(v, w)
	dists.sort()
	print(dists[-10:])

def draw_mesh(verts_filename, tris_filename):
	verts = read_from_json(verts_filename)
	verts = [rs.CreateVector(v) for v in verts]
	landscapes = read_from_json(tris_filename)
	upper_tris, lower_tris = landscapes
	rs.AddMesh( verts, upper_tris )
	rs.AddMesh( verts, lower_tris )

def draw_list_of_edges(filename):
	edges = read_from_json(filename)

	dists = []
	for edge in edges:
		edge = [transform_saul_to_henry_convention(rs.CreateVector(v)) for v in edge]
		v, w = edge
		dists.append(spherical_distance(v, w))
		add_arc_two_points(v, w)
	dists.sort()
	print(dists[-10:])

# draw_path('/Users/segerman/GitHub/Veering/scripts/Images/Cannon-Thurston/Spherical/cPcbbbiht_12_build_spherically_long_.1_float_[5564]_0_non_coast_factor_.5.json')
# 
# verts_filename = '/Users/segerman/GitHub/Veering/scripts/Images/Cannon-Thurston/Spherical/cPcbbbiht_12_build_spherically_long_.015_float_[1300026]_0.json'
# tris_filename = '/Users/segerman/GitHub/Veering/scripts/Images/Cannon-Thurston/Spherical/cPcbbbiht_12_build_spherically_long_.015_float_[1300026]_tris_0.json'
# draw_mesh(verts_filename, tris_filename)



verts_filename = '/Users/segerman/GitHub/Veering/scripts/Images/Cannon-Thurston/Spherical/gLLAQbecdfffhhnkqnc_120012_build_spherically_long_.03_float_[484218]_0.json'
tris_filename = '/Users/segerman/GitHub/Veering/scripts/Images/Cannon-Thurston/Spherical/gLLAQbecdfffhhnkqnc_120012_build_spherically_long_.03_float_[484218]_tris_0.json'
draw_mesh(verts_filename, tris_filename)

# draw_list_of_edges('/Users/segerman/Library/CloudStorage/Dropbox/Schleimer-Segerman/Figure_eight/Seifert_surface/bar_0.00265752_large_835774.json')



