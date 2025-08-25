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

# draw_path('/Users/segerman/GitHub/Veering/scripts/Images/Cannon-Thurston/Spherical/cPcbbbiht_12_build_spherically_long_.1_float_[20226]_0.json')
# 
draw_list_of_edges('/Users/segerman/Library/CloudStorage/Dropbox/Schleimer-Segerman/Figure_eight/Seifert_surface/foo_scale_5.json')