import rhinoscriptsyntax as rs
import json

def read_from_json(filename):                           
    with open(filename) as f:
  		data = json.load(f)
    return data

def add_arc_two_points(v, w):
	u = rs.VectorUnitize(v+w)
	rs.AddArc3Pt(v, w, u)

path = read_from_json('/Users/segerman/GitHub/Veering/scripts/Images/Cannon-Thurston/Spectrum/cPcbbbiht_12_build_spherically_long_.03_float_[73256]_0.json')
path = [rs.CreateVector(v) for v in path]

n = len(path)
for i in range(n):
	v = path[i]
	w = path[(i+1)%n]
	add_arc_two_points(v, w)