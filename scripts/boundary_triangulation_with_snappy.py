from snappy import Manifold

from veering.file_io import parse_data_file
from veering.taut import isosig_to_tri_angle

from boundary_triangulation import generate_boundary_triangulation

def draw_ladders_and_geometric_boundary_for_veering_isosig(sig, args = {}):
	if args == {}:
		args = {'draw_boundary_triangulation':True, 'draw_triangles_near_poles': False, 'ct_depth':-1, 'ct_epsilon':0.03, 'global_drawing_scale': 4, 'delta': 0.2, 'ladder_width': 10.0, 'ladder_height': 20.0, 'draw_labels': True}
	
	out_dir_ladders = 'Images/Ladders'
	out_dir_geometric = 'Images/Geometric'

	output_filename = sig + '.pdf'
	tri, angle = isosig_to_tri_angle(sig)

	M = Manifold(tri)
	tet_shapes = M.tetrahedra_shapes()
	tet_shapes = [complex(shape["rect"]) for shape in tet_shapes]
	args['tet_shapes'] = tet_shapes

	B = generate_boundary_triangulation(tri, angle, args = args, output_filename = output_filename)
	args_ladder = args.copy()
	args_ladder['style'] = 'ladders'
	output_filename_ladders = out_dir_ladders + '/' + output_filename
	B.draw(output_filename_ladders, args = args_ladder) 

	args_geometric = args.copy()
	args_geometric['style'] = 'geometric'
	output_filename_geometric = out_dir_geometric + '/' + output_filename
	B.draw(output_filename_geometric, args = args_geometric) 

def draw_ladders_and_geometric_boundary_for_veering_isosig_list(filename, args = {}):
	if args == {}:
		args = {'draw_boundary_triangulation':True, 'draw_triangles_near_poles': False, 'ct_depth':-1, 'ct_epsilon':0.03, 'global_drawing_scale': 4, 'delta': 0.2, 'ladder_width': 10.0, 'ladder_height': 20.0, 'draw_labels': True}
	sigs = parse_data_file(filename)
	for sig in sigs:
		draw_ladders_and_geometric_boundary_for_veering_isosig(sig, args = args)
	
