#
# veering_census_data.py
#

### Used to create or update the file "veering_census_with_data"

import regina
import snappy

from veering.file_io import parse_data_file, write_data_file
from veering.taut import isosig_to_tri_angle, apply_isom_to_angle_struct_list
from veering.veering_tri import is_veering
from veering.transverse_taut import is_transverse_taut, symmetry_group_size
from veering.edge_orientability import is_edge_orientable
from veering.taut_polytope import depth as tp_depth

from boundary_triangulation import generate_boundary_triangulation

def compute_census_data(filename_in, filename_out, functions, verbose = 0):
	### each function takes in data about the triangulation, returns a string
	census_data = parse_data_file(filename_in) 
	out = []
	for i, line in enumerate(census_data):
		line_data = line.split(' ') ## 0th is taut_sig, then may be other data we dont want to lose.
		taut_sig = line_data[0]
		regina_sig = taut_sig.split('_')[0]

		tri, angle = isosig_to_tri_angle(taut_sig)
		snappy_triang = snappy.Manifold(regina_sig)
		# snappy_triang = None
		triang_data = {'sig': taut_sig, 'angle': angle, 'tri': tri, 'snappy_triang': snappy_triang, 'old_data': line_data}

		line_out = []
		for func in functions:
			line_out.append(func(triang_data))

		out.append(line_out)

		if verbose > 0 and i % 1000 ==0:
			print(i, line_out)
			write_data_file(out, filename_out)
	write_data_file(out, filename_out)

def veering_isosig(triang_data):
	return triang_data['old_data'][0]

def num_cusps(triang_data):
	return str(triang_data['snappy_triang'].num_cusps())

def is_geometric(triang_data):
	if triang_data['snappy_triang'].verify_hyperbolicity()[0]:
		return 'G' # for geometric
	else:
		return 'N'

def num_toggles_red_blue(triang_data):
	tet_types = is_veering(triang_data['tri'], [int(a) for a in triang_data['angle']], return_type = 'tet_types')
	counts = [tet_types.count('toggle'), tet_types.count("red"), tet_types.count("blue")]
	return '[' + ','.join([str(i) for i in counts]) + ']'

def homology(triang_data):
	hom = triang_data['snappy_triang'].homology().__repr__()
	return hom.replace(' ','')

def other_names(triang_data):
	mflds = triang_data['snappy_triang'].identify()
	names = ["'" + mfld.name() + "'" for mfld in mflds]
	return '[' + ','.join(names) + ']'

def symmetry_group_size_data(triang_data):
	sig = triang_data['sig']
	count = symmetry_group_size(sig)
	return str(count)

def edge_orientable(triang_data):
	if is_edge_orientable(triang_data['sig']):
		return 'E'
	else:
		return 'N'

def ladder_counts(triang_data):
	b = generate_boundary_triangulation(triang_data['sig'], draw = False)
	counts = str(b.ladder_counts())
	return counts.replace(' ','')

def depth(triang_data):
	is_finite, cuts = tp_depth(triang_data['sig'])
	if is_finite:
		out = "F"
	else:
		out = "N"
	out = out + str(cuts)
	return out

# def LMN_from_old_data(triang_data): ## layered vs measurable vs non-measurable
# 	return triang_data['old_data'][1]

def depth_from_old_data(triang_data):
	return triang_data['old_data'][1]

def num_cusps_from_old_data(triang_data):
	return triang_data['old_data'][2]

def is_geometric_from_old_data(triang_data):
	return triang_data['old_data'][3]

def symmetry_group_size_from_old_data(triang_data):
	return triang_data['old_data'][4]

def edge_orientable_from_old_data(triang_data):
	return triang_data['old_data'][5]

def euler_class_from_old_data(triang_data):
	return triang_data['old_data'][6]

def ladder_counts_from_old_data(triang_data):
	return triang_data['old_data'][7]

def num_toggles_red_blue_from_old_data(triang_data):
	return triang_data['old_data'][8]

def homology_from_old_data(triang_data):
	return triang_data['old_data'][9]

def other_names_from_old_data(triang_data):
	return triang_data['old_data'][10]

functions_list = [veering_isosig, depth_from_old_data, num_cusps_from_old_data, is_geometric_from_old_data, symmetry_group_size_from_old_data, edge_orientable_from_old_data, euler_class_from_old_data, ladder_counts_from_old_data, num_toggles_red_blue_from_old_data, homology_from_old_data, other_names_from_old_data]

def recompute():
	compute_census_data('veering_census_with_data.txt', 'extra_data/veering_census_with_more_data.txt', functions_list, verbose = 1)

def search():
	census_data = parse_data_file('veering_census_with_data.txt') 
	out = []
	for i, line in enumerate(census_data):
		line_data = line.split(' ')
		if line_data[3] == 'N': # non geometric
			out.append([line_data[0]])
	write_data_file(out, 'extra_data/veering_non_geometric.txt')
