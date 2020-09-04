#
# veering_census_data.py
#

### Used to create or update the file "veering_census_with_data"


import regina
# import snappy

from file_io import parse_data_file, write_data_file
from taut import isosig_to_tri_angle, apply_isom_to_angle_struct_list
from veering import is_veering
from transverse_taut import is_transverse_taut, symmetry_group_size
from edge_orientability import is_edge_orientable

def compute_census_data(filename_in, filename_out, functions, verbose = 0):
	### each function takes in data about the triangulation, returns a string
	census_data = parse_data_file(filename_in) 
	out = []
	for i, line in enumerate(census_data):
		line_data = line.split(' ') ## 0th is taut_sig, then may be other data we dont want to lose.
		taut_sig = line_data[0]
		regina_sig = taut_sig.split('_')[0]

		tri, angle = isosig_to_tri_angle(taut_sig)
		# snappy_triang = snappy.Manifold(regina_sig)
		snappy_triang = None
		triang_data = {'sig': taut_sig, 'angle': angle, 'tri': tri, 'snappy_triang': snappy_triang, 'old_data': line_data}

		line_out = []
		for func in functions:
			line_out.append(func(triang_data))

		out.append(line_out)

		if verbose > 0 and i % 1000 ==0:
			print(i, line_out)
	write_data_file(out, filename_out)

def veering_isosig(triang_data):
	return triang_data['old_data'][0]

def num_cusps(triang_data):
	return str(triang_data['snappy_triang'].num_cusps())

def num_toggles_red_blue(triang_data):
	tet_types = is_veering(triang_data['tri'], [int(a) for a in triang_data['angle']], return_type = 'tet_types')
	counts = [tet_types.count('toggle'), tet_types.count('right'), tet_types.count('left')]
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

def LMN_from_old_data(triang_data): ## layered vs measurable vs non-measurable
	return triang_data['old_data'][1]

def num_cusps_from_old_data(triang_data):
	return triang_data['old_data'][2]

def symmetry_group_size_from_old_data(triang_data):
	return triang_data['old_data'][3]

def edge_orientable_from_old_data(triang_data):
	return triang_data['old_data'][4]

def euler_class_from_old_data(triang_data):
	return triang_data['old_data'][5]

def num_toggles_red_blue_from_old_data(triang_data):
	return triang_data['old_data'][6]

def homology_from_old_data(triang_data):
	return triang_data['old_data'][7]

def other_names_from_old_data(triang_data):
	return triang_data['old_data'][8]

functions_list = [veering_isosig, LMN_from_old_data, num_cusps_from_old_data, symmetry_group_size_from_old_data, edge_orientable_from_old_data, euler_class_from_old_data, num_toggles_red_blue_from_old_data, homology_from_old_data, other_names_from_old_data]

if __name__ == '__main__':
	compute_census_data('Data/veering_census_with_data.txt', 'Data/veering_census_with_more_data.txt', functions_list, verbose = 1)
