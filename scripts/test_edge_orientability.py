from veering.file_io import parse_data_file, read_from_pickle

import veering.taut
import veering.transverse_taut
from veering.edge_orientability import is_edge_orientable

veering_isosigs = parse_data_file("Data/veering_census.txt")

for sig in veering_isosigs[:87]:  # up to 12 tetrahedra
    print((sig + ' ' + str(is_edge_orientable(sig))))
