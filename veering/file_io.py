#
# file_io.py
#

import pickle

def parse_data_file(filename):
    """
    Parse a file into lines.
    """
    from pathlib import Path

    # 1. try relative path
    data_file = Path(filename)
    if not data_file.exists():
        # 2. try data_census path in the installation folder
        data_file = Path(__file__).parent / "data" / filename
        if not data_file.exists():
            raise ValueError('no data file {}'.format(filename))

    data_file = data_file.open()
    out = [line.strip() for line in data_file.readlines()]
    data_file.close()
    return out

def write_data_file(data, filename):
    """Compile lines into a file."""
    data_file = open(filename, 'w') # write
    for line in data:
        assert type(line) == type([])
        assert all(type(item) == type('') for item in line)
        data_line = ' '.join(line) + '\n' # join together and add a carriage return. 
        data_file.write(data_line)
    data_file.close()
    return None

def read_from_pickle(filename):                           
    f = open(filename, 'rb')
    data = pickle.load(f)
    f.close()
    return data

def output_to_pickle(out, filename):
    f = open(filename, 'wb')
    pickle.dump(out, f)
    f.close()
    return None

# a very common use case
def veering_census():
    return parse_data_file("veering_census.txt")
