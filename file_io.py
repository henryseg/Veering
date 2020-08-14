#
# file_io.py
#
import pickle

def parse_data_file(filename):
    """Parse a file into lines."""
    data_file = open(filename, 'r') # read mode
    out = []
    for line in data_file:
        data_line = line.split('\n')[0] # remove '\n' from the end of each line, if it is there
        out.append(data_line)
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
    f = open(filename, 'r')                                                                                                          
    data = pickle.load(f)                                                                                                           
    f.close()                                                                                                                        
    return data

def output_to_pickle(out, filename):
    f = open(filename, 'w')
    pickle.dump(out, f)
    f.close()
    return None