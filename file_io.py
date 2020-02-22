#
# file_io.py
#


def parse_data_file(filename):
    """
    Parse a file into lines.
    """
    data_file = open(filename, "r") # read mode
    out = []
    for line in data_file:
        data_line = line.split("\n")[0] # remove '\n' from the end of each line, if it is there
        out.append(data_line)
    data_file.close()
    return out


def write_data_file(data, filename):
    """
    Compile lines into a file.
    """
    data_file = open(filename, "w") # write mode
    for line in data:
        assert type(line) == type([])
        assert all(type(item) == type("") for item in line)
        data_line = ' '.join(line) + "\n" # join together and add a carriage return.
        data_file.write(data_line)
    data_file.close()
