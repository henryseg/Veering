from file_io import parse_data_file, write_data_file


lines = parse_data_file('Data/veering_census_with_data.txt')

out = []
for line in lines:
    foo = line.split(" ")
    sig = foo[0]
    lmn = foo[1]
    num_cusps = int(foo[2])
    if lmn == 'L' and num_cusps == 1:
        out.append([sig])

write_data_file(out, 'layered_one_cusp.txt')