from file_io import parse_data_file, write_data_file

lines = parse_data_file('Data/veering_census_with_data.txt')

out = []
for line in lines:
    data = line.split(" ")
    sig = data[0]
    edge_ori = data[4]
    euler_ord = data[5]
    ladder_counts = data[6]
    ladder_counts = ladder_counts[1:-1]
    ladder_counts = ladder_counts.split(",")
    ladder_counts = [int(a) for a in ladder_counts]
    ladders_even = all(a % 4 == 0 for a in ladder_counts)
    if edge_ori == 'N' and ladders_even:
        # print(sig, edge_ori, ladder_counts)
        out.append([sig])

write_data_file(out, 'foo.txt') 
