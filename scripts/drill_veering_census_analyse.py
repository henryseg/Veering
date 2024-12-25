
import json
import ast

class blob:
    def __init__(self, all_triangs, orig_triangs, orig_triang_indices):
    	self.all_triangs = all_triangs
    	self.orig_triangs = orig_triangs
    	self.orig_triang_indices = orig_triang_indices

    def isdisjoint(self, other):
    	return self.all_triangs.isdisjoint(other.all_triangs)

    def update(self, other):
    	self.all_triangs.update(other.all_triangs)
    	self.orig_triangs.extend(other.orig_triangs)
    	self.orig_triang_indices.extend(other.orig_triang_indices)
    	self.orig_triangs.sort()
    	self.orig_triang_indices.sort()

def parse_line(line):
	line = line.strip()
	first = line.index("_")
	second = line.index("_", first + 1) # start searching immediately after the first one
	undrilled_sig = line[:second]
	drill_data = line[second + 1:]
	drill_data = ast.literal_eval(drill_data)
	data = set([d[0] for d in drill_data])
	data.add(undrilled_sig)
	return data, undrilled_sig

def run(filename):
	with open(filename, 'r') as file:
		equiv_classes = []
		count = 0
		max_count = 2000

		while True:
			line = file.readline()
			
			if len(line) == 0 or count == max_count:
				break ### EOF
			all_triangs, undrilled_sig = parse_line(line)
			b = blob(all_triangs, [undrilled_sig], [count])
			overlaps = [e for e in equiv_classes if not e.isdisjoint(b)]
			
			for e in overlaps:
				b.update(e)
				equiv_classes.remove(e)
			equiv_classes.append(b)
			count += 1
			# print(count, len(equiv_classes))
			
		equiv_classes.sort(key=lambda x: (len(x.orig_triang_indices), -min(x.orig_triang_indices)))
		print('count', count, 'num equiv classes', len(equiv_classes))
		for e in equiv_classes:
			# if len(e.orig_triang_indices) > 1:
			print(len(e.all_triangs), len(e.orig_triang_indices), min(e.orig_triang_indices))
		




def main():
	# testline = "cPcbbbiht_12_[('gLLPQccdfeffhggaagb_201022', ((0, 0), (0, 5))), ('hLLPMkccdfeggghggaahah_2010221', ((0, 0), (0, 0), (0, 5))), ('iLLLQPcceegfghhhiimaimimi_10221212', ((0, 4), (1, 2))), ('iLLPMPcccdfeghhhhggaahhgb_20102211', ((0, 0), (0, 0), (0, 0), (0, 5))), ('jLLPMzQccdfeghiiihggaahhhah_201022112', ((0, 0), (0, 0), (0, 0), (0, 0), (0, 5))), ('jLLwQLQbeefgehiiixxxaaxxxcv_102211100', ((0, 0), (0, 0), (0, 5), (0, 0), (0, 5))), ('kLLLAzQkceegfihjjijiimaiigmiac_1022121021', ((0, 0), (0, 4), (1, 2))), ('kLLvMMQkcdfhfiijhjjhsawforeneq_1222011101', ((0, 0), (0, 0), (0, 5), (0, 5))), ('lLLLLPPQcbdefjghikkjkhruglokrmmpw_12000121111', ((0, 0), (0, 0), (0, 0), (0, 5), (0, 5))), ('mLLLAzPPQceegfijkillkliimaiimmaocog_102212102102', ((0, 0), (0, 0), (0, 4), (1, 2))), ('mLvLLLQQQcghgjillkjjkldeabbbaeacfef_120000112222', ((0, 0), (0, 4), (1, 0), (1, 2))), ('nLLvwLQPQkcegkijimilklmmiiiaccceftilwj_0112222011100', ((0, 0), (0, 4), (1, 2), (0, 5))), ('oLLLAzzPPQcceegfijimklnnmniimaiimaoogmioi_10221210210222', ((0, 0), (0, 0), (0, 0), (0, 4), (1, 2))), ('oLLvLQLMMQccegjiiljkkmnnmndqawwvqcildwvoo_12000112222210', ((0, 0), (0, 0), (0, 4), (1, 2), (0, 5))), ('oLvLLLLQQQccghgllkmnjknmlndeabaeaedfeddfc_12000011222211', ((0, 0), (0, 0), (0, 4), (1, 0), (1, 2))), ('ovLAMLvPQQccdfghgnmllmnlnmvoapaoppaappooo_01122000112211', ((0, 0), (0, 4), (1, 2), (0, 0), (0, 5))), ('pLLvvzLQQQQcehkmjjolnmomnnoiirkxkfptiadimww_011222200111100', ((0, 0), (0, 4), (1, 0), (1, 2), (0, 5))), ('qLLvwLALAQQkcegkijilmlponpoopiiiaccftmlokjcwwk_0112222011100012', ((0, 0), (0, 4), (1, 0), (1, 5), (1, 2))), ('tvLLvLQLvQQPQkegljmklmnlqopspsrrrsoaajojbfpouhdkoksbvg_0221110110000112021', ((0, 0), (0, 4), (1, 2), (0, 4), (1, 2)))]"
	run('/Users/segerman/Dropbox/Data/drillings_census_4_ladders_max_cycle_len_5.txt')
