import snappy
import regina

def are_cyclic_permutations(a, b):
	for i in range(len(a)):
		if a[i:] + a[:i] == b:
			return True
	return False

def flip(w):
	mytable = str.maketrans('LR', 'RL')
	return w.translate(mytable)

def is_in_list(word, word_list):
	# print('is in list', word, word_list)
	for w in word_list:
		if are_cyclic_permutations(word, w):
			# print('True')
			return True
	# print('False')
	return False

def torus_bundle_sigs_fixed_size(num_tets = 4):
	bundle_words = []
	for i in range(1, 2**num_tets - 1): ### don't do all L or all R
		b = bin(i)[2:]
		b = b.zfill(num_tets)
		# print(b)
		bundle_word = ""
		for letter in b:
			if letter == "0":
				bundle_word = bundle_word + "L"
			else:
				bundle_word = bundle_word + "R"
		if not is_in_list(bundle_word, bundle_words) and not is_in_list(flip(bundle_word), bundle_words):
			bundle_words.append(bundle_word)

	out = []
	for word in bundle_words:
		bundle_name = "b++" + word
		# print('bundle name', bundle_name)
		M = snappy.Manifold(bundle_name)
		T = regina.Triangulation3(M)
		out.append(T.isoSig())
	out.sort()
	return out

def torus_bundle_sigs(max_num_tets):
	out = []
	for i in range(2, max_num_tets + 1):
		out.extend(torus_bundle_sigs_fixed_size(num_tets = i))
	return out