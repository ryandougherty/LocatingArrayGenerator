import math

t = 2
v = 3
for d in range(1, 2+1):
	for k in range(10, 16+1):
		num_inters = math.comb(k,t) * v**t
		num_d_sets = math.comb(num_inters, d)
		num_pairs = math.comb(num_d_sets, 2)
		print(d, k, num_pairs)