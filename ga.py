import copy

import sys
import re
import random

def parse_dset(interaction):
	cols,vals = interaction.split('),')
	cols = list(map(int,cols.replace('(','').split(',')))
	vals = list(map(int,vals.replace('(','').replace(')', '').split(',')))
	return cols, vals

def read_file(d,t,k,v,l):
	filename = f'./outputs/d{d}_{v}^{k}-t{t}_l{l}.txt'
	with open(filename) as f:
		lines = f.readlines()[1:] # skip the runtime
	d_sets = []
	for line in lines:
		line = line.replace(',)', ')')
		s = line.split(' ')
		all_the_dsets = s[:-1]
		num_times_sep_already = int(s[2])
		parsed_dsets = []
		for i in range(0, len(all_the_dsets), 2):
			dset = [parse_dset(all_the_dsets[i]), parse_dset(all_the_dsets[i+1])]
			parsed_dsets.append(dset)
		d_sets.append((parsed_dsets,num_times_sep_already))
	return d_sets

def random_individual(N):
	return [[random.choice(range(v)) for _ in range(k)] for i in range(N)]

def rows_of_interaction(individual, interaction):
	cols,vals = interaction
	rows_it_appears = set()
	for idx, row in enumerate(individual):
		if all(row[col] == val for col,val in zip(cols,vals)):
			rows_it_appears.add(idx)
	return frozenset(rows_it_appears)

def fitness(individual, dsets):
	score = 0
	print(dsets[0])
	for dset1, dset2, num_times_sep_already in dsets:
		print(dset1, dset2, num_times_sep_already)
		requirement_to_reach = l - num_times_sep_already
		rows1 = set(rows_of_interaction(individual, I) for I in dset1)
		rows2 = set(rows_of_interaction(individual, I) for I in dset2)
		if len(rows1.symmetric_difference(rows2)) >= requirement_to_reach:
			score += 1
	return score

def try_N(N,dsets):
	pop = [random_individual(N) for _ in range(pop_size)]
	max_possible_fitness = len(dsets)
	for generation in range(num_generations): 
		fitnesses = list(sorted((I, fitness(I,dsets)) for I in pop))
		print(fitnesses)
		break
		best_ind, max_fitness = fitnesses[-1]
		if max_fitness == max_possible_fitness:
			return best_ind
	return None




def go(dsets):
	N = 10
	print('Going up...')
	while True:
		result = try_N(N,dsets)
		if N == 1 and result:
			return result
		if result:
			break
		N *= 2
		print(N)
	N_lo, N_hi = N // 2, N
	while N_lo < N_hi:
		N_mid = (N_lo + N_hi) // 2
		result2 = try_N(N_mid,dsets)
		if result2:
			N_hi = N_mid
			result = copy.deepcopy(result2)
		else:
			N_lo = N_mid + 1
		print(N_lo, N_hi)
	return result

num_generations = 100
pop_size = 100
d,t,k,v,l = map(int,sys.argv[1:])
dsets = read_file(d,t,k,v,l)
for dset in dsets:
	print(dset)
result = go(dsets)
print(result)
