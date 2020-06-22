#!/usr/bin/env python3

import argparse
import random
import sys
from grimoire.genome import Reader
from spliceomatic import isoformer

parser = argparse.ArgumentParser(
	description='Splice site usage estimator')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='gene region, just one sequence')
parser.add_argument('--gff', required=True, type=str,
	metavar='<path>', help='gff with RNASeq_splice data')
parser.add_argument('--transcripts', required=False, type=int, default=1000,
	metavar='<int>', help='number of transcripts to produce [%(default)i]')
parser.add_argument('--generations', required=False, type=int, default=100,
	metavar='<int>', help='number of generations to run [%(default)i]')
parser.add_argument('--population', required=False, type=int, default=100,
	metavar='<int>', help='population size [%(default)i]')
parser.add_argument('--breeders', required=False, type=int, default=50,
	metavar='<int>', help='number of breeders in population [%(default)i]')
parser.add_argument('--mut_rate', required=False, type=float, default=0.1,
	metavar='<float>', help='number of breeders in population [%(default).3f]')
parser.add_argument('--mut_mag', required=False, type=float, default=0.2,
	metavar='<float>', help='number of breeders in population [%(default).3f]')
parser.add_argument('--samples', required=False, type=int, default=1000,
	metavar='<float>', help='transcripts in fitness calculation [%(default)i]')
parser.add_argument('--seed', required=False, type=int,
	metavar='<int>', help='set random seed')
parser.add_argument('--verbose', action='store_true', help='show progress')
arg = parser.parse_args()

if arg.seed: random.seed(arg.seed)

def distance(p1, p2):
	sum = 0
	for i in p1:
		if i in p2: sum += abs(p1[i] - p2[i])
		else:       sum += p1[i]
	return sum

def fitness(ind, don, acc, obs, s):
	iso, iu, du, au = isoformer(don, acc, ind['pd'], ind['pa'], samples=s)
	d1 = distance(iu, obs)
	d2 = distance(obs, iu)
	return d1 + d2

def mutate(sites, f, v):
	for i in range(len(sites)):
		if random.random() < f:
			r = random.random() * v
			if random.random() > 0.5: sites[i] += r
			else:                     sites[i] -= r
			if   sites[i] < 0: sites[i] = 0
			elif sites[i] > 1: sites[i] = 1

def breed(p1, p2, f, v):

	# recombination
	child = {'pd':[], 'pa':[], 'fit':None}
	for i in range(len(p1['pd'])):
		if random.random() < 0.5: child['pd'].append(p1['pd'][i])
		else:                     child['pd'].append(p2['pd'][i])
	for i in range(len(p1['pa'])):
		if random.random() < 0.5: child['pa'].append(p1['pa'][i])
		else:                     child['pa'].append(p2['pa'][i])
	
	# mutation
	mutate(child['pd'], f, v)
	mutate(child['pa'], f, v)
	
	return child

def ga_stats(ga, t):
	min = ga[0]['fit']
	max = ga[0]['fit']
	sum = 0
	for i in ga:
		if i['fit'] < min: min = i['fit']
		if i['fit'] > max: max = i['fit']
		sum += i['fit']
	ave = sum / len(ga)
	return min, max, ave

##################################################
# Part 1: collect introns, donors, and acceptors #
##################################################

genome = Reader(gff=arg.gff, fasta=arg.fasta, source='wb.270')
chrom = genome.next()
intron = []
don = [] # positions of donor sites
acc = [] # positions of acceptor sites
for f in [f for f in chrom.ftable.features if f.source == 'RNASeq_splice']:
	intron.append({'beg':f.beg, 'end':f.end, 'score':f.score})
	if f.beg not in don: don.append(f.beg)
	if f.end not in acc: acc.append(f.end)
don.sort()
acc.sort()

obs = {}
total = 0
for i in intron: total += i['score']
for i in intron:
	i['prob'] = i['score'] / total
	obs[(i['beg'], i['end'])] = i['prob']

###############################################################
# Part 2: estimate don & acc values that fit observed introns #
###############################################################

pop = arg.population
gen = arg.generations
par = arg.breeders
mut = arg.mut_rate
mag = arg.mut_mag
sam = arg.samples

## intialize population
ga = []
for i in range(pop):
	pd = [random.random() for i in range(len(don))]
	pa = [random.random() for i in range(len(acc))]
	ind = {'pd':pd, 'pa':pa, 'fit':None}
	ind['fit'] = fitness(ind, don, acc, obs, sam)
	ga.append(ind)

## evolve population
for g in range(gen):
	
	# sort population by fitness (lower is better)
	ga = sorted(ga, key=lambda k: k['fit'])
	if arg.verbose:
		min, max, ave = ga_stats(ga, par)
		sys.stderr.write(f'{g} {min} {ave} {max}\n')	
	
	# replace less fit with children
	for i in range(par, pop):
		p1 = random.randint(0, par)
		p2 = random.randint(0, par)
		ga[i] = breed(ga[p1], ga[p2], mut, mag)
		ga[i]['fit'] = fitness(ga[i], don, acc, obs, sam)

ga = sorted(ga, key=lambda k: k['fit'])
best = ga[0]

##################
# Part 3: output #
##################

print('Donors')
for i in range(len(don)):
	print(don[i], best['pd'][i])
print('Acceptors')
for i in range(len(acc)):
	print(acc[i], best['pa'][i])


