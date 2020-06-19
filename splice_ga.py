#!/usr/bin/env python3

# estimates splice site probabilities to arrive at the observed pattern
# then generates splices at various depths

import argparse
import sys
import random
from grimoire.io import GFF_file
from grimoire.sequence import DNA
from grimoire.feature import Feature, mRNA, Gene, FeatureTable
from grimoire.genome import Reader

#random.seed(1)

parser = argparse.ArgumentParser(
	description='Splice site estimator and generator.')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='gene region, just one sequence')
parser.add_argument('--gff', required=True, type=str,
	metavar='<path>', help='gff with RNASeq_splice data')
parser.add_argument('--transcripts', required=False, type=int, default=1000,
	metavar='<int>', help='number of transcripts to produce [%(default)i]')
parser.add_argument('--min_intron', required=False, type=int, default=35,
	metavar='<int>', help='minimum intron size [%(default)i]')
parser.add_argument('--min_exon', required=False, type=int, default=10,
	metavar='<int>', help='minimum exon size [%(default)i]')
parser.add_argument('--generations', required=False, type=int, default=100,
	metavar='<int>', help='number of generations to run [%(default)i]')
parser.add_argument('--population', required=False, type=int, default=100,
	metavar='<int>', help='population size [%(default)i]')
parser.add_argument('--breeders', required=False, type=int, default=20,
	metavar='<int>', help='number of breeders in population [%(default)i]')
parser.add_argument('--mut_rate', required=False, type=float, default=0.1,
	metavar='<float>', help='number of breeders in population [%(default).3f]')
parser.add_argument('--mut_mag', required=False, type=float, default=0.2,
	metavar='<float>', help='number of breeders in population [%(default).3f]')
arg = parser.parse_args()

def splice_gen(ind, don, acc, imin, emin):
	pdon = ind['pd']
	pacc = ind['pa']
	
	i = 0
	splice = []
	while i < len(don):
		if random.random() < pdon[i]:
			j = 0
			while j < len(acc):
				d = acc[j] - don[i]
				if d > imin and random.random() < pacc[j]:
					splice.append((don[i], acc[j]))
					k = i
					while k < len(don):
						d = don[k] - acc[j]
						if d > emin:
							i = k
							break
						else:
							k += 1
							i = k
					break
				else:
					j += 1
		else:
			i += 1

	return splice

def fitness(ind, don, acc, obs, imin, emin):
	samples = 1000
	int_ct = {}
	for i in range(samples):
		introns = splice_gen(ind, don, acc, imin, emin)
		for i in introns:
			if i not in int_ct: int_ct[i] = 1
			else:               int_ct[i] += 1
	intf = {}
	for i in int_ct: intf[i] = int_ct[i] / samples
	
	# compare iso_freq with actual introns
	sum = 0
	for i in intf:
		if i in obs: sum += abs(intf[i] - obs[i])
		else: sum += 0.1 # default penalty...

	return sum

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
don = []
acc = []
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
mi  = arg.min_intron
me  = arg.min_exon

## intialize population
ga = []
for i in range(pop):
	pd = [random.random() for i in range(len(don))]
	pa = [random.random() for i in range(len(acc))]
	ind = {'pd':pd, 'pa':pa, 'fit':None}
	ind['fit'] = fitness(ind, don, acc, obs, mi, me)
	ga.append(ind)

## evolve population
for g in range(gen):
	
	# sort population by fitness (lower is better)
	ga = sorted(ga, key=lambda k: k['fit'])
	min, max, ave = ga_stats(ga, par)
	print(g, min, ave, max)	
	
	# replace less fit with children
	for i in range(par, pop):
		p1 = random.randint(0, par)
		p2 = random.randint(0, par)
		ga[i] = breed(ga[p1], ga[p2], mut, mag)
		ga[i]['fit'] = fitness(ga[i], don, acc, obs, mi, me)

ga = sorted(ga, key=lambda k: k['fit'])
best = ga[0]

#############################
# Part 3: generate isoforms #
#############################

count = {}
for i in range(arg.transcripts):
	isoforms = splice_gen(best, don, acc, mi, me)
	istr = str(isoforms)
	if istr not in count: count[istr] = 1
	else:                 count[istr] += 1

for istr in count:
	print(count[istr], istr)

