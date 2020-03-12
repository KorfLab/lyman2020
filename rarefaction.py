#!/usr/bin/env python3

import argparse
import sys
import json
import operator
import random
from grimoire.io import GFF_stream

## Command line stuff ##

parser = argparse.ArgumentParser(
	description='Estimates a rarefaction curve for introns.',
	formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--gff3', required=True, type=str,
	metavar='<path>', help='gff file, may be compressed')
parser.add_argument('--source', required=True, type=str,
	metavar='<str>', help='rule-based parsing based on gff source')
parser.add_argument('--maxcov', required=False, type=int, default = 10000,
	metavar='<int>', help='number of simulated reads [10000]')
parser.add_argument('--step', required=False, type=int, default = 100,
	metavar='<int>', help='coverage range increment [100]')	
arg = parser.parse_args()

# memorize coordinates of introns
gff = GFF_stream(arg.gff3)
wormbase = {}
splice = {}
for f in gff:
	if f.type != 'intron': continue
	stringy = '{}{}{}-{}'.format(f.chrom, f.strand, f.beg, f.end)
	if f.score == '.': f.score = 1.0
	else: f.score = float(f.score)
	if f.source == 'WormBase':
		wormbase[stringy] = f.score
	elif f.source == 'RNASeq_splice':
		splice[stringy] = f.score

def observations(ramp, cov):
	obs = set()
	for i in range(cov):
		p = random.random()
		for j in range(len(ramp)):
			if p <= ramp[j]:
				obs.add(j)
				break
	return(len(obs))

def mass_ramp(source):
	mass = 0
	for feature in source:
		mass += source[feature]
	pramp = []
	for fstring, fmass in sorted(source.items(), key = operator.itemgetter(1)):
		pramp.append(fmass / mass)
	cramp = []
	cramp.append(pramp[0])
	for i in range(1, len(pramp)):
		cramp.append(cramp[i - 1] + pramp[i])
	obs = []
	for i in range(0, arg.maxcov, arg.step):
		obs.append(observations(cramp, i))
	return(obs)	


	
anno_obs = mass_ramp(wormbase)
splice_obs = mass_ramp(splice)

print('cov\tWormBase\tRNASeq_splice')
for i in range(len(anno_obs)):
	print('{}\t{}\t{}'.format(i * arg.step, anno_obs[i], splice_obs[i]))
print('gff\t{}\t{}'.format(len(wormbase), len(splice)))


"""
for s in source:
	count = 0
	mass = 0
	for string in source[s]:
		count += 1
		mass += source[s][string]
	print(s, count, mass, mass/count)
"""