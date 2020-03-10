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
parser.add_argument('--coverage', required=False, type=int, default = 10000,
	metavar='<int>', help='number of simulated reads [10000]')
arg = parser.parse_args()

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

print('WormBase introns', len(wormbase))
print('RNASeq_splice introns', len(splice))


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
	obs = set()
	for i in range(arg.coverage):
		p = random.random()
		for j in range(len(cramp)):
			if p <= cramp[j]:
				obs.add(j)
				break
	return(len(obs))
	
anno_obs = mass_ramp(wormbase)
splice_obs = mass_ramp(splice)
print('observed WormBase introns', anno_obs)
print('observed RNASeq_splice introns', splice_obs)

"""
for s in source:
	count = 0
	mass = 0
	for string in source[s]:
		count += 1
		mass += source[s][string]
	print(s, count, mass, mass/count)
"""