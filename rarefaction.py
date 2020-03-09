#!/usr/bin/env python3

import argparse
import sys
import json
from grimoire.io import GFF_stream

## Command line stuff ##

parser = argparse.ArgumentParser(
	description='Estimates a rarefaction curve for introns.',
	formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--gff3', required=True, type=str,
	metavar='<path>', help='gff file, may be compressed')
parser.add_argument('--source', required=True, type=str,
	metavar='<str>', help='rule-based parsing based on gff source')
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

#print(json.dumps(wormbase, indent=4))

"""
for s in source:
	count = 0
	mass = 0
	for string in source[s]:
		count += 1
		mass += source[s][string]
	print(s, count, mass, mass/count)
"""