#!/usr/bin/env python3

import argparse
import sys
import operator
import os
import copy
import json
import math

from grimoire.io import GFF_file
from grimoire.sequence import DNA
from grimoire.feature import Feature, mRNA, Gene, FeatureTable
from grimoire.genome import Reader

## Command line stuff ##

extended_help = """
Creates a frequency distribution of intron phase variants.
"""

parser = argparse.ArgumentParser(
	description='Creates a frequency distribution of intron phase variants.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument('--gff', required=True, type=str,
	metavar='<str>', help='input GFF file')
parser.add_argument('--source', required=True, type=str,
	metavar='<str>', help='rule-based parsing based on gff source')
arg = parser.parse_args()

gff = GFF_file(file=arg.gff, source=arg.source)

shift_scores = {}
for i in range(3):
	shift_scores[i] = {'count':0, 'sum':0}

for chrom in gff.chroms:
	# memorize introns
	coords = {}
	for intron in gff.get(chrom=chrom, type='intron'):
		if intron.source == 'RNASeq_splice':
			if intron.beg not in coords:
				coords[intron.beg] = []
			if intron.end not in coords:
				coords[intron.end] = []
			coords[intron.beg].append(intron)
			coords[intron.end].append(intron)

	for coord in coords:
		# compare introns to each other
		# do not compare introns to themselves or to opposite-strand introns
		# do not compare introns to each other when the shared coordinate
		# is the end of one intron but the beginning of the other intron
		for i in range(len(coords[coord])):
			for j in range(i + 1, len(coords[coord])):
				if coords[coord][i] is coords[coord][j]: continue
				if coords[coord][i].strand != coords[coord][j].strand: continue
				shift_length = None
				if coords[coord][i].beg == coords[coord][j].beg:
					shift_length = abs(coords[coord][i].end - coords[coord][j].end)
# 				elif coords[coord][i].end == coords[coord][j].end:
# 					shift_length = abs(coords[coord][j].beg - coords[coord][i].beg)
				else: continue
				phase_shift = shift_length % 3
				score_fold = None
				if float(coords[coord][i].score) > float(coords[coord][j].score):
					score_fold = float(coords[coord][i].score) / float(coords[coord][j].score)
				else:
					score_fold = float(coords[coord][j].score) / float(coords[coord][i].score)
				shift_scores[phase_shift]['count'] += 1
				shift_scores[phase_shift]['sum'] += score_fold
				print('{}\t{:.3f}\t{}'.format(phase_shift, score_fold, shift_length))
	