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
parser.add_argument('--regions', required=True, type=str,
	metavar='<path>', help='directory containing region sub-directories')
parser.add_argument('--table', required=True, type=str,
	metavar='<path>', help='regions to examine')
parser.add_argument('--source', required=True, type=str,
	metavar='<str>', help='rule-based parsing based on gff source')
arg = parser.parse_args()

debug = 0
for region in os.listdir(arg.regions):
	debug += 1
	if debug == 200: sys.exit(0)
	
	prefix = arg.regions + '/' + region + '/' + region

	# read region status from JSON
	f = open(prefix + '.json')
	meta = json.loads(f.read())
	f.close()
	
	# skip multi-gene regions, non-coding regions, un-spliced and non-canonical genes
	if meta['pc_genes'] != 1 or meta['nc_genes'] > 0: continue
	if meta['introns'] == 0: continue
	if meta['pc_issues']: continue

	# memorize RNA-supported intron coordinates and count RNA-supported introns
	gff = GFF_file(file=prefix + '.gff', source=arg.source)
	genome = Reader(gff=prefix + '.gff', fasta=prefix + '.fa', source=arg.source)
	for chrom in genome:
		gene = chrom.ftable.build_genes()[0]
		txs = gene.transcripts()
		if len(txs) > 1: continue # for now, maybe address double-counting later
		tx = txs[0]
		if tx.strand == '-': continue # for now
		#rna_introns = []
		#for f in chrom.ftable.features:
		#	if f.type == 'intron' and f.source == 'RNASeq_splice' and f.strand == gene.strand:
		#		rna_introns.append(f)
		for i in range(len(tx.exons) -1):
			beg = tx.exons[i].beg
			end = tx.exons[i+1].end
			ibeg = tx.introns[i].beg
			iend = tx.introns[i].end
			count = 0
			stuff = []
			for f in chrom.ftable.fetch(beg, end):
				if f.type != 'intron': continue
				if f.source != 'RNASeq_splice': continue
				if f.strand != tx.strand: continue
				if f.beg < beg or f.end > end: continue
				stuff.append(f)
				if f.beg == ibeg and f.end == iend:
					count += 1
			if count == 0:
				print(beg, end, ibeg, iend)
				for f in stuff:
					print(f.gff())
				# same donor
				
				
				
		#		print(beg, end, f.gff())
		#print("")
		#sys.exit(0)
		#print("")
		
		#print(gene.id, len(rna_introns))

"""
genome = Reader(gff=arg.gff, fasta=arg.fasta, source=arg.source)
for chrom in genome:
	genes = chrom.ftable.build_genes()
	if len(genes) != 1 continue:
	gene = genes[0]
	if len(gene.introns) == 0: continue
	rna_introns = []
	for f in chrom.ftable.features:
		if f.type == 'intron' and f.source == 'RNASeq_splice':
			rna_introns.append
	print(gene.name, len(rna_introns))



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

"""