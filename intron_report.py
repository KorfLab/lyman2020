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
with open(arg.table) as table:
	for line in table.readlines():
		meta = line.split('\t')
		region = meta[0]
		if region == 'region': continue # skip header
		tx_count = int(float(meta[3]))
		if tx_count != 1: continue
		debug += 1
		if debug == 1000: sys.exit(0)	
		prefix = arg.regions + '/' + region + '/' + region
	
		# memorize RNA-supported intron coordinates and count RNA-supported introns
		genome = Reader(gff=prefix + '.gff', fasta=prefix + '.fa', source=arg.source)
		for chrom in genome:
			gene = chrom.ftable.build_genes()[0]
			tx = gene.transcripts()[0]
			for i in range(len(tx.cdss) -1):
				beg = tx.cdss[i].beg
				ibeg = tx.cdss[i].end + 1
				end = tx.cdss[i+1].end
				iend = tx.cdss[i+1].beg - 1
				iscore = None
				count = 0
				acc_pos = []
				don_pos = []
				acc_neg = []
				don_neg = []
				for f in chrom.ftable.fetch(beg, end):
					if f.type != 'intron': continue
					if f.source != 'RNASeq_splice': continue
					if f.strand != tx.strand: continue
					if f.beg < beg or f.end > end: continue
					if f.beg == ibeg and f.end == iend:
						count += 1
						iscore = f.score
					elif f.beg == ibeg:
						if gene.strand == '+':
							print('+ acceptor variant')
						else:
							print('- donor variant')
					elif f.end == iend:
						if gene.strand == '+':
							print('+ donor variant')
						else:
							print('- acceptor variant')
""""	
					if gene.strand == '+' and f.beg == ibeg:
						acc_pos.append(f) # acceptor sites differ, pos strand
					
					if gene.strand == '-' and f.beg == ibeg:
						don_neg.append(f) # donor sites differ, neg strand

					if gene.strand == '+' and f.end == iend:
						don_pos.append(f) # donor sites differ, pos strand

					if gene.strand == '-' and f.end == iend:
						acc_neg.append(f) # acceptor sites differ, neg strand

				if count != 1:
					raise Exception('number of RNASeq_splice features matching canonical intron is not one!')
				for f in acc_pos:
					print('acc_pos', abs(iend - f.end) % 3, abs(iend - f.end), f.score / iscore)
				for f in don_pos:
					print('don_pos', abs(ibeg - f.beg) % 3, abs(ibeg - f.beg), f.score / iscore)
				for f in acc_neg:
					print('acc_neg', abs(ibeg - f.beg) % 3, abs(ibeg - f.beg), f.score / iscore)
				for f in don_neg:
					print('don_neg', abs(iend - f.end) % 3, abs(iend - f.end), f.score / iscore)
				

				
			
		#		print(beg, end, f.gff())
		#print("")
		#sys.exit(0)
		#print("")
		
		#print(gene.id, len(rna_introns))


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