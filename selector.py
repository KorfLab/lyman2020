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
Summarizes and selects regions generated in setup.py for use in experiment.
"""

parser = argparse.ArgumentParser(
	description='Selects regions containing single protein-coding gene.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument('--regions', required=True, type=str,
	metavar='<path>', help='directory containing region sub-directories')
parser.add_argument('--source', required=True, type=str,
	metavar='<str>', help='rule-based parsing based on gff source')
parser.add_argument('--report', required=True, type=str,
	metavar='<path>', help='report output filename')
parser.add_argument('--out', required=True, type=str,
	metavar='<path>', help='qualified region list output filename')	
arg = parser.parse_args()

if not os.path.exists(arg.regions):
	raise Exception('invalid path to region directories')

def lkd(set1, set2):
	d = 0
	for beg in set1:
		for end in set1[beg]:
			for strand in set1[beg][end]:
				d += set1[beg][end][strand] * math.log(set1[beg][end][strand] / set2[beg][end][strand])

	return d

def convert_to_freq(count):
	total = 0
	for beg in count:
		for end in count[beg]:
			for strand in count[beg][end]:
				total += count[beg][end][strand]
	for beg in count:
		for end in count[beg]:
			for strand in count[beg][end]:
				count[beg][end][strand] /= total

def distance(region, gff, gene):

	# collect rnaseq introns
	rnaseq = {}
	rnamass = 0
	for f in gff.get(type='intron'):
		if f.source == 'RNASeq_splice':
			if f.beg not in rnaseq:
				rnaseq[f.beg] = {}
			if f.end not in rnaseq[f.beg]:
				rnaseq[f.beg][f.end] = {}
			if f.strand not in rnaseq[f.beg][f.end]:
				rnaseq[f.beg][f.end][f.strand] = float(f.score)
				rnamass += float(f.score)

	# collect annotated introns
	annotated = {}
	txs = []
	for tx in gene.transcripts():
		if len(tx.introns) > 0: txs.append(tx)
	txmass = rnamass / len(txs)
	for tx in txs:	
		iweight = txmass / (len(tx.introns))
		for f in tx.introns:
			f.score = iweight # not really used but maybe useful for debugging
			if f.beg not in annotated:
				annotated[f.beg] = {}
			if f.end not in annotated[f.beg]:
				annotated[f.beg][f.end] = {}
			if f.strand not in annotated[f.beg][f.end]:
				annotated[f.beg][f.end][f.strand] = 0
			annotated[f.beg][f.end][f.strand] += iweight

	# add missing introns from rnaseq with 1 pseudocount
	for beg in rnaseq:
		for end in rnaseq[beg]:
			for strand in rnaseq[beg][end]:
				if beg not in annotated:
					annotated[beg] = {}
				if end not in annotated[beg]:
					annotated[beg][end] = {}
				if strand not in annotated[beg][end]:
					annotated[beg][end][strand] = 1
	
	# convert to freqs and calculate distance
	convert_to_freq(rnaseq)
	convert_to_freq(annotated)
	return lkd(rnaseq, annotated), lkd(annotated, rnaseq)

coding = 0
isolated = 0
spliced = 0
canonical = 0
rna_supported = 0
training = []

o = open(arg.out, 'w+')
o.write('region\tlen\ttxs\texons\trna_introns\tlkd1\tlkd2\n')

for region in os.listdir(arg.regions):
	prefix = arg.regions + '/' + region + '/' + region

	# read region status from JSON
	f = open(prefix + '.json')
	meta = json.loads(f.read())
	f.close()
	# skip multi-gene regions, non-coding regions, un-spliced and non-canonical genes
	coding += meta['pc_genes']
	if meta['pc_genes'] != 1 or meta['nc_genes'] > 0: continue
	isolated += 1
	if meta['introns'] == 0: continue
	spliced += 1
	if meta['pc_issues']: continue
	canonical += 1
	#memorize RNA-seq supported intron coordinates and count RNA-supported introns
	gff = GFF_file(file=prefix + '.gff', source=arg.source)
	rna_introns = gff.get(type='intron')
	rna_count = 0
	intron_support = {}
	for intron in rna_introns:
		if intron.source != 'RNASeq_splice': continue
		rna_count += 1
		if intron.beg not in intron_support:
			intron_support[intron.beg] = {}
			intron_support[intron.beg][intron.end] = True
	#determine support for annotated introns
	genome = Reader(gff=prefix + '.gff', fasta=prefix + '.fa', source=arg.source)

	for chrom in genome:
		supported = 0
		tx_len = 0
		ex_count = 0
		gene = chrom.ftable.build_genes()[0]
		for tx in gene.transcripts():
			if tx.end - tx.beg > tx_len:
				tx_len = tx.end - tx.beg
				ex_count = len(tx.exons)
			for intron in tx.introns:
				if (intron.beg in intron_support and
					intron.end in intron_support[intron.beg]): supported += 1
		if meta['introns'] == supported:
			rna_supported += 1 
			training.append(region)
			glen = gene.end - gene.beg
			tx_count = len(gene.transcripts())
			# check how well annotation matches introns using new methods
			d1, d2 = distance(region, gff, gene)
			o.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(region, glen, tx_count,
				ex_count, rna_count, d1, d2))
				
o.close()

r = open(arg.report, 'w+')
r.write('coding:{}\nisolated:{}\nspliced:{}\ncanonical:{}\nsupported:{}'.format(coding,
	isolated, spliced, canonical, rna_supported))
r.close()