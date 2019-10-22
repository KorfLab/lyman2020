#!/usr/bin/env python3

import argparse
import sys
import operator
import os
import copy
import json

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


coding = 0
isolated = 0
spliced = 0
canonical = 0
rna_supported = 0
training = []

o = open(arg.out, 'w+')
o.write('region\tlen\ttxs\texons\trna_introns\tfit\n')

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
		genes = chrom.ftable.build_genes()
		for tx in genes[0].transcripts():
			if tx.end - tx.beg > tx_len:
				tx_len = tx.end - tx.beg
				ex_count = len(tx.exons)
			for intron in tx.introns:
				if (intron.beg in intron_support and
					intron.end in intron_support[intron.beg]): supported += 1
		if meta['introns'] == supported:
			rna_supported += 1 
			training.append(region)
			glen = genes[0].end - genes[0].beg
			tx_count = len(genes[0].transcripts())
			o.write('{}\t{}\t{}\t{}\t{}\n'.format(region, glen, tx_count, ex_count,
				rna_count))
				
o.close()

r = open(arg.report, 'w+')
r.write('coding:{}\nisolated:{}\nspliced:{}\ncanonical:{}\nsupported:{}'.format(coding,
	isolated, spliced, canonical, rna_supported))
r.close()


