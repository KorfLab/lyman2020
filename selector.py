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

def lkd (set1, set2):

	d = 0
	for f1 in set1:
		for f2 in set2:
			if f1.beg == f2.beg and f1.end == f2.end and f1.strand == f2.strand:
				d += f1.score * math.log(f1.score / f2.score)
	return d

def distance (region, gff, gene):

	# collect rnaseq introns
	rnaseq = []
	rtotal = 0
	for f in gff.get(type='intron'):
		if f.source == 'RNASeq_splice':
			rnaseq.append(f)
			rtotal += float(f.score)
	if len(rnaseq) == 0:
		print(region, "no rnaseq data?")
		return(0,0)
		sys.exit(1)

	# change to freq
	for f in rnaseq:
		f.score = float(f.score) / rtotal

	# collect annotated introns
	annotated = []
	txmass = rtotal / len(gene.transcripts())
	for tx in gene.transcripts():
		if len(tx.introns) > 0:
			iweight = txmass / (len(tx.introns))
			for f in tx.introns:
				f.score = iweight
				annotated.append(f)
	if len(annotated) == 0:
		print(region, "no annotation data?")
		return(0,0)
		sys.exit(1)
	
	# add missing introns from rnaseq with 1 pseudocount
	add = []
	for r in rnaseq:
		unique = True
		for f in annotated:
			if r.beg == f.beg and r.end == f.end and r.strand == f.strand:
				unique = False
				break
		if unique: add.append(Feature(f.dna, r.beg, r.end, r.strand, 'intron',
			score=1))
	annotated += add
	
	# change to frequency
	atotal = 0
	for f in annotated: atotal += f.score
	for f in annotated: f.score = f.score / atotal
	
	d1, d2 =  lkd(rnaseq, annotated), lkd(annotated, rnaseq)
	return d1, d2


coding = 0
isolated = 0
spliced = 0
canonical = 0
rna_supported = 0

training = []

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
	#memorize RNA-seq supported intron coordinates
	gff = GFF_file(file=prefix + '.gff', source=arg.source)
	rna_introns = gff.get(type='intron')
	intron_support = {}
	for intron in rna_introns:
		if intron.source != 'RNASeq_splice': continue
		if intron.beg not in intron_support:
			intron_support[intron.beg] = {}
			intron_support[intron.beg][intron.end] = True
	#determine support for annotated introns
	genome = Reader(gff=prefix + '.gff', fasta=prefix + '.fa', source=arg.source)
	supported = 0
	gene = None
	for chrom in genome:
		genes = chrom.ftable.build_genes()
		gene = genes[0]
		for tx in gene.transcripts():
			for intron in tx.introns:
				if intron.beg in intron_support and intron.end in intron_support[intron.beg]:
					supported += 1
	if meta['introns'] == supported:
		rna_supported += 1 
		training.append(region)
		# check how well annotation matches introns using new methods
		d1, d2 = distance(region, gff, gene)


r = open(arg.report, 'w+')
r.write('coding:{}\nisolated:{}\nspliced:{}\ncanonical:{}\nsupported:{}'.format(coding, isolated, spliced, canonical, rna_supported))
r.close()


o = open(arg.out, 'w+')
for region in training:
	gff = arg.regions + '/' + region + '/' + region + '.gff'
	fasta = arg.regions + '/' + region + '/' + region + '.fa'
	genome = Reader(gff=gff, fasta=fasta, source=arg.source)
	
	
o.close()
