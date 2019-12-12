#!/usr/bin/env python3

import argparse
import sys
import operator
import os
import copy
import json
import math
import matplotlib.pyplot as plot
import matplotlib.patches as patches
import matplotlib.figure as figure


from grimoire.io import GFF_file
from grimoire.sequence import DNA
from grimoire.feature import Feature, mRNA, Gene, FeatureTable
from grimoire.genome import Reader

## Command line stuff ##

extended_help = """
View the visible introns for a region and frequency threshold. Blue indicates a
valid intron, and red is invalid. Green lines are annotated beginnings (end of the first exon), and purple lines are known endings..
"""

parser = argparse.ArgumentParser(
	description='View the visible introns for a region and frequency threshold.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument('--regions', required=True, type=str,
	metavar='<path>', help='directory containing region sub-directories')
parser.add_argument('--region', required=True, type=str,
	metavar='<num>', help='region of interest')
parser.add_argument('--threshold', required=True, type=str,
	metavar='<str>', help='frequency threshold: 1e-4, 1e-6, or 1e-8')
parser.add_argument('--out', required=True, type=str,
	metavar='<path>', help='filename of output graphic (.pdf, .png, .svg ... )')
arg = parser.parse_args()

if not os.path.exists(arg.regions):
	raise Exception('invalid path to region directories')

prefix = arg.regions + '/' + arg.region + '/' + arg.region

f = open(prefix + '.json')
meta = json.loads(f.read())
f.close()

freq_thresh = float(arg.threshold)
dist_thresh = 30
# memorize RNA-supported intron coordinates and count RNA-supported introns
gff = GFF_file(file=prefix + '.gff', source='wb.270')
genome = Reader(gff=prefix + '.gff', fasta=prefix + '.fa', source='wb.270')
for chrom in genome:
	rna_count = 0
	rna_mass = 0
	supported = 0
	tx_len = 0
	ex_count = 0
	gene = chrom.ftable.build_genes()[0]
	rna_introns = []
	intron_support = {}
	for intron in chrom.ftable.features:
		if intron.type == 'intron' and intron.source == 'RNASeq_splice' and intron.strand == gene.strand:
			rna_count += 1
			rna_mass += int(float(intron.score))
			rna_introns.append(intron)
			if intron.beg not in intron_support:
				intron_support[intron.beg] = {}
				intron_support[intron.beg][intron.end] = True
	for intron in rna_introns:
		intron.score /= rna_mass

	for tx in gene.transcripts():
		if tx.end - tx.beg > tx_len:
			tx_len = tx.end - tx.beg
			ex_count = len(tx.exons)
		for intron in tx.introns:
			if (intron.beg in intron_support and
				intron.end in intron_support[intron.beg]): supported += 1
	if meta['introns'] == supported:
		glen = gene.end - gene.beg
		vis_introns = []
		mn = math.inf
		mx = 0
		for f in rna_introns:
			if(f.beg < mn): mn = f.beg
			if(f.end > mx): mx = f.end
			if f.score > freq_thresh and f.strand == gene.strand:
				vis_introns.append(f)
		vis_introns.sort(key = operator.attrgetter('beg'))

		# find known start/end points
		begs = set()
		ends = set()
		for tx in gene.transcripts():
			first = tx.exons[0]
			last = tx.exons[len(tx.exons)-1]
			begs.add(first.end + 1) # a valid intron starts 1 ahead of first exon
			ends.add(last.beg - 1)  # and ends 1 behind the start of last exon


		begs = sorted(list(begs))
		ends = sorted(list(ends))

		new_vis = []
		for end in ends:
			for i in range(len(vis_introns)-1):
				good = False
				for end in ends:
					good = not (vis_introns[i].end > end)
				if(good):
					new_vis.append(vis_introns[i])

		fig, axis = plot.subplots(1)
		axis.set_yticklabels([])
		axis.set(xlim=(mn-5, mx+5))
		axis.set(ylim=(0,len(vis_introns)*2))

		for i in range(len(vis_introns)):
			plot.plot([vis_introns[i].beg, vis_introns[i].end], [i+1, i+1],
				color=('blue' if vis_introns[i] in new_vis else 'red'))


		for b in begs:
			plot.plot([b, b], [0, len(vis_introns)*2], color='green', linewidth=0.1, linestyle='--')
		for e in ends:
			plot.plot([e, e], [0, len(vis_introns)*2], color='purple', linewidth=0.1, linestyle='--')

		plot.savefig(arg.out, dpi=400)
		
