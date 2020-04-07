#!/usr/bin/env python3

import argparse
import sys
import operator
import os
import copy
import json
import math
import itertools
import numpy as np

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
parser.add_argument('--write', '-w', action='store_true', required=False,
	help='if present, will write to .isoform files')
parser.add_argument('--verbose', '-v', action='store_true', required=False,
	help='if present, will print results to stdout')
parser.add_argument('--nruns', required=True,
	help='The number of times to generate paths per gene ("sample size")')
parser.add_argument('--nreport', required=True,
	help='The number of isoforms to report (the n most encountered in simulation)')
parser.add_argument('--target', required=False, type=str,
	metavar='<str>', help='run only on a specified region')

arg = parser.parse_args()

if not os.path.exists(arg.regions):
	raise Exception('invalid path to region directories')

def lkd(set1, set2):
	d = 0
	for beg in set1:
		for end in set1[beg]:
			d += set1[beg][end] * math.log(set1[beg][end] / set2[beg][end])

	return d

def convert_to_freq(count):
	total = 0
	for beg in count:
		for end in count[beg]:
			total += count[beg][end]
	for beg in count:
		for end in count[beg]:
			count[beg][end] /= total

def distance(region, gff, gene):

	# collect rnaseq introns
	rnaseq = {}
	rnamass = 0
	for f in gff.get(type='intron'):
		if f.strand != gene.strand: continue
		if f.source != 'RNASeq_splice': continue
		if f.beg not in rnaseq:
			rnaseq[f.beg] = {}
		if f.end not in rnaseq[f.beg]:
			rnaseq[f.beg][f.end] = float(f.score)
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
				annotated[f.beg][f.end] = 0
			annotated[f.beg][f.end] += iweight

	# add missing introns from rnaseq with 1 pseudocount
	for beg in rnaseq:
		for end in rnaseq[beg]:
			if beg not in annotated:
				annotated[beg] = {}
			if end not in annotated[beg]:
				annotated[beg][end] = 1

	# convert to freqs and calculate distance
	convert_to_freq(rnaseq)
	convert_to_freq(annotated)
	return lkd(rnaseq, annotated), lkd(annotated, rnaseq)

def fold(a, b):
	return a/b if a > b else b/a

def isoforms(gene, rna_introns, freq, n_gen, n_rep, dist_thr, fold_thr, file, region):
	vis_introns = []
	for f in rna_introns:
		if f.score > freq and f.strand == gene.strand:
			vis_introns.append(f)

	max_intr = max(vis_introns, key=operator.attrgetter('score'))

	vis_introns.sort(key = operator.attrgetter('beg')) # sort by most expressed

	begs = set()
	ends = set()

	for tx in gene.transcripts():
		first = tx.exons[0]
		last = tx.exons[len(tx.exons)-1]
		begs.add(first.end + 1)	
		ends.add(last.beg  - 1)

	begs = sorted(list(begs))
	ends = sorted(list(ends))

	maxgroup = [] # introns overlapping the most expressed

	expr_max = 0
	for j in range(0, len(vis_introns)):
		if max_intr.overlap(vis_introns[j]):
			expr_max += vis_introns[j].score
			maxgroup.append(vis_introns[j])	 
			
	paths = {}
	for _ in range(n_gen): # generate n_gen isoforms
		path = []
		used = []
		used.extend(maxgroup)
		prob = max_intr.score / expr_max
		if(np.random.choice([True, False], size=1, p=[prob, 1-prob])):
			path.append(max_intr)

		for i in range(len(vis_introns)):
			if(vis_introns[i] in used): continue
			prob = vis_introns[i].score / expr_max
			if(np.random.choice([True, False], size=1, p=[prob, 1-prob])):
				path.append(vis_introns[i])

			overlap = []
			for j in range(len(vis_introns)):
				if(vis_introns[i].overlap(vis_introns[j])):
					used.append(vis_introns[j])
			
		path.sort(key=operator.attrgetter('beg'))
		rep = encode_iso(path)
		paths.setdefault(rep, 0) # count the times we have seen this form
		paths[rep] += 1

	sort = sorted(paths, key=paths.get, reverse=True)
	tot = sum(paths.values())
	if(arg.verbose): print('\nregion '  + region + ' (' + gene.id + ')')
	for i in range(min(n_rep, len(sort))):
		if(arg.write):
			isofile.write(str(paths[sort[i]]/tot) +'\t'+(sort[i] if len(sort[i].strip()) > 0 else '[none] ') + '\n')
		if(arg.verbose):
			print(str(paths[sort[i]]/tot) +'\t'+(sort[i] if len(sort[i].strip()) > 0 else '[none]'))

def encode_iso(intrs):
	rep = '' #string representation of the path
	for intr in intrs:
		rep += str(intr.beg) + ',' + str(intr.end) + ' '

	return rep

coding = 0
isolated = 0
spliced = 0
canonical = 0
rna_supported = 0
training = []

o = open(arg.out, 'w+')
o.write('region\tgid\tlen\ttxs\texons\trna_introns\texp\tiso4\tiso6\tiso8\tlkd1\tlkd2\n')

for region in os.listdir(arg.regions):
	if(arg.target) and region != arg.target: continue
	prefix = arg.regions + '/' + region + '/' + region
	# read region status from JSON
	try:
		f = open(prefix + '.json')
	except:
		print('error opening json')
		continue
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
	# memorize RNA-supported intron coordinates and count RNA-supported introns
	gff = GFF_file(file=prefix + '.gff', source=arg.source)
	genome = Reader(gff=prefix + '.gff', fasta=prefix + '.fa', source=arg.source)

	dist_thr = 30 # minimum distance between introns to be valid
	fold_thr = 5  # introns must be within k-fold of the mean expression level
	nruns = int(arg.nruns)
	nreport = int(arg.nreport)
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
			rna_supported += 1
			training.append(region)
			glen = gene.end - gene.beg
			tx_count = len(gene.transcripts())
			# check how well annotation matches introns using new methods
			d1, d2 = distance(region, gff, gene)
			isofile = open(prefix + '.isoforms', 'w+')
			isoforms(gene, rna_introns, 1e-8, nruns, nreport, dist_thr, fold_thr, isofile, region)

			#o.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(region, gene.id, glen,
			#	tx_count, ex_count, rna_count, rna_mass, iso4, iso6, iso8, d1, d2))

o.close()

r = open(arg.report, 'w+')
r.write('coding:{}\nisolated:{}\nspliced:{}\ncanonical:{}\nsupported:{}'.format(coding,
	isolated, spliced, canonical, rna_supported))
r.close()
#!/usr/bin/env python3

