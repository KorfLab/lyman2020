#!/usr/bin/env python3

import argparse
import sys
import operator
import os
import copy
import json
import math
import itertools

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
parser.add_argument('--isomax', required=True, type=int,
	metavar='<path>', help='qualified region list output filename')
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

def write_isoforms(chunks, file):
	for iso in itertools.product(*chunks):
		for f in iso:
			file.write(str(f.beg) + "," + str(f.end) + ' ')
		file.write('\n')

def isoforms(gene, rna_introns, freq, dist_thr, fold_thr, file, region):
	vis_introns = []
	for f in rna_introns:
		if f.score > freq and f.strand == gene.strand:
			vis_introns.append(f)
	vis_introns.sort(key = operator.attrgetter('score'), reverse=True) # sort by most expressed
	
	chunks = [] #groups of overlapping introns
	used = []
	for i in range(len(vis_introns)):
		chunks.append([])
		curr = copy.copy(vis_introns[i])
		curr.beg = curr.beg - dist_thr # consider introns within a threshold part of the chunk
		curr.end = curr.end + dist_thr

		for j in range(0, len(vis_introns)):
			check = vis_introns[j]
			if(curr.overlap(check) and check not in used):
				for ch in chunks: # not sure about this
					for intr in ch:
						if(intr.overlap(check)):
							used.append(check)
							continue
				chunks[i].append(check)
				used.append(check)

	chunks = [ch for ch in chunks if len(ch) > 0]
	chunks.sort(key = lambda x: x[0].beg) # sort all chunks by position
		

	isoforms = 1 # estimate isoform count by multiplying number of introns in each overlapping chunk
	for ch in chunks:
		isoforms *= len(ch) if len(ch) > 1 else 1

	if(isoforms < arg.isomax):
		write_isoforms(chunks, file)
	else:
		print('skipping ' + str(region))
	return isoforms	


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
	for chrom in genome:
		rna_count = 0
		rna_mass = 0
		supported = 0
		tx_len = 0
		ex_count = 0
		gene = chrom.ftable.build_genes()[0]
		dist_thr = 30 # minimum distance between introns to be valid
		fold_thr = 5  # introns must be within k-fold of the mean expression level
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
			iso = open(prefix + '.isoforms', 'w+')
			iso.write('region: {} gid: {}\n'.format(region, gene.id))
			iso.write('\n1e-4 freq paths:\n')
			iso4 = isoforms(gene, rna_introns, 1e-4, dist_thr, fold_thr, iso, region)
			iso.write('\n1e-6 freq paths:\n')
			iso6 = isoforms(gene, rna_introns, 1e-6, dist_thr, fold_thr, iso, region)
			iso.write('\n1e-8 freq paths:\n')
			iso8 = isoforms(gene, rna_introns, 1e-8, dist_thr, fold_thr, iso, region)

			o.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(region, gene.id, glen,
				tx_count, ex_count, rna_count, rna_mass, iso4, iso6, iso8, d1, d2))

o.close()

r = open(arg.report, 'w+')
r.write('coding:{}\nisolated:{}\nspliced:{}\ncanonical:{}\nsupported:{}'.format(coding,
	isolated, spliced, canonical, rna_supported))
r.close()
#!/usr/bin/env python3

