#!/usr/bin/env python3

import argparse
import sys
import operator
import os
import copy
import json

from grimoire.sequence import DNA
from grimoire.feature import Feature, mRNA, Gene, FeatureTable
from grimoire.genome import Reader

## Command line stuff ##

extended_help = """
Custom version of haman used for Lyman2020
"""

parser = argparse.ArgumentParser(
	description='Segments chromosome sequences and features.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to input fasta file')
parser.add_argument('--gff', required=True, type=str,
	metavar='<path>', help='path to input GFF3 (or similar) file')
parser.add_argument('--out', required=True, type=str,
	metavar='<str>', help='output name (file or dir)')
parser.add_argument('--source', required=True, type=str,
	metavar='<str>', help='rule-based parsing based on gff source')
parser.add_argument('--padding', required=True, type=int,
	metavar='<int>', help='length of flanking sequence')
arg = parser.parse_args()

class HamanError(Exception):
	pass

def gene_regions(chrom, padding):

	# create gene-ish regions with padding added to each gene
	gene_ish = []
	for gene in [f for f in chrom.ftable.features if f.type == 'gene']:
		beg = gene.beg - padding
		end = gene.end + padding
		if beg < 1: beg = 1
		if end > len(chrom.seq): end = len(chrom.seq)
		gene_ish.append(Feature(chrom, beg, end, '.', 'gene-ish'))
	
	# create clusters of gene-ish regions
	gtable = FeatureTable(dna=chrom, features=gene_ish)
	gtable._sort()
	skip = {}
	for gene in gtable.features:
		if gene in skip: continue
		beg = gene.beg
		end = gene.end
		while True:
			got_all = True
			overlaps = gtable.fetch(beg, end)
			for ov in overlaps:
				skip[ov] = True
				if ov.end > end:
					end = ov.end
					got_all = False
			if got_all: break
		yield Feature(chrom, beg, end, '.', 'region', source='haman')

if __name__ == '__main__':	

	if not os.path.exists(arg.out):
		os.mkdir(arg.out)
	
	genome = Reader(gff=arg.gff, fasta=arg.fasta, source=arg.source)
	idx = 0
	for chrom in genome:
		for region in gene_regions(chrom, arg.padding):
			idx += 1
			dir = arg.out + '/' + str(idx)
			os.mkdir(dir)
			
			# fasta
			seq = chrom.seq[region.beg:region.end+1] # +1 or not?
			dna = DNA(name=str(idx), seq=seq,
				desc='chrom:{} beg:{} end:{}'.format(chrom.name,
				region.beg, region.end))
			ffp = open(dir + '/' + str(idx) + '.fa', 'w+')
			ffp.write(dna.fasta())
			ffp.close()
						
			# gff
			stuff = chrom.ftable.fetch(region.beg, region.end)
			keep = []
			for f in stuff:
				if f.beg <= region.beg or f.end >= region.end: continue
				keep.append(f)
			remap = []
			for f in keep:
				nbeg = f.beg - region.beg
				nend = f.end - region.beg
				remap.append(Feature(dna, nbeg, nend,
					f.strand, f.type, phase=f.phase,
					score=f.score, source=f.source, id=f.id, pid=f.pid))
			gfp = open(dir + '/' + str(idx) + '.gff', 'w+')
			for f in remap:
				dna.ftable.add_feature(f)
				gfp.write(f.gff())
				gfp.write('\n')
			
			# json metadata
			genes = dna.ftable.build_genes()
			names = []
			pc_genes = 0
			nc_genes = 0
			pc_issues = 0
			introns = 0
			strand = None
			splices = 0
			max_splices = 0
			reads = 0
			max_reads = 0
			
			for gene in genes:
				if strand == None:
					strand = gene.strand
				elif strand != gene.strand:
					strand = '.'
				if len(gene.transcripts()) == 0: nc_genes += 1
				else:
					pc_genes += 1
					if gene.issues: pc_issues += 1
					for tx in gene.transcripts():
						introns += len(tx.introns)
			for f in remap:
				if f.type == 'gene':
					names.append(f.id)
				elif f.source == 'RNASeq_splice':
					splices += 1
					if f.score > max_splices: max_splices = f.score
				elif f.source == 'RNASeq_reads':
					reads += 1
					if f.score > max_reads: max_reads = f.score
			
			meta = {
				'region': idx,
				'chrom': chrom.name,
				'beg': region.beg,
				'end': region.end,
				'length': len(dna.seq),
				'genes': len(names),
				'gene_names': names,
				'pc_genes': pc_genes,
				'nc_genes': nc_genes,
				'pc_issues': pc_issues,
				'introns': introns,
				'strand': strand,
				'RNASeq_splice': splices,
				'max_splices': max_splices,
				'RNASeq_reads': reads,
				'max_reads': max_reads,
			}
			mfp = open(dir + '/' + str(idx) + '.json', 'w+')
			mfp.write(json.dumps(meta, indent=4))
			mfp.close()

