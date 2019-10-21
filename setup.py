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

if __name__ == '__main__':	

	if not os.path.exists(arg.out):
		os.mkdir(arg.out)
	
	genome = Reader(gff=arg.gff, fasta=arg.fasta, source=arg.source)
	idx = 0
	for chrom in genome:
	
		# create regions of overlapping genes
		chrom.ftable._sort()
		genes = []
		for f in chrom.ftable.features:
			if f.type == 'gene': genes.append(f)
		regions = []
		skip = 0
		for i in range(len(genes)):
			if skip > i: continue
			beg = genes[i].beg - arg.padding
			end = genes[i].end + arg.padding
			if beg < 1: beg = 1
			if end > len(chrom.seq): end = len(chrom.seq)
			region = Feature(chrom, beg, end, '.', 'region', source='haman')
			for j in range(i + 1, len(genes)):
				if genes[j].overlap(region):
					if genes[j].end > region.end:
						region.end = genes[j].end
				else:
					skip = j
					break
			regions.append(region)
		
		# remap sequence and features smaller and export
		for region in regions:
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
			pc_genes = 0
			nc_genes = 0
			pc_issues = 0
			introns = 0
			strand = None
			rnaseq = 0
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
				if f.type == 'intron' and f.source == 'RNASeq_splice':
					rnaseq += 1
			meta = {
				'pc_genes': pc_genes,
				'nc_genes': nc_genes,
				'pc_issues': pc_issues,
				'introns': introns,
				'strand': strand,
				'rnaseq': rnaseq,
			}
			mfp = open(dir + '/' + str(idx) + '.json', 'w+')
			mfp.write(json.dumps(meta, indent=4))
			mfp.close()



