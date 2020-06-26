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

for region in os.listdir('region'):

	prefix = 'region/' + region + '/' + region
	f = open(prefix + '.json')
	meta = json.loads(f.read())
	f.close()

	# skips
	if meta['genes'] != 1: continue
	if meta['pc_genes'] != 1: continue
	if meta['introns'] == 0: continue
	if meta['pc_issues']: continue
	
	# process genes
	gfile = prefix + '.gff'
	genome = Reader(gff=prefix + '.gff', fasta=prefix + '.fa', source='wb.270')
	for chrom in genome:
		gene = chrom.ftable.build_genes()[0] # there is only ever one gene
		
		intron_support = {}
		for intron in chrom.ftable.features:
			if intron.source == 'RNASeq_splice' and intron.strand == gene.strand:
				intron_support[(intron.beg, intron.end)] = True

		isoforms = 0
		supported = 0
		nosupport = 0
		for tx in gene.transcripts():
			isoforms += 1
			for intron in tx.introns:
				it = (intron.beg, intron.end)
				if it in intron_support: supported += 1
				else:                    nosupport += 1
		if nosupport > 0: continue
		print(region)

