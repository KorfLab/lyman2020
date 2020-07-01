#!/usr/bin/env python3


import json
import os
import sys
import statistics

from grimoire.genome import Reader
from grimoire.toolbox import translate_str

def write_file(name, seqs):
	with open(name, 'w') as fp:
		for seq in seqs:
			fp.write(seq)
			fp.write('\n')

gfp = open('gen.fa', 'w')
tfp = open('tx.fa', 'w')
cfp = open('cds.fa', 'w')
pfp = open('pro.fa', 'w')

for region in os.listdir('favorites'):

	prefix = 'favorites/' + region + '/' + region
	f = open(prefix + '.json')
	meta = json.loads(f.read())
	f.close()

	genome = Reader(fasta=prefix+'.fa', gff=prefix+'.gff')
	chrom = genome.next()
	gene = chrom.ftable.build_genes()[0]
	
	tx = gene.transcripts()[0]
	cds = tx.cds_str()
	pro = translate_str(cds.upper())

	gfp.write(f'>{region}\n{chrom.seq}\n')
	tfp.write(f'>{region}\n{tx.seq_str()}\n')
	cfp.write(f'>{region}\n{cds}\n')
	pfp.write(f'>{region}\n{pro}\n')
	
