#!/usr/bin/env python3

import sys
import os
import json

from grimoire.io import GFF_file
from grimoire.sequence import DNA
from grimoire.feature import Feature, mRNA, Gene, FeatureTable
from grimoire.genome import Reader
from grimoire.toolbox import translate_str

EXP_LO = 1e3
EXP_HI = 1e6
PRO_LO = 100
PRO_HI = 1000
ALEN = 70
PCTID = 70

complex = 0
non_coding = 0
no_intron = 0
issues = 0
support = 0
hi_exp = 0
lo_exp = 0
mult_tx = 0
pro_short = 0
pro_long = 0

seqs = []
for region in os.listdir('region'):

	prefix = 'region/' + region + '/' + region
	f = open(prefix + '.json')
	meta = json.loads(f.read())
	f.close()

	# skips
	if meta['genes'] != 1:
		complex += 1
		continue
	if meta['pc_genes'] != 1:
		non_coding += 1
		continue
	if meta['introns'] == 0:
		no_intron += 1
		continue
	if meta['pc_issues']:
		issues += 1
		continue
	if meta['max_splices'] < EXP_LO:
		lo_exp += 1
		continue
	if meta['max_splices'] > EXP_HI:
		hi_exp += 1
		continue
	
	# process genes
	gfile = prefix + '.gff'
	genome = Reader(gff=prefix + '.gff', fasta=prefix + '.fa', source='wb.270')
	for chrom in genome:
		gene = chrom.ftable.build_genes()[0] # there is only ever one gene
		
		isoforms = len(gene.transcripts())
		if isoforms != 1:
			mult_tx += 1
			continue
		
		intron_support = {}
		for intron in chrom.ftable.features:
			if intron.source == 'RNASeq_splice' and intron.strand == gene.strand:
				intron_support[(intron.beg, intron.end)] = True

		supported = 0
		nosupport = 0
		for tx in gene.transcripts():
			isoforms += 1
			for intron in tx.introns:
				it = (intron.beg, intron.end)
				if it in intron_support: supported += 1
				else:                    nosupport += 1
		if nosupport > 0:
			support += 1
			continue

		tx = gene.transcripts()[0]
		cds = tx.cds_str()
		pro = translate_str(cds.upper())
		
		if len(pro) < PRO_LO:
			pro_short += 1
			continue
		if len(pro) > PRO_HI:
			pro_long += 1
			continue
		
		seqs.append((region, pro))

# BLASTP analysis - single linkage clustering to remove near duplicates

with open('temp.fa', 'w') as fp:
	for reg, pro in seqs:
		fp.write(f'>{reg}\n{pro}\n')

os.system('xdformat -p temp.fa')
os.system('blastp temp.fa temp.fa W=5 mformat=3 > temp.blastp')

cluster = {}
used = {}
with open('temp.blastp') as fp:
	for line in fp.readlines():
		if line.startswith('#'): continue
		f = line.split()
		if f[0] == f[1]: continue # self match	
		qid = f[0]
		sid = f[1]
		alen = int(f[6])
		pct_id = float(f[10])
		if alen >= ALEN and pct_id > PCTID:
			if sid in used: continue
			if qid not in cluster:
				cluster[qid] = []
				used[qid] = True
			cluster[qid].append(sid)
			used[sid] = True

remove = {}
clustered = 0
for parent in cluster:
	for child in cluster[parent]:
		remove[child] = True
		clustered += 1

# Done, output section

with open('isoset.txt', 'w') as fp:
	for reg, pro in seqs:
		if reg not in remove:
			fp.write(f'{reg}\n')


with open('skipped.txt', 'w') as fp:
	fp.write(f'Complex Locus: {complex}\n')
	fp.write(f'Non-coding: {non_coding}\n')
	fp.write(f'No introns: {no_intron}\n')
	fp.write(f'Non-canonical issues: {issues}\n')
	fp.write(f'Low expression: {lo_exp}\n')
	fp.write(f'High expression: {hi_exp}\n')
	fp.write(f'Multiple isoforms: {mult_tx}\n')
	fp.write(f'Lacking RNA support: {support}\n')
	fp.write(f'Short protein: {pro_short}\n')
	fp.write(f'Long protein: {pro_long}\n')
	fp.write(f'Too similar to others: {clustered}\n')
	