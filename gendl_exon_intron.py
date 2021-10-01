#!/usr/bin/env python3

import random
from grimoire.genome import Reader
from grimoire.toolbox import revcomp_str

def write_file(name, stub, seqs):
	with open(name, 'w') as fp:
		for i, seq in enumerate(seqs):
			fp.write(f'>{stub}:{i}\n{seq}\n')

cds = []
with open('isoset.txt') as fp:
	for line in fp.readlines():
		reg = line.rstrip()
		prefix = f'favorites/{reg}/{reg}'	
		genome = Reader(fasta=prefix+'.fa', gff=prefix+'.gff')
		chrom = genome.next()
		for gene in chrom.ftable.build_genes():
			tx = gene.transcripts()[0]
			s = tx.cds_str()
			#s2 = ''
			#for i in range(len(s)):
			#	if i % 3 == 0:
			#		s2 += s[i]
			#	else:
			#		s2 += s[i].lower()
			cds.append(s)

random.seed(1)

SIZE = 51

# file 1: initial exon
atg = [s[0:SIZE] for s in cds]

# file 2: phase 0 cds
SAMPLES = 5
ph0 = []
for s in cds:
	for n in range(SAMPLES):
		limit = (len(s) - SIZE) / 3
		r = random.randint(0, limit)
		r *= 3
		ph0.append(s[r:r+SIZE])

# file 3 and 4: any phase cds (and reverse-complement)
phx = []
phy = []
for s in cds:
	for n in range(SAMPLES):
		limit = len(s) - SIZE
		r = random.randint(0, limit)
		plus = s[r:r+SIZE]
		anti = revcomp_str(plus)
		phx.append(plus)
		phy.append(anti)

write_file('cdsi.fa', 'cdsi', atg)
write_file('cds0.fa', 'cds0', ph0)
write_file('cds+.fa', 'cds+', phx)
write_file('cds-.fa', 'cds-', phy)
