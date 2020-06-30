#!/usr/bin/env python3

import json
import sys
import statistics

from grimoire.genome import Reader

def write_file(name, seqs):
	with open(name, 'w') as fp:
		for seq in seqs:
			fp.write(seq)
			fp.write('\n')

MIN_EXP = 10000
MAX_EXP = 1000000
FLANK = 20
WINDOW = FLANK * 2 + 2
CANONICAL = True
CUT_HI = 0.01

sd_hi = {}   # splice donors above threshold
sd_lo = {}   # splice donors below threshold
sd_fake = {} # fake splice donors
sa_hi = {}   # acceptors
sa_lo = {}   # acceptors
sa_fake = {} # acceptors

with open('isoset.txt') as fp:
	for line in fp.readlines():
		reg = line.rstrip()
		prefix = f'region/{reg}/{reg}'
		f = open(prefix + '.json')
		m = json.loads(f.read())
		f.close()
		exp = m['max_splices']
		
		if exp < MIN_EXP: continue
		if exp > MAX_EXP: continue		
		genome = Reader(fasta=prefix+'.fa', gff=prefix+'.gff')
		chrom = genome.next()

		# hi and lo (which contain some double-counts if they share sites)
		for f in chrom.ftable.features:
			if f.source == 'RNASeq_splice':
				s = f.seq_str(off5=-20, off3=20).upper()
				d = s[0:WINDOW]
				a = s[-WINDOW:]
				p = f.score / exp
				
				# filters
				if len(d) != WINDOW or len(a) != WINDOW:
					#sys.stderr.write(f'skipping edge {f} \n')
					continue
				gt = d[FLANK:FLANK+2]
				ag = a[FLANK:FLANK+2]
				if CANONICAL and (gt != 'GT' or ag != 'AG'):
					#sys.stderr.write(f'skipping non-canonical {f}\n')
					continue
				
				if p >= CUT_HI:
					sd_hi[d] = True
					sa_hi[a] = True
				else:
					sd_lo[d] = True
					sa_lo[a] = True
		
		# fake
		for f in chrom.ftable.features:
			if f.type == 'gene':
				s = f.seq_str().upper()
				for i in range(FLANK, len(s) - FLANK +1):
					n2 = s[i:i+2]
					if n2 == 'GT':
						d = s[i-FLANK:i+2+FLANK]
						if d not in sd_hi and d not in sd_lo:
							sd_fake[d] = True
					if n2 == 'AG':
						a = s[i-FLANK:i+2+FLANK]
						if a not in sa_hi and a not in sa_lo:
							sa_fake[a] = True

# prune hi from lo (double-counted introns)
for d in sd_hi:
	if d in sd_lo: del sd_lo[d]
for a in sa_hi:
	if a in sa_lo: del sa_lo[a]

# outputs
write_file('don.lo.true.txt', sd_lo.keys())
write_file('don.hi.true.txt', sd_hi.keys())
write_file('acc.lo.true.txt', sa_lo.keys())
write_file('acc.hi.true.txt', sa_hi.keys())
write_file('don.fake.txt', sd_fake.keys())
write_file('acc.fake.txt', sa_fake.keys())
