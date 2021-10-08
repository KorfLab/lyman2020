#!/usr/bin/env python3

from grimoire.genome import Reader

cds = []
with open('isoset.txt') as fp:
	for line in fp.readlines():
		reg = line.rstrip()
		prefix = f'favorites/{reg}/{reg}'	
		genome = Reader(fasta=prefix+'.fa', gff=prefix+'.gff')
		chrom = genome.next()
		for gene in chrom.ftable.build_genes():
			tx = gene.transcripts()[0]
		lseq = ['g'] * len(chrom.seq)
		for exon in tx.exons:
			for i in range(exon.beg, exon.end+1):
				lseq[i-1] = 'e'
		for intron in tx.introns:
			for i in range(intron.beg, intron.end+1):
				lseq[i-1] = 'i'
		print(chrom.name, chrom.seq, ''.join(lseq))

		