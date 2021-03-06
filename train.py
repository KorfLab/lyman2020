#!/usr/bin/env python3

import argparse
import sys
import os
import json
#import copy

import grimoire.toolbox as toolbox
import grimoire.hmm as hmm
from grimoire.sequence import DNA
from grimoire.feature import Feature, FeatureTable
from grimoire.genome import Reader

## Command line stuff ##

extended_help = """
%(prog)s is used for training States and State Arrays that gene components
such as exons or introns. The results are structured as JSON files and
stored in the specified directory.
"""

parser = argparse.ArgumentParser(
	description='State and State Array trainer.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to fasta file')
parser.add_argument('--gff', required=True, type=str,
	metavar='<path>', help='path to GFF (or related) file')
parser.add_argument('--dir', required=True, type=str,
	metavar='<path>', help='path to directory containing JSON output files')
arg = parser.parse_args()

###############
## Functions ##
###############

class MilwaError(Exception):
	pass


#############
## Globals ##
#############

## single states ##

int_ctx_max = 5
exon_ctx_max = 5
u5_ctx_max = 3
u3_ctx_max = 3

gen_ctx_max = 5
null_ctx_max = 5
cds_ctx_max = 2

## state arrays ##
sparam = {
	'DON': {'o5':[0, 2], 'o3':[2, 6], 'ctx': [0, 1]},
	'ACC': {'o5':[2, 8], 'o3':[0, 2], 'ctx': [0, 1]},
	'KOZ': {'o5':[10],   'o3':[3],    'ctx':[0]},
	'TER': {'o5':[3],    'o3':[10],   'ctx':[0]},
}

def train_genomic(ctx, name):
	state = hmm.null_state_factory(file=arg.fasta, context=ctx)
	state.name = name
	path = '{}/{}-genomic-{}.json'.format(arg.dir, name, ctx)
	with open(path, 'w+') as file:
		file.write(state.to_json())

def train_state(genes, type, ctx):
	seqs = []
	for gene in genes:
		for tx in gene.transcripts():
			features = None
			if   type == 'INT':  features = tx.introns
			elif type == 'EXON': features = tx.exons
			elif type == 'UTR5': features = tx.utr5s
			elif type == 'UTR3': features = tx.utr3s
			for f in features:
				seqs.append({'seq' : f.seq_str(), 'weight' : 1})
					
	em = hmm.train_emission(seqs, context=ctx)
	state = hmm.State(name=type, context=ctx, emits=em)
	path = '{}/{}-gene_models-{}.json'.format(arg.dir, type, ctx)
	with open(path, 'w+') as file:
		file.write(state.to_json())

def train_cds_state(genes,  ctx):
	seqs = []
	for gene in genes:
		for tx in gene.transcripts():
			seqs.append({'seq' : tx.cds_str()[3:-3], 'weight' : 1})
			
	em = hmm.train_cds(seqs, context=ctx)
	state = hmm.State(name='CDS', context=ctx, emits=em)
	path = '{}/{}-gene_models-{}.json'.format(arg.dir, 'CDS', ctx)
	with open(path, 'w+') as file:
		file.write(state.to_json())

def train_state_array(genes, type, o5, o3, ctx):
	seqs = []
	for gene in genes:
		for tx in gene.transcripts():
			if type == 'DON':
				for intron in tx.introns:
					ilen = o5 + o3
					iseq = intron.seq_str(off5=-o5)[0:ilen]
					seqs.append({'seq':iseq, 'weight':1})
			elif type == 'ACC':
				for intron in tx.introns:
					ilen = o5 + o3
					iseq = intron.seq_str(off3=o3)[-ilen:]
					seqs.append({'seq':iseq, 'weight':1})
			elif type == 'KOZ':
				cds = None
				if tx.strand == '+': cds = tx.cdss[0]
				else:                cds = tx.cdss[-1]
				if o3 > cds.end - cds.beg +1:
					continue # split start codon
				seqs.append({'seq':cds.seq_str(off5=-o5)[0:o5+o3], 'weight':1})
			elif type == 'TER':
				cds = None
				if tx.strand == '+': cds = tx.cdss[-1]
				else:                cds = tx.cdss[0]
				if o5 > cds.end - cds.beg +1:
					continue # split stop codon
				s = cds.seq_str(off3=o3)[-(o5+o3):]
				seqs.append({'seq':s, 'weight':1})

	em = hmm.train_emissions(seqs, context=ctx)
	states = hmm.state_factory(type, em)
	
	path = '{}/{}-gene_models-{}-{}-{}.json'.format(arg.dir, type, o5, o3, ctx)
	with open(path, 'w+') as file:
		file.write(json.dumps(states, indent=4, cls=hmm.HMMdecoder))

if __name__ == '__main__':

	if not os.path.exists(arg.dir):
		os.mkdir(arg.dir)

	# stats
	stats = {
		'genomic_length'   : 0,
		'chromosome_count' : 0,
		'gene_count'       : 0,
		'gene_length'      : 0,
		'intron_length'    : 0,
		'intron_count'     : 0,
		'exon_length'      : 0,
		'exon_count'       : 0,
		'utr5_length'      : 0,
		'utr5_count'       : 0,
		'utr3_length'      : 0,
		'utr3_count'       : 0,
		'cds_length'       : 0,
		'cds_count'        : 0,
	}

	# genomic state
	for ctx in range(gen_ctx_max +1): train_genomic(ctx, 'GEN')
	
	# null state
	for ctx in range(null_ctx_max +1): train_genomic(ctx, 'NULL')
	
	# setup
	genes = []
	genome = Reader(fasta=arg.fasta, gff=arg.gff)
	for chr in genome:
		genes += chr.ftable.build_genes()
		stats['genomic_length'] += len(chr.seq)
		stats['chromosome_count'] += 1

	# simple states
	for ctx in range(int_ctx_max +1):  train_state(genes, 'INT', ctx)
	for ctx in range(exon_ctx_max +1): train_state(genes, 'EXON', ctx)
	for ctx in range(u5_ctx_max +1):   train_state(genes, 'UTR5', ctx)
	for ctx in range(u3_ctx_max +1):   train_state(genes, 'UTR3', ctx)
	
	# coding states
	for ctx in range(cds_ctx_max +1): train_cds_state(genes, ctx)	
	
	# state arrays
	for name in sparam:
		for o5 in sparam[name]['o5']:
			for o3 in sparam[name]['o3']:
				for ctx in sparam[name]['ctx']:
					train_state_array(genes, name, o5, o3, ctx)
	
	# stats
	for gene in genes:
		stats['gene_count'] += 1
		stats['gene_length'] += gene.end - gene.beg + 1
		for tx in gene.transcripts():
			for f in tx.introns:
				stats['intron_length'] += f.end - f.beg + 1
				stats['intron_count'] += 1
			for f in tx.exons:
				stats['exon_length'] += f.end - f.beg + 1
				stats['exon_count'] += 1
			for f in tx.utr5s:
				stats['utr5_length'] += f.end - f.beg + 1
				stats['utr5_count'] += 1
			for f in tx.utr3s:
				stats['utr3_length'] += f.end - f.beg + 1
				stats['utr3_count'] += 1
			stats['cds_count'] += 1
			stats['cds_length'] += len(tx.cds_str())
	path = '{}/{}'.format(arg.dir, 'stats.json')
	with open(path, 'w+') as file:
		file.write(json.dumps(stats, indent=4))
	
# are sequences getting reverse-complimented before counting? how smart were we?
	


"""	
## 
	
	
	# donors





	gen_seqs.append({'seq' : chr.seq, 'weight' : 1})
	for gene in genes:
		if gene.issues and arg.canonical: continue
		if not gene.transcripts(): continue
		for tx in gene.transcripts():
			if len(tx.exons) < 2: continue
			mRNAs += 1
			for exon in tx.exons:
				exon_seq = exon.seq_str()
				exon_len += len(exon_seq)
				exon_seqs.append({'seq' : exon_seq, 'weight' : 1})
			for intron in tx.introns:
				iseq = intron.seq_str()
				continue
				don_seq = iseq[0:arg.don_len]
				int_seq = iseq[arg.don_len:-arg.acc_len]
				acc_seq = iseq[-arg.acc_len:len(iseq)]
				don_seqs.append({'seq' : don_seq, 'weight' : 1})
				int_seqs.append({'seq' : int_seq, 'weight' : 1})
				acc_seqs.append({'seq' : acc_seq, 'weight' : 1})
				intron_len += len(int_seq)
				splices += 1




	
	gen_emits  = hmm.train_emission(gen_seqs,   context=arg.gen_ctx)
	exon_emits = hmm.train_emission(exon_seqs,  context=arg.exon_ctx)
	don_emits  = hmm.train_emissions(don_seqs,  context=arg.don_ctx)
	int_emits  = hmm.train_emission(int_seqs,   context=arg.int_ctx)
	acc_emits  = hmm.train_emissions(acc_seqs,  context=arg.acc_ctx)
	
	gen_state = hmm.State(name='GEN', context=arg.gen_ctx, emits=gen_emits)
	exon_state = hmm.State(name='EXON', context=arg.exon_ctx, emits=exon_emits)
	don_states = hmm.state_factory('DON', don_emits)
	int_state = hmm.State(name='INT', context=arg.int_ctx, emits=int_emits)
	acc_states = hmm.state_factory('ACC', acc_emits)
	
	gen_state.init = 1
	gen_state.term = 1
	
	hmm.connect2(gen_state, gen_state, 0.999)
	hmm.connect2(gen_state, exon_state, 0.001)
	hmm.connect2(exon_state, exon_state, 0.95)
	hmm.connect2(exon_state, gen_state, 0.01)
	hmm.connect2(exon_state, don_states[0], 0.04)
	hmm.connect_all(don_states)
	hmm.connect2(don_states[-1], int_state, 1)
	hmm.connect2(int_state, int_state, 0.99)
	hmm.connect2(int_state, acc_states[0], 0.01)
	hmm.connect_all(acc_states)
	hmm.connect2(acc_states[-1], exon_state, 1)
	
	null_state = hmm.null_state_factory(file=arg.fasta, context=arg.null_ctx)
	
	model = hmm.HMM(name=arg.hmm, null=null_state,
		states=[gen_state] + [exon_state] + don_states + [int_state] + acc_states)

"""

