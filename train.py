#!/usr/bin/env python3

import argparse
import sys
import json
import copy

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
parser.add_argument('--weight', required=False, type=str,
	help='method for weighting sequence of gff objects')
arg = parser.parse_args()

###############
## Functions ##
###############

class MilwaError(Exception):
	pass


#############
## Globals ##
#############

## state arrays ##
don5 = [0, 2]    # add more later
don3 = [2, 6]    # add more later
don_ctx_max = 2  # maximum donor emission order
acc5 = [2, 8]    # add more later
acc3 = [0, 2]    # add more later
acc_ctx_max = 2  # maximum acceptor emission order
flank5 = [5, 10, 15]
flank3 = [5, 10, 15]

## single states ##
gen_ctx_max = 5
int_ctx_max = 5
exon_ctx_max = 5
cds_ctx_max = 5
u5_ctx_max = 3
u3_ctx_max = 3


gen_seqs = []
int_seqs = []
don_seqs = []
acc_seqs = []
exon_seqs = []
cds_seqs = []
up_seqs = []
down_seqs = []
utr5_seqs = []
utr3_seqs = []
gen_len = 0
exon_len = 0
splices = 0
intron_len = 0
mRNAs = 0

genome = Reader(fasta=arg.fasta, gff=arg.gff)
for chr in genome:
	genes = chr.ftable.build_genes()



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


"""		
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

