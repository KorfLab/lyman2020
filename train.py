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
parser.add_argument('--source', required=False, type=str,
	metavar='<name>', help='output sources to <name>.fasta, <name>.gff')
parser.add_argument('--weighted', required=False, type=str,
	help='method for weighting sequence of gff objects')

arg = parser.parse_args()

###############
## Functions ##
###############

class MilwaError(Exception):
	pass

def make_source(name, comps):
	seq = ''
	features = []
	for comp in comps: seq += comp['seq']
	dna = DNA(seq=seq, name=name)
	beg, end = 1, None
	for comp in comps:
		end = beg + len(comp['seq']) - 1
		features.append(Feature(dna, beg, end, '+', comp['label']))
		beg += len(comp['seq'])
	return {'dna':dna, 'features':features}

#############
## Globals ##
#############

sources = [] # used for outputting sources and performing replicant analysis
gen_seqs = []
int_seqs = []
don_seqs = []
acc_seqs = []
ex_seqs = []
cds_seqs = []
up_seqs = []
down_seqs = []
utr5_seqs = []
utr3_seqs = []
exon_len = 0
splices = 0
intron_len = 0

genome = Reader(fasta=arg.fasta, gff=arg.gff)


for chr in genome:
	genes = chr.ftable.build_genes()	
	for gene in genes:
		if gene.issues: continue
		if not gene.transcripts(): continue
			txs = []
			txs = gene.transcripts()
			for tx in txs:
				if len(tx.exons) < 2: continue
				rcount = 0
				for i in range(len(tx.exons) - 1):
					exp_seq = tx.exons[i].seq_str()
					iseq = tx.introns[i].seq_str()
					don_seq = iseq[0:arg.don_len]
					int_seq = iseq[arg.don_len:-arg.acc_len]
					exp_seqs.append({'seq' : exp_seq, 'weight' : 1})
					don_seqs.append({'seq' : don_seq, 'weight' : 1})
					int_seqs.append({'seq' : int_seq, 'weight' : 1})
					intron_len += len(int_seq)
					exon_len += len(exp_seq)
					splices += 1
					if arg.source:
						rcount += 1
						sources.append(make_source(tx.id + '-region-' + str(rcount),
							[{'label':'EXON', 'seq':exp_seq},
							{'label':'DON', 'seq':don_seq},
							{'label':'INT', 'seq':int_seq}]))
			
	exp_emits = hmm.train_emission(exp_seqs, context=arg.exon_ctx)
	don_emits = hmm.train_emissions(don_seqs, context=arg.don_ctx)
	int_emits = hmm.train_emission(int_seqs, context=arg.int_ctx)
	
	exp_state = hmm.State(name='EXON', context=arg.exon_ctx, emits=exp_emits)
	exp_state.init = 1
	don_states = hmm.state_factory('DON', don_emits)
	int_state = hmm.State(name='INT', context=arg.int_ctx, emits=int_emits)
	int_state.term = 1
	
	hmm.connect2(exp_state, exp_state, 1 - splices/exon_len)
	hmm.connect2(exp_state, don_states[0], splices/exon_len)
	hmm.connect_all(don_states)
	hmm.connect2(don_states[arg.don_len-1], int_state, 1)
	hmm.connect2(int_state, int_state, 1)
	
	null_state = hmm.null_state_factory(file=arg.fasta, context=arg.null_ctx)
	
	model = hmm.HMM(name=arg.hmm, null=null_state,
		states=[exp_state] + don_states + [int_state])

elif arg.model == 'acc':
	genome = Reader(fasta=arg.fasta, gff=arg.gff)
	exn_seqs = []
	acc_seqs = []
	int_seqs = []
	exon_len = 0
	splices = 0
	intron_len = 0
	for chr in genome:
		genes = chr.ftable.build_genes()	
		for gene in genes:
			if gene.issues and arg.canonical: continue
			if not gene.transcripts(): continue
			txs = []
			if arg.first: txs.append(gene.transcripts()[0])
			else: txs = gene.transcripts()
			for tx in txs:
				if len(tx.exons) < 2: continue
				rcount = 0
				for i in range(len(tx.exons) - 1):
					exn_seq = tx.exons[i + 1].seq_str()
					iseq = tx.introns[i].seq_str()
					int_seq = iseq[arg.don_len:-arg.acc_len]
					acc_seq = iseq[-arg.acc_len:len(iseq)]
					int_seqs.append({'seq' : int_seq, 'weight' : 1})
					acc_seqs.append({'seq' : acc_seq, 'weight' : 1})
					exn_seqs.append({'seq' : exn_seq, 'weight' : 1})
					intron_len += len(int_seq)
					exon_len += len(exn_seq)
					splices += 1
					if arg.source:
						rcount += 1
						sources.append(make_source(tx.id + '-region-' + str(rcount),
							[{'label':'INT', 'seq':int_seq},
							{'label':'ACC', 'seq':acc_seq},
							{'label':'EXON', 'seq':exn_seq}]))
	
	int_emits = hmm.train_emission(int_seqs, context=arg.int_ctx)
	acc_emits = hmm.train_emissions(acc_seqs, context=arg.acc_ctx)
	exn_emits = hmm.train_emission(exn_seqs, context=arg.exon_ctx)
	
	int_state = hmm.State(name='INT', context=arg.int_ctx, emits=int_emits)
	int_state.init = 1
	acc_states = hmm.state_factory('ACC', acc_emits)
	exn_state = hmm.State(name='EXON', context=arg.exon_ctx, emits=exn_emits)
	exn_state.term = 1
	
	hmm.connect2(int_state, int_state, 1 - splices/intron_len)
	hmm.connect2(int_state, acc_states[0], splices/intron_len)
	hmm.connect_all(acc_states)
	hmm.connect2(acc_states[arg.acc_len-1], exn_state, 1)
	hmm.connect2(exn_state, exn_state, 1)
	
	null_state = hmm.null_state_factory(file=arg.fasta, context=arg.null_ctx)
	
	model = hmm.HMM(name=arg.hmm, null=null_state,
		states=[int_state] + acc_states + [exn_state])

elif arg.model == 'exon':
	genome = Reader(fasta=arg.fasta, gff=arg.gff)
	intron_seqs = []
	acc_seqs = []
	don_seqs = []
	exon_seqs = []
	exon_len = 0
	splices = 0
	inta_len = 0
	
	for chr in genome:
		genes = chr.ftable.build_genes()	
		for gene in genes:
			if gene.issues and arg.canonical: continue
			if not gene.transcripts(): continue
			txs = []
			if arg.first: txs.append(gene.transcripts()[0])
			else: txs = gene.transcripts()
			for tx in txs:
				if len(tx.exons) < 2: continue
				rcount = 0
				for i in range(1, len(tx.exons) -1):
					iprev = tx.introns[i-1].seq_str()
					inext = tx.introns[i].seq_str()
					acc_seq = iprev[-arg.acc_len:len(iprev)]
					exon_seq = tx.exons[i].seq_str()
					don_seq = inext[0:arg.don_len]
					ip_seq = iprev[0:-len(acc_seq)]
					in_seq = inext[len(don_seq):]
					acc_seqs.append({'seq' : acc_seq, 'weight' : 1})
					exon_seqs.append({'seq' : exon_seq, 'weight' : 1})
					don_seqs.append({'seq' : don_seq, 'weight' : 1})
					intron_seqs.append({'seq' : ip_seq, 'weight' : 1})
					intron_seqs.append({'seq' : in_seq, 'weight' : 1})
					splices += 1
					exon_len += tx.exons[i].end - tx.exons[i].beg + 1
					inta_len += len(in_seq)
					if arg.source:
						rcount += 1
						sources.append(make_source(tx.id + '-region-' + str(rcount), [
							{'label':'INTa', 'seq':ip_seq},
							{'label':'ACC',  'seq':acc_seq},
							{'label':'EXON', 'seq':exon_seq},
							{'label':'DON',  'seq':don_seq},
							{'label':'INTb', 'seq':in_seq}]))

	acc_emits = hmm.train_emissions(acc_seqs, context=arg.acc_ctx)
	don_emits = hmm.train_emissions(don_seqs, context=arg.don_ctx)
	exon_emits = hmm.train_emission(exon_seqs, context=arg.exon_ctx)
	inta_emits = hmm.train_emission(intron_seqs, context=arg.int_ctx)
	intb_emits = hmm.train_emission(intron_seqs, context=arg.int_ctx)
	
	acc_states = hmm.state_factory('ACC', acc_emits)
	don_states = hmm.state_factory('DON', don_emits)
	exon_state = hmm.State(name='EXON', context=arg.exon_ctx, emits=exon_emits)
	inta_state = hmm.State(name='INTa', context=arg.int_ctx, emits=inta_emits)
	intb_state = hmm.State(name='INTb', context=arg.int_ctx, emits=intb_emits)

	inta_state.init = 1
	intb_state.term = 1
	
	hmm.connect2(inta_state, inta_state, 1 - splices/inta_len)
	hmm.connect2(inta_state, acc_states[0], splices/inta_len)
	hmm.connect_all(acc_states)
	hmm.connect2(acc_states[-1], exon_state, 1)
	hmm.connect2(exon_state, exon_state, 1 - splices/exon_len)
	hmm.connect2(exon_state, don_states[0], splices/exon_len)
	hmm.connect_all(don_states)
	hmm.connect2(don_states[-1], intb_state, 1)
	hmm.connect2(intb_state, intb_state, 1)
	
	null_state = hmm.null_state_factory(file=arg.fasta, context=arg.null_ctx)
	
	model = hmm.HMM(name=arg.hmm, states=[inta_state] + acc_states + [exon_state]
		+ don_states + [intb_state], null = null_state)

elif arg.model == 'intron':
	genome = Reader(fasta=arg.fasta, gff=arg.gff)
	exp_seqs = []
	exn_seqs = []
	don_seqs = []
	acc_seqs = []
	int_seqs = []
	exon_len = 0
	splices = 0
	intron_len = 0
	
	for chr in genome:
		genes = chr.ftable.build_genes()	
		for gene in genes:
			if gene.issues and arg.canonical: continue
			if not gene.transcripts(): continue
			txs = []
			if arg.first: txs.append(gene.transcripts()[0])
			else: txs = gene.transcripts()
			for tx in txs:
				if len(tx.exons) < 2: continue
				rcount = 0
				for i in range(len(tx.exons) - 1):
					exp_seq = tx.exons[i].seq_str()
					exn_seq = tx.exons[i + 1].seq_str()
					iseq = tx.introns[i].seq_str()
					don_seq = iseq[0:arg.don_len]
					int_seq = iseq[arg.don_len:-arg.acc_len]
					acc_seq = iseq[-arg.acc_len:len(iseq)]
					exp_seqs.append({'seq' : exp_seq, 'weight' : 1})
					don_seqs.append({'seq' : edon_seq, 'weight' : 1})
					int_seqs.append({'seq' : eint_seq, 'weight' : 1})
					acc_seqs.append({'seq' : eacc_seq, 'weight' : 1})
					exn_seqs.append({'seq' : eexn_seq, 'weight' : 1})
					intron_len += len(int_seq)
					exon_len += len(exp_seq)
					splices += 1
					if arg.source:
						rcount += 1
						sources.append(make_source(tx.id + '-region-' + str(rcount),
							[{'label':'EXP', 'seq':exp_seq},
							{'label':'DON', 'seq':don_seq},
							{'label':'INT', 'seq':int_seq},
							{'label':'ACC', 'seq':acc_seq},
							{'label':'EXN', 'seq':exn_seq}]))
	
	exp_emits = hmm.train_emission(exp_seqs, context=arg.exon_ctx)
	don_emits = hmm.train_emissions(don_seqs, context=arg.don_ctx)
	int_emits = hmm.train_emission(int_seqs, context=arg.int_ctx)
	acc_emits = hmm.train_emissions(acc_seqs, context=arg.acc_ctx)
	exn_emits = hmm.train_emission(exn_seqs, context=arg.exon_ctx)
	
	exp_state = hmm.State(name='EXP', context=arg.exon_ctx, emits=exp_emits)
	exp_state.init = 1
	don_states = hmm.state_factory('DON', don_emits)
	int_state = hmm.State(name='INT', context=arg.int_ctx, emits=int_emits)
	acc_states = hmm.state_factory('ACC', acc_emits)
	exn_state = hmm.State(name='EXN', context=arg.exon_ctx, emits=exn_emits)
	exn_state.term = 1
	
	hmm.connect2(exp_state, exp_state, 1 - splices/exon_len)
	hmm.connect2(exp_state, don_states[0], splices/exon_len)
	hmm.connect_all(don_states)
	hmm.connect2(don_states[arg.don_len-1], int_state, 1)
	hmm.connect2(int_state, int_state, 1 - splices/intron_len)
	hmm.connect2(int_state, acc_states[0], splices/intron_len)
	hmm.connect_all(acc_states)
	hmm.connect2(acc_states[arg.acc_len-1], exn_state, 1)
	hmm.connect2(exn_state, exn_state, 1)
	
	null_state = hmm.null_state_factory(file=arg.fasta, context=arg.null_ctx)
	
	model = hmm.HMM(name=arg.hmm, null=null_state,
		states=[exp_state] + don_states + [int_state] + acc_states + [exn_state])

elif arg.model == 'gene':
	if arg.source:
		raise MilwaError('source not available in gene')
		
	genome = Reader(fasta=arg.fasta, gff=arg.gff)

	gen_seqs = []
	exon_seqs = []
	don_seqs = []
	int_seqs = []
	acc_seqs = []
	
	gen_len = 0
	exon_len = 0
	splices = 0
	intron_len = 0
	mRNAs = 0
	
	for chr in genome:
		genes = chr.ftable.build_genes()	
		gen_seqs.append({'seq' : echr.seq, 'weight' : 1})
		for gene in genes:
			if gene.issues and arg.canonical: continue
			if not gene.transcripts(): continue
			txs = []
			if arg.first: txs.append(gene.transcripts()[0])
			else: txs = gene.transcripts()
			for tx in txs:
				if len(tx.exons) < 2: continue
				mRNAs += 1
				for exon in tx.exons:
					exon_seq = exon.seq_str()
					exon_len += len(exon_seq)
					exon_seqs.append({'seq' : eexon_seq, 'weight' : 1})
				for intron in tx.introns:
					iseq = intron.seq_str()
					don_seq = iseq[0:arg.don_len]
					int_seq = iseq[arg.don_len:-arg.acc_len]
					acc_seq = iseq[-arg.acc_len:len(iseq)]
					don_seqs.append({'seq' : edon_seq, 'weight' : 1})
					int_seqs.append({'seq' : eint_seq, 'weight' : 1})
					acc_seqs.append({'seq' : eacc_seq, 'weight' : 1})
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

elif arg.model == 'mRNA':
	genome = Reader(fasta=arg.fasta, gff=arg.gff)
	gen_seqs = []
	ut5_seqs = []
	koz_seqs = []
	atg_seqs = []
	cds_seqs = []
	stp_seqs = []
	ut3_seqs = []
	
	ut5_len = 0
	mRNAs = 0
	cds_len = 0
	ut3_len = []
	
	for chr in genome:
		genes = chr.ftable.build_genes()	
		gen_seqs.append({'seq' : chr.seq, 'weight' : 1})
		for gene in genes:
			if gene.issues and arg.canonical: continue
			if not gene.transcripts(): continue
			txs = []
			if arg.first: txs.append(gene.transcripts()[0])
			else: txs = gene.transcripts()
			for tx in txs:
				cds = tx.cds_str()
				ptx = tx.tx_str()
				beg = ptx.find(cds)
				end = beg + len(cds)
				if beg < arg.u5_ctx + arg.koz_len + 1: continue
				if end >= len(ptx): continue
				ut5_seq = ptx[0 : beg - arg.koz_len]
				koz_seq = ptx[beg-arg.koz_len:beg]
				atg_seq = cds[0:3]
				cds_seq = cds[3:-3]
				stp_seq = cds[-3:]
				ut3_seq = ptx[end:len(ptx)]
				ut5_seqs.append({'seq' : chr.sequt5_seq, 'weight' : 1})
				koz_seqs.append({'seq' : chr.seqkoz_seq, 'weight' : 1})
				atg_seqs.append({'seq' : chr.seqatg_seq, 'weight' : 1})
				cds_seqs.append({'seq' : chr.seqcds_seq, 'weight' : 1})
				stp_seqs.append({'seq' : chr.seqstp_seq, 'weight' : 1})
				ut3_seqs.append({'seq' : chr.sequt3_seq, 'weight' : 1})
				ut5_len += len(ut5_seq)
				cds_len += len(cds_seq)
				ut3_len == len(ut3_seq)
				mRNAs += 1
				if arg.source:
					sources.append(make_source(tx.id + '-cDNA',
						[{'label':'UTR5', 'seq':ut5_seq},
						{'label':'KOZ', 'seq':koz_seq},
						{'label':'ATG', 'seq':atg_seq},
						{'label':'CDS', 'seq':cds_seq},
						{'label':'STOP', 'seq':stp_seq},
						{'label':'UTR3', 'seq':ut3_seq}]))
		
	ut5_emits = hmm.train_emission(ut5_seqs, context=arg.u5_ctx)
	koz_emits = hmm.train_emissions(koz_seqs, context=arg.koz_ctx)
	atg_emits = hmm.train_emissions(atg_seqs, context=arg.atg_ctx)
	cds_emits = hmm.train_cds(cds_seqs, context=arg.cds_ctx)
	ut3_emits = hmm.train_emission(ut3_seqs, context=arg.u3_ctx)
	stp_emits = hmm.train_emissions(stp_seqs, context=arg.stop_ctx)
	
	u5_state = hmm.State(name='UTR5', context=arg.u5_ctx, emits=ut5_emits)
	koz_states = hmm.state_factory('KOZ', koz_emits)
	atg_states = hmm.state_factory('ATG', atg_emits)
	cds_states = hmm.state_factory('CDS', cds_emits)
	u3_state = hmm.State(name='UTR3', context=arg.u3_ctx, emits=ut3_emits)
	stp_states = hmm.state_factory('STOP', stp_emits)
	u5_state.init = 1
	u3_state.term = 1

	hmm.connect2(u5_state, u5_state, 1 - mRNAs/ut5_len)
	hmm.connect2(u5_state, koz_states[0], mRNAs/ut5_len)
	hmm.connect_all(koz_states)
	hmm.connect2(koz_states[-1], atg_states[0], 1)
	hmm.connect_all(atg_states)
	hmm.connect2(atg_states[-1], cds_states[0], 1)
	hmm.connect_all(cds_states)
	hmm.connect2(cds_states[2], cds_states[0], 1 - mRNAs / (cds_len/3))
	hmm.connect2(cds_states[2], stp_states[0], mRNAs / (cds_len/3))
	hmm.connect_all(stp_states)
	hmm.connect2(stp_states[-1], u3_state, 1)
	hmm.connect2(u3_state, u3_state, 1)
	
	null_state = hmm.null_state_factory(file=arg.fasta, context=arg.null_ctx)
	
	model = hmm.HMM(name=arg.hmm, states=[u5_state] + koz_states + atg_states
		+ cds_states + [u3_state] + stp_states, null=null_state)

else:
	raise MilwaError('unknown model type: ' + arg.model)

########################
## Post-build Actions ##
########################

if arg.hmm:
	model.write(arg.hmm)

if arg.source:
	with open(arg.source + '.fa', 'w+') as file:
		for source in sources:
			file.write(source['dna'].fasta())
	with open(arg.source + '.gff', 'w+') as file:
		for source in sources:
			for f in source['features']:
				file.write(f.gff() + '\n')

"""
if arg.replicant:
	model.convert2log()
	for source in sources:
		v = hmm.Viterbi(model=model, dna=source['dna'])
		p = v.generate_path()
		truth = FeatureTable(dna=source['dna'], features=source['features'])
		preds = FeatureTable(dna=source['dna'], features=p.features())
		-- some kind of standard comparision needed --
"""

