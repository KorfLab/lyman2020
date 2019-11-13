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
%(prog)s is used for training States and State Arrays representing gene
components such as exons or introns from RNA-Seq data. The results are
structured as JSON files and stored in the specified directory. This program
does not produce States or State Arrays trained on gene models, which must
be generated separately. Currently, only training for introns and splice
sites is implemented.

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
parser.add_argument('--weight', required=True, type=str,
	metavar='<str>', help='method for weighting sequence of gff objects')
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

## state arrays ##
sparam = {
	'DON': {'o5':[0, 2], 'o3':[2, 6], 'ctx': [0, 1]},
	'ACC': {'o5':[2, 8], 'o3':[0, 2], 'ctx': [0, 1]},
}

def train_state(features, type, ctx):
	seqs = []
	for f in features:
		if type == 'INT' and arg.weight == 'wb.270':
			if f.type == 'intron' and f.source == 'RNASeq_splice':
				seqs.append({'seq' : f.seq_str(), 'weight' : f.score})
		else:
			raise NotImplementedError('feature or weight type not yet supported')			

	em = hmm.train_emission(seqs, context=ctx)
	state = hmm.State(name=type, context=ctx, emits=em)
	path = '{}/{}-{}-{}.json'.format(arg.dir, type, ctx, arg.weight)
	with open(path, 'w+') as file:
		file.write(state.to_json())

def train_state_array(features, type, o5, o3, ctx):
	seqs = []
	for f in features:
		if type == 'DON' and arg.weight == 'wb.270':
			if f.type == 'intron' and f.source == 'RNASeq_splice':
					iseq = f.seq_str(off5=-o5)[0:o5+o3]
					seqs.append({'seq':iseq, 'weight':f.score})
		elif type == 'ACC' and arg.weight == 'wb.270':
			if f.type == 'intron' and f.source == 'RNASeq_splice':
				iseq = f.seq_str(off3=o3)[-(o5 + o3):]
				seqs.append({'seq':iseq, 'weight':f.score})
		else:
			raise NotImplementedError('feature or weight type not yet supported')
	
	em = hmm.train_emissions(seqs, context=ctx)
	states = hmm.state_factory(type, em)
	
	path = '{}/{}-{}-{}-{}-{}.json'.format(arg.dir, type, o5, o3, ctx, arg.weight)
	with open(path, 'w+') as file:
		file.write(json.dumps(states, indent=4, cls=hmm.HMMdecoder))

if __name__ == '__main__':

	if not os.path.exists(arg.dir):
		os.mkdir(arg.dir)

	# stats
	stats = {
		'intron_length'    : 0,
		'intron_count'     : 0,
	}

	genome = Reader(fasta=arg.fasta, gff=arg.gff)
	flist = []
	for chrom in genome:
		flist += chrom.ftable.features
		
	# simple states
	for ctx in range(int_ctx_max +1):  train_state(flist, 'INT', ctx)
	
	# state arrays
	for name in sparam:
		for o5 in sparam[name]['o5']:
			for o3 in sparam[name]['o3']:
				for ctx in sparam[name]['ctx']:
					train_state_array(flist, name, o5, o3, ctx)
	
	# stats
	for f in flist:
		if arg.weight == 'wb.270':
			if f.type == 'intron' and f.source == 'RNASeq_splice':
				stats['intron_length'] += f.end - f.beg + 1
				stats['intron_count']  += f.score
	
	path = '{}/{}-{}.json'.format(arg.dir, 'stats', arg.weight)
	with open(path, 'w+') as file:
		file.write(json.dumps(stats, indent=4))
	