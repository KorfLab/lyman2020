#!/usr/bin/env python3

import argparse
import sys
import json
import os

import grimoire.hmm as hmm

## Command line stuff ##

extended_help = """
%(prog)s is used for assemblingHMMs from previously trained States and State
Arrays. The HMMs are structured as JSON files.
"""

parser = argparse.ArgumentParser(
	description='HMM builder genes and their components.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument('--dir', required=True, type=str,
	metavar='<path>', help='path to directory containing trained states')
parser.add_argument('--model', required=True, type=str,
	metavar='<model>', help='acc|don|exon|intron|mRNA1|mRNA2|gene')
parser.add_argument('--hmm', required=True, type=str,
	metavar='<path>', help='output HMM file')
parser.add_argument('--null_ctx', required=False, type=int, default=None,
	metavar='<int>', help='null model context [%(default)d]')
parser.add_argument('--acc_off5', required=False, type=int, default=None,
	metavar='<int>', help='acceptor 5 prime offset  [%(default)d]')
parser.add_argument('--acc_off3', required=False, type=int, default=None,
	metavar='<int>', help='acceptor 3 prime offset [%(default)d]')
parser.add_argument('--acc_ctx', required=False, type=int, default=None,
	metavar='<int>', help='acceptor context [%(default)d]')
parser.add_argument('--don_off5', required=False, type=int, default=None,
	metavar='<int>', help='donor 5 prime offset [%(default)d]')
parser.add_argument('--don_off3', required=False, type=int, default=None,
	metavar='<int>', help='don 3 prime offset [%(default)d]')
parser.add_argument('--don_ctx', required=False, type=int, default=None,
	metavar='<int>', help='donor context [%(default)d]')
parser.add_argument('--exon_ctx', required=False, type=int, default=None,
	metavar='<int>', help='exon context [%(default)d]')
parser.add_argument('--gen_ctx', required=False, type=int, default=None,
	metavar='<int>', help='genomic context [%(default)d]')
parser.add_argument('--int_ctx', required=False, type=int, default=None,
	metavar='<int>', help='intron context [%(default)d]')
parser.add_argument('--u5_ctx', required=False, type=int, default=None,
	metavar='<int>', help='UTR5 context [%(default)d]')
parser.add_argument('--u3_ctx', required=False, type=int, default=None,
	metavar='<int>', help='UTR3 context [%(default)d]')
parser.add_argument('--koz_off5', required=False, type=int, default=None,
	metavar='<int>', help='Kozak 5 prime offset [%(default)d]')
parser.add_argument('--koz_off3', required=False, type=int, default=None,
	metavar='<int>', help='Kozak 3 prime offset [%(default)d]')
parser.add_argument('--koz_ctx', required=False, type=int, default=None,
	metavar='<int>', help='Kozak context [%(default)d]')
parser.add_argument('--cds_ctx', required=False, type=int, default=None,
	metavar='<int>', help='CDS context [%(default)d]')
parser.add_argument('--ter_off5', required=False, type=int, default=None,
	metavar='<int>', help='transcription termination 5 prime offset [%(default)d]')
parser.add_argument('--ter_off3', required=False, type=int, default=None,
	metavar='<int>', help='transcription termination 3 prime offset [%(default)d]')
parser.add_argument('--ter_ctx', required=False, type=int, default=None,
	metavar='<int>', help='transcription termination context [%(default)d]')
parser.add_argument('--weight_method', required=False, type=str, default=None,
	metavar='<str>', help=' [%(default)d]')
#parser.add_argument('--weight_feature', required=False, type=str, default=None,
#	metavar='<str>', action='append', help='adds a feature type to be weighted. user may use this option to add any combination of DON, ACC, and INT. All features will be weighted using the same weight method.')
parser.add_argument('--weight_features', required=False, type=str,
	default=None, nargs='+', metavar='<str>',
	help='space-delimited list of features to weight. DON, ACC, and INT are currently implemented.')
arg = parser.parse_args()

###############
## Functions ##
###############

class MilwaError(Exception):
	pass

def get_parameters(type):
	(off5, off3, ctx) = (None, None, None)
	if   type == 'NULL': ctx = arg.null_ctx
	elif type == 'GEN':  ctx = arg.gen_ctx
	elif type == 'UTR5': ctx  = arg.u5_ctx
	elif type == 'EXON': ctx = arg.exon_ctx
	elif type == 'CDS':  ctx = arg.cds_ctx
	elif type == 'INT':  ctx = arg.int_ctx
	elif type == 'UTR3': ctx = arg.u3_ctx
	elif type == 'KOZ':
		off5 = arg.koz_off5
		off3 = arg.koz_off3
		ctx  = arg.koz_ctx
	elif type == 'DON':
		off5 = arg.don_off5
		off3 = arg.don_off3
		ctx  = arg.don_ctx
	elif type == 'ACC':
		off5 = arg.acc_off5
		off3 = arg.acc_off3
		ctx  = arg.acc_ctx
	elif type == 'TER':
		off5 = arg.ter_off5
		off3 = arg.ter_off3
		ctx  = arg.ter_ctx
	else:
		raise MilwaError('unrecognized state name', state)
	return (off5, off3, ctx)

def get_path(type):
	path = None
	(off5, off3, ctx) = get_parameters(type)
	if off5 != None and off3 != None:
		if arg.weight_method and arg.weight_features and type in arg.weight_features:
			path = '{}/{}-{}-{}-{}-{}.json'.format(arg.dir, type, off5, off3, ctx, arg.weight_method)
		else:
			path = '{}/{}-{}-{}-{}.json'.format(arg.dir, type, off5, off3, ctx)
	elif ctx != None:
		if arg.weight_method and arg.weight_features and type in arg.weight_features:
			path = '{}/{}-{}-{}.json'.format(arg.dir, type, ctx, arg.weight_method)
		else:
			path = '{}/{}-{}.json'.format(arg.dir, type, ctx)
	else:
		raise MilwaError('missing parameters for {} state(s)'.format(type))
	return path
	
def check_training(states):
	errors = []
	for state in states:
		path = get_path(state)
		if not os.path.exists(path):
			errors.append('parameter combination not trained for {} state(s)'.format(state))
	if errors:
		raise MilwaError('\n'.join(errors))

def read_state(type):
	path = get_path(type)
	state = None
	with open(path) as file:
		state = hmm.State.from_json(file.read())# not working
	return state

def read_states(type):
	path = get_path(type)
	states = []
	with open(path) as file:
		for jstr in json.loads(file.read()):
			states.append(hmm.State.from_json(json.dumps(jstr)))# also not working
	return states

############################
## Model Building Section ##
############################
if not os.path.exists(arg.dir + '/stats.json'):
	raise MilwaError('{}/stats.json missing'.format(arg.dir))
stats = None
with open(arg.dir + '/stats.json') as file:
	stats = json.loads(file.read())
file.close()

if arg.weight_method and not os.path.exists('{}/stats-{}.json'.format(arg.dir, arg.weight_method)):
		raise MilwaError('{}/stats-{}.json missing'.format(arg.dir, arg.weight_method))
wstats = None
if os.path.exists('{}/stats-{}.json'.format(arg.dir, arg.weight_method)):
	with open('{}/stats-{}.json'.format(arg.dir, arg.weight_method)) as file:
		wstats = json.loads(file.read())
	file.close()

model = None

if arg.model == 'don':
	states = ['INT', 'DON']
	check_training(states)
#	exon_state = read_state('EXON')
	don_states = read_states('DON')
	for state in don_states:
		print(state.name)
	sys.exit(1)
	int_state = read_state('INT')
#	null_state = read_state('NULL')
#	hmm.connect2(exon_state, exon_state, 1 - stats['exon_count']/stats['exon_length'])
#	hmm.connect2(exon_state, don_states[0], stats['exon_count']/stats['exon_length'])
	hmm.connect_all(don_states)
	hmm.connect2(don_states[-1], int_state, 1)
	hmm.connect2(int_state, int_state, 1)
	model = hmm.HMM(name=arg.hmm, null=None,
		states= don_states + [int_state])
	
else:
	raise MilwaError('unknown model type: ' + arg.model)


""""
elif arg.model == 'acc':
	
	hmm.connect2(int_state, int_state, 1 - splices/intron_len)
	hmm.connect2(int_state, acc_states[0], splices/intron_len)
	hmm.connect_all(acc_states)
	hmm.connect2(acc_states[arg.acc_len-1], exn_state, 1)
	hmm.connect2(exn_state, exn_state, 1)
	
	null_state = hmm.null_state_factory(file=arg.fasta, context=arg.null_ctx)
	
	model = hmm.HMM(name=arg.hmm, null=null_state,
		states=[int_state] + acc_states + [exn_state])

elif arg.model == 'exon':
	
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
"""

model.write(arg.hmm)
