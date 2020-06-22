#!/usr/bin/env python3

import argparse
import random
from spliceomatic import isoformer
from collections import OrderedDict

parser = argparse.ArgumentParser(
	description='Isoform generator')
parser.add_argument('--sites', required=True, type=str,
	metavar='<path>', help='file of donor and acceptor sites')
parser.add_argument('--transcripts', required=False, type=int, default=1000,
	metavar='<int>', help='number of transcripts to produce [%(default)i]')
parser.add_argument('--seed', required=False, type=int,
	metavar='<int>', help='set random seed')
arg = parser.parse_args()

if arg.seed: random.seed(arg.seed)

don = {}
acc = {}
with open(arg.sites) as fp:
	for line in fp.readlines():
		(type, pos, prob) = line.split()
		if   type == 'don': don[int(pos)] = float(prob)
		elif type == 'acc': acc[int(pos)] = float(prob)
		else: raise Exception('unknown type', type)

don = OrderedDict(sorted(don.items()))
acc = OrderedDict(sorted(acc.items()))

iso, iu, du, au = isoformer(don, acc, samples=arg.transcripts)

for k, v in sorted(iso.items(), key=lambda item: item[1], reverse=True):
	print(v, k)

