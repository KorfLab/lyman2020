#!/usr/bin/env python3

import json
import sys
import statistics

lengths = []
introns = []
splices = []
expression = []
set = []
with open('isoset.txt') as fp:
	for line in fp.readlines():
		reg = line.rstrip()
		f = open('region' + '/' + reg + '/' + reg + '.json')
		m = json.loads(f.read())
		f.close()
		lengths.append(int(m['length']))
		introns.append(int(m['introns']))
		splices.append(int(m['RNASeq_splice']))
		expression.append(int(m['max_splices']))
		set.append(m)

for reg in set: reg['score'] = 0

lmean = statistics.mean(lengths)
lsdev = statistics.stdev(lengths)
for reg in set: reg['score'] += abs((reg['length'] - lmean) / lsdev)

imean = statistics.mean(introns)
isdev = statistics.stdev(introns)
for reg in set: reg['score'] += abs((reg['introns'] - imean) / isdev)

smean = statistics.mean(splices)
ssdev = statistics.stdev(splices)
for reg in set: reg['score'] += abs((reg['RNASeq_splice'] - smean) / ssdev)

emean = statistics.mean(expression)
esdev = statistics.stdev(expression)
for reg in set: reg['score'] += abs((reg['max_splices'] - emean) / esdev)

print(f'region\tgene\tlength\tintrons\tsplices\texpression')
for reg in sorted(set, key=lambda r: r['score']):
	print(f"{reg['region']}\t{reg['gene_names'][0]}\t{reg['length']}\t{reg['introns']}\t{reg['RNASeq_splice']}\t{reg['max_splices']}\t{reg['score']}")

