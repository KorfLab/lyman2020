import matplotlib.pyplot as plot
import matplotlib.patches as patches
import matplotlib.figure as figure
import math
import random
import argparse
import numpy as np
import os
from grimoire.io import GFF_file

extended_help = """
Creates a wire diagram of isoforms of a specified gene.
"""

parser = argparse.ArgumentParser(
	description='Diagrams isoforms of a gene.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument('--regions', required=True, type=str,
	metavar='<path>', help='directory containing region sub-directories')
parser.add_argument('--region', required=True, type=str,
	metavar='<num>', help='region containing desired gene')
parser.add_argument('--threshold', required=True, type=str,
	metavar='<str>', help='frequency threshold: 1e-4, 1e-6, or 1e-8')
parser.add_argument('--out', required=True, type=str,
	metavar='<path>', help='filename of output graphic (.pdf, .png, .svg ... )')
parser.add_argument('--txt', required=False, type=str,
	metavar='<path>', help='text file of isoform counts')
arg = parser.parse_args()

if not os.path.exists(arg.regions):
	raise Exception('invalid path to region directories')

prefix = arg.regions + '/' + arg.region + '/' + arg.region
ipath = prefix + '.isoforms'
fapath = prefix + '.fa'
color = (1, .498, .01)
y_level = 10
pad = 60
height = 3

freq_target = arg.threshold + ' freq paths:'

isoforms = []
with open(ipath) as ifile:
	gff = GFF_file(file=prefix + '.gff')
	start_intrs = False
	mn = math.inf
	mx = 0
	for line in ifile:
		line = line[:-1] # remove \n
		if(start_intrs):
			if('e-' in line):
				break
			introns = [x.split(',') for x in line.split(' ')]
			min_freq = math.inf
			for i in range(len(introns)):
				introns[i] = [int(x) for x in introns[i]]
				results = gff.get(type='intron', beg=introns[i][0], end=introns[i][1])
				results = [x for x in results if (x.beg == introns[i][0] and x.end == introns[i][1] and x.score != '.')]
				curr_freq = int(float(results[0].score))
				min_freq = curr_freq if curr_freq < min_freq else min_freq
			isoforms.append((introns, min_freq))
			first, last = introns[0][0], introns[len(introns) - 1][1]
			mn = first if first < mn else mn
			mx = last if last > mx else mx
		if(line == freq_target):
			start_intrs = True


isoforms.sort(key=lambda x:x[1])
ymax = len(isoforms) * (height + 5) + y_level
fig, axis = plot.subplots(1)
axis.set( ylim=(0, ymax))
axis.set_yticklabels([])
axis.set_yticks([])
axis.set( xlim=(gff.get(type='exon', end=mn-1)[0].beg-pad, gff.get(type='exon', beg=mx+1)[0].end+pad))
freq_x = plot.gca().get_xlim()[1] + 10

if arg.txt:
	txtfh = open(arg.txt, 'w+')

for i in range(len(isoforms)):
	introns = isoforms[i][0]
	freq = isoforms[i][1]
	length = len(introns)
	rects = []
	boost = (height + 5) * i
	start = gff.get(type='exon', end=introns[0][0]-1)[0]
	end = gff.get(type='exon', beg=introns[length-1][1]+1)[0]
	
	rects.append(patches.Rectangle((start.beg, y_level + boost), start.end - start.beg, height, facecolor=color))
	for j in range(len(introns)-1):
		beg = introns[j][1]
		rects.append(patches.Rectangle((beg, y_level + boost), introns[j+1][0] - beg, height, facecolor=color))
	rects.append(patches.Rectangle((introns[length-1][1], y_level + boost), end.end - end.beg, height, facecolor=color))
	for j in range(1, len(rects)):
		left_x = rects[j-1].get_x() + rects[j-1].get_width()
		left_y = rects[j-1].get_y() + rects[j-1].get_height()
		right_x = rects[j].get_x()
		right_y = rects[j].get_y() + rects[j].get_height()
		mid_x = (right_x + left_x) / 2
		mid_y = right_y + 2
		axis.add_patch(rects[j-1])
		plot.plot([left_x, mid_x],[left_y, mid_y],'k-', linewidth=0.5)
		plot.plot([mid_x, right_x],[mid_y, right_y],'k-', linewidth=0.5)
	last = rects[len(rects)-1]
	axis.add_patch(last)
	freq_y = last.get_y()
	plot.text(freq_x, freq_y, str(isoforms[i][1]), fontsize=6)
	if arg.txt:
		txtfh.write(str(isoforms[i][1]) + '\n')
	plot.plot([last.get_x() + last.get_width(), freq_x], [last.get_y() + last.get_height()/2, last.get_y() + last.get_height()/2],
	linestyle='dotted', color='grey')

#plot.show()

plot.savefig(arg.out, dpi=400)
