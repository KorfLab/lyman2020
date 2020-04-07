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
parser.add_argument('--out', required=False, type=str,
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

isoforms = []
with open(ipath) as ifile:
	gff = GFF_file(file=prefix + '.gff')
	start_intrs = False
	mn = math.inf
	mx = 0
	for line in ifile:
		line = line[:-2] # remove \n
		percent = 100 * float(line.split('\t')[0])
		introns = [x.split(',') for x in line.split('\t')[1].split(' ')]
		min_freq = math.inf
		if introns[0][0] == '[none]':
			isoforms.append((introns, percent))
			continue
		for i in range(len(introns)):
			introns[i] = [int(x) for x in introns[i]]
		isoforms.append((introns, percent))
		first, last = introns[0][0], introns[len(introns) - 1][1]
		mn = first if first < mn else mn
		mx = last if last > mx else mx


isoforms.sort(key=lambda x:x[1])
ymax = len(isoforms) * (height + 5) + y_level
fig, axis = plot.subplots(1)
axis.set( ylim=(0, ymax))
axis.set_yticklabels([])
axis.set_yticks([])
axis.set( xlim=(gff.get(type='exon', end=mn-1)[0].beg-pad, gff.get(type='exon', beg=mx+1)[0].end+pad))
percent_x = plot.gca().get_xlim()[1] + 10

if arg.txt:
	txtfh = open(arg.txt, 'w+')

for i in range(len(isoforms)):
	introns = isoforms[i][0]
	percent = isoforms[i][1]
	length = len(introns)
	rects = []
	boost = (height + 5) * i
	if(introns[0][0] == '[none]'):
		start = gff.get(type='exon', end=mn-1)[0]
		end = gff.get(type='exon', beg=mx+1)[0]
	else:
		start = gff.get(type='exon', end=introns[0][0]-1)[0]
		end = gff.get(type='exon', beg=introns[length-1][1]+1)[0]
	
	if arg.txt:
		txtfh.write(str(isoforms[i][1]) + '\n')
	
	if arg.out:	
		if(introns[0][0] == '[none]'):
			last = patches.Rectangle((start.beg, y_level + boost), end.end - start.beg, height, facecolor=color)
			axis.add_patch(last)
			percent_y = last.get_y()
		else:	
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
			percent_y = last.get_y()
		
		plot.text(percent_x, percent_y, '{:.3f}%'.format(isoforms[i][1]), fontsize=6)
		plot.plot([last.get_x() + last.get_width(), percent_x], [last.get_y() + last.get_height()/2, last.get_y() + last.get_height()/2],
		linestyle='dotted', color='grey')

if arg.out:
	plot.savefig(arg.out, dpi=400)
