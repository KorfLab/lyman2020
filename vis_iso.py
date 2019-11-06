import matplotlib.pyplot as plot
import matplotlib.patches as patches
import math
import argparse
import os

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
	metavar='<path>', help='region containing desired gene')
parser.add_argument('--gid', required=True, type=str,
	metavar='<path>', help='the id of the desired gene')
parser.add_argument('--threshold', required=True, type=str,
	metavar='<path>', help='frequency threshold: 1e-4, 1e-6, or 1e-8')
parser.add_argument('--out', required=True, type=str,
	metavar='<path>', help='filename of output graphic')
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

ftr_target = 'region: ' + arg.region + ' gid: ' + arg.gid
freq_target = '@' + arg.threshold + ' freq paths:'

isoforms = []
with open(ipath) as ifile:
    start_ftr = False
    start_intrs = False
    mn = math.inf
    mx = 0
    for line in ifile:
        line = line[:-1] # remove \n
        if(start_intrs):
            if(line[0] == '@'):
                break
            introns = [x.split(',') for x in line.split(' ')]
            for i in range(len(introns)):
                introns[i] = [int(x) for x in introns[i]]
            isoforms.append(introns)
            first, last = introns[0][0], introns[len(introns) - 1][1]
            mn = first if first < mn else mn
            mx = last if last > mx else mx
        if(start_ftr and line == freq_target):
            start_intrs = True
        if(line == ftr_target):
            start_ftr = True
# todo: coordinates
with open(fapath) as fafile:
	first = fafile.readline()
	start = int(first.split(':')[2].split(' ')[0])
	end = int(first.split(':')[3])

start = mn - pad
end = mx + pad
ymax = len(isoforms) * (height + 5) + y_level
fig, axis = plot.subplots(1)
axis.set(xlim=(start-pad, end+pad), ylim=(0, ymax))
axis.set_yticklabels([])
axis.set_yticks([])
for i in range(len(isoforms)):
    introns = isoforms[i]
    length = len(introns)
    rects = []
    boost = (height + 5) * i
    rects.append(patches.Rectangle((start, y_level + boost), introns[0][0]-start, height, facecolor=color))
    for j in range(len(introns)-1):
        beg = introns[j][1]
        rects.append(patches.Rectangle((beg, y_level + boost), introns[j+1][0] - beg, height, facecolor=color))
    rects.append(patches.Rectangle((introns[length-1][1], y_level + boost), end - introns[length-1][1], height, facecolor=color))
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
    axis.add_patch(rects[len(rects)-1])

plot.savefig(arg.out, dpi=400)