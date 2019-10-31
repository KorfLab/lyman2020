#!/bin/bash

python3 setup.py --fasta c_elegans.PRJNA13758.WS270.genomic.fa --gff c_elegans.PRJNA13758.WS270.annotations.gff3 --out region --padding 100 --source wb.270

python3 selector.py --regions region --source wb.270 --report stats.txt --out table.txt --isomax 50
