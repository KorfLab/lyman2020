#!/bin/bash

python3 setup.py --fasta c_elegans.PRJNA13758.WS270.genomic.fa --gff c_elegans.PRJNA13758.WS270.annotations.gff3 --out r2 --padding 100 --source wb.270

python3 selector.py --regions region --source wb.270 --report stats.txt --out table.txt

cat table.txt | awk '{if ($4 == 1) print}' | head -n1000 > training.txt

cut -f2 training.txt > training_gids.txt

grep -vf training_gids.txt  table.txt > testing.txt

rm training_gids.txt

