#!/bin/bash

# run setup and selector

SIZE=`ls region | wc -w`
if [ $SIZE != 27356 ]; then
	python3 setup.py --fasta c_elegans.PRJNA13758.WS270.genomic.fa --gff c_elegans.PRJNA13758.WS270.annotations.gff3 --out region --padding 100 --source wb.270
	python3 selector.py --regions region --source wb.270 --report stats.txt --out table.txt --isomax 50
fi

# build testing and training sets (add conditional later)

cat table.txt | awk '{if ($4 == 1) print}' | head -n1000 > training.txt
cut -f2 training.txt > training_gids.txt
grep -vf training_gids.txt table.txt | grep -v ^region > testing.txt
rm training_gids.txt
rm -f training.gff training.fa testing.gff testing.fa
for reg in `cut -f1 training.txt`; do
	cat region/$reg/$reg.gff >> training.gff
	cat region/$reg/$reg.fa >> training.fa
done
for reg in `cut -f1 testing.txt`; do
	cat region/$reg/$reg.gff >> testing.gff
	cat region/$reg/$reg.fa >> testing.fa
done

# make figures or whatever
