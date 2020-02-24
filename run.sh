#!/bin/bash

# run setup and selector

SIZE=`ls region | wc -w`
if [ $SIZE != 25571 ]; then
	python3 setup.py --fasta c_elegans.PRJNA13758.WS274.genomic.fa --gff c_elegans.PRJNA13758.WS274.annotations.gff3 --out region --padding 100 --source wb.270
	python3 selector.py --regions region --source wb.270 --report stats.txt --out table.txt --isomax 500
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

for reg in `cat table.txt | awk '{if ($4 == 1) print}' | cut -f 1`; do python3 vis_iso.py --regions region --region $reg --threshold 1e-4 --txt region/$reg/$reg.iscount --out region/$reg/$reg.pdf; done

# make figures or whatever
