#!/bin/bash

# run setup 

SIZE=`ls region | wc -w`
if [ $SIZE != 21056 ]; then
	python3 regioner.py --fasta ce.fa.gz --gff ce.gff3.gz
fi

if [ ! -e isoset.txt ]; then
	python3 ikselect.py > isoset.txt
fi
