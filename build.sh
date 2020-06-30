
while read p; do
	cp -r region/$p favorites
done < isoset.txt
