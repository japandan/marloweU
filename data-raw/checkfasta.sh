#!/bin/bash
Taxon=$(zgrep -m1 ">" $1|awk -F= '{print $2}')
Taxon=${Taxon::-3}
TaxonID=$(zgrep -m1 ">" $1|awk -F= '{print $3}')
TaxonID=${TaxonID::-3}
TaxonCount=$(zgrep  "OX=$TaxonID" $1|wc -l)
proteins=$(zgrep ">" $1|wc -l)
#echo "Filename,Taxon,TaxonID,Proteins,TaxonCount"
echo "$1,Taxon:$Taxon,TaxonID:$TaxonID,$proteins,$TaxonCount"
if [ $TaxonCount -ne $proteins ]; then
	echo "Error: Identical TaxonID not present for file, $1"
	mv $1 errors/
fi
