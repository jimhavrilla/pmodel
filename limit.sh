bash query.sh $1/lolli.db '%' $2 $3 $4| python formatvar.py | sort -k6,6 -k8,8 > $1/allint2maf.$2.$3-$4.bed
python maketable.py $1/allint2maf.$2.$3-$4.bed > $1/foo.bed
python mergetable.py -f $1/foo.bed $1/human_pfam.counts $1/sumlist.bed > $1/dtablemaf.$2.$3-$4.txt; rm $1/foo.bed