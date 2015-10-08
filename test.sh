DATA=$1


check(){
	obs=$1
	exp=$2

if [[ "$exp" != "$obs" ]]; then
	echo "FAIL"
else
	echo "    OK"
fi
}

echo "pmodel.t1... test that number of uniq genes before filtering on maximum length is equal to the number in transcripts.txt"
exp=$(cut -f 2 $DATA/appris_data.principal.txt | sort -u | wc -l | cut -d " " -f 1)
obs=$(wc -l < $DATA/transcripts.txt)
check $obs $exp


echo "pmodel.t2... test that the nodoms don't intersect the doms"

obs=$(bedtools intersect -a $DATA/doms-any-coverage.bed -b $DATA/nodoms-any-coverage.bed -wo \
	| python test-different-transcripts.py)
exp=""
check $obs $exp


echo "pmodel.t3... test that exons-(doms+nodoms) gives empty results"

obs=$(bedtools subtract -a <(cat $DATA/doms-any-coverage.bed $DATA/nodoms-any-coverage.bed|cut -f 1-22) -b $DATA/exons_sans_utrs.bed)
exp=""
check $obs $exp

echo "pmodel.t4... test that (doms+nodoms)-exons gives empty results"
obs=$(bedtools subtract -b <(cat $DATA/doms-any-coverage.bed $DATA/nodoms-any-coverage.bed | cut -f 1-22) -a <(cut -f 2 $DATA/transcripts.txt | grep -wFf - $DATA/exons_sans_utrs.bed))
exp=""
check $obs $exp
