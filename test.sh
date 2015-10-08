DATA=$1

exp=$(cut -f 2 $DATA/appris_data.principal.txt | sort -u | wc -l | cut -d " " -f 1)
obs=$(wc -l < $DATA/transcripts.txt)

echo "pmodel.t1... test that number of uniq genes before filtering on maximum length is equal to the number in transcripts.txt"
if [[ $exp != $obs ]]; then
	echo "FAIL the number of lines in transcripts should equal the number of uniq genes from appris"
else
	echo "    OK"
fi


echo "pmodel.t2... test that the nodoms don't intersect the doms"

obs=$(bedtools intersect -a $DATA/doms-any-coverage.bed -b $DATA/nodoms-any-coverage.bed -wo \
	| python test-different-transcripts.py)
exp=""

if [[ "$exp" != "$obs" ]]; then
	echo "FAIL"
else
	echo "    OK"
fi
