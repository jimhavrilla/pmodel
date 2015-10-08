DATA=$1

exp=$(cut -f 2 $DATA/appris_data.principal.txt | sort -u | wc -l | cut -d " " -f 1)
obs=$(wc -l < $DATA/transcripts.txt)

echo "pmodel.t1... test that number of uniq genes is equal to uniq appris transcripts"
if [[ $exp != $obs ]]; then
	echo "FAIL the number of lines in transcripts should equal the number of uniq genes from appris"
else
	echo "    OK"
fi
