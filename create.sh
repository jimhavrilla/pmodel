#create.sh - generate nodoms and domains, also make filtered transcript file
set -e

COV=$1 # coverage factor to filter on

# remove utrs and introns

bedtools intersect -a $DATA/pfam.bed -b $DATA/exons_sans_utrs.bed -wb | awk {'print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$17'} FS='\t' OFS='\t' | perl uniq.pl | perl -pe 's/"|;//g' > $DATA/foo.bed

# get number one appris transcripts by length or randomness

python select-transcript.py <(grep -v ALTERNATIVE $DATA/appris_data.principal.txt) <(sed '1d' $DATA/transcriptlengths.txt) > $DATA/transcripts.txt

# use appris to remove non-canonical transcripts; sort by ENSL gene_id and Pfam autoreg to merge uniqids into single domain

awk 'NR==FNR{a[$2]}$29 in a{print $0}' $DATA/transcripts.txt $DATA/foo.bed | tr -s " " "\t" | sort -k15,15 -k25,25 -k1,1 -k2,2n -k3,3n > $DATA/blah; mv $DATA/blah $DATA/foo.bed

# domain coverage and rearranging: // based on histograms, used 5x as a filter

python make-doms.py <(awk '{if ($15==$29) print}' $DATA/foo.bed | tr -s " " "\t" | cut -f -45) | sort -k1,1 -k2,2n > $DATA/doms-any-coverage.bed

bedtools intersect -a <(cat $DATA/doms-any-coverage.bed | tr -s "\t" " " | cut -d " " -f 1-45 | tr -s " " "\t") -b <(awk -v cov=$COV '{if ($4>=cov) print}' $DATA/coverage.bed) -wa -wb -sorted \
| awk '{ct[$1 $2 $3 $45]++; len[$1 $2 $3 $45]=$4; row[$1 $2 $3 $45]=$0} END {for (i in ct) print row[i],(ct[i]==0 ? ct[i]=0: ct[i]), (len[i]==0 ? len[i]=1: len[i]),ct[i]/(len[i]==0 ? len[i]=1: len[i])}' \
| tr -s " " "\t" | cut -f 1-45,50- > $DATA/doms.bed

# nodoms
perl -pe 's/tag\s*\S*?(?=\n|\s)//g' $DATA/exons_sans_utrs.bed | perl -pe 's/ccds_id\s*\S*?(?=\n|\s)//g' | perl -pe 's/"|;//g' | awk 'NR==FNR{a[$2]}$10 in a{print}' $DATA/transcripts.txt - | tr -s " " "\t" | sort -k10,10 > $DATA/foo2
awk '{split($15,trans,","); for (i in trans) print $1,$2,$3,$13,$25,trans[i],$43}' $DATA/doms-any-coverage.bed | sort -k6,6 | tr -s " " "\t"  > $DATA/doms-split.bed
python nodom.py $DATA/foo2 $DATA/doms-split.bed > $DATA/nodoms-any-coverage.bed; 
rm $DATA/foo2 $DATA/foo3

bedtools intersect -a <(sort -k1,1 -k2,2n $DATA/nodoms-any-coverage.bed) -b <(awk '{if ($4>=5) print}' $DATA/coverage.bed) -wa -wb -sorted \
| awk '{ct[$1 $2 $3 $8 $24]++; len[$1 $2 $3 $8 $24]=$4; row[$1 $2 $3 $8 $24]=$0} END {for (i in row) print row[i],(ct[i]==0 ? ct[i]=0: ct[i]),(len[i]==0 ? len[i]=1: len[i]),ct[i]/(len[i]==0 ? len[i]=1: len[i])}' \
| tr -s " " "\t" | cut -f -24,29- | sort -k1,1 -k2,2n > $DATA/nodoms.bed
