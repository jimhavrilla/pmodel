#tables.sh - make regions tables, and soon, gene tables

# sort domain occurrence count from bill, remove weird carriage return characters that screw things up

sort -k2,2 $DATA/human_pfam.counts | perl -pe 's/\r$//g'  > $DATA/blah; mv $DATA/blah $DATA/human_pfam.counts

# by uniqid a table

sed '1d' $DATA/allintfilter.bed | python table.py > $DATA/regionstable.txt

#finds non-intersecting regions

cat <(bedtools intersect -a <(awk '{$13=$33; print $0}' OFS="\t" $DATA/alluniq.bed | cut -f 1,2,3,11,13,25,27,43,45,47,48 | python lencount.py | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -v -sorted | cut -f 1,2,3,4,5,6,7,8,9,10,14,15,18) \
<(bedtools intersect -a <(cut -f 1,2,3,12,24,26,27 $DATA/nodom.bed | awk '{t=$8;$8=$7;$7=t;print}' | tr -s " " "\t" | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -v -sorted | cut -f 1,2,3,4,5,6,7,8,12,13,16 | awk '{print $1,$2,$3,".",$4,$5,$5,$6,$7,$8}') > $DATA/nointregions.txt

# adds noint regions and pfam occurrences to file

python noints.py $DATA/regionstable.txt $DATA/nointregions.txt
awk 'NR==FNR{a[$2]=$3} NR!=FNR{if ($1 in a) print $0"\t"a[$1]; else print $0"\t"1}' $DATA/human_pfam.counts $DATA/regionstable.txt > $DATA/blah; mv $DATA/blah $DATA/regionstable.txt

# make MAF spanning file

cat <(gawk 'NR==FNR{a[$5 $6][$1 $2 $3]=$1 " " $2 " " $3 " " $4} NR!=FNR{if ($3 $2 in a) {for (i in a[$3 $2]) print a[$3 $2][i],$0}}' <(cut -f 1,2,3,15,25,33,48 $DATA/alluniq.bed) <(cut -d " " -f 1,2,3,7,8,9,10,11 $DATA/regionstable.txt)) \
<(gawk 'NR==FNR{a[$5 $6][$1 $2 $3]=$1 " " $2 " " $3 " " $4} NR!=FNR{if ($3 $2 in a) {for (i in a[$3 $2]) print a[$3 $2][i],$0}}' <(cut -f 1,2,3,8,12,24,27 $DATA/nodom.bed | awk '{t=$6; $6=$5; $5=t; print}' OFS='\t') <(cut -d " " -f 1,2,3,7,8,9,10,11 $DATA/regionstable.txt)) \
| tr -s " " "\t" | sort -k1,1 -k2,2n > $DATA/regioncoordsdnds.bed

sort -k1,1 -k2,2n $DATA/allintfilter.bed > $DATA/blah; mv $DATA/blah $DATA/allintfilter.bed

#also adds noint regions

cat <(bedtools intersect -a $DATA/regioncoordsdnds.bed -b $DATA/allintfilter.bed -wa -wb -sorted | awk '{if ($7==$23) print}' | cut -f -16,28,29,30) <(bedtools intersect -a $DATA/regioncoordsdnds.bed -b $DATA/allintfilter.bed -v -sorted | awk '{print $0"\t.\t.\t."}') \
| bedtools groupby -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 -c 17,18,19 -o collapse,collapse,collapse | sort -k1,1 -k6,6 -k5,5 -k7,7 \
| bedtools groupby -i - -g 1,6,5,7 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,3 -o distinct,min,max,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,collapse,collapse,collapse,collapse,collapse | cut -f 5-26 \
| python recompute.py | cat <(printf "#chr\tstart\tend\ttranscript\tdomain\tgene\tautoregs(uniqid for non-domain regions)\tcov_ratio\tlength\tdn\tds\tna\tdn/ds\tdensity\tfvrv\tprevalence\tmafs\timpacts\ttype\tstarts\tends\tmaf_modifier\n") - | bgzip > $DATA/regionsmafsdnds.bed.gz 