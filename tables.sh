#!/bin/bash
#create.sh - generate nodoms and domains, also make filtered transcript file
set -e


## BEGIN SCRIPT
usage()
{
    cat << EOF

usage: $0 OPTIONS

OPTIONS can be:
    -d      Data directory

EOF
}

# Show usage when there are no arguments.
if test -z "$1"
then
    usage
    exit
fi

DATA=

# Check options passed in.
while getopts "h d:c:" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        d)
            DATA=$OPTARG
            ;;
        c)
            COV=$OPTARG
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

if [ -z "$DATA" ]
then
    echo "Data directory not set"
    usage
    exit
fi

if [ ! -d "$DATA" ]
then
    echo "Data directory does not exist"
    usage
    exit
fi

#MAF=$1

# sort domain occurrence count from bill, remove weird carriage return characters that screw things up

sort -k2,2 $DATA/human_pfam.counts | perl -pe 's/\r$//g'  > $DATA/blah; mv $DATA/blah $DATA/human_pfam.counts

# by uniqid a table

cat $DATA/doms.fix.sort.bed.NI.VEPEXAC3filter.vcf.gz.bed \
| python append_fields.py \
    -f chr,start,end,POS,REF,ALT,transcript_id,pfamA_id,pfamA_auto_reg,uniq_id,len,ratio_of_cover_to_total,gene_name,maf,impact,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM \
> $DATA/doms.fix.sort.bed.NI.VEPEXAC3filter.vcf.gz.bed.xtra_cols.bed

cat $DATA/nodoms.fix.sort.bed.NI.VEPEXAC3filter.vcf.gz.bed \
| python append_fields.py -f chr,start,end,POS,REF,ALT,transcript_id,pfamA_id,pfamA_auto_reg,uniq_id,len,ratio_of_cover_to_total,gene_name,maf,impact,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM \
> $DATA/nodoms.fix.sort.bed.NI.VEPEXAC3filter.vcf.gz.bed.xtra_cols.bed

(cat $DATA/allintfilter.bed; \
 tail -n+2 $DATA/doms.fix.sort.bed.NI.VEPEXAC3filter.vcf.gz.bed.xtra_cols.bed;
 tail -n+2 $DATA/nodoms.fix.sort.bed.NI.VEPEXAC3filter.vcf.gz.bed.xtra_cols.bed) \
| bedtools sort -header \
> $DATA/merged_allintfilter_doms.NI_nodoms.NI.bed

cat $DATA/merged_allintfilter_doms.NI_nodoms.NI.bed \
| python table.py \
    -p $DATA/human_pfam.counts \
> $DATA/doms_nodoms_exac.bed

(head -n 1 $DATA/doms_nodoms_exac.bed; \
 tail -n+2 $DATA/doms_nodoms_exac.bed | sort -k 4,4 -k 2,2n) \
| gzip \
> $DATA/doms_nodoms_exac.bed.gz

#sed '1d' $DATA/allintfilter.bed | python table.py > $DATA/regionstable.txt
#
##finds non-intersecting regions
# 
#cat <(bedtools intersect -a <(cut -f 1,2,3,11,15,25,27,33,43,45,47,48 $DATA/doms.bed | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -v -sorted | cut -f 1,2,3,4,6,8,10,11,12) \
#<(bedtools intersect -a <(cut -f 1,2,3,12,24,26,27 $DATA/nodoms.bed | awk '{t=$4;$4=$5;$5=t;print}' | tr -s " " "\t" | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -v -sorted | cut -f 1,2,3,4,5,6,7,8,12,13,16 | awk '{print $1,$2,$3,".",$4,$5,$4,$6,$7,$8}') > $DATA/nointregions.txt
#
## adds noint regions and pfam occurrences to file
#
#python noints.py $DATA/regionstable.txt $DATA/nointregions.txt
#awk 'NR==FNR{a[$2]=$3} NR!=FNR{if ($1 in a) print $0"\t"a[$1]; else print $0"\t"1}' $DATA/human_pfam.counts $DATA/regionstable.txt > $DATA/blah; mv $DATA/blah $DATA/regionstable.txt
#
## make MAF spanning file
#
#cat <(gawk 'NR==FNR{a[$5 $6][$1 $2 $3]=$1 " " $2 " " $3 " " $4} NR!=FNR{if ($3 $2 in a) {for (i in a[$3 $2]) print a[$3 $2][i],$0}}' <(cut -f 1,2,3,15,25,33,48 $DATA/doms.bed) <(cut -d " " -f 1,2,3,7,8,9,10,11 $DATA/regionstable.txt)) \
#<(gawk 'NR==FNR{a[$5 $6][$1 $2 $3]=$1 " " $2 " " $3 " " $4} NR!=FNR{if ($3 $2 in a) {for (i in a[$3 $2]) print a[$3 $2][i],$0}}' <(cut -f 1,2,3,8,12,24,27 $DATA/nodoms.bed | awk '{t=$6; $6=$5; $5=t; print}' OFS='\t') <(cut -d " " -f 1,2,3,7,8,9,10,11 $DATA/regionstable.txt)) \
#| tr -s " " "\t" | sort -k1,1 -k2,2n > $DATA/regioncoordsdnds.bed
#
#sort -k1,1 -k2,2n $DATA/allintfilter.bed > $DATA/blah; mv $DATA/blah $DATA/allintfilter.bed
#
##also adds noint regions
#
#cat <(bedtools intersect -a $DATA/regioncoordsdnds.bed -b $DATA/allintfilter.bed -wa -wb -sorted | awk '{if ($7==$23) print}' | cut -f -16,28,29,30) <(bedtools intersect -a $DATA/regioncoordsdnds.bed -b $DATA/allintfilter.bed -v -sorted | awk '{print $0"\t.\t.\t."}') \
#| bedtools groupby -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 -c 17,18,19 -o collapse,collapse,collapse | sort -k1,1 -k6,6 -k5,5 -k7,7 \
#| bedtools groupby -i - -g 1,6,5,7 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,3 -o distinct,min,max,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,collapse,collapse,collapse,collapse,collapse | cut -f 5-26 \
#| python recompute.py | cat <(printf "#chr\tstart\tend\ttranscript\tdomain\tgene\tautoregs(uniqid for non-domain regions)\tcov_ratio\tlength\tdn\tds\tna\tdn/ds\tdensity\tfvrv\tprevalence\tmafs\timpacts\ttype\tstarts\tends\tmaf_modifier\n") - | bgzip > $DATA/regionsmafsdnds.bed.gz 
#
#python mafcalc.py -f $DATA/regionsmafsdnds.bed.gz -m $MAF > $DATA/regions.$MAF.bed # only greater than at the moment...
