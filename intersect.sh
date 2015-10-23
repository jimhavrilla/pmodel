#intersect.sh - runs intersections

# do intersections for variants, merge lengths for autoreg+gene combos (lencount.py), get MAF, variant type (FG), gene, domain, chr, start, end, impact, other info for analysis

# filter if AN_Adj >= 97129.6:

if [[ ! -f $DATA/VEPEXAC3filter.vcf.gz ]]; then
cat <(zgrep "^#" $DATA/VEPEXAC3.vcf.gz) <(zgrep -v "^#" $DATA/VEPEXAC3.vcf.gz \
	| awk -F ';' '{if ($1 ~ /X|Y/) {t=$17} else t=$16} {sub(/\w+=/,"",t); if (t>=0.8*60706*2) print}' \
	| grep -w PASS \
	| sort -k1,1 -k2,2n) \
	| bgzip -c > $DATA/VEPEXAC3filter.vcf.gz
fi

if [[ ! -f  $DATA/doms.fix.sort.bed ]]
then
    cat $DATA/doms.bed  \
    | python fix_doms.py  \
    | bedtools sort -i stdin -header \
    > $DATA/doms.fix.sort.bed
fi

if [[ ! -f  $DATA/nodoms.fix.sort.bed ]]
then
    cat $DATA/nodoms.bed \
    | python fix_nodoms.py  \
    | bedtools sort -i stdin -header \
    > $DATA/nodoms.fix.sort.bed
fi

#intersect

if [[ ! -f  $DATA/doms.fix.sort.bed.I.VEPEXAC3filter.vcf.gz.bed ]]
then
    (paste \
        <(head -n1 $DATA/doms.fix.sort.bed) \
        <(bcftools view -h $DATA/VEPEXAC3filter.vcf.gz | tail -n1 ); \
    bedtools intersect \
        -sorted -wa -wb -header \
        -a $DATA/doms.fix.sort.bed \
        -b $DATA/VEPEXAC3filter.vcf.gz \
        | bedtools sort -i stdin) \
    > $DATA/doms.fix.sort.bed.I.VEPEXAC3filter.vcf.gz.bed
fi

if [[ ! -f  $DATA/nodoms.fix.sort.bed.I.VEPEXAC3filter.vcf.gz.bed ]]
then
    (paste \
        <(head -n1 $DATA/nodoms.fix.sort.bed) \
        <(bcftools view -h $DATA/VEPEXAC3filter.vcf.gz | tail -n1 ); \
    bedtools intersect \
        -sorted -wa -wb -header \
        -a $DATA/nodoms.fix.sort.bed \
        -b $DATA/VEPEXAC3filter.vcf.gz \
    | bedtools sort -i stdin) \
    > $DATA/nodoms.fix.sort.bed.I.VEPEXAC3filter.vcf.gz.bed
fi

if [[ ! -f $DATA/doms.fix.sort.bed.NI.VEPEXAC3filter.vcf.gz.bed ]] 
then
    bedtools intersect \
        -v -sorted -header\
        -a $DATA/doms.fix.sort.bed \
        -b $DATA/doms.fix.sort.bed.I.VEPEXAC3filter.vcf.gz.bed \
        -sorted \
        > $DATA/doms.fix.sort.bed.NI.VEPEXAC3filter.vcf.gz.bed
fi

if [[ ! -f $DATA/nodoms.fix.sort.bed.NI.VEPEXAC3filter.vcf.gz.bed ]] 
then
    bedtools intersect \
        -v -sorted -header\
        -a $DATA/nodoms.fix.sort.bed \
        -b $DATA/nodoms.fix.sort.bed.I.VEPEXAC3filter.vcf.gz.bed \
        > $DATA/nodoms.fix.sort.bed.NI.VEPEXAC3filter.vcf.gz.bed
fi

if [[ ! -f  $DATA/nodomintfilter.bed ]]
then
    cat $DATA/nodoms.fix.sort.bed.I.VEPEXAC3filter.vcf.gz.bed \
        | python filter_doms.py \
        > $DATA/nodomintfilter.bed
fi

if [[ ! -f $DATA/domintfilter.bed ]]
then
    cat $DATA/doms.fix.sort.bed.I.VEPEXAC3filter.vcf.gz.bed \
        | python filter_doms.py \
        > $DATA/domintfilter.bed
fi

if [[ ! -f $DATA/allintfilter.bed ]]
then
    (cat $DATA/domintfilter.bed; tail -n+2 $DATA/nodomintfilter.bed) \
        | bedtools sort -i stdin -header  \
        > $DATA/allintfilter.bed
fi
