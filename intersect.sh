#intersect.sh - runs intersections

# do intersections for variants, merge lengths for autoreg+gene combos (lencount.py), get MAF, variant type (FG), gene, domain, chr, start, end, impact, other info for analysis

# filter if AN_Adj >= 97129.6:
cat <(zgrep "^#" $DATA/VEPEXAC3.vcf.gz) <(zgrep -v "^#" $DATA/VEPEXAC3.vcf.gz \
	| awk -F ';' '{if ($1 ~ /X|Y/) {t=$17} else t=$16} {sub(/\w+=/,"",t); if (t>=0.8*60706*2) print}' | grep -w PASS \
	| sort -k1,1 -k2,2n) \
	| bgzip -c $DATA/VEPEXAC3filter.vcf.gz

#intersect

bedtools intersect -a <(cut -f 1,2,3,11,15,25,27,33,43,45,47,48 $DATA/doms.bed \
	| sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -sorted -wb \
    | cut -f 1,2,3,4,5,6,8,10,11,12,16,17,20 \
	| python var.py -d | sort -k1,1 -k2,2n > $DATA/domintfilter.bed

bedtools intersect -a <(cut -f 1,2,3,8,12,24,26,27 $DATA/nodoms.bed \
    | tr -s " " "\t" | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -sorted -wb \
    | cut -f 1,2,3,4,5,6,7,8,12,13,16 \
    | python var.py -n | sort -k1,1 -k2,2n > $DATA/nodomintfilter.bed

# assign types to impacts

cat $DATA/domintfilter.bed $DATA/nodomintfilter.bed | cat <(printf "#chr,start,end,ref,alt,pfamA_id,autoreg,uniqid,length_of_region,covratio,gene_symbol,maf,impact,type,codons,amino_acids,gene_id_csq,gene_symbol_csq,transcript_id_csq,exon_number_csq,polyphen,sift,protein_position,biotype\n") <(awk 'NR==FNR{a[$2]}$19 in a{print $0}' $DATA/transcripts.txt - | sort -k11,11 -k7,7) > $DATA/allintfilter.bed
