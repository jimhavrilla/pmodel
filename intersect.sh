#intersect.sh - runs intersections

# do intersections for variants, merge lengths for autoreg+gene combos (lencount.py), get MAF, variant type (FG), gene, domain, chr, start, end, impact, other info for analysis

# filter if AN_Adj >= 97129.6:
cat <(zgrep "^#" $DATA/VEPEXAC3.vcf.gz) <(zgrep -v "^#" $DATA/VEPEXAC3.vcf.gz \
	| awk -F ';' '{if ($1 ~ /X|Y/) {t=$17} else t=$16} {sub(/\w+=/,"",t); if (t>=0.8*60706*2) print}' | grep -w PASS \
	| sort -k1,1 -k2,2n) \
	| bgzip -c $DATA/VEPEXAC3filter.vcf.gz

#intersect

bedtools intersect -a <(cut -f 1,2,3,11,25,27,33,43,45,47,48 $DATA/doms.bed | python lencount.py | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -sorted -wb | cut -f 1,2,3,4,5,6,7,8,9,13,14,17 | python var.py -d > $DATA/uniqintfilter.bed
bedtools intersect -a <(cut -f 1,2,3,12,24,26,27 $DATA/nodoms.bed | awk '{t=$7;$7=$6;$6=t;print}' | tr -s " " "\t" | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -sorted -wb | cut -f 1,2,3,4,5,6,7,11,12,15 | python var.py -n > $DATA/nodomintfilter.bed

# assign types to impacts

cat $DATA/uniqintfilter.bed $DATA/nodomintfilter.bed | cat <(printf "#chr,start,end,ref,alt,pfamA_id,autoreg,uniqid,length_of_region,covratio,gene_symbol,maf,impact,type,codons,amino_acids,gene_id_csq,gene_symbol_csq,transcript_id_csq,exon_number_csq,polyphen,sift,protein_position,biotype\n") <(awk 'NR==FNR{a[$2]}$19 in a{print $0}' $DATA/transcripts.txt - | sort -k11,11 -k7,7) > $DATA/allintfilter.bed
