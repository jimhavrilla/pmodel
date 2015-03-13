DATA=$(perl -pe 's/(.*\/)\w*\.\w*/$1/g' <<< $1)
#$1=database,$2=omimfile,$3=greater/less/range,$4=leftbound,$5=rightbound
# pick one impact per variant and an maf limit for variants

bash query.sh $1 '%' $3 $4 $5 | python formatvar.py | sort -k6,6 -k8,8 > $DATA/allint2maf.$3.$4-$5.bed

# make table by gene/domain pairing

python dg.py -f $DATA/allint2maf.$3.$4-$5.bed $DATA/human_pfam.counts $DATA/sumlist.bed $DATA/dngpair.txt
grep -v -w "\." $DATA/dngpair.txt > $DATA/dgpair.txt

## to compare with RVIS and other scores for divergence

awk '{if ($6>0) print $0}' $DATA/dgpair.txt | python diverge.py -p | awk '{if ($3>3) print $0}' | awk '{gene[$1]=$12; row[$1]=$0} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]+=1}} END {for (i in rank) print row[i],rank[i]/NR}' | tr -s " " "\t" | sort -k11,11nr  > $DATA/diverge.pair.txt

#turn XLS supplement from RVIS paper into tab-delimited $DATA/rvis.txt with RVIS scores/percentiles
#then use R to pull out relevant rows and columns

Rscript divtable.r $DATA $3.$4-$5

#in unix use awk NR==FNR to get rows with RVIS percentiles for each gene and omim phenotype counts

OMIM=$2 

awk 'NR==FNR{a[$1]=$3;next}$1 in a{print $0,a[$1]}' OFS='\t' $DATA/rvis.txt $DATA/dps.pair.$3.$4-$5.txt | sort -k8,8nr > RVISint.$3.$4-$5.txt

awk 'NR==FNR{a[$1]=$3"\t"$4"\t"$5"\t"$6"\t"$7;next}$1 in a{print $0,a[$1]}' FS='\t' OFS='\t' $OMIM RVISint.$3.$4-$5.txt | sort -k10,10r | cat <(printf "gene, domains, num_of_domains, dnds_list, z-scores, z-score_range, divergence_metric, div_rank, rvis_rank, omim_count, de_novo, dominant_negative, haploinsufficient, recessive\n") - > ordivtable.$3.$4-$5.txt