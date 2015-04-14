DATA=$2; OMIM=$3; MOD=$4; MAF1=$5; MAF2=$6;
#$2=database - pmodeldata,$3=omimfile - $DATA/omim.txt,$4=greater/less/range,$5=leftbound,$6=rightbound

case "$1" in
	-u|--uniqid)
		#coverage,var density,nc,sct per uniqid

		bash query.sh $DATA '%' $MOD $MAF1 $MAF2 | sort -t '|' -k7,7 -k1,1 -k2,2 -k3,3 | python formatvar.py > $DATA/allint2uniq.$MOD.$MAF1-$MAF2.bed

		DATA=$(perl -pe 's/(.*\/)\w*\.\w*/$1/g' <<< $DATA)

		#calculates median maf, counts
		gawk 'function median(v) {c=asort(v,j); if (c % 2) return j[(c+1)/2]; else return (j[c/2+1]+j[c/2])/2} {{if ($14=="ds") {sct[$7]++; ct[$7]++} else if ($14=="dn") nct[$7]++; ct[$7]++} len[$7]=$9; maf[$7][$1 $2 $3 $14]=$12; row[$7]=$6 " " $11 " " $7 " " $8 " " $9 " " $10} END {for (i in ct) print row[i],(nct[i]==0 ? nct[i]=0: nct[i]),(sct[i]==0 ? sct[i]=0: sct[i]),(nct[i]+sct[i]),nct[i]/(sct[i]==0 ? sct[i]+1: sct[i]),(nct[i]+sct[i])/len[i],median(maf[i])}' $DATA/allint2uniq.$MOD.$MAF1-$MAF2.bed > $DATA/uniqtable.$MOD.$MAF1-$MAF2.txt

		awk 'NR==FNR{a[$2]=$3} NR!=FNR{if ($1 in a) print $0,a[$1]; else print $0,0}' $DATA/human_pfam.counts $DATA/uniqtable.$MOD.$MAF1-$MAF2.txt > $DATA/foo; mv $DATA/foo $DATA/uniqtable.$MOD.$MAF1-$MAF2.txt
		
		grep -v -w "\." $DATA/uniqtable.$MOD.$MAF1-$MAF2.txt | awk '{if ($6>=.9) print $0}' | sort -k1,1 -k2,2 | python dg.py -u | sort -k2,2 -k3,3 | python diverge.py -u | awk '{if ($3>3) print $0}'| awk '{uniq[$1 $2]=$1; row[$1 $2]=$0; gene[$1]=$12} END {for (i in gene) ct++} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]++}} END {for (j in uniq) for (i in gene) {if (i==uniq[j]) print row[j],rank[i]/ct}}' | tr -s " " "\t" | sort -k16,16nr  > $DATA/diverge.pair.txt

		awk 'NR==FNR{a[$1]=$3;next}$1 in a{print $0,a[$1]}' OFS='\t' $DATA/rvis.txt $DATA/diverge.pair.txt > $DATA/foo

		awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0,"n"}' $OMIM/denovo.txt $DATA/foo \
		| awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0,"n"}' $OMIM/domneg.txt - \
		| awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0, "n"}' $OMIM/haplo.txt - \
		| awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0,"n"}' $OMIM/recessive.txt - | tr -s " " "\t" > $DATA/divtable.txt
		;;

		(grep ^# <(gzcat $DATA/clinvar_20150305.vcf.gz); grep 'CLNSIG=5' <(gzcat $DATA/clinvar_20150305.vcf.gz)) | bedtools intersect -a $DATA/alluniq.bed -b stdin | sort -k45,45 | bedtools groupby -i stdin -g 45 -c 45 -o count > $DATA/clinvaruniq.txt

		awk 'FNR==NR{uniq[$1]=$2} {u[$2]=$0} END {for (i in uniq) {if (i in u) print u[i],uniq[i]}}' $DATA/clinvaruniq.txt $DATA/divtable.txt | tr -s " " "\t" | sort -k20,20r -k16,16r > $DATA/truedivtable.txt

	-d|--domain)

		# pick one impact per variant and an maf limit for variants

		bash query.sh $DATA '%' $MOD $MAF1 $MAF2 | python formatvar.py | sort -k6,6 -k8,8 > $DATA/allint2dom.$MOD.$MAF1-$MAF2.bed

		DATA=$(perl -pe 's/(.*\/)\w*\.\w*/$1/g' <<< $DATA)

		# make table by gene/domain pairing

		python dg.py -f $DATA/allint2maf.$MOD.$MAF1-$MAF2.bed $DATA/human_pfam.counts $DATA/sumlist.bed $DATA/dngpair.txt
		grep -v -w "\." $DATA/dngpair.txt > $DATA/dgpair.txt

		## to compare with RVIS and other scores for divergence

		awk '{if ($6>0) print $0}' $DATA/dgpair.txt | python diverge.py -p | awk '{if ($3>3) print $0}' | awk '{gene[$1]=$12; row[$1]=$0} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]+=1}} END {for (i in rank) print row[i],rank[i]/NR}' | tr -s " " "\t" | sort -k11,11nr  > $DATA/diverge.pair.txt

		#turn XLS supplement from RVIS paper into tab-delimited $DATA/rvis.txt with RVIS scores/percentiles
		#then use R to pull out relevant rows and columns

		Rscript divtable.r $DATA $MOD.$MAF1-$MAF2

		#in unix use awk NR==FNR to get rows with RVIS percentiles for each gene and omim phenotype counts

		awk 'NR==FNR{a[$1]=$3;next}$1 in a{print $0,a[$1]}' OFS='\t' $DATA/rvis.txt $DATA/dps.pair.$MOD.$MAF1-$MAF2.txt | sort -k8,8nr > RVISint.$MOD.$MAF1-$MAF2.txt

		awk 'NR==FNR{a[$1]=$3"\t"$4"\t"$5"\t"$6"\t"$7;next}$1 in a{print $0,a[$1]}' FS='\t' OFS='\t' $OMIM RVISint.$MOD.$MAF1-$MAF2.txt | sort -k10,10r | cat <(printf "gene, domains, num_of_domains, dnds_list, z-scores, z-score_range, divergence_metric, div_rank, rvis_rank, omim_count, de_novo, dominant_negative, haploinsufficient, recessive\n") - > ordivtable.$MOD.$MAF1-$MAF2.txt
		;;
esac