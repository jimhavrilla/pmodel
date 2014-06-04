#for initial p53 pass:
#(grep ^# $DATA/ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf; grep -v ^# -h $DATA/ESP*.vcf) > $DATA/foo.vcf; mv $DATA/foo.vcf $DATA/ESPALL.vcf
#python extract.py
#hgTables is pfam from UCSC GenomeBrowser
#cat $DATA/hgTables | awk '{t=$1;$1=$2;$2=t;t=$2;$2=$3;$3=t;t=$3;$3=$4;$4=t;print;}' | tr -s " " | tr -s " " "\t" | tail -n +2 | sort -k1,1 -k2,2n > $DATA/hgTables.bed
#bedtools intersect -a $DATA/hgTables.bed -b $DATA/gene.bed -sorted | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,5 -o distinct > $DATA/p53domain.bed
#for bill's gtfs:
export DATA=~/work/data/pmodeldata
for chrom in {1..22}
do
python rearrange.py $DATA/Homo_sapiens.chr$chrom.pfam.gtf | sort -k1,1 -k2,2n > $DATA/chr$chrom.bed
bedtools intersect -a $DATA/chr$chrom.bed -b $DATA/ESP6500SI-V2-SSA137.updatedRsIds.chr$chrom.snps_indels.vcf -sorted > $DATA/intchr$chrom.bed
bedtools intersect -a <(cat $DATA/chr$chrom.bed | tr -s " " "\t" | cut -f 1,2,3,5,6,7,8,9,10,11 ) -b $DATA/ESP6500SI-V2-SSA137.updatedRsIds.chr$chrom.snps_indels.vcf -sorted > $DATA/chrcount$chrom.bed
done
cat $DATA/chr*.bed | sort -k1,1 -k2,2n > $DATA/all.bed
cat $DATA/intchr*.bed | sort -k1,1 -k2,2n > $DATA/allint.bed
cat $DATA/chrcount*.bed | sort -k10,10 | bedtools groupby -g 10 -c 10 -o count | sort -k1,1  > $DATA/allcount.bed
rm $DATA/chr* $DATA/intchr* rm $DATA/chrcount*
sort -k2,2nr $DATA/allcount.bed | head > top10.txt
sort -k2,2nr $DATA/allcount.bed | tail > bottom10.txt
cat $DATA/all.bed | tr -s " " "\t" | cut -f 4,11 | sort -k2,2 | less
cat $DATA/all.bed | tr -s " " "\t" | cut -f 4,11 | awk '{arr[$2]+=$1} END {for (i in arr) {print i,arr[i]}}' | sort -k1,1 > $DATA/sumlist.bed
python divide.py $DATA/allcount.bed $DATA/sumlist.bed > $DATA/allratio.bed
sort -k2,2n $DATA/allratio.bed | head > top10ratio.txt
sort -k2,2n $DATA/allratio.bed | tail > bottom10ratio.txt