#for initial p53 pass:
#(grep ^# ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf; grep -v ^# -h ESP*.vcf) > foo.vcf; mv foo.vcf ESPALL.vcf
#python extract.py
#hgTables is pfam from UCSC GenomeBrowser
#cat hgTables | awk '{t=$1;$1=$2;$2=t;t=$2;$2=$3;$3=t;t=$3;$3=$4;$4=t;print;}' | tr -s " " | tr -s " " "\t" | tail -n +2 | sort -k1,1 -k2,2n > hgTables.bed
#bedtools intersect -a hgTables.bed -b gene.bed -sorted | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,5 -o distinct > p53domain.bed
#for bill's gtfs:
export DATA=~/work/data/pmodeldata
for chrom in {1..22}
do
cat $DATA/Homo_sapiens.chr$chrom.pfam.gtf | awk -F$'\t' '{OFS="\t"} {t=$2;$2=$4;$4=t;t=$3;$3=$5;$5=t;print;}' | awk -F$'\t' '{OFS="\t"} {if (($7 == "-") && ($2 > $3)) {t=$2; $2=$3; $3=t;} print $0}' | sort -k1,1 -k2,2n | tr -s " " "\t" > chr$chrom.bed
bedtools intersect -a <(cut -f 1-24 chr$chrom.bed) -b $DATA/ESP6500SI-V2-SSA137.updatedRsIds.chr$chrom.snps_indels.vcf -sorted > intchr$chrom.bed
done
(cat chr1.bed; cat chr*.bed) > all.bed
(cat intchr1.bed; cat intchr*.bed) > allint.bed
rm chr* intchr*
sort allint.bed -k22,22 | bedtools groupby -g 22 -c 22 -o count | sort -k2,2nr | head > top10.txt
sort allint.bed -k22,22 | bedtools groupby -g 22 -c 22 -o count | sort -k2,2nr | tail > bottom10.txt
sum=$(sort allint.bed -k22,22 | bedtools groupby -g 22 -c 22 -o count | awk '{ sum += $2 } END {print sum}')
cat top10.txt | awk '{ratio=$2/sum} {print $1"\t"ratio}' sum="$sum" > top10ratio.txt
cat bottom10.txt | awk '{ratio=$2/sum} {print $1"\t"ratio}' sum="$sum" > bottom10ratio.txt