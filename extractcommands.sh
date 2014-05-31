#for initial p53 pass:
(grep ^# ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf; grep -v ^# -h ESP*.vcf) > foo.vcf; mv foo.vcf ESPALL.vcf
python extract.py
#hgTables is pfam from UCSC GenomeBrowser
cat hgTables | awk '{t=$1;$1=$2;$2=t;t=$2;$2=$3;$3=t;t=$3;$3=$4;$4=t;print;}' | tr -s " " | tr -s " " "\t" | tail -n +2 | sort -k1,1 -k2,2n > hgTables.bed
bedtools intersect -a hgTables.bed -b gene.bed -sorted | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,5 -o distinct > p53domain.bed
#for bill's gtf:
cat Homo_sapiens.chr17.pfam.gtf | awk -F$'\t' '{OFS="\t"} {t=$2;$2=$4;$4=t;t=$3;$3=$5;$5=t;print;}' | awk -F$'\t' '{OFS="\t"} {if (($7 == "-") && ($2 > $3)) {t=$2; $2=$3; $3=t;} print "chr"$0}' | sort -k1,1 -k2,2n | tr -s " " "\t" > chr17.bed