(grep ^# ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf; grep -v ^# -h ESP*.vcf) > foo.vcf; mv foo.vcf ESPALL.vcf
python extract.py
#hgTables is pfam from UCSC GenomeBrowser
cat hgTables | awk '{t=$1;$1=$2;$2=t;t=$2;$2=$3;$3=t;t=$3;$3=$4;$4=t;print;}' > hgTables.bed
cat hgTables.bed | tr -s " " | tr -s " " "\t" > hg.bed
cp hg.bed hgTables.bed
tail -n +2 hgTables.bed > hg.bed
cp hg.bed hgTables.bed
sort -k1,1 -k2,2n hgTables.bed > hg.bed
mv hg.bed hgTables.bed
bedtools intersect -a hgTables.bed -b gene.bed -sorted > pdomaintp53.bed
bedtools merge -i <(sort -k1,1 -k2,2n pdomaintp53.bed) -c 4,5 -o distinct | less