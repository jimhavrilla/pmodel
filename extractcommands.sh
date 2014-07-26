#for initial p53 pass:
#python extract.py
#hgTables is pfam from UCSC GenomeBrowser
#cat $DATA/hgTables | awk '{t=$1;$1=$2;$2=t;t=$2;$2=$3;$3=t;t=$3;$3=$4;$4=t;print;}' | tr -s " " | tr -s " " "\t" | tail -n +2 | sort -k1,1 -k2,2n > $DATA/hgTables.bed
#bedtools intersect -a $DATA/hgTables.bed -b $DATA/gene.bed -sorted | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,5 -o distinct > $DATA/p53domain.bed
#bill's domain count file:
mysql -N --raw -h wrpxdb.its.virginia.edu -u web_user -pfasta_www pfam27 < count_human_pfam.sql > human_pfam.counts
#for bill's gtfs:
grep ^# $DATA/ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf; grep -v ^# -h $DATA/ESP*.vcf) > $DATA/foo.vcf; mv $DATA/foo.vcf $DATA/ESPALL.vcf
export DATA=~/work/data/pmodeldata
gzcat $DATA/Homo_sapiens.GRCh37.75.gtf.gz | grep 'protein_coding      exon' > $DATA/GRCh37.gtf
#grep -v -E "C|G" $DATA/mart_export.bed | sort -k1,1 -k2,2n | sed '395352d' > foo.txt; mv foo.txt $DATA/mart_export.bed; # human genes from Ensembl
for chrom in {1..22} X Y ; do wget -P $DATA http://fastademo.bioch.virginia.edu/pfam_dna/Homo_sapiens.chr$chrom.pfam.gtf; done
for chrom in {1..22} X Y
do
sort -k2,2n $DATA/Homo_sapiens.chr$chrom.pfam.gtf | python rearrange.py | sort -k1,1 -k2,2n  > $DATA/chr$chrom.bed
done
cat $DATA/chr*.bed | sort -k1,1 -k2,2n > $DATA/all.bed
rm $DATA/chr*
# remove utrs and introns
bedtools intersect -a $DATA/all.bed -b $DATA/GRCh37.gtf -wb | awk {'print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$19'} FS='\t' OFS='\t' | sort -k1,1 -k2,2n | python rearrange2.py > $DATA/foo.bed; mv $DATA/foo.bed $DATA/all.bed
# get a count for variants per domain
bedtools intersect -a <(sort -k1,1 -k2,2n -k3,3n $DATA/all.bed | tr -s " " "\t" | cut -f 1,2,3,8,11,13 ) -b $DATA/ESPALL.vcf -sorted -wb | cut -f 1,2,3,4,5,6,14,15 | tee $DATA/foo.txt | sort -k5,5 | bedtools groupby -g 5 -c 5 -o count | sort -k1,1 > $DATA/allcount.bed
# if you want to check type count: awk -F";" '{print $1,$2,$7,$20}' foo.txt | sed 's/FG=.*://g' | sed 's/MAF=.*,//g' | tr -s " " "\t" | cut -f 8 | sort -k1,1 | uniq -c
# get MAF and convert from percent to fraction, variant type (FG), gene, domain, chr, start, end for analysis
awk -F";" '{print $1,$2,$7,$20}' $DATA/foo.txt | sed 's/FG=.*://g' | sed 's/MAF=.*,//g' | awk '{print $1,$2,$3,$4,$5,$6,$7/100,$8}' | tr -s " " "\t" | sort -k5,5 > $DATA/allint.bed; rm foo.txt
#sort counts from bill by domain
sort -k2,2 $DATA/human_pfam.counts > $DATA/foo.txt; mv $DATA/foo.txt $DATA/human_pfam.counts
#get total bp for each domain
cat $DATA/all.bed | tr -s " " "\t" | cut -f 4,11 | awk '{arr[$2]+=$1} END {for (i in arr) {print i,arr[i]}}' | sort -k1,1 > $DATA/sumlist.bed
# make table of counts, non-syn, syn, total var, total bp/exome per domain
python maketable.py $DATA/allint.bed > $DATA/foo.txt
python mergetable.py $DATA/foo.txt $DATA/human_pfam.counts $DATA/sumlist.bed > $DATA/table.txt; rm $DATA/foo.txt

R commands:

library(calibrate)
table <- read.delim("~/work/data/pmodeldata/table.txt")
pointmod <- function(table,l=c(2,5,10,15,30,50,80,120,300,500,1000,2000,4000,8000),w=sturges()){
	cex <- 1 - (table$domcount<l[1])*.2 - (table$domcount<l[2])*.35 - (table$domcount<l[3])*.5 + (table$domcount>l[4])*.1 + (table$domcount>l[5])*.05 + (table$domcount>l[6])*.07 + (table$domcount>l[7])*.08 + (table$domcount>l[8])*.1 + (table$domcount>l[9])*.2 + (table$domcount>l[10])*.25 + (table$domcount>l[11])*.3 + (table$domcount>l[12])*.4 + (table$domcount>l[13])*.5 + (table$domcount>l[14])*.7
	m <- seq(56,56+length(w)-1)
	
	col <- 56 - (table$mmaf<m[1])*1 - (table$mmaf<m[2])*1 - (table$mmaf<m[3])*1 + (table$mmaf<m[4])*1 + (table$mmaf<m[5])*1 + (table$mmaf<m[6])*1 + (table$mmaf<m[7])*1 + (table$mmaf>m[8])*1 + (table$mmaf>m[9])*1 + (table$mmaf>m[10])*1 + (table$mmaf>m[11])*1 + (table$mmaf>m[12])*1 + (table$mmaf>m[13])*1 + (table$mmaf>m[14])*1
	t<-data.frame(cex,col)
	return(t)
}
sturges <- function(m=table$mmaf){
	w <- (range(m)[2]-range(m)[1])/(1+3.322*log(length(m),base=10))
	tot <- min(m)
	widthlist <- c(tot)
	while (tot < max(m)){
		widthlist[length(widthlist)+1] <- tot+w
		tot <- tot+w
	}
	if (tot > max(m)){
		return(widthlist)
	}
}
plot(table$var.bp.ratio,log(table$totalbp,base=10),xlab="variants/total length of domain",ylab="total length of domain (log)") #xlim=c(0,0.02),ylim=c(10,14)
textxy(table$var.bp.ratio,log(table$totalbp,base=10),labs=table$domain)
plot(table$dn.ds,table$domcount,xlab="dn/ds",ylab="domain occurrences in exome")
textxy(table$dn.ds,table$domcount,labs=table$domain)
