#for initial p53 pass:
#python extract.py
#hgTables is pfam from UCSC GenomeBrowser
#cat $DATA/hgTables | awk '{t=$1;$1=$2;$2=t;t=$2;$2=$3;$3=t;t=$3;$3=$4;$4=t;print;}' | tr -s " " | tr -s " " "\t" | tail -n +2 | sort -k1,1 -k2,2n > $DATA/hgTables.bed
#bedtools intersect -a $DATA/hgTables.bed -b $DATA/gene.bed -sorted | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,5 -o distinct > $DATA/p53domain.bed
#bill's domain count file:
mysql -N --raw -h wrpxdb.its.virginia.edu -u web_user -pfasta_www pfam27 < count_human_pfam.sql > human_pfam.counts
#for bill's gtfs:
grep ^# $DATA/ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf; grep -v ^# -h $DATA/ESP*.vcf) > $DATA/foo.vcf; mv $DATA/foo.vcf $DATA/ESPALL.vcf
perl variant_effect_predictor.pl -i $DATA/ESPALL.vcf -o $DATA/VEPESPALL.vcf --offline --format vcf --vcf --force_overwrite;
export DATA=~/work/data/pmodeldata
gzcat $DATA/Homo_sapiens.GRCh37.75.gtf.gz | grep 'protein_coding\texon' | sed 's/protein_coding\texon\t//g' > $DATA/GRCh37.bed
#grep -v -E "C|G" $DATA/mart_export.bed | sort -k1,1 -k2,2n | sed '395352d' > foo.bed; mv foo.bed $DATA/mart_export.bed; # human genes from Ensembl
for chrom in {1..22} X Y ; do wget -P $DATA http://fastademo.bioch.virginia.edu/pfam_dna/Homo_sapiens.chr$chrom.pfam.gtf; done
for chrom in {1..22} X Y
do
sort -k2,2n $DATA/Homo_sapiens.chr$chrom.pfam.gtf | python rearrange.py | sort -k1,1 -k2,2n  > $DATA/chr$chrom.bed
done
cat $DATA/chr*.bed | sort -k1,1 -k2,2n > $DATA/all.bed
rm $DATA/chr*
# remove utrs and introns
bedtools intersect -a $DATA/all.bed -b $DATA/GRCh37.bed -wb | awk {'print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$17'} FS='\t' OFS='\t' | perl uniq.pl | tr -s " " "\t" | sort -k45,45 -k1,1 -k2,2n > $DATA/foo.bed
python rearrange2.py $DATA/foo.bed $DATA/alluniq.bed
sort -k11,11 -k1,1 -k2,2n $DATA/foo.bed | python rearrange3.py $DATA/alldom.bed
1	25817901	25818075	.	+	.	gene^C
bedtools intersect -a <(perl -pe 's/tag\s*\S*?(?=\n|\s)//g' $DATA/GRCh37.bed | perl -pe 's/ccds_id\s*\S*?(?=\n|\s)//g' | tr -s " " "\t" | sort -k1,1 -k2,2n) -b <(sort -k1,1 -k2,2n $DATA/alluniq.bed) -wa -wb | awk '{ if ($8==$51 && $24==$67) print $0}' | tr -s " " "\t" | cut -f -35 | python nodom.py $DATA/nodom.bed
# get a count for variants per domain
bedtools intersect -a <(sort -k1,1 -k2,2n -k3,3n $DATA/alldom.bed | tr -s " " "\t" | cut -f 1,2,3,8,11,13,45 ) -b $DATA/VEPESPALL.vcf -sorted -wb | cut -f 1,2,3,4,5,6,7,15 | tee $DATA/foo.bed | sort -k6,6 | bedtools groupby -g 6 -c 6 -o count | sort -k1,1 > $DATA/allcount.bed
bedtools intersect -a <(sort -k1,1 -k2,2n -k3,3n $DATA/allnodom.bed | tr -s " " "\t" | cut -f 1,2,3,9,17 ) -b $DATA/VEPESPALL.vcf -sorted -wb | cut -f 1,2,3,4,5,13 | awk -F";" '{print $1,$2,$7,$20}' | sed 's/FG=.*://g' | sed 's/MAF=.*,//g' > $DATA/nodomint.bed
# if you want to check type count: awk -F";" '{print $1,$2,$7,$20}' foo.bed | sed 's/FG=.*://g' | sed 's/MAF=.*,//g' | tr -s " " "\t" | cut -f 8 | sort -k1,1 | uniq -c
# get MAF and convert from percent to fraction, variant type (FG), gene, domain, chr, start, end for analysis
awk -F";" '{print $1,$2,$3,$8,$21}' $DATA/foo.bed | sed 's/FG=.*://g' | sed 's/MAF=.*,//g' | awk '{print $1,$2,$3,$4,$5,$6,$7,$8/100,$9}' | tr -s " " "\t" | sort -k6,6 > $DATA/allint.bed; rm $DATA/foo.bed
awk -F";" '{print $1,$2,$3,$8,$21,$23}' $DATA/foo.bed | perl -pe 's/[XN]M_[A-Z0-9.]*?:intron//g' | perl -pe 's/((\S*\s*){7})MAF=.*,(.*)\s*FG=(.*:).*\n/$1$3\n/g'
#sort counts from bill by domain
sort -k2,2 $DATA/human_pfam.counts > $DATA/foo.bed; mv $DATA/foo.bed $DATA/human_pfam.counts
#get total bp for each domain
cat $DATA/all.bed | tr -s " " "\t" | cut -f 4,11 | awk '{arr[$2]+=$1} END {for (i in arr) {print i,arr[i]}}' | sort -k1,1 > $DATA/sumlist.bed
# make table of counts, non-syn, syn, total var, total bp/exome per domain
python maketable.py $DATA/allint.bed > $DATA/foo.bed
python mergetable.py $DATA/foo.bed $DATA/human_pfam.counts $DATA/sumlist.bed > $DATA/dtable.txt; rm $DATA/foo.bed

R commands:

library(ggplot2)
library(plotrix)
library(RColorBrewer)
dtable <- read.delim("~/work/data/pmodeldata/dtable.txt")
clans <- read.delim("~/work/data/pmodeldata/count_human_pfam_clan.tab")
m<-merge(dtable,clans,by.x="domain",by.y="pfamA_id")
write.table(m,paste(DATA,"/clans.txt",sep=""),sep="\t")
echo -e "clan_acc\tclan_id\tnonsynct\tsynct\ttotalvarct\tdomcount\ttotalbp\tmmaf\tdn.ds\tvar.bp.ratio" > $DATA/ctable.txt
sed '1d' $DATA/clans.txt | sort -k13,13 | bedtools groupby -g 13,14 -c 4,5,6,7,8,9,10,11 -o sum,sum,sum,sum,sum,sum,sum,sum | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"; if ($4!=0) printf $3/$4; else printf $3/(1+$4); printf "\t"$5/$7"\n"}' >> $DATA/ctable.txt
ctable <- read.delim("~/work/data/pmodeldata/ctable.txt")
shuffle <-function(list){
	l<-c()
	i=1
	while (i < length(list)/2){
		l[length(l)+1]<-list[i]
		l[length(l)+1]<-list[length(list)-i+1]
		i=i+1
	}
	return(l)
}
pointmod <- function(i=dtable,ll=0,ul=14000,w=sturges(t(dtable["mmaf"])),nbin=6,pchcolumn="totalvarct",colorcolumn="mmaf",pchbool=TRUE,colorbool=TRUE){
	if (pchbool==FALSE){
		pch=NULL
	} else{
		pch<-c()
		plist<-seq(0,nbin-1)
		pchleg<-cut_number(t(i[pchcolumn]),nbin)
		for (k in 1:length(t(i[pchcolumn]))){
			for (j in 1:length(levels(pchleg))){
				if (toString(pchleg[k])==pchleg[j]){
					pch[length(pch)+1] <- plist[j]
					j_=j
					break
				}
			}
		}
	}
	if (colorbool==FALSE){
		col<-1
	} else{
		col <- c()
		pal <- shuffle(colorRampPalette(brewer.pal(11,"Spectral"))(length(w)+1))
		for (k in 1:length(t(i[colorcolumn]))){
			for (j in 1:length(w)){
				if (t(i[colorcolumn])[k]<=w[j]){
					col[length(col)+1] <- pal[j]
					j_=j
					break
				}
			}
			if (t(i[colorcolumn])[k]>w[j_]){
				col[length(col)+1]<-pal[j_+1]
			}
		}
	}
	pchleg<-levels(pchleg)
	l<-list(pch,col,pchleg)
	return(l)
}
sturges <- function(m=t(dtable["mmaf"]),column="mmaf"){
	w <- (range(m)[2]-range(m)[1])/(1+3.322*log(length(m),base=10))
	tot <- min(m)
	widthlist <- c(tot)
	while (tot < max(m)){
		widthlist[length(widthlist)+1] <- tot+w
		tot <- tot+w
	}
	if (tot => max(m)){
		return(widthlist)
	}
}
subplot <- function(i=dtable,yll=0,yul=3,pchcolumn="totalvarct",colorcolumn="mmaf",nbin=6,x="dn.ds",y="domcount",xlab="dn/ds",ylab="domain occurrences in exome",ylimiter="domcount",label="domain",labelbool=TRUE,colorbool=TRUE,pchbool=TRUE,limitbool=TRUE){
	if (limitbool==TRUE){
		i=i[yll<i[ylimiter] & i[ylimiter]<=yul,]
	}
	if (colorbool==TRUE){
		w<-round(sturges(m=t(i[colorcolumn]),column=colorcolumn),5)
		c<-shuffle(colorRampPalette(brewer.pal(11,"Spectral"))(length(w)+1))
	} else{
		w=NULL
	}
	l<-pointmod(i=i,pchcolumn=pchcolumn,colorcolumn=colorcolumn,nbin=nbin,w=w,pchbool=pchbool,colorbool=colorbool)
	pch<-l[1][[1]]
	col<-l[2][[1]]
	pchleg<-l[3][[1]]
	# par(fig=c(0.1,0.9,0,0.9),new=TRUE)
	plot(t(i[x]),t(i[y]),xlab=xlab,ylab=ylab,pch=pch,col=col)
	if (labelbool == TRUE){
		text(t(i[x]),t(i[y]),labels=t(i[label]),adj=c(.5,-.4),cex=0.5,offset=0.5)
	}
	# par(fig=c(0,.1,.3,.5),xpd=NA)
	color.legend((max(i[x])-min(i[x]))*.985+min(i[x]),(max(i[y])-min(i[y]))*.75+min(i[y]),max(i[x]),max(i[y]),legend=w,rect.col=c,gradient="y")
	legend(par()$usr[1],mean(par()$usr[3:4]),border=rep("white",nbin),legend=pchleg,pch=seq(0,nbin-1),bty="n",xpd=TRUE,xjust=0,yjust=0.5,cex=0.65) #inset=c(0,0) #min(i[x]),(max(i[y])-min(i[y]))*.75+min(i[y]),(max(i[x])-min(i[x]))*.05+min(i[x]),max(i[y])
}
plot(table$var.bp.ratio,log(table$totalbp,base=10),xlab="variants/total length of domain",ylab="total length of domain (log)",cex=cex,col=col)
color.legend(.10,3,.14,5,legend=w,rect.col=c,gradient="y")
