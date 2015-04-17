#for initial p53 pass:

#python extract.py
#hgTables is pfam from UCSC GenomeBrowser
#cat $DATA/hgTables | awk '{t=$1;$1=$2;$2=t;t=$2;$2=$3;$3=t;t=$3;$3=$4;$4=t;print;}' | tr -s " " | tr -s " " "\t" | tail -n +2 | sort -k1,1 -k2,2n > $DATA/hgTables.bed
#bedtools intersect -a $DATA/hgTables.bed -b $DATA/gene.bed -sorted | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,5 -o distinct > $DATA/p53domain.bed

#bill's domain count file:

export DATA=~/work/data/pmodeldata
export SOFTWARE=~/software
mysql -N --raw -h wrpxdb.its.virginia.edu -u web_user -pfasta_www pfam27 < count_human_pfam.sql > human_pfam.counts

#for vcfs:

#grep ^# $DATA/ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf; grep -v ^# -h $DATA/ESP*.vcf) > $DATA/foo.vcf; mv $DATA/foo.vcf $DATA/ESPALL.vcf
#perl $SOFTWARE/ensembl-tools-release-76/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/ESPALL.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --total_length -o $DATA/VEPESPALL.vcf --vcf --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --force_overwrite --port 3337 --dir ~/.vep/ --compress "gunzip -c"
sudo perl $SOFTWARE/ensembl-tools-release-76/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/ExAC.r0.3.sites.vep.vcf.gz --cache --sift b --polyphen b --symbol --numbers --biotype --total_length -o $DATA/VEPEXAC3.vcf --vcf --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --force_overwrite --port 3337 --dir ~/.vep/ --compress "gunzip -c" --offline
gzcat $DATA/Homo_sapiens.GRCh37.75.gtf.gz | grep 'protein_coding\texon' | perl -pe 's/protein_coding\texon\t//g' | grep -e '^\d*[0-9X-Y]\t'| perl gtf2bed.pl | sort -k1,1 -k2,2n > $DATA/GRCh37.bed
#grep -v -E "C|G" $DATA/mart_export.bed | sort -k1,1 -k2,2n | sed '395352d' > foo.bed; mv foo.bed $DATA/mart_export.bed; # human genes from Ensembl

#for vcf coverage:
for chrom in {1..22} X Y ; do wget -P $DATA ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/coverage/Panel.chr$chrom.coverage.txt.gz; done
gzcat $DATA/Panel.chr*.coverage.txt.gz | awk '/^#/ {sub(/#.*/,"");getline;}1 {print $1,$2-1,$2,$3}' | tr -s " " "\t" > $DATA/coverage.bed
#bedtools intersect -a $DATA/VEPEXAC.vcf -b $DATA/coverage.bed -wa -wb | awk '{if ($12 >= 5) print $1,$2,$3,$4,$5,$6,$7,$8,$12}' OFS="\t" > $DATA/covVEPEXAC.vcf

#bills gtfs:

for chrom in {1..22} X Y ; do wget -P $DATA http://fastademo.bioch.virginia.edu/pfam_dna/Homo_sapiens.chr$chrom.pfam.gtf; done
for chrom in {1..22} X Y
do
sort -k2,2n $DATA/Homo_sapiens.chr$chrom.pfam.gtf | python rearrange.py | sort -k1,1 -k2,2n  > $DATA/chr$chrom.bed
done
cat $DATA/chr*.bed | sort -k1,1 -k2,2n > $DATA/all.bed
rm $DATA/chr*

# remove utrs and introns

bedtools intersect -a $DATA/all.bed -b $DATA/GRCh37.bed -wb | awk {'print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$17'} FS='\t' OFS='\t' | perl uniq.pl | perl -pe 's/"|;//g' > $DATA/foo.bed

#use appris to remove non-canonical transcripts; sort by ENSL gene_id and Pfam autoreg to merge uniqids into single domain

awk 'NR==FNR{a[$3];next}$15 in a{print $0}' /Users/jmh2tt/work/data/pmodeldata/appris_data.principal.txt $DATA/foo.bed | tr -s " " "\t" | sort -k27,27 -k25,25 -k1,1 -k2,2 > $DATA/blah.txt; mv $DATA/blah.txt $DATA/foo.bed

# domain coverage and rearranging:  // based on histograms, used 5x as a filter

python rearrange2.py -u $DATA/foo.bed $DATA/bar
bedtools intersect -a <(cat $DATA/bar | tr -s "\t" " " | cut -d " " -f 1-45 | tr -s " " "\t") -b <(awk '{if ($4>=5) print}' $DATA/coverage.bed) -wa -wb \
| awk '{ct[$1 $2 $3 $45]++; len[$1 $2 $3 $45]=$4; row[$1 $2 $3 $45]=$0} END {for (i in ct) print row[i],ct[i],len[i],ct[i]/len[i]}' \
| tr -s " " "\t" | cut -f 1-45,50- > $DATA/alluniq.bed; rm $DATA/bar
sort -k11,11 -k1,1 -k2,2n $DATA/foo.bed | python rearrange2.py -d $DATA/alldom.bed
bedtools intersect -a <(perl -pe 's/tag\s*\S*?(?=\n|\s)//g' $DATA/GRCh37.bed | perl -pe 's/ccds_id\s*\S*?(?=\n|\s)//g' | tr -s " " "\t" | sort -k1,1 -k2,2n) -b <(sort -k1,1 -k2,2n $DATA/alluniq.bed) -wa -wb \
| perl -pe 's/"|;//g' | awk '{ if ($8==$51 && $24==$67) print $0}' | tr -s " " "\t" | cut -f -35 \
| cat - <(bedtools intersect -v -a <(perl -pe 's/tag\s*\S*?(?=\n|\s)//g' $DATA/GRCh37.bed | perl -pe 's/ccds_id\s*\S*?(?=\n|\s)//g' | tr -s " " "\t" | sort -k1,1 -k2,2n) -b $DATA/alluniq.bed | sort -k1,1 -k2,2n) \
| sort -k1,1 -k2,2n | python nodom.py $DATA/foo
bedtools intersect -a $DATA/foo -b <(awk '{if ($4>=5) print}' $DATA/coverage.bed) -wa -wb \
| awk '{ct[$1 $2 $3 $24]++; len[$1 $2 $3 $24]=$4; row[$1 $2 $3 $24]=$0} END {for (i in ct) print row[i],ct[i],len[i],ct[i]/len[i]}' \
| tr -s " " "\t" | cut -f -24,29- > $DATA/nodom.bed; rm $DATA/foo
cat $DATA/alluniq.bed $DATA/nodom.bed | sort -k1,1 -k2,2n | cat <(printf "#header for nodoms:\n#chr,start,end,length,info\n#info field contains gene_id, transcript_id, exon_number, gene_name, gene_source, gene_biotype, transcript_name, transcript_source, exon_id, near_pfamA_id (if applicable, describes what domain it is near), uniq_id\n#header for domains separated by exon:\n#chr,start,end,length,pfam_database_ver,type_of_sequence,blank_field,strand,blank_field2,info\n#info field contains pfamA_id, gene_name, transcript_id, protein_id, pfamseq_acc, pfamseq_id, pfamA_acc, pfamA_auto_reg, gene_id, matched_transcript_id, expn_number, matched_gene_name, gene_source, gene_biotype, transcript_name, transcript_source, exon_id, uniq_id, ccds_id (if applicable)\n") - > $DATA/allregions.bed

# get MAF and convert from percent to fraction, variant type (FG), gene, domain, chr, start, end, impact, other info for analysis

cat <(grep "^#" $DATA/VEPEXAC.vcf) <(grep -v "^#" $DATA/VEPEXAC.vcf | sort -k1,1 -k2,2n) > $DATA/foo; mv $DATA/foo $DATA/VEPEXAC.vcf
bedtools intersect -a <(sort -k1,1 -k2,2n -k3,3n $DATA/alldom.bed | tr -s " " "\t" | cut -f 1,2,3,11,13,45 ) -b $DATA/VEPEXAC.vcf -sorted -wb | cut -f 1,2,3,4,5,6,10,11,14 | python var.py -d > $DATA/domint.bed
bedtools intersect -a <(sort -k1,1 -k2,2n -k3,3n $DATA/alluniq.bed | tr -s " " "\t" | cut -f 1,2,3,11,13,45,46,47,48 ) -b $DATA/VEPEXAC.vcf -sorted -wb | cut -f 1,2,3,4,5,6,7,8,9,13,14,17 | python var.py -d > $DATA/uniqint.bed
bedtools intersect -a <(sort -k1,1 -k2,2n -k3,3n $DATA/nodom.bed | tr -s " " "\t" | cut -f 1,2,3,12,24,25,26,27 ) -b $DATA/VEPEXAC.vcf -sorted -wb | cut -f 1,2,3,4,5,6,7,8,12,13,16 | python var.py -n > $DATA/nodomint.bed

####NEW####bedtools intersect -a <(sort -k1,1 -k2,2n -k3,3n $DATA/alluniq.bed | tr -s " " "\t" | cut -f 1,2,3,11,13,45 ) -b <(gzcat ~/Downloads/ExAC.r0.1.sites.vep.vcf.gz) -sorted -wb | cut -f 1,2,3,4,5,6,10,11,14 | python var.py > foodom

cat $DATA/uniqint.bed $DATA/nodomint.bed | sort -k1,1 -k2,2n | cat <(printf "#chr,start,end,ref,alt,pfamA_id,uniqid,covct,length_of_region,covratio,gene_symbol,maf,impact,codons,amino_acids,gene_id_csq,gene_symbol_csq,transcript_id_csq,exon_number_csq,polyphen,sift,protein_position,biotype\n") - > $DATA/allint.bed

# sort domain occurrence count from bill

sort -k2,2 $DATA/human_pfam.counts > $DATA/blah.bed; mv $DATA/blah.bed $DATA/human_pfam.counts

# get total bp for each domain, counts

cat $DATA/alldom.bed | tr -s " " "\t" | cut -f 4,11 | awk '{arr[$2]+=$1} END {for (i in arr) {print i,arr[i]}}' | sort -k1,1 > $DATA/sumlist.bed

# make db for queries AND filter variants by canonical transcripts

bash makedb.sh $DATA/var.db <(awk 'NR==FNR{a[$3];next}$18 in a{print $0}' $DATA/appris_data.principal.txt <(sed '1d' $DATA/allint.bed)) 

# pick one impact per variant for domains

bash query.sh $DATA/var.db '%' g 0 | sort -k6,6 -k8,8 | python formatvar.py > $DATA/allint2dom.bed

# for uniqids

bash query.sh $DATA/var.db '%' g 0 | sort -t '|' -k7,7 -k1,1 -k2,2 -k3,3 | python formatvar.py > $DATA/allint2uniq.bed

# do variant analysis by gene and maf (filters out entries that are "na," i.e, not dn or ds)

GENE="FLG"; MOD=g; MAF1=0.01 MAF2=''
bash lollipop.sh $DATA $OUT $GENE $MOD $MAF1 $MAF2

# make table of counts, non-syn, syn, total var, total bp/exome per domain

python maketable.py $DATA/allint2.bed > $DATA/foo.bed
python mergetable.py -f $DATA/foo.bed $DATA/human_pfam.counts $DATA/sumlist.bed > $DATA/dtable.txt; rm $DATA/foo.bed

# by uniqid a table

awk '{{if ($14=="ds") {sct[$7]++; ct[$7]++} else if ($14=="dn") nct[$7]++; ct[$7]++} len[$7]=$9; row[$7]=$6 " " $11 " " $7 " " $8 " " $9 " " $10} END {for (i in ct) print row[i],(nct[i]==0 ? nct[i]=0: nct[i]),(sct[i]==0 ? sct[i]=0: sct[i]),ct[i],nct[i]/(sct[i]==0 ? sct[i]+1: sct[i]),ct[i]/len[i]}' $DATA/allint2uniq.bed > $DATA/uniqtable.txt

# divergence by domain table

sed '1d' $DATA/dtablemaf.g.01-.txt | awk '{if ($6>0) print $0}' | python diverge.py -d | awk '{if ($3>3) print $0}' | awk '{gene[$1]=$10; row[$1]=$0} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]+=1}} END {for (i in rank) print row[i],rank[i]/NR}' | tr -s " " "\t" | sort -k11,11nr  > $DATA/diverge.01.txt

# creates gene-protein pair tables and compares z score/mad metric to rvis and omim
# uses omim genemap of all genes and phenotypes for those genes as well as a list of genes that hit certain keywords: recessive, haploinsufficient, de novo, dominant negative, autosomal dominant/heterozygous mut 

OMIM=~/work/omim

#python omim.py $OMIM/genemap2.txt $DATA/foo.txt

# awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/denovo.txt)) $DATA/foo.txt \
# | awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/domneg.txt)) - \
# | awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0 "n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/haplo.txt)) - \
# | awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/recessive.txt)) - \
# | awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/autodom.txt)) - >$DATA/omim.txt


bash div.sh $DATA/var.db $DATA/omim.txt g 0.00001 

# to make limited dtable-based dn/ds distributions and files use:

bash limit.sh $DATA g 0.001

# alldomint<-read.delim(paste(DATA,"/domint.bed",sep=""),header=FALSE)
# alldomint$V6=gsub(";","",alldomint$V6)
# m<-merge(alldomint,clans,by.x="V6",by.y="pfamA_id",all.x=TRUE)
# write.table(m,paste(DATA,"/alldomsclansint.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
# cat <(printf "domain\tchr\tstart\tend\tref\talt\tuniqid\tgene\teamaf\taamaf\tmaf\timpact\tcodon\taminoacid\tgeneid\tgenecsq\ttransid\texonno\tpolyphen\tsift\tproteinpos\tbiotype\tpfamacc\tclanacc\tclanid\tclandomcount\n") <(sed '1d' $DATA/alldomsclansint.txt) > $DATA/foo.txt; mv $DATA/foo.txt $DATA/alldomsclansint.txt

# perl -pe 's/ /\t/g' $DATA/alldom.bed | cut -f -45 > $DATA/foo.bed
# alldom<-read.delim(paste(DATA,"/foo.bed",sep=""),header=FALSE)
# alldom$V11=gsub(";","",alldom$V11)
# m<-merge(alldom,clans,by.x="V11",by.y="pfamA_id",all.x=TRUE)
# m<-m[,c(2:11,1,12:49)]
# write.table(m,paste(DATA,"/alldomsclans.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
# cat <(printf "chr\tstart\tend\tbplength\tdb\tseqtype\tblank\tstrand\tblank\tpfamA_id\tgene_name\ttranscript_id\tprotein_id\tpfamseq_acc\tpfamseq_id\tpfamA_acc\tpfamA_auto_reg\tgene_id\ttranscript_id2\texon_number\tgene_name2\tgene_source\tgene_biotype\ttranscript_name\ttranscript_source\texon_id\tuniq_id\tpfamA_acc\tclan_acc\tclan_id\tdom_cnt\n") <(sed '1d' $DATA/alldomsclans.txt) > $DATA/foo.txt; mv $DATA/foo.txt $DATA/alldomsclans.txt

R commands:

#generate dn/ds distro plots
uniqtable <- read.table("~/work/data/pmodeldata/uniqtable.txt", quote="\"")
u<-uniqtable[uniqtable["V1"]!='.',]
library(dplyr)
u<-arrange(u,desc(V10))
plot(u$V10,ylim=c(0,6),ylab="dn/ds",main="Dn/Ds distribution across domain regions")
abline(h=mean(u$V10))
abline(h=median(u$V10))
text(x=c(length(u$V10)*3/4,length(u$V10)*3/4),cex=0.7,y=c(mean(u$V10),median(u$V10)),pos=3,labels=c('mean','median'))

library(ggplot2)
library(plotrix)
library(RColorBrewer)
DATA<-"~/work/data/pmodeldata"
dtable <- read.delim("~/work/data/pmodeldata/dtable.txt")
clans <- read.delim("~/work/data/pmodeldata/count_human_pfam_clan.tab")
m<-merge(dtable,clans,by.x="domain",by.y="pfamA_id")
write.table(m,paste(DATA,"/clans.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
echo -e "clan_acc\tclan_id\tnonsynct\tsynct\ttotalvarct\tdomcount\ttotalbp\tmmaf\tdn.ds\tvar.bp.ratio" > $DATA/ctable.txt
sed '1d' $DATA/clans.txt | sort -k12,12 | grep -v -w "NA" | bedtools groupby -g 12,13 -c 3,4,5,6,7,8,9,10 -o sum,sum,sum,sum,sum,sum,sum,sum | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"; if ($4!=0) printf $3/$4; else printf $3/(1+$4); printf "\t"$5/$7"\n"}' >> $DATA/ctable.txt
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
	if (tot >= max(m)){
		return(widthlist)
	}
}
subplot <- function(i=dtable,xll=0,xul=6,yll=0,yul=3,pchcolumn="totalvarct",colorcolumn="mmaf",nbin=6,x="dn.ds",y="domcount",xlab="dn/ds",ylab="domain occurrences in exome",xlimiter="dn.ds",ylimiter="domcount",label="domain",labelbool=TRUE,colorbool=TRUE,pchbool=TRUE,xlimitbool=TRUE,ylimitbool=TRUE){
	if (xlimitbool==TRUE){
		i=i[xll<i[xlimiter] & i[xlimiter]<=xul,]
	}
	if (ylimitbool==TRUE){
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
	if (labelbool==TRUE){
		text(t(i[x]),t(i[y]),labels=t(i[label]),adj=c(.5,-.4),cex=0.5,offset=0.5)
	}
	# par(fig=c(0,.1,.3,.5),xpd=NA)
	if (colorbool==TRUE){
		color.legend((max(i[x])-min(i[x]))*.985+min(i[x]),(max(i[y])-min(i[y]))*.75+min(i[y]),max(i[x]),max(i[y]),legend=w,rect.col=c,gradient="y")
	}
	if (pchbool==TRUE){
		legend(par()$usr[1],mean(par()$usr[3:4]),border=rep("white",nbin),legend=pchleg,pch=seq(0,nbin-1),bty="n",xpd=TRUE,xjust=0,yjust=0.5,cex=0.65) #inset=c(0,0) #min(i[x]),(max(i[y])-min(i[y]))*.75+min(i[y]),(max(i[x])-min(i[x]))*.05+min(i[x]),max(i[y])
	}
}
plot(dtable$var.bp.ratio,log(dtable$totalbp,base=10),xlab="variants/total length of domain",ylab="total length of domain (log)",cex=cex,col=col)
color.legend(.10,3,.14,5,legend=w,rect.col=c,gradient="y")

f<-dtable[order(-dtable$dn.ds),]
l<-f[f$dn.ds<6]
plot(l$dn.ds,ylab="dn/ds")
text(l$dn.ds,labels=l$domain,cex=.5,pos=4,srt=45)

f<-ctable[order(-ctable$dn.ds),]
l<-f[f$dn.ds<6,]
plot(l$dn.ds,ylab="dn/ds")
text(l$dn.ds,labels=l$clan_id,cex=.5,pos=4,srt=45)
