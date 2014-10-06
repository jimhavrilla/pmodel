#for initial p53 pass:
#python extract.py
#hgTables is pfam from UCSC GenomeBrowser
#cat $DATA/hgTables | awk '{t=$1;$1=$2;$2=t;t=$2;$2=$3;$3=t;t=$3;$3=$4;$4=t;print;}' | tr -s " " | tr -s " " "\t" | tail -n +2 | sort -k1,1 -k2,2n > $DATA/hgTables.bed
#bedtools intersect -a $DATA/hgTables.bed -b $DATA/gene.bed -sorted | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,5 -o distinct > $DATA/p53domain.bed
#bill's domain count file:
export DATA=~/work/data/pmodeldata
mysql -N --raw -h wrpxdb.its.virginia.edu -u web_user -pfasta_www pfam27 < count_human_pfam.sql > human_pfam.counts
#for bill's gtfs:
grep ^# $DATA/ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf; grep -v ^# -h $DATA/ESP*.vcf) > $DATA/foo.vcf; mv $DATA/foo.vcf $DATA/ESPALL.vcf
perl variant_effect_predictor.pl -i $DATA/ESPALL.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --total_length -o $DATA/VEPESPALL.vcf --vcf --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --force_overwrite --port 3337 --dir ~/.vep/ --compress "gunzip -c"
gzcat $DATA/Homo_sapiens.GRCh37.75.gtf.gz | grep 'protein_coding\texon' | perl -pe 's/protein_coding\texon\t//g' | grep -e '^\d*[0-9X-Y]\t'| perl gtf2bed.pl > $DATA/GRCh37.bed
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
bedtools intersect -a <(perl -pe 's/tag\s*\S*?(?=\n|\s)//g' $DATA/GRCh37.bed | perl -pe 's/ccds_id\s*\S*?(?=\n|\s)//g' | tr -s " " "\t" | sort -k1,1 -k2,2n) -b <(sort -k1,1 -k2,2n $DATA/alluniq.bed) -wa -wb | awk '{ if ($8==$51 && $24==$67) print $0}' | tr -s " " "\t" | cut -f -35 | cat - <(bedtools intersect -v -a <(perl -pe 's/tag\s*\S*?(?=\n|\s)//g' $DATA/GRCh37.bed | perl -pe 's/ccds_id\s*\S*?(?=\n|\s)//g' | tr -s " " "\t" | sort -k1,1 -k2,2n) -b $DATA/alluniq.bed | sort -k1,1 -k2,2n) | sort -k1,1 -k2,2n | python nodom.py $DATA/nodom.bed
cat $DATA/alluniq.bed $DATA/nodom.bed | sort -k1,1 -k2,2n | cat <(printf "#header for nodoms:\n#chr,start,end,length,info\n#info field contains gene_id, transcript_id, exon_number, gene_name, gene_source, gene_biotype, transcript_name, transcript_source, exon_id, near_pfamA_id (if applicable, describes what domain it is near), uniq_id\n#header for domains separated by exon:\n#chr,start,end,length,pfam_database_ver,type_of_sequence,blank_field,strand,blank_field2,info\n#info field contains pfamA_id, gene_name, transcript_id, protein_id, pfamseq_acc, pfamseq_id, pfamA_acc, pfamA_auto_reg, gene_id, matched_transcript_id, expn_number, matched_gene_name, gene_source, gene_biotype, transcript_name, transcript_source, exon_id, uniq_id, ccds_id (if applicable)\n") - > $DATA/allregions.bed
# get MAF and convert from percent to fraction, variant type (FG), gene, domain, chr, start, end, impact, other info for analysis
bedtools intersect -a <(sort -k1,1 -k2,2n -k3,3n $DATA/alldom.bed | tr -s " " "\t" | cut -f 1,2,3,11,13,45 ) -b $DATA/VEPESPALL.vcf -sorted -wb | cut -f 1,2,3,4,5,6,10,11,14 | python var.py -d > $DATA/domint.bed
bedtools intersect -a <(sort -k1,1 -k2,2n -k3,3n $DATA/alluniq.bed | tr -s " " "\t" | cut -f 1,2,3,11,13,45 ) -b $DATA/VEPESPALL.vcf -sorted -wb | cut -f 1,2,3,4,5,6,10,11,14 | python var.py -d > $DATA/uniqint.bed
bedtools intersect -a <(sort -k1,1 -k2,2n -k3,3n $DATA/allnodom.bed | tr -s " " "\t" | cut -f 1,2,3,9,17 ) -b $DATA/VEPESPALL.vcf -sorted -wb | cut -f 1,2,3,4,5,9,10,13 | python var.py -n > $DATA/nodomint.bed
cat $DATA/uniqint.bed $DATA/nodomint.bed | sort -k1,1 -k2,2n | cat <(printf "#chr,start,end,ref,alt,pfamA_id,uniqid,gene_symbol,ea_maf,aa_maf,all_maf,impact,codons,amino_acids,gene_id_csq,gene_symbol_csq,transcript_id_csq,exon_number_csq,polyphen,sift,protein_position,biotype\n") - > $DATA/allint.bed
#sort domain occurrence count from bill
sort -k2,2 $DATA/human_pfam.counts > $DATA/foo.bed; mv $DATA/foo.bed $DATA/human_pfam.counts
#get total bp for each domain, counts
cat $DATA/alldom.bed | tr -s " " "\t" | cut -f 4,11 | awk '{arr[$2]+=$1} END {for (i in arr) {print i,arr[i]}}' | sort -k1,1 > $DATA/sumlist.bed
# pick one impact per variant
sed '1d' $DATA/allint.bed | python formatvar.py > $DATA/allint2.bed
# make table of counts, non-syn, syn, total var, total bp/exome per domain
sort -k6,6 $DATA/allint2.bed | python maketable.py > $DATA/foo.bed
python mergetable.py $DATA/foo.bed $DATA/human_pfam.counts $DATA/sumlist.bed > $DATA/dtable.txt; rm $DATA/foo.bed

alldomint<-read.delim(paste(DATA,"/domint.bed",sep=""),header=FALSE)
alldomint$V6=gsub(";","",alldomint$V6)
m<-merge(alldomint,clans,by.x="V6",by.y="pfamA_id",all.x=TRUE)
write.table(m,paste(DATA,"/alldomsclans.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

perl -pe 's/ /\t/g' $DATA/alldom.bed | cut -f -45 > $DATA/foo.bed
alldom<-read.delim(paste(DATA,"/foo.bed",sep=""),header=FALSE)
alldom$V11=gsub(";","",alldom$V11)
m<-merge(alldom,clans,by.x="V11",by.y="pfamA_id",all.x=TRUE)
write.table(m,paste(DATA,"/alldomsclans.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

R commands:

library(ggplot2)
library(plotrix)
library(RColorBrewer)
DATA<-"~/work/data/pmodeldata"
dtable <- read.delim("~/work/data/pmodeldata/dtable.txt")
clans <- read.delim("~/work/data/pmodeldata/count_human_pfam_clan.tab")
m<-merge(dtable,clans,by.x="domain",by.y="pfamA_id")
write.table(m,paste(DATA,"/clans.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
echo -e "clan_acc\tclan_id\tnonsynct\tsynct\ttotalvarct\tdomcount\ttotalbp\tmmaf\tdn.ds\tvar.bp.ratio" > $DATA/ctable.txt
sed '1d' $DATA/clans.txt | sort -k12,12 | bedtools groupby -g 12,13 -c 3,4,5,6,7,8,9,10 -o sum,sum,sum,sum,sum,sum,sum,sum | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"; if ($4!=0) printf $3/$4; else printf $3/(1+$4); printf "\t"$5/$7"\n"}' >> $DATA/ctable.txt
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
