#graphics.r

#generate dn/ds distro plots

library(dplyr)
uniqtablefilter <- read.table("~/work/data/pmodeldata/uniqtablefilter.g.0.0001-.txt", quote="\"")
u<-uniqtablefilter
u<-arrange(u,desc(V11))
plot(u$V11,ylim=c(0,6),ylab="dn/ds",main="dN/dS distribution across domain regions") # dnds dist
abline(h=mean(u$V11))
abline(h=median(u$V11))
text(x=c(length(u$V11)*3/4,length(u$V11)*3/4),cex=0.7,y=c(mean(u$V11),median(u$V11)),pos=3,labels=c(paste('mean=',mean(u$V11),sep=''),paste('median=',median(u$V11),sep='')))
text(x=c(length(n$V11)*.3),y=5,labels=c("MAF>0.001%")) ## for common
hist(u$V11,breaks=100,xlab="dN/dS",main="dN/dS distribution for domains",xlim=c(0,5)) # histogram
legend('topright','MAF>0.001%',bty='n') ## for common
plot(u$V7,u$V11,xlab="length",ylab="dn/ds",main="dN/dS vs length of region for domain regions") # length vs dnds
legend('topright','MAF>0.001%',bty='n') ## for common
plot(u$V14,u$V11,xlab="domain count across genome",ylab="dn/ds",main="dN/dS vs domain prevalence") # domain counts vs dnds
legend('topright','MAF>0.001%',bty='n') ## for common
plot(u$V12,u$V11,xlab="density",ylab="dn/ds",main="dN/dS vs density of region for domain regions")
legend('topright','MAF>0.001%',bty='n') ## for common

nodomtablefilter <- read.table("~/work/data/pmodeldata/nodomtablefilter.g.0.0001-.txt", quote="\"")
n<-nodomtablefilter
n<-arrange(n,desc(V11))
plot(n$V11,ylim=c(0,6),ylab="dn/ds",main="dN/dS distribution across non-domain regions")
abline(h=mean(n$V11))
abline(h=median(n$V11))
text(x=c(length(n$V11)*3/4,length(n$V11)*3/4),cex=0.7,y=c(mean(n$V11),median(n$V11)),pos=3,labels=c(paste('mean=',mean(n$V11),sep=''),paste('median=',median(n$V11),sep='')))
text(x=c(length(n$V11)*.3),y=5,labels=c("MAF>0.001%")) ## for common
hist(n$V11,breaks=100,xlab="dN/dS",main="dN/dS distribution for non-domain regions",xlim=c(0,5))
legend('topright','MAF>0.001%',bty='n') ## for common
plot(n$V7,n$V11,xlab="length",ylab="dn/ds",main="dN/dS vs length of region for non-domain regions")
legend('topright','MAF>0.001%',bty='n') ## for common
plot(n$V12,n$V11,xlab="density",ylab="dn/ds",main="dN/dS vs density of region for non-domain regions")
legend('topright','MAF>0.001%',bty='n') ## for common

genetablefilter <- read.table("~/work/data/pmodeldata/genetablefilter.g.0.0001-.txt", quote="\"")
g<-genetablefilter
g<-arrange(g,desc(V6))
png("~/Documents/workimages/maf0.01/genednds.png",width=1603,height=800,units="px")
plot(g$V6,ylim=c(0,6),ylab="dn/ds",main="dN/dS distribution across genes")
abline(h=mean(g$V6))
abline(h=median(g$V6))
text(x=c(length(g$V6)*3/4,length(g$V6)*3/4),cex=0.7,y=c(mean(g$V6),median(g$V6)),pos=3,labels=c(paste('mean=',mean(g$V6),sep=''),paste('median=',median(g$V6),sep='')))
text(x=c(length(g$V6)*.3),y=5,labels=c("MAF>0.001%")) ## for common
dev.off()
png("~/Documents/workimages/maf0.01/genedndshist.png",width=1603,height=800,units="px")
hist(g$V6,breaks=100,xlab="dN/dS",main="dN/dS distribution for genes",xlim=c(0,5))
legend('topright','MAF>0.001%',bty='n') ## for common
dev.off()
png("~/Documents/workimages/maf0.01/genedndslength.png",width=1603,height=800,units="px")
plot(g$V2,g$V6,xlab="length",ylab="dn/ds",main="dN/dS vs length of region for genes")
legend('topright','MAF>0.001%',bty='n') ## for common
dev.off()
png("~/Documents/workimages/maf0.01/genedndsdensity.png",width=1603,height=800,units="px")
plot(g$V7,g$V6,xlab="density",ylab="dn/ds",main="dN/dS vs density of region for genes")
legend('topright','MAF>0.001%',bty='n') ## for common
dev.off()

uniqgenetablefilter <- read.table("~/work/data/pmodeldata/uniqgenetablefilter.g.0.0001-.txt", quote="\"")
gu<-uniqgenetablefilter
gu<-arrange(gu,desc(V6))
plot(gu$V6,ylim=c(0,6),ylab="dn/ds",main="dN/dS distribution across genes without nodoms")
abline(h=mean(gu$V6))
abline(h=median(gu$V6))
text(x=c(length(gu$V6)*3/4,length(gu$V6)*3/4),cex=0.7,y=c(mean(gu$V6),median(gu$V6)),pos=3,labels=c(paste('mean=',mean(gu$V6),sep=''),paste('median=',median(gu$V6),sep='')))
text(x=c(length(gu$V6)*.3),y=5,labels=c("MAF>0.001%")) ## for common
hist(gu$V6,breaks=100,xlab="dN/dS",main="dN/dS distribution for genes without nodoms",xlim=c(0,5))
legend('topright','MAF>0.001%',bty='n') ## for common
plot(gu$V2,gu$V6,xlab="length",ylab="dn/ds",main="dN/dS vs length of genes without nodoms")
legend('topright','MAF>0.001%',bty='n') ## for common
plot(gu$V7,gu$V6,xlab="density",ylab="dn/ds",main="dN/dS vs density of region for genes without nodoms")
legend('topright','MAF>0.001%',bty='n') ## for common

regionstable <- read.table("~/work/data/pmodeldata/regionstable.g.0.0001-.txt", quote="\"")
r<-regionstable
r<-arrange(r,desc(V11))
plot(r$V11,ylim=c(0,6),ylab="dn/ds",main="dN/dS distribution across all protein regions")
abline(h=mean(r$V11))
abline(h=median(r$V11))
text(x=c(length(r$V11)*3/4,length(r$V11)*3/4),cex=0.7,y=c(mean(r$V11),median(r$V11)),pos=3,labels=c(paste('mean=',mean(r$V11),sep=''),paste('median=',median(r$V11),sep='')))
text(x=c(length(r$V11)*.3),y=5,labels=c("MAF>0.001%")) ## for common
hist(r$V11,breaks=100,xlab="dN/dS",main="dN/dS distribution for all protein regions",xlim=c(0,5))
legend('topright','MAF>0.001%',bty='n') ## for common
plot(r$V7,r$V11,xlab="length",ylab="dn/ds",main="dN/dS vs length of region for all protein regions")
legend('topright','MAF>0.001%',bty='n') ## for common
plot(r$V12,r$V11,xlab="density",ylab="dn/ds",main="dN/dS vs density of region for all protein regions")
legend('topright','MAF>0.001%',bty='n') ## for common


#for density plots and other plots

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


#limiting the maf and then testing distribution of domcount
dtable <- read.delim("~/work/data/pmodeldata/dtablemaf.g.0.001-.txt")
par(mfrow=c(2,4))
hist(t(dtable[dtable["domcount"]==1,]["dn.ds"]),main="domcount=1",xlab="dn.ds",breaks=30,xlim=c(0,ceiling(max(dtable["dn.ds"]))))
text(.85*par('usr')[2],.9*par('usr')[4],labels=paste("mean:",round(mean(t(dtable[dtable["domcount"]==1,]["dn.ds"])),2)))
text(.85*par('usr')[2],.86*par('usr')[4],labels=paste("median:",round(median(t(dtable[dtable["domcount"]==1,]["dn.ds"])),2)))
hist(t(dtable[dtable["domcount"]>1 & dtable["domcount"]<10,]["dn.ds"]),main="10>domcount>1",xlab="dn.ds",breaks=30,xlim=c(0,ceiling(max(dtable["dn.ds"])))) 
text(.85*par('usr')[2],.9*par('usr')[4],labels=paste("mean:",round(mean(t(dtable[dtable["domcount"]>1 & dtable["domcount"]<10,]["dn.ds"])),2)))
text(.85*par('usr')[2],.86*par('usr')[4],labels=paste("median:",round(median(t(dtable[dtable["domcount"]>1 & dtable["domcount"]<10,]["dn.ds"])),2)))
hist(t(dtable[dtable["domcount"]>10 & dtable["domcount"]<25,]["dn.ds"]),main="25>domcount>10",xlab="dn.ds",breaks=30,xlim=c(0,ceiling(max(dtable["dn.ds"]))))
text(.85*par('usr')[2],.9*par('usr')[4],labels=paste("mean:",round(mean(t(dtable[dtable["domcount"]>10 & dtable["domcount"]<25,]["dn.ds"])),2)))
text(.85*par('usr')[2],.86*par('usr')[4],labels=paste("median:",round(median(t(dtable[dtable["domcount"]>10 & dtable["domcount"]<25,]["dn.ds"])),2)))
hist(t(dtable[dtable["domcount"]>25 & dtable["domcount"]<50,]["dn.ds"]),main="50>domcount>25",xlab="dn.ds",breaks=30,xlim=c(0,ceiling(max(dtable["dn.ds"]))))
text(.85*par('usr')[2],.9*par('usr')[4],labels=paste("mean:",round(mean(t(dtable[dtable["domcount"]>25 & dtable["domcount"]<50,]["dn.ds"])),2)))
text(.85*par('usr')[2],.86*par('usr')[4],labels=paste("median:",round(median(t(dtable[dtable["domcount"]>25 & dtable["domcount"]<50,]["dn.ds"])),2)))
hist(t(dtable[dtable["domcount"]>50 & dtable["domcount"]<75,]["dn.ds"]),main="75>domcount>50",xlab="dn.ds",breaks=30,xlim=c(0,ceiling(max(dtable["dn.ds"]))))
text(.85*par('usr')[2],.9*par('usr')[4],labels=paste("mean:",round(mean(t(dtable[dtable["domcount"]>50 & dtable["domcount"]<75,]["dn.ds"])),2)))
text(.85*par('usr')[2],.86*par('usr')[4],labels=paste("median:",round(median(t(dtable[dtable["domcount"]>50 & dtable["domcount"]<75,]["dn.ds"])),2)))
hist(t(dtable[dtable["domcount"]>75 & dtable["domcount"]<100,]["dn.ds"]),main="100>domcount>75",xlab="dn.ds",breaks=30,xlim=c(0,ceiling(max(dtable["dn.ds"]))))
text(.85*par('usr')[2],.9*par('usr')[4],labels=paste("mean:",round(mean(t(dtable[dtable["domcount"]>75 & dtable["domcount"]<100,]["dn.ds"])),2)))
text(.85*par('usr')[2],.86*par('usr')[4],labels=paste("median:",round(median(t(dtable[dtable["domcount"]>75 & dtable["domcount"]<100,]["dn.ds"])),2)))
hist(t(dtable[dtable["domcount"]>100,]["dn.ds"]),main="domcount>100",xlab="dn.ds",breaks=30,xlim=c(0,ceiling(max(dtable["dn.ds"]))))
text(.85*par('usr')[2],.9*par('usr')[4],labels=paste("mean:",round(mean(t(dtable[dtable["domcount"]>100,]["dn.ds"])),2)))
text(.85*par('usr')[2],.86*par('usr')[4],labels=paste("median:",round(median(t(dtable[dtable["domcount"]>100,]["dn.ds"])),2)))
mtext("maf>0.001", side = 3, line = -2, outer = TRUE, col=2)