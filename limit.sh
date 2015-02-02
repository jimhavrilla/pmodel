bash query.sh $1/lolli.db '%' $2 $3 $4| python formatvar.py | sort -k6,6 -k8,8 > $1/allint2maf.$2.$3-$4.bed
python maketable.py $1/allint2maf.$2.$3-$4.bed > $1/foo.bed
python mergetable.py -f $1/foo.bed $1/human_pfam.counts $1/sumlist.bed > $1/dtablemaf.$2.$3-$4.txt; rm $1/foo.bed

# in R:
dtable <- read.delim("~/work/data/pmodeldata/dtablemaf.g.0.01-.txt")
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
mtext("maf<0.0001", side = 3, line = -2, outer = TRUE, col=2)