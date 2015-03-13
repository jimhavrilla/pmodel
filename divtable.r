args <- commandArgs(TRUE)
DATA<-args[1]
NAME<-args[2]
diverge <- read.delim(paste(DATA,"/diverge.pair.txt",sep=""),header=FALSE)
y=diverge[diverge["V13"]>=.8,][,c("V1","V2","V3","V6","V10","V11","V12","V13")]
write.table(y,paste(DATA,"/dps.pair.",NAME,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)