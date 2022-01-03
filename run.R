#--- install necessary packages

if(!require(optparse)){
	install.packages("optparse")
    library(optparse)
}


option_list = list(
make_option(c("-f", "--file"), type="character", default=NA,help="input file", metavar="character"),
make_option(c("--output"), type="character", default=NA,help="output folder", metavar="character"),
make_option(c("--script_folder"), type ="character", default=NA,help="folder which contains cd8 ... markers ", metavar="character" )
);
parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)

if(is.na(opt$file)) {
	print('Missing input file; run Rscript this.R --help');
} else if(is.na(opt$output)) {
	print('Missing output folder; run Rscript this.R --help');
} else if(is.na(opt$script_folder)) {
	print('Missing folder which contains cd8 etc makers file; run Rscript this.R --help');
} else {
inputfilename = opt$file;
outputfolder = opt$output;
scriptfolder = opt$script_folder;
system(paste("mkdir ",outputfolder,sep=''))
#--- uploading libraries
if(!require(org.Hs.eg.db)){
	if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	BiocManager::install("org.Hs.eg.db")
    library(org.Hs.eg.db)
}
if(!require(annotate)){
	if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	BiocManager::install("annotate")
    suppressMessages(library(annotate))
}

if(!require(data.table)){
	install.packages("data.table")
    library(data.table)
}

if(!require(ggpubr)){
	if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	BiocManager::install("ggpubr")
    library(ggpubr)
}

#---- data input
#file = "data_RNA_Seq_v2_mRNA_median_Zscores.txt"
file = inputfilename
data <- read.table(file,header=T,sep="\t")
getSYMBOL(as.character(data[,2]), data='org.Hs.eg') -> data[,1]
#-----
mesen <- read.table(paste(scriptfolder,'/mesenchymal.txt',sep=""))
epithelial <- read.table(paste(scriptfolder,'/epithelial.txt',sep=""))
cd8 <- read.table(paste(scriptfolder,'/cd8.txt',sep=""))
checkpoint <- read.table(paste(scriptfolder,'/immunecheckpoint.txt',sep=""))
inflam <- read.table(paste(scriptfolder,'/immuno-inflamatory.txt',sep=""))
#----
cd8=intersect(cd8$V1,data[,1])
epithelial=intersect(epithelial$V1,data[,1])
mesen=intersect(mesen$V1,data[,1])
checkpoint=intersect(checkpoint$V1,data[,1])
inflam=intersect(inflam$V1,data[,1])
#----
data[c(which(data[,1] %in% mesen)),] -> d_mesen
rownames(d_mesen) <- d_mesen[,1]
data[c(which(data[,1] %in% epithelial)),] -> d_epi
rownames(d_epi) <- d_epi[,1]
data[c(which(data[,1] %in% cd8)),] -> d_cd8
rownames(d_cd8) <- d_cd8[,1]
data[c(which(data[,1] %in% checkpoint)),] -> d_checkpoint
rownames(d_checkpoint) <- d_checkpoint[,1]
data[c(which(data[,1] %in% inflam)),] -> d_inflam
rownames(d_inflam) <- d_inflam[,1]


#----- Using ratio
round(colMeans(d_mesen[,-c(1,2)]),2) -> y
round(colMeans(d_epi[,-c(1,2)]),2) -> z
colnames(d_mesen[,-c(1,2)]) -> x
emt <- data.frame(name=x,mesen=y,epi=z)
emt$score <- emt$mesen-emt$epi
ifelse (emt$score <= quantile(emt$score,0.25), "epi",
	ifelse (emt$score >= quantile(emt$score,0.75),"mesen","intermediate")
 ) -> emt$class

nrow(emt) -> total_samples
nrow(emt[emt$class=="epi",]) -> total_epi
nrow(emt[emt$class=="mesen",]) -> total_mesen

write.table(emt,paste(outputfolder,'/classification.txt',sep=''),row.names=F,quote=F,sep="\t")
#---- using clustering
set.seed(786)
rbind(d_mesen,d_epi) -> values
dist(t(values[,-c(1,2)])) -> d
prcomp(d) -> pr
k <- kmeans(t(values[,-c(1,2)]),4)
#----
jpeg(paste(outputfolder,'/pca.jpeg',sep=''),unit='in',res=300,height=3*4,width=3*5)
rbPal <- colorRampPalette(c('cyan','darkblue'))
opar<-par(mfrow=c(4,5))
plot(pr$rotation[,1],pr$rotation[,2],col=c('yellow','orange','darkblue','red4')[as.factor(emt$class)],pch=20,main='clustering: EMT score',xlab='eigenvector 1',ylab='eigenvector 2')

plot(pr$rotation[,1],pr$rotation[,2],col=c('yellow','orange','darkblue','red4')[as.vector(k$cluster)],pch=20,main='clustering: Kclust',xlab='eigenvector 1',ylab='eigenvector 2')

write.table(pr$rotation,"pca.txt",quote=F,sep="\t")

for(i in c(1:16)) 
	{
		x   <- as.numeric(values[i,-c(1,2)])
		Col <- rbPal(10)[as.numeric(cut(x,breaks = 100))]
		Col[!complete.cases(Col)] <- 'black'
		plot(pr$rotation[,1],pr$rotation[,2],col=Col,pch=20,main=values[i,1],xlab='eigenvector 1',ylab='eigenvector 2')
		}
dev.off()

#-----
c('mesen','epi') -> classes

print('working on inflamatory markers')
t(d_inflam[,-c(1,2)]) -> value_inflam
cbind(value_inflam,emt) -> value_inflam
xx=value_inflam
n <- xx[xx$class %in% classes,]
droplevels(as.factor(n$class)) -> n$class

jpeg(paste(outputfolder,'/inflam.jpeg',sep=''),unit='in',res=300,width=2*2,height=3*9)

output=matrix(nrow=ncol(n)-5,ncol=6)
colnames(output) <- c('gene','epithelial','mesenchymal','p-value','number_epithelial','number_mesenchymal')

opar<-par(mfrow=c(9,2))
for(i in c(1: (ncol(n)-5))){
	xvalue = n[n$class=='mesen',i]
	yvalue = n[n$class=='epi',i]
	l_xvalue=length(xvalue[complete.cases(xvalue)])
	l_yvalue=length(yvalue[complete.cases(yvalue)])
	output[i,5] <- paste(l_yvalue,total_epi,total_samples,sep="/"); output[i,6] <- paste(l_xvalue,total_mesen,total_samples,sep="/")
	mean_xvalue=mean(xvalue); sd_xvalue=sd(xvalue); xx <- paste(round(mean_xvalue,2),round(sd_xvalue,2),sep=" ± ")
	mean_yvalue=mean(yvalue); sd_yvalue=sd(yvalue); yy <- paste(round(mean_yvalue,2),round(sd_yvalue,2),sep=" ± ")
	output[i,1] <- colnames(n)[i]
	output[i,2] <- yy; output[i,3] <- xx
	zvalue = ks.test(xvalue,yvalue)
	output[i,4] <- formatC(zvalue$p.value, format = "e", digits = 2)
	l=paste('p-value',formatC(zvalue$p.value, format = "e", digits = 2),sep='=')
	boxplot(n[,i]~n[,ncol(n)],main=paste(colnames(n)[i],"\n",l,sep=""))
}
dev.off() 
write.table(output,paste(outputfolder,"/inflam.stats.txt",sep=''),sep="\t",quote=F,row.names=F)
write.table(n,paste(outputfolder,"/inflam.txt",sep=''),sep="\t",quote=F,row.names=F)
#---
print('working on checkpoint markers')
t(d_checkpoint[,-c(1,2)]) -> value_checkpoint
cbind(value_checkpoint,emt) -> value_checkpoint
xx=value_checkpoint
n <- xx[xx$class %in% classes,]
droplevels(as.factor(n$class)) -> n$class

#jpeg('checkpoint.jpeg',unit='in',res=300,width=2*(ncol(n) - 5),height=3)
jpeg(paste(outputfolder,'/checkpoint.jpeg',sep=''),unit='in',res=300,width=2*2,height=3*4)
output=matrix(nrow=ncol(n)-5,ncol=6)
colnames(output) <- c('gene','epithelial','mesenchymal','p-value','number_epithelial','number_mesenchymal')
opar<-par(mfrow=c(4,2))
for(i in c(1: (ncol(n)-5))){
	xvalue = n[n$class=='mesen',i]
	yvalue = n[n$class=='epi',i]
	l_xvalue=length(xvalue[complete.cases(xvalue)])
	l_yvalue=length(yvalue[complete.cases(yvalue)])
	output[i,5] <- paste(l_yvalue,total_epi,total_samples,sep="/"); output[i,6] <- paste(l_xvalue,total_mesen,total_samples,sep="/")
	mean_xvalue=mean(xvalue); sd_xvalue=sd(xvalue); xx <- paste(round(mean_xvalue,2),round(sd_xvalue,2),sep=" ± ")
	mean_yvalue=mean(yvalue); sd_yvalue=sd(yvalue); yy <- paste(round(mean_yvalue,2),round(sd_yvalue,2),sep=" ± ")
	output[i,1] <- colnames(n)[i]
	output[i,2] <- yy; output[i,3] <- xx
	zvalue = ks.test(xvalue,yvalue)
	output[i,4] <- formatC(zvalue$p.value, format = "e", digits = 2)
	l=paste('p-value',formatC(zvalue$p.value, format = "e", digits = 2),sep='=')
	boxplot(n[,i]~n[,ncol(n)],main=paste(colnames(n)[i],"\n",l,sep=""))
}
dev.off() 
write.table(output,paste(outputfolder,"/checkpoint.stats.txt",sep=''),sep="\t",quote=F,row.names=F)
write.table(n,paste(outputfolder,"/checkpoint.txt",sep=''),sep="\t",quote=F,row.names=F)
#---
print('working on cd8 markers')
t(d_cd8[,-c(1,2)]) -> value_cd8
cbind(value_cd8,emt) -> value_cd8
xx=value_cd8
n <- xx[xx$class %in% classes,]
droplevels(as.factor(n$class)) -> n$class

jpeg(paste(outputfolder,'/cd8.jpeg',sep=''),unit='in',res=300,width=2*(ncol(n) - 5),height=3)
output=matrix(nrow=ncol(n)-5,ncol=6)
colnames(output) <- c('gene','epithelial','mesenchymal','p-value','number_epithelial','number_mesenchymal')
opar<-par(mfrow=c(1,ncol(n) - 5))
for(i in c(1: (ncol(n)-5))){
	xvalue = n[n$class=='mesen',i]
	yvalue = n[n$class=='epi',i]
	l_xvalue=length(xvalue[complete.cases(xvalue)])
	l_yvalue=length(yvalue[complete.cases(yvalue)])
	output[i,5] <- paste(l_yvalue,total_epi,total_samples,sep="/"); output[i,6] <- paste(l_xvalue,total_mesen,total_samples,sep="/")
	mean_xvalue=mean(xvalue); sd_xvalue=sd(xvalue); xx <- paste(round(mean_xvalue,2),round(sd_xvalue,2),sep=" ± ")
	mean_yvalue=mean(yvalue); sd_yvalue=sd(yvalue); yy <- paste(round(mean_yvalue,2),round(sd_yvalue,2),sep=" ± ")
	output[i,1] <- colnames(n)[i]
	output[i,2] <- yy; output[i,3] <- xx
	zvalue = ks.test(xvalue,yvalue)
	output[i,4] <- formatC(zvalue$p.value, format = "e", digits = 2)
	l=paste('p-value',formatC(zvalue$p.value, format = "e", digits = 2),sep='=')
	boxplot(n[,i]~n[,ncol(n)],main=paste(colnames(n)[i],"\n",l,sep=""))
}
dev.off() 
write.table(output,paste(outputfolder,"/cd8.stats.txt",sep=""),sep="\t",quote=F,row.names=F)
write.table(n,paste(outputfolder,"/cd8.txt",sep=""),sep="\t",quote=F,row.names=F)
}
