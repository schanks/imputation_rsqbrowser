#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(data.table))

percent=function(x){
   x=x[which(!is.na(x))]
   perc=length(which(x>=0.8))/length(x)
   return(perc)
}

optionList <- list(
  make_option(c("-r", "--results"), type="character", help="Merged result file"),
  make_option(c("--one"), type="character", help="Site list- 100"),
  make_option(c("--two"), type="character", help="Site list- 200"),
  make_option(c("--three"), type="character", help="Site list- 300"),
  make_option(c("--four"), type="character", help="Site list- 400"),
  make_option(c("--five"), type="character", help="Site list- 500"),
  make_option(c("--six"), type="character", help="Site list- 1000"),
  make_option(c("--seven"), type="character", help="Site list- 2000"))

parser <- OptionParser(
  option_list=optionList
)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

results = fread(opt$results, stringsAsFactors = FALSE)
results[is.na(results)]=0
results$MAF=as.numeric(results$AF)
results$MAF[which(results$MAF>0.5)]<-(1-results$MAF[which(results$MAF>0.5)])
breaks=c(seq(0, 0.00175, 0.00025), seq(0.002, 0.50014, 0.001))
results$bin=cut(results$MAF, breaks=breaks, labels=FALSE)

s100=fread(opt$one, header=FALSE)
s200=fread(opt$two, header=FALSE)
s300=fread(opt$three, header=FALSE)
s400=fread(opt$four, header=FALSE)
s500=fread(opt$five, header=FALSE)
s1000=fread(opt$six, header=FALSE)
s2000=fread(opt$seven, header=FALSE)

colnames(s100)=c("CHR","POS","REF","ALT","maf100")
colnames(s200)=c("CHR","POS","REF","ALT","maf200")
colnames(s300)=c("CHR","POS","REF","ALT","maf300")
colnames(s400)=c("CHR","POS","REF","ALT","maf400")
colnames(s500)=c("CHR","POS","REF","ALT","maf500")
colnames(s1000)=c("CHR","POS","REF","ALT","maf1000")
colnames(s2000)=c("CHR","POS","REF","ALT","maf2000")

s100$maf100=as.numeric(s100$maf100)
s200$maf200=as.numeric(s200$maf200)
s300$maf300=as.numeric(s300$maf300)
s400$maf400=as.numeric(s400$maf400)
s500$maf500=as.numeric(s500$maf500)
s1000$maf1000=as.numeric(s1000$maf1000)
s2000$maf2000=as.numeric(s2000$maf2000)

s100$maf100[which(s100$maf100>0.5)]<-(1-s100$maf100[which(s100$maf100>0.5)])
s200$maf200[which(s200$maf1200>0.5)]<-(1-s200$maf200[which(s200$maf200>0.5)])
s300$maf300[which(s300$maf300>0.5)]<-(1-s300$maf300[which(s300$maf300>0.5)])
s400$maf400[which(s400$maf400>0.5)]<-(1-s400$maf400[which(s400$maf400>0.5)])
s500$maf500[which(s500$maf500>0.5)]<-(1-s500$maf500[which(s500$maf500>0.5)])
s1000$maf1000[which(s1000$maf1000>0.5)]<-(1-s1000$maf1000[which(s1000$maf1000>0.5)])
s2000$maf2000[which(s2000$maf2000>0.5)]<-(1-s2000$maf2000[which(s2000$maf2000>0.5)])

results=merge(results, s100, by=c())
results=merge(results, s200, by=c())
results=merge(results, s300, by=c())
results=merge(results, s400, by=c())
results=merge(results, s500, by=c())
results=merge(results, s1000, by=c())
results=merge(results, s2000, by=c())

print(colnames(results))

results$n100[which(results$maf100==0)]=NA
results$n200[which(results$maf200==0)]=NA
results$n300[which(results$maf300==0)]=NA
results$n400[which(results$maf400==0)]=NA
results$n500[which(results$maf500==0)]=NA
results$n1000[which(results$maf1000==0)]=NA
results$n2000[which(results$maf2000==0)]=NA

#Aggregate
resultsag=aggregate(results[,c("n100", "n200", "n300", "n400", "n500", "n1000", "n2000")],by=list(results$bin), FUN=percent)
resultsag$MAF=aggregate(results[,c("MAF")], by=list(results$bin), FUN=mean, na.rm=TRUE)$MAF
write.table(resultsag[,-1], "metsim.all_perc8.tab", sep="\t", quote=FALSE, row.names=FALSE)

