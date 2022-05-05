#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
optionList <- list(
  make_option(c("-i", "--imputed"), type="character", help="imputed file"),
  make_option(c("-w", "--wgs"), type="character", help="wgs file"),
  make_option(c("-m", "--maf"), type="numeric", help="max maf between 0 and 0.5"),
  make_option(c("-n", "--min"), type="numeric", help="min maf between 0 and 0.5"),
  make_option(c("-t", "--threshold"), type="numeric", help="r2 threshold between 0 and 1"),
  make_option(c("-c", "--centromere"), type="character", help="centromere file"),
  make_option(c("-o", "--out"), type="character", help="output filename")
)

parser <- OptionParser(
  usage="%prog -i imputed -w wgs -m maf -t threshold -c centromere -o out",
  option_list=optionList
)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#Read in wgs
wgs = fread(opt$wgs)
colnames(wgs)=c("CHR", "POS", "REF", "ALT", "AF")
print("WGS file read")
chr=as.character(wgs[1,"CHR"])
wgs$POS=as.numeric(wgs$POS)

results=fread(opt$imputed)
results$POS=as.numeric(results$POS)
results$Panel=rep(1, nrow(results))
print("Imputed file read")
results=merge(wgs[,-1], results[,-7], by=c("POS","REF","ALT"), all.x=TRUE)
results$Rsq=as.numeric(results$Rsq)
results$Rsq[which(is.na(results$Rsq))]=0
results$MAF=as.numeric(results$AF)
results$MAF[which(results$MAF>0.5)]<-(1-results$MAF[which(results$MAF>0.5)])
results=results[which(results$MAF>0),]
results = results[order(results$POS),]
print("Results merged")

maf=opt$maf
min=opt$min
rsq=opt$threshold
results = results[which(results$MAF<maf & results$MAF>min),]

centro=fread(opt$centromere, header=TRUE)
centro=centro[which(centro$chrom==chr),]
centstart=min(centro$chromStart)
centend=max(centro$chromEnd)

arm1=results[which(results$POS<centstart)]
arm2=results[which(results$POS>centend)]

switches = as.data.frame(matrix(c(0,0,0,0,0), nrow=1, ncol=5))
colnames(switches)=c("start", "end", "length", "nvar", "pass")

results=arm1
switchrow=1
i=1
while(i<length(results$POS)+1){
  start = as.numeric(results$POS[i])
  pass = (results$Rsq[i] >= rsq)
  nvar = 0
  if (pass){
    while(results$Rsq[i] >= rsq) {
      nvar = nvar + 1
      i = i+1
      if (i>length(results$POS)) {break}
    }
  }
  else{
    while(results$Rsq[i] < rsq){
      nvar = nvar + 1
      i = i+1
      if (i>length(results$POS)) {break}
    }
  }
  end = as.numeric(results$POS[i-1])
  switches[switchrow,]=c(start, end, end-start+1, nvar, pass)
  switchrow = switchrow + 1
}
arm1=switches
print("Arm 1 complete")
results=arm2
switchrow=1
i=1
while(i<length(results$POS)+1){
  start = as.numeric(results$POS[i])
  pass = (results$Rsq[i] >= rsq)
  nvar = 0
  if (pass){
    while(results$Rsq[i] >= rsq) {
      nvar = nvar + 1
      i = i+1
      if (i>length(results$POS)) {break}
    }
  }
  else{
    while(results$Rsq[i] < rsq){
      nvar = nvar + 1
      i = i+1
      if (i>length(results$POS)) {break}
    }
  }
  end = as.numeric(results$POS[i-1])
  switches[switchrow,]=c(start, end, end-start+1, nvar, pass)
  switchrow = switchrow + 1
}
arm2=switches

switches=rbind(arm1,arm2)

out=opt$out
write.table(switches, out, sep="\t", quote = FALSE, row.names = FALSE)

