#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(data.table))

optionList <- list(
  make_option(c("-r", "--results"), type="character", help="Merged result file with repeat column"),
  make_option(c("-c", "--recombination"), type="character", help="Recombination file"),
  make_option(c("-w", "--window"), type="character", help="Window size"),
  make_option(c("-m", "--chromosome"), type="character", help="Chromosome (numeric)"),
  make_option(c("-o", "--output"), type="character", help="Output file")
)

parser <- OptionParser(
  usage="%prog -r results -c conservation -o output",
  option_list=optionList
)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

print("Beginning results file read")
results = fread(opt$results)
print("Results file read")
chr=opt$chromosome
colnames(results)=c("CHR","POS","REF","ALT","AF")
results=results[which(results$CHR==chr),]
results$MAF=as.numeric(results$AF)
results$MAF[which(results$MAF>0.5)]<-(1-results$MAF[which(results$MAF>0.5)])
print("MAF calculated")
window=as.numeric(opt$window)
results$START=results$POS-window
results$END=results$POS+window

recomb = fread(opt$recombination, fill=TRUE)
print("Recombination file read")
colnames(recomb)=c("CHR","POS","RECOMB_POS","POS_CM")
recomb$POS=as.numeric(recomb$POS)
recomb$RECOMB_POS=as.numeric(recomb$RECOMB_POS)
results=merge(results, recomb[,c("POS","RECOMB_POS")], by="POS", all.x=TRUE)
results$RECOMB_REG=rep(0, length(results$POS))
for (i in 1:length(results$POS)){
	rec = recomb[which(recomb$POS > results$START[i] & recomb$POS < results$END[i]),]
	results$RECOMB_REG[i]=mean(rec$RECOMB_POS, na.rm=TRUE)
	if (i %% 1000==0) {
		Sys.sleep(0.01)
		print(round(i*100/nrow(results), digits=2))
	}
}

write.table(results, opt$output, quote=FALSE, sep="\t", row.names=FALSE)
