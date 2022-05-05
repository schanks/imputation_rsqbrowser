#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

optionList <- list(
  make_option(c("-r", "--results"), type="character", help="Merged result file with repeat column"),
  make_option(c("-s", "--segdup"), type="character", help="Segmental duplications file"),
  make_option(c("-w", "--window"), type="character", help="Window size"),
  make_option(c("-c", "--chromosome"), type="character", help="Chromosome (numeric)"),
  make_option(c("-o", "--output"), type="character", help="Output file")
)

parser <- OptionParser(
  usage="%prog -r results -s segdup -o output",
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

seg = fread(opt$segdup, fill=TRUE)
print("Segmental duplications file read")
colnames(seg)=c("CHR","START","END")
seg$START=as.numeric(seg$START)
seg$END=as.numeric(seg$END)

results$SEG=rep(0,length(results$POS))
for (i in 1:length(results$POS)){
	segs = seg[which((seg$START>results$START[i] & seg$START<results$END[i]) | seg$END> results$START[i] & seg$END<results$END[i]),]
	results$SEG[i]=min(1, nrow(segs))
	if (i %% 1000==0) {
		Sys.sleep(0.01)
		print(round(i*100/nrow(results), digits=2))
	}
}

write.table(results, opt$output, quote=FALSE, sep="\t", row.names=FALSE)
