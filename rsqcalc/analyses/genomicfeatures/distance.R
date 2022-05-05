#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

optionList <- list(
  make_option(c("-r", "--results"), type="character", help="Merged result file"),
  make_option(c("-a", "--array"), type="character", help="Array site list"),
  make_option(c("-o", "--output"), type="character", help="Output file"),
  make_option(c("-c", "--chromosome"), type="character", help="Chromosome (numeric)")
)

parser <- OptionParser(
  usage="%prog -r results -a array -w window -o output -c chromosome",
  option_list=optionList
)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

print("Beginning sites file read")
results = fread(opt$results)
print("Results file read")
chr=opt$chromosome
colnames(results)=c("CHR","POS","REF","ALT","AF")
results=results[which(results$CHR==chr),]
results$MAF=as.numeric(results$AF)
results$MAF[which(results$MAF>0.5)]<-(1-results$MAF[which(results$MAF>0.5)])
print("MAF calculated")

arr = fread(opt$array)
print("Array file read")
colnames(arr)=c("CHR", "POS1", "POS")

results$DIST=rep(0,length(results$DIST))
for (i in 1:length(results$POS)){
	results$DIST[i]=min(abs(results$POS[i]-arr$POS))
	if (i %% 1000==0) {
		Sys.sleep(0.01)
		print(round(i*100/nrow(results), digits=2))
	}
}

write.table(results, opt$output, quote=FALSE, sep="\t", row.names=FALSE)
