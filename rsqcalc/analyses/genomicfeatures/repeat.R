#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

optionList <- list(
  make_option(c("-r", "--results"), type="character", help="Merged result file"),
  make_option(c("-p", "--repeats"), type="character", help="Repeats file"),
  make_option(c("-w", "--window"), type="character", help="Window size"),
  make_option(c("-o", "--output"), type="character", help="Output file"),
  make_option(c("-c", "--chromosome"), type="character", help="Chromosome (numeric)")
)

parser <- OptionParser(
  usage="%prog -r results -p repeats -w window -o output -c chromosome",
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
window=as.numeric(opt$window)
results$START=results$POS-window
results$END=results$POS+window


repeats = fread(opt$repeats)
print("Repeat file read")
colnames(repeats)=c("CHR", "START", "END", "CLASS")

results$REP=rep(0,length(results$POS))
for (i in 1:length(results$POS)){
	reps = repeats[which((repeats$START>results$START[i] & repeats$START<results$END[i])|(
				repeats$END>results$START[i] & repeats$END<results$END[i])),]
	results$REP[i]=nrow(reps)
	if (i %% 1000==0) {
		Sys.sleep(0.01)
		print(round(i*100/nrow(results), digits=2))
	}
}

write.table(results, opt$output, quote=FALSE, sep="\t", row.names=FALSE)
