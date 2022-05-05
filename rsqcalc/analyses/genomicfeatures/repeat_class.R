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

results$DNA=rep(0,length(results$POS))
results$SINE=rep(0,length(results$POS))
results$LINE=rep(0,length(results$POS))
results$LTR=rep(0,length(results$POS))
results$Simple=rep(0,length(results$POS))
results$LowComplex=rep(0,length(results$POS))
results$Satellite=rep(0,length(results$POS))
results$RNA=rep(0,length(results$POS))
results$RC=rep(0,length(results$POS))
results$Unknown=rep(0,length(results$POS))


for (i in 1:length(results$POS)){
	reps = repeats[which((repeats$START>results$START[i] & repeats$START<results$END[i])|(
				repeats$END>results$START[i] & repeats$END<results$END[i])),]
	results$DNA[i]=length(which(reps$CLASS=="DNA" | reps$CLASS=="DNA?"))
	results$SINE[i]=length(which(reps$CLASS=="SINE"))
	results$LINE[i]=length(which(reps$CLASS=="LINE"))
	results$LTR[i]=length(which(reps$CLASS=="LTR" | reps$CLASS=="LTR?" | reps$CLASS=="Retroposon"))
	results$Simple[i]=length(which(reps$CLASS=="Simple_repeat"))
	results$LowComplex[i]=length(which(reps$CLASS=="Low_complexity"))
	results$Satellite[i]=length(which(reps$CLASS=="Satellite"))
	results$RNA[i]=length(which(reps$CLASS=="RNA" | reps$CLASS=="rRNA" | reps$CLASS=="scRNA" | reps$CLASS=="snRNA" | reps$CLASS=="srpRNA" | reps$CLASS=="tRNA"))
	results$RC[i]=length(which(reps$CLASS=="RC" | reps$CLASS=="RC?"))
	results$Unknown[i]=length(which(reps$CLASS=="Unknown"))
	if (i %% 1000==0) {
		Sys.sleep(0.01)
		print(round(i*100/nrow(results), digits=2))
	}
}

write.table(results, opt$output, quote=FALSE, sep="\t", row.names=FALSE)
