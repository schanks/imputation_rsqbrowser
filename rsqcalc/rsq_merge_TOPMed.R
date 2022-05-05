#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

optionList <- list(
  make_option(c("-c", "--CoreTOP"), type="character", help="Core TOPMed file"),
  make_option(c("-g", "--ExpressTOP"), type="character", help="Express TOP file"),
  make_option(c("-l", "--OmniTOP"), type="character", help="Omni TOPMed file"),
  make_option(c("-p", "--GlobalTOP"), type="character", help="MEGA TOPMed file"),
  make_option(c("-z", "--out"), type="character", help="Output file name"),
  make_option(c("-n", "--wgs"), type="character", help="WGS variant file")
)

parser <- OptionParser(
  usage="%prog -c CoreTOP -g ExpressTOP -l OmniTOP -p GlobalTOP -n wgs -z out",
  option_list=optionList
)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

wgs = fread(opt$wgs, header = FALSE, stringsAsFactors = FALSE, fill=TRUE)
colnames(wgs)=c("CHR", "POS", "REF", "ALT", "AF")
wgs = wgs[which(wgs$AF>0 & wgs$AF<1),]
wgs = wgs[which(wgs$REF!="REF"),]
wgs$CHR=as.factor(wgs$CHR)

#Core TOP
results=read.delim(opt$CoreTOP, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="Core_TOPMed"
merged=merge(wgs, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Express TOP
results=read.delim(opt$ExpressTOP, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="OmniExp_TOPMed"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Omni TOP
results=read.delim(opt$OmniTOP, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="Omni25_TOPMed"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Global TOP
results=read.delim(opt$GlobalTOP, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="MEGA_TOPMed"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

outfile=opt$out
write.table(merged, paste(outfile, ".tab", sep=""), quote = FALSE, sep="\t", row.names = FALSE)


