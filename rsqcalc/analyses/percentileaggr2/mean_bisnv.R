#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

optionList <- list(
  make_option(c("-r", "--results"), type="character", help="Merged result file"),
  make_option(c("-o", "--out"), type="character", help="Output filename")
)

parser <- OptionParser(
  usage="%prog -r results -o out",
  option_list=optionList
)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

print("Beginning results file read")
results = fread(opt$results, stringsAsFactors = FALSE)
print("Results file read")
print("Results subset")
results[is.na(results)]=0
print("NAs set to 0")
results$MAF=as.numeric(results$AF)
results$MAF[which(results$MAF>0.5)]<-(1-results$MAF[which(results$MAF>0.5)])
print("MAF calculated")
results=results[which(results$MAF>0),]
breaks=c(seq(0, 0.00175, 0.00025), seq(0.002, 0.50014, 0.001))
results$bin=cut(results$MAF, breaks=breaks, labels=FALSE)
print("Results file formatted")

print(colnames(results))
#Aggregate
resultsag=aggregate(results[,c("Core_1000G", "Core_HRC", "Core_TOPMed", "OmniExp_1000G","OmniExp_HRC", "OmniExp_TOPMed", "Omni25_1000G", "Omni25_HRC", "Omni25_TOPMed","MEGA_1000G", "MEGA_HRC", "MEGA_TOPMed")],by=list(results$bin), FUN=mean)
print("Percentiles calculated")
resultsag$MAF=aggregate(results[,c("MAF")], by=list(results$bin), FUN=mean, na.rm=TRUE)$MAF
print("Rsq aggregated")
outfile=opt$out
write.table(resultsag[,-1], outfile, sep="\t", quote=FALSE, row.names=FALSE)
