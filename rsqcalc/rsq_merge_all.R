#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

optionList <- list(
  make_option(c("-a", "--Core1000G"), type="character", help="Core 1000G file"),
  make_option(c("-b", "--CoreHRC"), type="character", help="Core HRC file"),
  make_option(c("-c", "--CoreTOP"), type="character", help="Core TOPMed file"),
  make_option(c("-e", "--Express1000G"), type="character", help="Express 1000G file"),
  make_option(c("-f", "--ExpressHRC"), type="character", help="Express HRC file"),
  make_option(c("-g", "--ExpressTOP"), type="character", help="Express TOP file"),
  make_option(c("-j", "--Omni1000G"), type="character", help="Omni 1000G file"),
  make_option(c("-k", "--OmniHRC"), type="character", help="Omni HRC file"),
  make_option(c("-l", "--OmniTOP"), type="character", help="Omni TOPMed file"),
  make_option(c("-m", "--Global1000G"), type="character", help="MEGA 1000G file"),
  make_option(c("-o", "--GlobalHRC"), type="character", help="MEGA HRC file"),
  make_option(c("-p", "--GlobalTOP"), type="character", help="MEGA TOPMed file"),
  make_option(c("-n", "--wgs"), type="character", help="WGS variant file"),
  make_option(c("-z", "--fname"), type="character", help="Output file name")
)

parser <- OptionParser(
  usage="%prog -a Core1000G -b CoreHRC -c CoreTOP -e Express1000G -f ExpressHRC -g ExpressTOP -j Omni1000G -k OmniHRC -l OmniTOP -m Global1000G -o GlobalHRC -p GlobalTOP -n wgs -z fname",
  option_list=optionList
)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

wgs = fread(opt$wgs, header = FALSE, stringsAsFactors = FALSE, fill=TRUE)
colnames(wgs)=c("CHR", "POS", "REF", "ALT", "AF")
wgs = wgs[which(wgs$AF>0 & wgs$AF<1),]
wgs = wgs[which(wgs$REF!="REF"),]
wgs$POS=as.numeric(wgs$POS)

#Core 1000G
results=fread(opt$Core1000G, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="Core_1000G"
merged=merge(wgs, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Core HRC
results=fread(opt$CoreHRC, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="Core_HRC"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Core TOP
results=fread(opt$CoreTOP, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="Core_TOPMed"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Express 1000G
results=fread(opt$Express1000G, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="OmniExp_1000G"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Express HRC
results=fread(opt$ExpressHRC, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="OmniExp_HRC"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Core TOP
results=fread(opt$ExpressTOP, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="OmniExp_TOPMed"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Omni 1000G
results=fread(opt$Omni1000G, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="Omni25_1000G"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Omni HRC
results=fread(opt$OmniHRC, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="Omni25_HRC"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Omni TOP
results=fread(opt$OmniTOP, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="Omni25_TOPMed"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Global 1000G
results=fread(opt$Global1000G, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="MEGA_1000G"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Global HRC
results=fread(opt$GlobalHRC, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="MEGA_HRC"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

#Global TOP
results=fread(opt$GlobalTOP, header = TRUE, stringsAsFactors = FALSE)
colnames(results)[9]="MEGA_TOPMed"
merged=merge(merged, results[,c(1,2,3,4,9)], by=c("CHR", "POS", "REF", "ALT"), all.x=TRUE)

outfile=opt$fname
write.table(merged, paste(outfile, ".tab", sep=""), quote = FALSE, sep="\t", row.names = FALSE)


