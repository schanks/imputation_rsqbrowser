library(data.table)
library(ggplot2)

inpsyght=fread("../inpsyght.bi.snv.tab")
inpsyght[is.na(inpsyght)]=0
inpsyght$MAF=as.numeric(inpsyght$AF)
inpsyght$MAF[which(inpsyght$MAF>0.5)]=1-inpsyght$MAF[which(inpsyght$MAF>0.5)]
inpsyght$POS=inpsyght$POS/1000000
inpsyght=inpsyght[which(inpsyght$MAF>0.05),]
inpsyght$Ancestry=rep("African", nrow(inpsyght))

biome=fread("../biome.bi.snv.tab")
biome[is.na(biome)]=0
biome$MAF=as.numeric(biome$AF)
biome$MAF[which(biome$MAF>0.5)]=1-biome$MAF[which(biome$MAF>0.5)]
biome$POS=biome$POS/1000000
biome=biome[which(biome$MAF>0.05),]
biome$Ancestry=rep("Hispanic/Latino", nrow(biome))

mlof=fread("../mlof.bi.snv.tab")
mlof[is.na(mlof)]=0
mlof$MAF=as.numeric(mlof$AF)
mlof$MAF[which(mlof$MAF>0.5)]=1-mlof$MAF[which(mlof$MAF>0.5)]
mlof$POS=mlof$POS/1000000
mlof=mlof[which(mlof$MAF>0.05),]
mlof$Ancestry=rep("European", nrow(mlof))

metsim=fread("../metsim.bi.snv.tab")
metsim[is.na(metsim)]=0
metsim$MAF=as.numeric(metsim$AF)
metsim$MAF[which(metsim$MAF>0.5)]=1-metsim$MAF[which(metsim$MAF>0.5)]
metsim$POS=metsim$POS/1000000
metsim=metsim[which(metsim$MAF>0.05),]
metsim$Ancestry=rep("Finnish", nrow(metsim))

example=rbind(inpsyght, biome, mlof, metsim)
example=example[,c("CHR","POS","Core_TOPMed","OmniExp_TOPMed","MEGA_TOPMed","Omni25_TOPMed","Ancestry")]
example=melt(example, id.vars=c("CHR","POS","Ancestry"))
example$Array=rep("Core", nrow(example))
example$Array[which(example$variable=="OmniExp_TOPMed")]="Omni Express"
example$Array[which(example$variable=="MEGA_TOPMed")]="MEGA"
example$Array[which(example$variable=="Omni25_TOPMed")]="Omni 2.5M"
example$Array=factor(example$Array, levels=c("Omni 2.5M","MEGA","Omni Express","Core"))
example$Ancestry=factor(example$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))

for (i in 1:22){
	chr=paste("chr",i,sep="")
	filename=paste("~/FigS9_",i,".tiff",sep="")
	tiff(file=filename, height=85, width=174, units="mm", res=300)
	ggplot(example[which(example$CHR==chr)], aes(x=POS, y=value))+ggtitle(paste("Chromosome",i,sep=" "))+xlab("Position (Mb)")+ylab(expression(Observed~imputation~r^2))+geom_point(aes(col=Ancestry),size=0.2)+geom_hline(yintercept=0.8)+theme_bw()+facet_grid(Ancestry~Array)+theme(legend.position="none",legend.text=element_text(size=7), legend.title=element_text(size=8),strip.background=element_rect(color=NA, fill="white"), strip.text=element_text(size=6), plot.title=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=7))+scale_color_manual(values=c("#e41a1c","#4daf4a","#377eb8","#984ea3"))
	dev.off()
}







