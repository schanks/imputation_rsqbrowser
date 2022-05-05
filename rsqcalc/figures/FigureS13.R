library(data.table)
library(ggplot2)

inpsyght=fread("../inpsyght.bi.snv.tab")
biome=fread("../biome.bi.snv.tab")
mlof=fread("../mlof.bi.snv.tab")
metsim=fread("../metsim.bi.snv.tab")
inpsyght=inpsyght[,c("CHR","POS","REF","ALT","AF")]
biome=biome[,c("CHR","POS","REF","ALT","AF")]
mlof=mlof[,c("CHR","POS","REF","ALT","AF")]
metsim=metsim[,c("CHR","POS","REF","ALT","AF")]

inpsyght$MAF=as.numeric(inpsyght$AF)
inpsyght$MAF[which(inpsyght$MAF>0.5)]=1-inpsyght$MAF[which(inpsyght$MAF>0.5)]
biome$MAF=as.numeric(biome$AF)
biome$MAF[which(biome$MAF>0.5)]=1-biome$MAF[which(biome$MAF>0.5)]
mlof$MAF=as.numeric(mlof$AF)
mlof$MAF[which(mlof$MAF>0.5)]=1-mlof$MAF[which(mlof$MAF>0.5)]
metsim$MAF=as.numeric(metsim$AF)
metsim$MAF[which(metsim$MAF>0.5)]=1-metsim$MAF[which(metsim$MAF>0.5)]

inpsyght$Ancestry=rep("African", nrow(inpsyght))
biome$Ancestry=rep("Hispanic/Latino", nrow(biome))
mlof$Ancestry=rep("European", nrow(mlof))
metsim$Ancestry=rep("Finnish", nrow(metsim))

all=rbind(inpsyght, biome, mlof, metsim)
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))
all=all[which(all$MAF>0),]

breaks=c(0,0.0005,0.005,0.05,0.5)
all$bin=cut(all$MAF, breaks=breaks, labels=FALSE)

allag=aggregate(all$MAF, by=list(all$bin, all$Ancestry), FUN=length)
allag$`WGS MAF`=rep("Very rare: <0.005%", nrow(allag))
allag$`WGS MAF`[which(allag$Group.1==2)]="Rare: 0.05%-0.5%"
allag$`WGS MAF`[which(allag$Group.1==3)]="Low frequency: 0.5%-5%"
allag$`WGS MAF`[which(allag$Group.1==4)]="Common: >5%"
allag$`WGS MAF`=factor(allag$`WGS MAF`, levels=c("Very rare: <0.005%","Rare: 0.05%-0.5%","Low frequency: 0.5%-5%","Common: >5%"))
colnames(allag)[2]="Ancestry"
colnames(allag)[3]="Count"
afrsum=sum(allag$Count[which(allag$Ancestry=="African")])
hlsum=sum(allag$Count[which(allag$Ancestry=="Hispanic/Latino")])
mlsum=sum(allag$Count[which(allag$Ancestry=="European")])
metsum=sum(allag$Count[which(allag$Ancestry=="Finnish")])
allag$Proportion=allag$Count/afrsum
allag$Proportion[which(allag$Ancestry=="Hispanic/Latino")]=allag$Count[which(allag$Ancestry=="Hispanic/Latino")]/hlsum
allag$Proportion[which(allag$Ancestry=="European")]=allag$Count[which(allag$Ancestry=="European")]/mlsum
allag$Proportion[which(allag$Ancestry=="Finnish")]=allag$Count[which(allag$Ancestry=="Finnish")]/metsum

colnames(allag)[3]="Count of SNVs in WGS dataset"
colnames(allag)[5]="Proportion of SNVs in WGS dataset"

allmelt=melt(allag[,-1], id.vars=c("Ancestry","WGS MAF"))

tiff(file="~/FigS13.tiff", height=114, width=174, units="mm", res=300)
ggplot(allmelt, aes(x=Ancestry, fill=Ancestry, alpha=`WGS MAF`, y=value))+geom_bar(stat="identity", position="dodge")+facet_wrap(variable~., scales="free",strip.position="left",nrow=2)+theme_bw()+ylab("Biallelic SNVs")+theme(strip.text=element_text(size=8), axis.title=element_text(size=8), legend.title=element_text(size=8), legend.text=element_text(size=7), axis.text=element_text(size=7),strip.background=element_rect(color=NA, fill="white"), strip.placement="outside")+scale_fill_manual(values=c("#cb181d","#238b45","#2171b5","#6a51a3"))+scale_alpha_manual(values=c(1.0,0.8,0.6,0.4))
dev.off()

