library(data.table)
library(ggplot2)
library(reshape2)
library(gridExtra)

btc=fread("../analyses/runs/biome.all.common.stretch")
colnames(btc)=c("Start","End","length","nvar","pass")
btc=btc[which(btc$pass==1),]
btc$MAF=rep("Common", nrow(btc))
btc$nvar=as.numeric(btc$nvar)
btl=fread("../analyses/runs/biome.all.lowfreq.stretch")
colnames(btl)=c("Start","End","length","nvar","pass")
btl=btl[which(btl$pass==1),]
btl$nvar=as.numeric(btl$nvar)
btl$MAF=rep("Lowfreq",nrow(btl))
btr=fread("../analyses/runs/biome.all.rare.stretch")
colnames(btr)=c("Start","End","length","nvar","pass")
btr=btr[which(btr$pass==1),]
btr$MAF=rep("Rare", nrow(btr))
btr$nvar=as.numeric(btr$nvar)
btop=rbind(btc, btl, btr)
btop$`Reference Panel`=rep("TOPMed", nrow(btop))
btop$Ancestry=rep("Hispanic/Latino", nrow(btop))

itc=fread("../analyses/runs/inpsyght.all.common.stretch")
colnames(itc)=c("Start","End","length","nvar","pass")
itc=itc[which(itc$pass==1),]
itc$MAF=rep("Common", nrow(itc))
itc$nvar=as.numeric(itc$nvar)
itl=fread("../analyses/runs/inpsyght.all.lowfreq.stretch")
colnames(itl)=c("Start","End","length","nvar","pass")
itl=itl[which(itl$pass==1),]
itl$MAF=rep("Lowfreq",nrow(itl))
itl$nvar=as.numeric(itl$nvar)
itr=fread("../analyses/runs/inpsyght.all.rare.stretch")
colnames(itr)=c("Start","End","length","nvar","pass")
itr=itr[which(itr$pass==1),]
itr$MAF=rep("Rare", nrow(itr))
itr$nvar=as.numeric(itr$nvar)
itop=rbind(itc, itl, itr)
itop$`Reference Panel`=rep("TOPMed", nrow(itop))
itop$Ancestry=rep("African", nrow(itop))

mtc=fread("../analyses/runs/mlof.all.common.stretch")
colnames(mtc)=c("Start","End","length","nvar","pass")
mtc=mtc[which(mtc$pass==1),]
mtc$MAF=rep("Common", nrow(mtc))
mtc$nvar=as.numeric(mtc$nvar)
mtl=fread("../analyses/runs/mlof.all.lowfreq.stretch")
colnames(mtl)=c("Start","End","length","nvar","pass")
mtl=mtl[which(mtl$pass==1),]
mtl$MAF=rep("Lowfreq",nrow(mtl))
mtl$nvar=as.numeric(mtl$nvar)
mtr=fread("../analyses/runs/mlof.all.rare.stretch")
colnames(mtr)=c("Start","End","length","nvar","pass")
mtr=mtr[which(mtr$pass==1),]
mtr$MAF=rep("Rare", nrow(mtr))
mtr$nvar=as.numeric(mtr$nvar)
mtop=rbind(mtc, mtl, mtr)
mtop$`Reference Panel`=rep("TOPMed", nrow(mtop))
mtop$Ancestry=rep("European", nrow(mtop))

mttc=fread("../analyses/runs/metsim.all.common.stretch")
colnames(mttc)=c("Start","End","length","nvar","pass")
mttc=mttc[which(mttc$pass==1),]
mttc$MAF=rep("Common", nrow(mttc))
mttc$nvar=as.numeric(mttc$nvar)
mttl=fread("../analyses/runs/metsim.all.lowfreq.stretch")
colnames(mttl)=c("Start","End","length","nvar","pass")
mttl=mttl[which(mttl$pass==1),]
mttl$MAF=rep("Lowfreq",nrow(mttl))
mttl$nvar=as.numeric(mttl$nvar)
mttr=fread("../analyses/runs/metsim.all.rare.stretch")
colnames(mttr)=c("Start","End","length","nvar","pass")
mttr=mttr[which(mttr$pass==1),]
mttr$MAF=rep("Rare", nrow(mttr))
mttr$nvar=as.numeric(mttr$nvar)
mttop=rbind(mttc, mttl, mttr)
mttop$`Reference Panel`=rep("TOPMed", nrow(mttop))
mttop$Ancestry=rep("Finnish", nrow(mttop))

all=rbind(btop, mtop, itop, mttop)
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))
all$MAF[which(all$MAF=="Lowfreq")]="Low frequency"
a=ggplot(all, aes(x=nvar,group=Ancestry,color=Ancestry,linetype=MAF))+stat_ecdf(size=1)+facet_wrap(MAF~., nrow=3,scales="free")+scale_x_continuous(expand=c(0,0))+scale_color_manual(values=c("#e41a1c","#4daf4a","#377eb8","#984ea3"))+scale_linetype_manual(values=c("solid","longdash","dotted"))+xlab("Number of consecutive well-imputed biallelic SNVs")+ylab("Cumulative distribution")+theme_bw()+theme(legend.key.width=unit(1.2, "cm"),strip.background=element_rect(fill="white",color=NA),strip.text=element_text(size=8),axis.text=element_text(size=7), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8))


inpsyght=fread("../inpsyght.bi.snv.tab")
inpsyght=inpsyght[which(inpsyght$CHR=="chr20"),c("POS","AF","Omni25_TOPMed","Omni25_TOPMed_in")]
inpsyght[is.na(inpsyght)]=0
inpsyght$MAF=as.numeric(inpsyght$AF)
inpsyght$MAF[which(inpsyght$MAF>0.5)]=1-inpsyght$MAF[which(inpsyght$MAF>0.5)]
inpsyght$POS=inpsyght$POS/1000000
inpsyght=inpsyght[which(inpsyght$MAF>0.05),]
inpsyght$Ancestry=rep("African", nrow(inpsyght))

biome=fread("../biome.bi.snv.tab")
biome=biome[which(biome$CHR=="chr20"),c("POS","AF","Omni25_TOPMed")]
biome[is.na(biome)]=0
biome$MAF=as.numeric(biome$AF)
biome$MAF[which(biome$MAF>0.5)]=1-biome$MAF[which(biome$MAF>0.5)]
biome$POS=biome$POS/1000000
biome=biome[which(biome$MAF>0.05),]
biome$Ancestry=rep("Hispanic/Latino", nrow(biome))

mlof=fread("../mlof.bi.snv.tab")
mlof=mlof[which(mlof$CHR=="chr20"),c("POS","AF","Omni25_TOPMed")]
mlof[is.na(mlof)]=0
mlof$MAF=as.numeric(mlof$AF)
mlof$MAF[which(mlof$MAF>0.5)]=1-mlof$MAF[which(mlof$MAF>0.5)]
mlof$POS=mlof$POS/1000000
mlof=mlof[which(mlof$MAF>0.05),]
mlof$Ancestry=rep("European", nrow(mlof))

metsim=fread("../metsim.bi.snv.tab")
metsim=metsim[which(metsim$CHR=="chr20"),c("POS","AF","Omni25_TOPMed")]
metsim[is.na(metsim)]=0
metsim$MAF=as.numeric(metsim$AF)
metsim$MAF[which(metsim$MAF>0.5)]=1-metsim$MAF[which(metsim$MAF>0.5)]
metsim$POS=metsim$POS/1000000
metsim=metsim[which(metsim$MAF>0.05),]
metsim$Ancestry=rep("Finnish", nrow(metsim))

example=rbind(inpsyght, biome, mlof, metsim)
example$Ancestry=factor(example$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))

b=ggplot(example, aes(x=POS, y=Omni25_TOPMed))+xlab("Chromosome 20 Position (Mb)")+ylab(expression(Observed~imputation~r^2))+geom_point(aes(col=Ancestry), size=0.2)+geom_hline(yintercept=0.8)+theme_bw()+facet_grid(Ancestry~.)+theme(legend.position="none",strip.background=element_rect(color=NA, fill="white"), strip.text=element_text(size=8), axis.title=element_text(size=8),axis.text=element_text(size=7))+scale_color_manual(values=c("#e41a1c","#4daf4a","#377eb8","#984ea3"))


tiff(file="~/Fig3.tiff", height=100, width=174, units="mm", res=300)
ab=grid.arrange(b,a, nrow=1, widths=c(0.8,1))
dev.off()




