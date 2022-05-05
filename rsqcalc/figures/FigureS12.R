library(data.table)
library(ggplot2)
library(gridExtra)
library(reshape2)

biome=fread("../analyses/percentileaggr2/biome.bi.snv.top.perc8.tab")
biome=biome[c(which(abs(biome$MAF-0.0002)==min(abs(biome$MAF-0.0002))),which(abs(biome$MAF-0.0004)==min(abs(biome$MAF-0.0004))),which(abs(biome$MAF-0.0006)==min(abs(biome$MAF-0.0006))),which(abs(biome$MAF-0.001)==min(abs(biome$MAF-0.001))),which(abs(biome$MAF-0.0032)==min(abs(biome$MAF-0.0032))),which(abs(biome$MAF-0.01)==min(abs(biome$MAF-0.01))),which(abs(biome$MAF-0.0316)==min(abs(biome$MAF-0.0316))),which(abs(biome$MAF-0.1)==min(abs(biome$MAF-0.1))),which(abs(biome$MAF-0.5)==min(abs(biome$MAF-0.5)))),]
biome=melt(biome, id.vars="MAF")
biome$Type=rep("Biallelic SNVs", nrow(biome))
bbsnv=biome
biome=fread("../analyses/percentileaggr2/biome.bi.indel.top.perc8.tab")
biome=biome[c(which(abs(biome$MAF-0.0002)==min(abs(biome$MAF-0.0002))),which(abs(biome$MAF-0.0004)==min(abs(biome$MAF-0.0004))),which(abs(biome$MAF-0.0006)==min(abs(biome$MAF-0.0006))),which(abs(biome$MAF-0.001)==min(abs(biome$MAF-0.001))),which(abs(biome$MAF-0.0032)==min(abs(biome$MAF-0.0032))),which(abs(biome$MAF-0.01)==min(abs(biome$MAF-0.01))),which(abs(biome$MAF-0.0316)==min(abs(biome$MAF-0.0316))),which(abs(biome$MAF-0.1)==min(abs(biome$MAF-0.1))),which(abs(biome$MAF-0.5)==min(abs(biome$MAF-0.5)))),]
biome=melt(biome, id.vars="MAF")
biome$Type=rep("Biallelic indels",nrow(biome))
bbindel=biome
biome=fread("../analyses/percentileaggr2/biome.multi.snv.top.perc8.tab")
biome=biome[c(which(abs(biome$MAF-0.0002)==min(abs(biome$MAF-0.0002))),which(abs(biome$MAF-0.0004)==min(abs(biome$MAF-0.0004))),which(abs(biome$MAF-0.0006)==min(abs(biome$MAF-0.0006))),which(abs(biome$MAF-0.001)==min(abs(biome$MAF-0.001))),which(abs(biome$MAF-0.0032)==min(abs(biome$MAF-0.0032))),which(abs(biome$MAF-0.01)==min(abs(biome$MAF-0.01))),which(abs(biome$MAF-0.0316)==min(abs(biome$MAF-0.0316))),which(abs(biome$MAF-0.1)==min(abs(biome$MAF-0.1))),which(abs(biome$MAF-0.5)==min(abs(biome$MAF-0.5)))),]
biome=melt(biome, id.vars="MAF")
biome$Type=rep("Multiallelic SNVs",nrow(biome))
bmsnv=biome
biome=fread("../analyses/percentileaggr2/biome.multi.indel.top.perc8.tab")
biome=biome[c(which(abs(biome$MAF-0.0002)==min(abs(biome$MAF-0.0002))),which(abs(biome$MAF-0.0004)==min(abs(biome$MAF-0.0004))),which(abs(biome$MAF-0.0006)==min(abs(biome$MAF-0.0006))),which(abs(biome$MAF-0.001)==min(abs(biome$MAF-0.001))),which(abs(biome$MAF-0.0032)==min(abs(biome$MAF-0.0032))),which(abs(biome$MAF-0.01)==min(abs(biome$MAF-0.01))),which(abs(biome$MAF-0.0316)==min(abs(biome$MAF-0.0316))),which(abs(biome$MAF-0.1)==min(abs(biome$MAF-0.1))),which(abs(biome$MAF-0.5)==min(abs(biome$MAF-0.5)))),]
biome=melt(biome, id.vars="MAF")
biome$Type=rep("Multiallelic indels",nrow(biome))
bmindel=biome
biome=rbind(bbsnv, bbindel, bmsnv, bmindel)
biome$Ancestry=rep("Hispanic/Latino", nrow(biome))

mlof=fread("../analyses/percentileaggr2/mlof.bi.snv.top.perc8.tab")
mlof=mlof[c(which(abs(mlof$MAF-0.0002)==min(abs(mlof$MAF-0.0002))),which(abs(mlof$MAF-0.0004)==min(abs(mlof$MAF-0.0004))),which(abs(mlof$MAF-0.0006)==min(abs(mlof$MAF-0.0006))),which(abs(mlof$MAF-0.001)==min(abs(mlof$MAF-0.001))),which(abs(mlof$MAF-0.0032)==min(abs(mlof$MAF-0.0032))),which(abs(mlof$MAF-0.01)==min(abs(mlof$MAF-0.01))),which(abs(mlof$MAF-0.0316)==min(abs(mlof$MAF-0.0316))),which(abs(mlof$MAF-0.1)==min(abs(mlof$MAF-0.1))),which(abs(mlof$MAF-0.5)==min(abs(mlof$MAF-0.5)))),]
mlof=melt(mlof, id.vars="MAF")
mlof$Type=rep("Biallelic SNVs", nrow(mlof))
mbsnv=mlof
mlof=fread("../analyses/percentileaggr2/mlof.bi.indel.top.perc8.tab")
mlof=mlof[c(which(abs(mlof$MAF-0.0002)==min(abs(mlof$MAF-0.0002))),which(abs(mlof$MAF-0.0004)==min(abs(mlof$MAF-0.0004))),which(abs(mlof$MAF-0.0006)==min(abs(mlof$MAF-0.0006))),which(abs(mlof$MAF-0.001)==min(abs(mlof$MAF-0.001))),which(abs(mlof$MAF-0.0032)==min(abs(mlof$MAF-0.0032))),which(abs(mlof$MAF-0.01)==min(abs(mlof$MAF-0.01))),which(abs(mlof$MAF-0.0316)==min(abs(mlof$MAF-0.0316))),which(abs(mlof$MAF-0.1)==min(abs(mlof$MAF-0.1))),which(abs(mlof$MAF-0.5)==min(abs(mlof$MAF-0.5)))),]
mlof=melt(mlof, id.vars="MAF")
mlof$Type=rep("Biallelic indels", nrow(mlof))
mbindel=mlof
mlof=fread("../analyses/percentileaggr2/mlof.multi.snv.top.perc8.tab")
mlof=mlof[c(which(abs(mlof$MAF-0.0002)==min(abs(mlof$MAF-0.0002))),which(abs(mlof$MAF-0.0004)==min(abs(mlof$MAF-0.0004))),which(abs(mlof$MAF-0.0006)==min(abs(mlof$MAF-0.0006))),which(abs(mlof$MAF-0.001)==min(abs(mlof$MAF-0.001))),which(abs(mlof$MAF-0.0032)==min(abs(mlof$MAF-0.0032))),which(abs(mlof$MAF-0.01)==min(abs(mlof$MAF-0.01))),which(abs(mlof$MAF-0.0316)==min(abs(mlof$MAF-0.0316))),which(abs(mlof$MAF-0.1)==min(abs(mlof$MAF-0.1))),which(abs(mlof$MAF-0.5)==min(abs(mlof$MAF-0.5)))),]
mlof=melt(mlof, id.vars="MAF")
mlof$Type=rep("Multiallelic SNVs", nrow(mlof))
mmsnv=mlof
mlof=fread("../analyses/percentileaggr2/mlof.multi.indel.top.perc8.tab")
mlof=mlof[c(which(abs(mlof$MAF-0.0002)==min(abs(mlof$MAF-0.0002))),which(abs(mlof$MAF-0.0004)==min(abs(mlof$MAF-0.0004))),which(abs(mlof$MAF-0.0006)==min(abs(mlof$MAF-0.0006))),which(abs(mlof$MAF-0.001)==min(abs(mlof$MAF-0.001))),which(abs(mlof$MAF-0.0032)==min(abs(mlof$MAF-0.0032))),which(abs(mlof$MAF-0.01)==min(abs(mlof$MAF-0.01))),which(abs(mlof$MAF-0.0316)==min(abs(mlof$MAF-0.0316))),which(abs(mlof$MAF-0.1)==min(abs(mlof$MAF-0.1))),which(abs(mlof$MAF-0.5)==min(abs(mlof$MAF-0.5)))),]
mlof=melt(mlof, id.vars="MAF")
mlof$Type=rep("Multiallelic indels", nrow(mlof))
mmindel=mlof
mlof=rbind(mbsnv, mbindel, mmsnv, mmindel)
mlof$Ancestry=rep("European", nrow(mlof))

inpsyght=fread("../analyses/percentileaggr2/inpsyght.bi.snv.top.perc8.tab")
inpsyght=inpsyght[c(which(abs(inpsyght$MAF-0.0002)==min(abs(inpsyght$MAF-0.0002))),which(abs(inpsyght$MAF-0.0004)==min(abs(inpsyght$MAF-0.0004))),which(abs(inpsyght$MAF-0.0006)==min(abs(inpsyght$MAF-0.0006))),which(abs(inpsyght$MAF-0.001)==min(abs(inpsyght$MAF-0.001))),which(abs(inpsyght$MAF-0.0032)==min(abs(inpsyght$MAF-0.0032))),which(abs(inpsyght$MAF-0.01)==min(abs(inpsyght$MAF-0.01))),which(abs(inpsyght$MAF-0.0316)==min(abs(inpsyght$MAF-0.0316))),which(abs(inpsyght$MAF-0.1)==min(abs(inpsyght$MAF-0.1))),which(abs(inpsyght$MAF-0.5)==min(abs(inpsyght$MAF-0.5)))),]
inpsyght=melt(inpsyght, id.vars="MAF")
inpsyght$Type=rep("Biallelic SNVs", nrow(inpsyght))
ibsnv=inpsyght
inpsyght=fread("../analyses/percentileaggr2/inpsyght.bi.indel.top.perc8.tab")
inpsyght=inpsyght[c(which(abs(inpsyght$MAF-0.0002)==min(abs(inpsyght$MAF-0.0002))),which(abs(inpsyght$MAF-0.0004)==min(abs(inpsyght$MAF-0.0004))),which(abs(inpsyght$MAF-0.0006)==min(abs(inpsyght$MAF-0.0006))),which(abs(inpsyght$MAF-0.001)==min(abs(inpsyght$MAF-0.001))),which(abs(inpsyght$MAF-0.0032)==min(abs(inpsyght$MAF-0.0032))),which(abs(inpsyght$MAF-0.01)==min(abs(inpsyght$MAF-0.01))),which(abs(inpsyght$MAF-0.0316)==min(abs(inpsyght$MAF-0.0316))),which(abs(inpsyght$MAF-0.1)==min(abs(inpsyght$MAF-0.1))),which(abs(inpsyght$MAF-0.5)==min(abs(inpsyght$MAF-0.5)))),]
inpsyght=melt(inpsyght, id.vars="MAF")
inpsyght$Type=rep("Biallelic indels", nrow(inpsyght))
ibindel=inpsyght
inpsyght=fread("../analyses/percentileaggr2/inpsyght.multi.snv.top.perc8.tab")
inpsyght=inpsyght[c(which(abs(inpsyght$MAF-0.0002)==min(abs(inpsyght$MAF-0.0002))),which(abs(inpsyght$MAF-0.0004)==min(abs(inpsyght$MAF-0.0004))),which(abs(inpsyght$MAF-0.0006)==min(abs(inpsyght$MAF-0.0006))),which(abs(inpsyght$MAF-0.001)==min(abs(inpsyght$MAF-0.001))),which(abs(inpsyght$MAF-0.0032)==min(abs(inpsyght$MAF-0.0032))),which(abs(inpsyght$MAF-0.01)==min(abs(inpsyght$MAF-0.01))),which(abs(inpsyght$MAF-0.0316)==min(abs(inpsyght$MAF-0.0316))),which(abs(inpsyght$MAF-0.1)==min(abs(inpsyght$MAF-0.1))),which(abs(inpsyght$MAF-0.5)==min(abs(inpsyght$MAF-0.5)))),]
inpsyght=melt(inpsyght, id.vars="MAF")
inpsyght$Type=rep("Multiallelic SNVs", nrow(inpsyght))
imsnv=inpsyght
inpsyght=fread("../analyses/percentileaggr2/inpsyght.multi.indel.top.perc8.tab")
inpsyght=inpsyght[c(which(abs(inpsyght$MAF-0.0002)==min(abs(inpsyght$MAF-0.0002))),which(abs(inpsyght$MAF-0.0004)==min(abs(inpsyght$MAF-0.0004))),which(abs(inpsyght$MAF-0.0006)==min(abs(inpsyght$MAF-0.0006))),which(abs(inpsyght$MAF-0.001)==min(abs(inpsyght$MAF-0.001))),which(abs(inpsyght$MAF-0.0032)==min(abs(inpsyght$MAF-0.0032))),which(abs(inpsyght$MAF-0.01)==min(abs(inpsyght$MAF-0.01))),which(abs(inpsyght$MAF-0.0316)==min(abs(inpsyght$MAF-0.0316))),which(abs(inpsyght$MAF-0.1)==min(abs(inpsyght$MAF-0.1))),which(abs(inpsyght$MAF-0.5)==min(abs(inpsyght$MAF-0.5)))),]
inpsyght=melt(inpsyght, id.vars="MAF")
inpsyght$Type=rep("Multiallelic indels", nrow(inpsyght))
imindel=inpsyght
inpsyght=rbind(ibsnv, ibindel, imsnv, imindel)
inpsyght$Ancestry=rep("African", nrow(inpsyght))

metsim=fread("../analyses/percentileaggr2/metsim.bi.snv.top.perc8.tab")
metsim=metsim[c(which(abs(metsim$MAF-0.0002)==min(abs(metsim$MAF-0.0002))),which(abs(metsim$MAF-0.0004)==min(abs(metsim$MAF-0.0004))),which(abs(metsim$MAF-0.0006)==min(abs(metsim$MAF-0.0006))),which(abs(metsim$MAF-0.001)==min(abs(metsim$MAF-0.001))),which(abs(metsim$MAF-0.0032)==min(abs(metsim$MAF-0.0032))),which(abs(metsim$MAF-0.01)==min(abs(metsim$MAF-0.01))),which(abs(metsim$MAF-0.0316)==min(abs(metsim$MAF-0.0316))),which(abs(metsim$MAF-0.1)==min(abs(metsim$MAF-0.1))),which(abs(metsim$MAF-0.5)==min(abs(metsim$MAF-0.5)))),]
metsim=melt(metsim, id.vars="MAF")
metsim$Type=rep("Biallelic SNVs", nrow(metsim))
mtbsnv=metsim
metsim=fread("../analyses/percentileaggr2/metsim.bi.indel.top.perc8.tab")
metsim=metsim[c(which(abs(metsim$MAF-0.0002)==min(abs(metsim$MAF-0.0002))),which(abs(metsim$MAF-0.0004)==min(abs(metsim$MAF-0.0004))),which(abs(metsim$MAF-0.0006)==min(abs(metsim$MAF-0.0006))),which(abs(metsim$MAF-0.001)==min(abs(metsim$MAF-0.001))),which(abs(metsim$MAF-0.0032)==min(abs(metsim$MAF-0.0032))),which(abs(metsim$MAF-0.01)==min(abs(metsim$MAF-0.01))),which(abs(metsim$MAF-0.0316)==min(abs(metsim$MAF-0.0316))),which(abs(metsim$MAF-0.1)==min(abs(metsim$MAF-0.1))),which(abs(metsim$MAF-0.5)==min(abs(metsim$MAF-0.5)))),]
metsim=melt(metsim, id.vars="MAF")
metsim$Type=rep("Biallelic indels", nrow(metsim))
mtbindel=metsim
metsim=fread("../analyses/percentileaggr2/metsim.multi.snv.top.perc8.tab")
metsim=metsim[c(which(abs(metsim$MAF-0.0002)==min(abs(metsim$MAF-0.0002))),which(abs(metsim$MAF-0.0004)==min(abs(metsim$MAF-0.0004))),which(abs(metsim$MAF-0.0006)==min(abs(metsim$MAF-0.0006))),which(abs(metsim$MAF-0.001)==min(abs(metsim$MAF-0.001))),which(abs(metsim$MAF-0.0032)==min(abs(metsim$MAF-0.0032))),which(abs(metsim$MAF-0.01)==min(abs(metsim$MAF-0.01))),which(abs(metsim$MAF-0.0316)==min(abs(metsim$MAF-0.0316))),which(abs(metsim$MAF-0.1)==min(abs(metsim$MAF-0.1))),which(abs(metsim$MAF-0.5)==min(abs(metsim$MAF-0.5)))),]
metsim=melt(metsim, id.vars="MAF")
metsim$Type=rep("Multiallelic SNVs", nrow(metsim))
mtmsnv=metsim
metsim=fread("../analyses/percentileaggr2/metsim.multi.indel.top.perc8.tab")
metsim=metsim[c(which(abs(metsim$MAF-0.0002)==min(abs(metsim$MAF-0.0002))),which(abs(metsim$MAF-0.0004)==min(abs(metsim$MAF-0.0004))),which(abs(metsim$MAF-0.0006)==min(abs(metsim$MAF-0.0006))),which(abs(metsim$MAF-0.001)==min(abs(metsim$MAF-0.001))),which(abs(metsim$MAF-0.0032)==min(abs(metsim$MAF-0.0032))),which(abs(metsim$MAF-0.01)==min(abs(metsim$MAF-0.01))),which(abs(metsim$MAF-0.0316)==min(abs(metsim$MAF-0.0316))),which(abs(metsim$MAF-0.1)==min(abs(metsim$MAF-0.1))),which(abs(metsim$MAF-0.5)==min(abs(metsim$MAF-0.5)))),]
metsim=melt(metsim, id.vars="MAF")
metsim$Type=rep("Multiallelic indels", nrow(metsim))
mtmindel=metsim
metsim=rbind(mtbsnv, mtbindel, mtmsnv, mtmindel)
metsim$Ancestry=rep("Finnish", nrow(metsim))

all=rbind(biome, mlof, inpsyght,metsim)
all$Type=factor(all$Type, levels=c("Biallelic SNVs","Biallelic indels","Multiallelic SNVs","Multiallelic indels"))
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))
all$Array=rep(0, nrow(all))
all$Array[which(all$variable=="Core_TOPMed")]="Core (0.3M)"
all$Array[which(all$variable=="OmniExp_TOPMed")]="Omni Express (0.7M)"
all$Array[which(all$variable=="Omni25_TOPMed")]="Omni 2.5M (2.4M)"
all$Array[which(all$variable=="MEGA_TOPMed")]="MEGA (1.7M)"
all$Array=factor(all$Array, levels=c("Omni 2.5M (2.4M)","MEGA (1.7M)","Omni Express (0.7M)","Core (0.3M)"))

ybrks<-c(0,0.2,0.4,0.6,0.8,1)
brks<-c(0.5, 0.1, 0.01,0.001,0.0001)
labs=c(expression(0.5),expression(10^-1), expression(10^-2), expression(10^-3),expression(10^-4))

tiff(file="~/FigS12.tiff", height=174, width=174, units="mm", res=300)
ggplot(all, aes(x=MAF,y=value, color=Type))+scale_x_log10(breaks=brks, labels=labs, limits=c(1e-4,0.5),expand=c(0,0))+scale_y_continuous(breaks=ybrks,limits=c(0,1), expand=c(0,0))+theme_bw()+facet_grid(Array~Ancestry)+theme(legend.position="bottom",legend.text=element_text(size=7),panel.grid.major=element_line(colour="#e7e7e7"),panel.spacing=unit(0.7, "lines"), strip.text=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=7),panel.border=element_rect(color="black"),strip.background=element_rect(fill="white",color=NA), legend.title=element_text(size=8), plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"),legend.key.width=unit(3,"lines"))+geom_line(aes(group=Type), size=0.7)+xlab("WGS Minor Allele Frequency")+ylab(expression(atop(Proportion~of,well-imputed~variants)))+guides(color=guide_legend(nrow=2))+scale_color_manual(values=c("#1b9e77", "#e6ab02","#7570b3","#e7298a"))
dev.off()

