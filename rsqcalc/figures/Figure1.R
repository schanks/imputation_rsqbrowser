library(data.table)
library(ggplot2)
library(reshape2)
library(gridExtra)

biome=fread("../analyses/percentileaggr2/biome.bi.snv.perc8.tab")
biome=biome[c(which(abs(biome$MAF-0.0002)==min(abs(biome$MAF-0.0002))),which(abs(biome$MAF-0.0004)==min(abs(biome$MAF-0.0004))),which(abs(biome$MAF-0.0006)==min(abs(biome$MAF-0.0006))),which(abs(biome$MAF-0.001)==min(abs(biome$MAF-0.001))),which(abs(biome$MAF-0.0032)==min(abs(biome$MAF-0.0032))),which(abs(biome$MAF-0.01)==min(abs(biome$MAF-0.01))),which(abs(biome$MAF-0.0316)==min(abs(biome$MAF-0.0316))),which(abs(biome$MAF-0.1)==min(abs(biome$MAF-0.1))),which(abs(biome$MAF-0.316)==min(abs(biome$MAF-0.316))),which(abs(biome$MAF-0.5)==min(abs(biome$MAF-0.5)))),]
biome=melt(biome, id="MAF")
biome$Array=rep(0, nrow(biome))
biome$Array[which(biome$variable=="Core_1000G" | biome$variable=="Core_HRC" | biome$variable=="Core_TOPMed")]="Core (0.3M)"
biome$Array[which(biome$variable=="OmniExp_1000G" | biome$variable=="OmniExp_HRC" | biome$variable=="OmniExp_TOPMed")]="Omni Express (0.7M)"
biome$Array[which(biome$variable=="Omni25_1000G" | biome$variable=="Omni25_HRC" | biome$variable=="Omni25_TOPMed")]="Omni 2.5M (2.4M)"
biome$Array[which(biome$variable=="MEGA_1000G" | biome$variable=="MEGA_HRC" | biome$variable=="MEGA_TOPMed")]="MEGA (1.7M)"
biome$Array=as.factor(biome$Array)
biome$`Reference Panel`=rep(0, nrow(biome))
biome$`Reference Panel`[which(biome$variable=="Core_1000G" | biome$variable=="OmniExp_1000G" | biome$variable=="Omni25_1000G" | biome$variable=="MEGA_1000G")]="1000G"
biome$`Reference Panel`[which(biome$variable=="Core_HRC" | biome$variable=="OmniExp_HRC" | biome$variable=="Omni25_HRC"| biome$variable=="MEGA_HRC")]="HRC"
biome$`Reference Panel`[which(biome$variable=="Core_TOPMed" | biome$variable=="OmniExp_TOPMed" | biome$variable=="Omni25_TOPMed" | biome$variable=="MEGA_TOPMed")]="TOPMed"
biome$`Reference Panel`=as.factor(biome$`Reference Panel`)
biome$Ancestry=rep("Hispanic/Latino", nrow(biome))

mlof=fread("../analyses/percentileaggr2/mlof.bi.snv.perc8.tab")
mlof=mlof[c(which(abs(mlof$MAF-0.0002)==min(abs(mlof$MAF-0.0002))),which(abs(mlof$MAF-0.0004)==min(abs(mlof$MAF-0.0004))),which(abs(mlof$MAF-0.0006)==min(abs(mlof$MAF-0.0006))),which(abs(mlof$MAF-0.001)==min(abs(mlof$MAF-0.001))),which(abs(mlof$MAF-0.0032)==min(abs(mlof$MAF-0.0032))),which(abs(mlof$MAF-0.01)==min(abs(mlof$MAF-0.01))),which(abs(mlof$MAF-0.0316)==min(abs(mlof$MAF-0.0316))),which(abs(mlof$MAF-0.1)==min(abs(mlof$MAF-0.1))),which(abs(mlof$MAF-0.316)==min(abs(mlof$MAF-0.316))),which(abs(mlof$MAF-0.5)==min(abs(mlof$MAF-0.5)))),]
mlof=melt(mlof, id="MAF")
mlof$Array=rep(0, nrow(mlof))
mlof$Array[which(mlof$variable=="Core_1000G" | mlof$variable=="Core_HRC" | mlof$variable=="Core_TOPMed")]="Core (0.3M)"
mlof$Array[which(mlof$variable=="OmniExp_1000G" | mlof$variable=="OmniExp_HRC" | mlof$variable=="OmniExp_TOPMed")]="Omni Express (0.7M)"
mlof$Array[which(mlof$variable=="Omni25_1000G" | mlof$variable=="Omni25_HRC" | mlof$variable=="Omni25_TOPMed")]="Omni 2.5M (2.4M)"
mlof$Array[which(mlof$variable=="MEGA_1000G" | mlof$variable=="MEGA_HRC" | mlof$variable=="MEGA_TOPMed")]="MEGA (1.7M)"
mlof$Array=as.factor(mlof$Array)
mlof$`Reference Panel`=rep(0, nrow(mlof))
mlof$`Reference Panel`[which(mlof$variable=="Core_1000G" | mlof$variable=="OmniExp_1000G" | mlof$variable=="Omni25_1000G" | mlof$variable=="MEGA_1000G")]="1000G"
mlof$`Reference Panel`[which(mlof$variable=="Core_HRC" | mlof$variable=="OmniExp_HRC" | mlof$variable=="Omni25_HRC"| mlof$variable=="MEGA_HRC")]="HRC"
mlof$`Reference Panel`[which(mlof$variable=="Core_TOPMed" | mlof$variable=="OmniExp_TOPMed" | mlof$variable=="Omni25_TOPMed" | mlof$variable=="MEGA_TOPMed")]="TOPMed"
mlof$`Reference Panel`=as.factor(mlof$`Reference Panel`)
mlof$Ancestry=rep("European", nrow(mlof))

inpsyght=fread("../analyses/percentileaggr2/inpsyght.bi.snv.perc8.tab")
inpsyght=inpsyght[c(which(abs(inpsyght$MAF-0.0002)==min(abs(inpsyght$MAF-0.0002))),which(abs(inpsyght$MAF-0.0004)==min(abs(inpsyght$MAF-0.0004))),which(abs(inpsyght$MAF-0.0006)==min(abs(inpsyght$MAF-0.0006))),which(abs(inpsyght$MAF-0.001)==min(abs(inpsyght$MAF-0.001))),which(abs(inpsyght$MAF-0.0032)==min(abs(inpsyght$MAF-0.0032))),which(abs(inpsyght$MAF-0.01)==min(abs(inpsyght$MAF-0.01))),which(abs(inpsyght$MAF-0.0316)==min(abs(inpsyght$MAF-0.0316))),which(abs(inpsyght$MAF-0.1)==min(abs(inpsyght$MAF-0.1))),which(abs(inpsyght$MAF-0.316)==min(abs(inpsyght$MAF-0.316))),which(abs(inpsyght$MAF-0.5)==min(abs(inpsyght$MAF-0.5)))),]
inpsyght=melt(inpsyght, id="MAF")
inpsyght$Array=rep(0, nrow(inpsyght))
inpsyght$Array[which(inpsyght$variable=="Core_1000G" | inpsyght$variable=="Core_HRC" | inpsyght$variable=="Core_TOPMed")]="Core (0.3M)"
inpsyght$Array[which(inpsyght$variable=="OmniExp_1000G" | inpsyght$variable=="OmniExp_HRC" | inpsyght$variable=="OmniExp_TOPMed")]="Omni Express (0.7M)"
inpsyght$Array[which(inpsyght$variable=="Omni25_1000G" | inpsyght$variable=="Omni25_HRC" | inpsyght$variable=="Omni25_TOPMed")]="Omni 2.5M (2.4M)"
inpsyght$Array[which(inpsyght$variable=="MEGA_1000G" | inpsyght$variable=="MEGA_HRC" | inpsyght$variable=="MEGA_TOPMed")]="MEGA (1.7M)"
inpsyght$Array=as.factor(inpsyght$Array)
inpsyght$`Reference Panel`=rep(0, nrow(inpsyght))
inpsyght$`Reference Panel`[which(inpsyght$variable=="Core_1000G" | inpsyght$variable=="OmniExp_1000G" | inpsyght$variable=="Omni25_1000G" | inpsyght$variable=="MEGA_1000G")]="1000G"
inpsyght$`Reference Panel`[which(inpsyght$variable=="Core_HRC" | inpsyght$variable=="OmniExp_HRC" | inpsyght$variable=="Omni25_HRC"| inpsyght$variable=="MEGA_HRC")]="HRC"
inpsyght$`Reference Panel`[which(inpsyght$variable=="Core_TOPMed" | inpsyght$variable=="OmniExp_TOPMed" | inpsyght$variable=="Omni25_TOPMed" | inpsyght$variable=="MEGA_TOPMed")]="TOPMed"
inpsyght$`Reference Panel`=as.factor(inpsyght$`Reference Panel`)
inpsyght$Ancestry=rep("African", nrow(inpsyght))

metsim=fread("../analyses/percentileaggr2/metsim.bi.snv.perc8.tab")
metsim=metsim[c(which(abs(metsim$MAF-0.0002)==min(abs(metsim$MAF-0.0002))),which(abs(metsim$MAF-0.0004)==min(abs(metsim$MAF-0.0004))),which(abs(metsim$MAF-0.0006)==min(abs(metsim$MAF-0.0006))),which(abs(metsim$MAF-0.001)==min(abs(metsim$MAF-0.001))),which(abs(metsim$MAF-0.0032)==min(abs(metsim$MAF-0.0032))),which(abs(metsim$MAF-0.01)==min(abs(metsim$MAF-0.01))),which(abs(metsim$MAF-0.0316)==min(abs(metsim$MAF-0.0316))),which(abs(metsim$MAF-0.1)==min(abs(metsim$MAF-0.1))),which(abs(metsim$MAF-0.316)==min(abs(metsim$MAF-0.316))),which(abs(metsim$MAF-0.5)==min(abs(metsim$MAF-0.5)))),]
metsim=melt(metsim, id="MAF")
metsim$Array=rep(0, nrow(metsim))
metsim$Array[which(metsim$variable=="Core_1000G" | metsim$variable=="Core_HRC" | metsim$variable=="Core_TOPMed")]="Core (0.3M)"
metsim$Array[which(metsim$variable=="OmniExp_1000G" | metsim$variable=="OmniExp_HRC" | metsim$variable=="OmniExp_TOPMed")]="Omni Express (0.7M)"
metsim$Array[which(metsim$variable=="Omni25_1000G" | metsim$variable=="Omni25_HRC" | metsim$variable=="Omni25_TOPMed")]="Omni 2.5M (2.4M)"
metsim$Array[which(metsim$variable=="MEGA_1000G" | metsim$variable=="MEGA_HRC" | metsim$variable=="MEGA_TOPMed")]="MEGA (1.7M)"
metsim$Array=as.factor(metsim$Array)
metsim$`Reference Panel`=rep(0, nrow(metsim))
metsim$`Reference Panel`[which(metsim$variable=="Core_1000G" | metsim$variable=="OmniExp_1000G" | metsim$variable=="Omni25_1000G" | metsim$variable=="MEGA_1000G")]="1000G"
metsim$`Reference Panel`[which(metsim$variable=="Core_HRC" | metsim$variable=="OmniExp_HRC" | metsim$variable=="Omni25_HRC"| metsim$variable=="MEGA_HRC")]="HRC"
metsim$`Reference Panel`[which(metsim$variable=="Core_TOPMed" | metsim$variable=="OmniExp_TOPMed" | metsim$variable=="Omni25_TOPMed" | metsim$variable=="MEGA_TOPMed")]="TOPMed"
metsim$`Reference Panel`=as.factor(metsim$`Reference Panel`)
metsim$Ancestry=rep("Finnish", nrow(metsim))

all=rbind(biome, mlof, inpsyght,metsim)
all$Array=factor(all$Array, levels=c("Omni 2.5M (2.4M)","MEGA (1.7M)","Omni Express (0.7M)","Core (0.3M)"))
all$`Reference Panel`=factor(all$`Reference Panel`, levels=c("TOPMed","HRC","1000G"))
ybrks<-c(0,0.2,0.4,0.6,0.8,1)
brks<-c(0.5, 0.1, 0.01,0.001,0.0001)
labs=c(expression(0.5),expression(10^-1), expression(10^-2), expression(10^-3),expression(10^-4))
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))
all=all[order(all$Array),]


all$`Ancestry: Reference Panel`=rep("Hispanic/Latino: TOPMed", nrow(all))
all$`Ancestry: Reference Panel`[which(all$Ancestry=="Hispanic/Latino" & all$`Reference Panel`=="HRC")]="Hispanic/Latino: HRC"
all$`Ancestry: Reference Panel`[which(all$Ancestry=="Hispanic/Latino" & all$`Reference Panel`=="1000G")]="Hispanic/Latino: 1000G"
all$`Ancestry: Reference Panel`[which(all$Ancestry=="African" & all$`Reference Panel`=="TOPMed")]="African: TOPMed"
all$`Ancestry: Reference Panel`[which(all$Ancestry=="African" & all$`Reference Panel`=="HRC")]="African: HRC"
all$`Ancestry: Reference Panel`[which(all$Ancestry=="African" & all$`Reference Panel`=="1000G")]="African: 1000G"
all$`Ancestry: Reference Panel`[which(all$Ancestry=="European" & all$`Reference Panel`=="TOPMed")]="European: TOPMed"
all$`Ancestry: Reference Panel`[which(all$Ancestry=="European" & all$`Reference Panel`=="HRC")]="European: HRC"
all$`Ancestry: Reference Panel`[which(all$Ancestry=="European" & all$`Reference Panel`=="1000G")]="European: 1000G"
all$`Ancestry: Reference Panel`[which(all$Ancestry=="Finnish" & all$`Reference Panel`=="TOPMed")]="Finnish: TOPMed"
all$`Ancestry: Reference Panel`[which(all$Ancestry=="Finnish" & all$`Reference Panel`=="HRC")]="Finnish: HRC"
all$`Ancestry: Reference Panel`[which(all$Ancestry=="Finnish" & all$`Reference Panel`=="1000G")]="Finnish: 1000G"
all$`Ancestry: Reference Panel`=factor(all$`Ancestry: Reference Panel`, levels=c("African: TOPMed","African: HRC","African: 1000G","Hispanic/Latino: TOPMed","Hispanic/Latino: HRC","Hispanic/Latino: 1000G","European: TOPMed","European: HRC","European: 1000G","Finnish: TOPMed","Finnish: HRC","Finnish: 1000G"))


omni=all[which(all$Array=="Omni 2.5M (2.4M)"),]
omni$Array=rep("Omni 2.5M", nrow(omni))

a=ggplot(omni, aes(x=MAF,y=value, color=`Ancestry: Reference Panel`, size=`Reference Panel`))+scale_x_log10(breaks=brks, labels=labs, limits=c(1e-4,0.5),expand=c(0,0))+scale_y_continuous(breaks=ybrks,limits=c(0,1), expand=c(0,0))+theme_bw()+facet_grid(Array~Ancestry)+theme(legend.position="none",panel.grid.major=element_line(colour="#e7e7e7"),panel.grid.minor=element_blank(),panel.spacing.y=unit(1, "lines"), strip.text=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=7),panel.border=element_rect(color="black"),strip.background=element_rect(fill="white",color=NA), panel.spacing.x=unit(0.8, "lines"),plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"), legend.key.width=unit(3, "lines"),legend.text=element_text(size=8),legend.title=element_text(size=8))+geom_line(aes(group=`Reference Panel`))+scale_size_manual(values=c(0.7,0.6,0.5))+scale_color_manual(values=c("#cb181d","#ef3b2c","#fb6a4a","#238b45","#41ab5d","#74c476","#2171b5","#4292c6","#6baed6","#6a51a3","#807dba","#9e9ac8"))+xlab("WGS Minor Allele Frequency")+ylab(expression(atop(Proportion~of,well-imputed~SNVs)))+guides(color=guide_legend(nrow=3))

b=ggplot(omni, aes(x=MAF,y=value, color=`Ancestry: Reference Panel`, size=`Reference Panel`))+scale_x_log10(breaks=brks, labels=labs, limits=c(1e-4,0.5),expand=c(0,0))+scale_y_continuous(breaks=ybrks,limits=c(0,1), expand=c(0,0))+theme_bw()+facet_grid(`Reference Panel`~Array)+theme(legend.position="none",panel.grid.major=element_line(colour="#e7e7e7"),panel.grid.minor=element_blank(),panel.spacing.y=unit(1, "lines"), strip.text=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=7),panel.border=element_rect(color="black"),strip.background=element_rect(fill="white",color=NA), panel.spacing=unit(0.5, "lines"),plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"), legend.key.width=unit(3, "lines"),legend.text=element_text(size=8),legend.title=element_text(size=8))+geom_line(aes(group=Ancestry))+scale_size_manual(values=c(0.7,0.6,0.5))+scale_color_manual(values=c("#cb181d","#ef3b2c","#fb6a4a","#238b45","#41ab5d","#74c476","#2171b5","#4292c6","#6baed6","#6a51a3","#807dba","#9e9ac8"))+xlab("WGS Minor Allele Frequency")+ylab(expression(atop(Proportion~of,well-imputed~SNVs)))+guides(color=guide_legend(nrow=3))

c=ggplot(all, aes(x=MAF,y=value, color=`Ancestry: Reference Panel`, linetype=Array, size=`Ancestry: Reference Panel`))+scale_x_log10(breaks=brks, labels=labs, limits=c(1e-4,0.5),expand=c(0,0))+scale_y_continuous(breaks=ybrks,limits=c(0,1), expand=c(0,0))+theme_bw()+facet_grid(`Reference Panel`~Ancestry)+theme(legend.position="bottom",panel.grid.major=element_line(colour="#e7e7e7"),panel.grid.minor=element_blank(),panel.spacing.y=unit(1, "lines"), strip.text=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=7),panel.border=element_rect(color="black"),strip.background=element_rect(fill="white",color=NA), panel.spacing.x=unit(0.8, "lines"),plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"), legend.key.width=unit(2, "lines"), legend.key.height=unit(0.5, "lines"),legend.text=element_text(size=7, margin=margin(0,0,0,0)),legend.spacing=unit(0.05, "lines"),legend.title=element_text(size=8), legend.box="vertical")+geom_line(aes(group=Array))+scale_size_manual(values=c(rep(c(0.7,0.6,0.5),4)))+scale_color_manual(values=c("#cb181d","#ef3b2c","#fb6a4a","#238b45","#41ab5d","#74c476","#2171b5","#4292c6","#6baed6","#6a51a3","#807dba","#9e9ac8"))+scale_linetype_manual(values=c("solid","F1","dashed","dotted"))+xlab("WGS Minor Allele Frequency")+ylab(expression(atop(Proportion~of,well-imputed~SNVs)))+guides(color=guide_legend(title.position="top", title.hjust=0.5), linetype=guide_legend(title.position="top", keywidth=3, title.hjust=0.5), size=guide_legend(title.position="top", title.hjust=0.5))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend=get_legend(c)
c=c+theme(legend.position="none")

tiff(file="~/Fig1.tiff", height=174, width=174, units="mm", res=300)
grid.arrange(ggplotGrob(a),ggplotGrob(b),ggplotGrob(c),legend, widths=c(1.1,3), heights=c(1,2.3,0.7), layout_matrix=rbind(c(NA,1),c(2,3),c(5)))
dev.off()
