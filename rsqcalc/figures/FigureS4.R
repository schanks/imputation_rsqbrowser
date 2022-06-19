library(data.table)
library(ggplot2)
library(gridExtra)

gatk=fread("../analyses/varcaller/metsim.gatk.bi.snv.perc8.tab")
gatk=gatk[,c(1:13)]
gatk=gatk[c(which(abs(gatk$MAF-0.0002)==min(abs(gatk$MAF-0.0002))),which(abs(gatk$MAF-0.0004)==min(abs(gatk$MAF-0.0004))),which(abs(gatk$MAF-0.0006)==min(abs(gatk$MAF-0.0006))),which(abs(gatk$MAF-0.001)==min(abs(gatk$MAF-0.001))),which(abs(gatk$MAF-0.0032)==min(abs(gatk$MAF-0.0032))),which(abs(gatk$MAF-0.01)==min(abs(gatk$MAF-0.01))),which(abs(gatk$MAF-0.0316)==min(abs(gatk$MAF-0.0316))),which(abs(gatk$MAF-0.1)==min(abs(gatk$MAF-0.1))),which(abs(gatk$MAF-0.5)==min(abs(gatk$MAF-0.5)))),]
gatk=melt(gatk, id="MAF")
gatk$Array=rep(0, nrow(gatk))
gatk$Array[which(gatk$variable=="Core_1000G" | gatk$variable=="Core_HRC" | gatk$variable=="Core_TOPMed")]="Core"
gatk$Array[which(gatk$variable=="OmniExp_1000G" | gatk$variable=="OmniExp_HRC" | gatk$variable=="OmniExp_TOPMed")]="Omni Express"
gatk$Array[which(gatk$variable=="Omni25_1000G" | gatk$variable=="Omni25_HRC" | gatk$variable=="Omni25_TOPMed")]="Omni 2.5M"
gatk$Array[which(gatk$variable=="MEGA_1000G" | gatk$variable=="MEGA_HRC" | gatk$variable=="MEGA_TOPMed")]="MEGA"
gatk$Array=as.factor(gatk$Array)
gatk$`Reference Panel`=rep(0, nrow(gatk))
gatk$`Reference Panel`[which(gatk$variable=="Core_1000G" | gatk$variable=="OmniExp_1000G" | gatk$variable=="Omni25_1000G" | gatk$variable=="MEGA_1000G")]="1000G"
gatk$`Reference Panel`[which(gatk$variable=="Core_HRC" | gatk$variable=="OmniExp_HRC" | gatk$variable=="Omni25_HRC"| gatk$variable=="MEGA_HRC")]="HRC"
gatk$`Reference Panel`[which(gatk$variable=="Core_TOPMed" | gatk$variable=="OmniExp_TOPMed" | gatk$variable=="Omni25_TOPMed" | gatk$variable=="MEGA_TOPMed")]="TOPMed"
gatk$`Reference Panel`=as.factor(gatk$`Reference Panel`)
gatk$VariantCaller=rep("GATK", nrow(gatk))

metsim=fread("../analyses/percentileaggr2/metsim.bi.snv.perc8.tab")
metsim=metsim[,c(1:13)]
metsim=metsim[c(which(abs(metsim$MAF-0.0002)==min(abs(metsim$MAF-0.0002))),which(abs(metsim$MAF-0.0004)==min(abs(metsim$MAF-0.0004))),which(abs(metsim$MAF-0.0006)==min(abs(metsim$MAF-0.0006))),which(abs(metsim$MAF-0.001)==min(abs(metsim$MAF-0.001))),which(abs(metsim$MAF-0.0032)==min(abs(metsim$MAF-0.0032))),which(abs(metsim$MAF-0.01)==min(abs(metsim$MAF-0.01))),which(abs(metsim$MAF-0.0316)==min(abs(metsim$MAF-0.0316))),which(abs(metsim$MAF-0.1)==min(abs(metsim$MAF-0.1))),which(abs(metsim$MAF-0.5)==min(abs(metsim$MAF-0.5)))),]
metsim=melt(metsim, id="MAF")
metsim$Array=rep(0, nrow(metsim))
metsim$Array[which(metsim$variable=="Core_1000G" | metsim$variable=="Core_HRC" | metsim$variable=="Core_TOPMed")]="Core"
metsim$Array[which(metsim$variable=="OmniExp_1000G" | metsim$variable=="OmniExp_HRC" | metsim$variable=="OmniExp_TOPMed")]="Omni Express"
metsim$Array[which(metsim$variable=="Omni25_1000G" | metsim$variable=="Omni25_HRC" | metsim$variable=="Omni25_TOPMed")]="Omni 2.5M"
metsim$Array[which(metsim$variable=="MEGA_1000G" | metsim$variable=="MEGA_HRC" | metsim$variable=="MEGA_TOPMed")]="MEGA"
metsim$Array=as.factor(metsim$Array)
metsim$`Reference Panel`=rep(0, nrow(metsim))
metsim$`Reference Panel`[which(metsim$variable=="Core_1000G" | metsim$variable=="OmniExp_1000G" | metsim$variable=="Omni25_1000G" | metsim$variable=="MEGA_1000G")]="1000G"
metsim$`Reference Panel`[which(metsim$variable=="Core_HRC" | metsim$variable=="OmniExp_HRC" | metsim$variable=="Omni25_HRC"| metsim$variable=="MEGA_HRC")]="HRC"
metsim$`Reference Panel`[which(metsim$variable=="Core_TOPMed" | metsim$variable=="OmniExp_TOPMed" | metsim$variable=="Omni25_TOPMed" | metsim$variable=="MEGA_TOPMed")]="TOPMed"
metsim$`Reference Panel`=as.factor(metsim$`Reference Panel`)
metsim$VariantCaller=rep("GotCloud", nrow(metsim))

prop=rbind(gatk,metsim)
prop$Measure=rep("proportion", nrow(prop))

gatk=fread("../analyses/varcaller/metsim.gatk.bi.snv.mean.tab")
gatk=fread("../analyses/varcaller/metsim.gatk.bi.snv.perc8.tab")
gatk=gatk[,c(1:13)]
gatk=gatk[c(which(abs(gatk$MAF-0.0002)==min(abs(gatk$MAF-0.0002))),which(abs(gatk$MAF-0.0004)==min(abs(gatk$MAF-0.0004))),which(abs(gatk$MAF-0.0006)==min(abs(gatk$MAF-0.0006))),which(abs(gatk$MAF-0.001)==min(abs(gatk$MAF-0.001))),which(abs(gatk$MAF-0.0032)==min(abs(gatk$MAF-0.0032))),which(abs(gatk$MAF-0.01)==min(abs(gatk$MAF-0.01))),which(abs(gatk$MAF-0.0316)==min(abs(gatk$MAF-0.0316))),which(abs(gatk$MAF-0.1)==min(abs(gatk$MAF-0.1))),which(abs(gatk$MAF-0.5)==min(abs(gatk$MAF-0.5)))),]
gatk=melt(gatk, id="MAF")
gatk$Array=rep(0, nrow(gatk))
gatk$Array[which(gatk$variable=="Core_1000G" | gatk$variable=="Core_HRC" | gatk$variable=="Core_TOPMed")]="Core"
gatk$Array[which(gatk$variable=="OmniExp_1000G" | gatk$variable=="OmniExp_HRC" | gatk$variable=="OmniExp_TOPMed")]="Omni Express"
gatk$Array[which(gatk$variable=="Omni25_1000G" | gatk$variable=="Omni25_HRC" | gatk$variable=="Omni25_TOPMed")]="Omni 2.5M"
gatk$Array[which(gatk$variable=="MEGA_1000G" | gatk$variable=="MEGA_HRC" | gatk$variable=="MEGA_TOPMed")]="MEGA"
gatk$Array=as.factor(gatk$Array)
gatk$`Reference Panel`=rep(0, nrow(gatk))
gatk$`Reference Panel`[which(gatk$variable=="Core_1000G" | gatk$variable=="OmniExp_1000G" | gatk$variable=="Omni25_1000G" | gatk$variable=="MEGA_1000G")]="1000G"
gatk$`Reference Panel`[which(gatk$variable=="Core_HRC" | gatk$variable=="OmniExp_HRC" | gatk$variable=="Omni25_HRC"| gatk$variable=="MEGA_HRC")]="HRC"
gatk$`Reference Panel`[which(gatk$variable=="Core_TOPMed" | gatk$variable=="OmniExp_TOPMed" | gatk$variable=="Omni25_TOPMed" | gatk$variable=="MEGA_TOPMed")]="TOPMed"
gatk$`Reference Panel`=as.factor(gatk$`Reference Panel`)
gatk$VariantCaller=rep("GATK", nrow(gatk))

metsim=fread("../analyses/percentileaggr2/metsim.bi.snv.mean.tab")
metsim=metsim[,c(1:13)]
metsim=metsim[c(which(abs(metsim$MAF-0.0002)==min(abs(metsim$MAF-0.0002))),which(abs(metsim$MAF-0.0004)==min(abs(metsim$MAF-0.0004))),which(abs(metsim$MAF-0.0006)==min(abs(metsim$MAF-0.0006))),which(abs(metsim$MAF-0.001)==min(abs(metsim$MAF-0.001))),which(abs(metsim$MAF-0.0032)==min(abs(metsim$MAF-0.0032))),which(abs(metsim$MAF-0.01)==min(abs(metsim$MAF-0.01))),which(abs(metsim$MAF-0.0316)==min(abs(metsim$MAF-0.0316))),which(abs(metsim$MAF-0.1)==min(abs(metsim$MAF-0.1))),which(abs(metsim$MAF-0.5)==min(abs(metsim$MAF-0.5)))),]
metsim=melt(metsim, id="MAF")
metsim$Array=rep(0, nrow(metsim))
metsim$Array[which(metsim$variable=="Core_1000G" | metsim$variable=="Core_HRC" | metsim$variable=="Core_TOPMed")]="Core"
metsim$Array[which(metsim$variable=="OmniExp_1000G" | metsim$variable=="OmniExp_HRC" | metsim$variable=="OmniExp_TOPMed")]="Omni Express"
metsim$Array[which(metsim$variable=="Omni25_1000G" | metsim$variable=="Omni25_HRC" | metsim$variable=="Omni25_TOPMed")]="Omni 2.5M"
metsim$Array[which(metsim$variable=="MEGA_1000G" | metsim$variable=="MEGA_HRC" | metsim$variable=="MEGA_TOPMed")]="MEGA"
metsim$Array=as.factor(metsim$Array)
metsim$`Reference Panel`=rep(0, nrow(metsim))
metsim$`Reference Panel`[which(metsim$variable=="Core_1000G" | metsim$variable=="OmniExp_1000G" | metsim$variable=="Omni25_1000G" | metsim$variable=="MEGA_1000G")]="1000G"
metsim$`Reference Panel`[which(metsim$variable=="Core_HRC" | metsim$variable=="OmniExp_HRC" | metsim$variable=="Omni25_HRC"| metsim$variable=="MEGA_HRC")]="HRC"
metsim$`Reference Panel`[which(metsim$variable=="Core_TOPMed" | metsim$variable=="OmniExp_TOPMed" | metsim$variable=="Omni25_TOPMed" | metsim$variable=="MEGA_TOPMed")]="TOPMed"
metsim$`Reference Panel`=as.factor(metsim$`Reference Panel`)
metsim$VariantCaller=rep("GotCloud", nrow(metsim))

mean=rbind(gatk,metsim)
mean$Measure=rep("mean", nrow(mean))

##Plotting
all=rbind(prop,mean)
all$Array=factor(all$Array, levels=c("Omni 2.5M","MEGA","Omni Express","Core"))
all$`Reference Panel`=factor(all$`Reference Panel`, levels=c("TOPMed","HRC","1000G"))
all$VariantCaller=factor(all$VariantCaller, levels=c("GotCloud","GATK"))
ybrks<-c(0,0.2,0.4,0.6,0.8,1)
brks<-c(0.5, 0.1, 0.01,0.001,0.0001)
labs=c(expression(0.5),expression(10^-1), expression(10^-2), expression(10^-3),expression(10^-4))

a=ggplot(all[which(all$Measure=="proportion"),], aes(x=MAF,y=value, color=VariantCaller))+scale_x_log10(breaks=brks, labels=labs, limits=c(1e-4,0.5),expand=c(0,0))+scale_y_continuous(breaks=ybrks,limits=c(0,1), expand=c(0,0))+theme_bw()+facet_grid(`Reference Panel`~Array)+theme(legend.position="none",panel.grid.major=element_line(colour="#e7e7e7"),panel.spacing=unit(0.7, "lines"), strip.text=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=7),panel.border=element_rect(color="black"),strip.background=element_rect(fill="white",color=NA), plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"))+geom_line(aes(group=VariantCaller), size=0.7)+xlab("WGS Minor Allele Frequency")+ylab(expression(atop(Mean~observed,imputation~r^2)))

b=ggplot(all[which(all$Measure=="mean"),], aes(x=MAF,y=value, color=VariantCaller))+scale_x_log10(breaks=brks, labels=labs, limits=c(1e-4,0.5),expand=c(0,0))+scale_y_continuous(breaks=ybrks,limits=c(0,1), expand=c(0,0))+theme_bw()+facet_grid(`Reference Panel`~Array)+theme(legend.position="bottom",legend.text=element_text(size=7),panel.grid.major=element_line(colour="#e7e7e7"),panel.spacing=unit(0.7, "lines"), strip.text=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=7),panel.border=element_rect(color="black"),strip.background=element_rect(fill="white",color=NA), legend.title=element_text(size=8), plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"),legend.key.width=unit(2,"lines"))+geom_line(aes(group=VariantCaller), size=0.7)+xlab("WGS Minor Allele Frequency")+ylab(expression(atop(Mean~observed,imputation~r^2)))+guides(color=guide_legend(nrow=2))

tiff(file="~/FigS4.tiff", height=174, width=174, units="mm", res=300)
grid.arrange(a,b,ncol=1, heights=c(1,1.3))
dev.off()

