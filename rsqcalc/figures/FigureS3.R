library(data.table)
library(ggplot2)
library(gridExtra)

ho=fread("../analyses/omniexp/metsim.humanomni.bi.snv.perc8.tab")
ho=ho[c(which(abs(ho$MAF-0.0002)==min(abs(ho$MAF-0.0002))),which(abs(ho$MAF-0.0004)==min(abs(ho$MAF-0.0004))),which(abs(ho$MAF-0.0006)==min(abs(ho$MAF-0.0006))),which(abs(ho$MAF-0.001)==min(abs(ho$MAF-0.001))),which(abs(ho$MAF-0.0032)==min(abs(ho$MAF-0.0032))),which(abs(ho$MAF-0.01)==min(abs(ho$MAF-0.01))),which(abs(ho$MAF-0.0316)==min(abs(ho$MAF-0.0316))),which(abs(ho$MAF-0.1)==min(abs(ho$MAF-0.1))),which(abs(ho$MAF-0.5)==min(abs(ho$MAF-0.5)))),]
ho=melt(ho, id="MAF")
ho$`Reference Panel`=rep(0, nrow(ho))
ho$`Reference Panel`[which(ho$variable=="HO_1000G")]="1000G"
ho$`Reference Panel`[which(ho$variable=="HO_HRC")]="HRC"
ho$`Reference Panel`[which(ho$variable=="HO_TOPMed")]="TOPMed"
ho$`Reference Panel`=as.factor(ho$`Reference Panel`)
ho$Array=rep("Real array", nrow(ho))

metsim=fread("../analyses/percentileaggr2/metsim.bi.snv.perc8.tab")
metsim=metsim[,c("MAF","OmniExp_1000G","OmniExp_HRC","OmniExp_TOPMed")]
metsim=metsim[c(which(abs(metsim$MAF-0.0002)==min(abs(metsim$MAF-0.0002))),which(abs(metsim$MAF-0.0004)==min(abs(metsim$MAF-0.0004))),which(abs(metsim$MAF-0.0006)==min(abs(metsim$MAF-0.0006))),which(abs(metsim$MAF-0.001)==min(abs(metsim$MAF-0.001))),which(abs(metsim$MAF-0.0032)==min(abs(metsim$MAF-0.0032))),which(abs(metsim$MAF-0.01)==min(abs(metsim$MAF-0.01))),which(abs(metsim$MAF-0.0316)==min(abs(metsim$MAF-0.0316))),which(abs(metsim$MAF-0.1)==min(abs(metsim$MAF-0.1))),which(abs(metsim$MAF-0.5)==min(abs(metsim$MAF-0.5)))),]
metsim=melt(metsim, id="MAF")
metsim$`Reference Panel`=rep(0, nrow(metsim))
metsim$`Reference Panel`[which(metsim$variable=="Core_1000G" | metsim$variable=="OmniExp_1000G" | metsim$variable=="Omni25_1000G" | metsim$variable=="MEGA_1000G")]="1000G"
metsim$`Reference Panel`[which(metsim$variable=="Core_HRC" | metsim$variable=="OmniExp_HRC" | metsim$variable=="Omni25_HRC"| metsim$variable=="MEGA_HRC")]="HRC"
metsim$`Reference Panel`[which(metsim$variable=="Core_TOPMed" | metsim$variable=="OmniExp_TOPMed" | metsim$variable=="Omni25_TOPMed" | metsim$variable=="MEGA_TOPMed")]="TOPMed"
metsim$`Reference Panel`=as.factor(metsim$`Reference Panel`)
metsim$Array=rep("WGS subset", nrow(metsim))

prop=rbind(ho,metsim)
prop$Measure=rep("proportion", nrow(prop))

ho=fread("../analyses/omniexp/metsim.humanomni.bi.snv.mean.tab")
ho=ho[c(which(abs(ho$MAF-0.0002)==min(abs(ho$MAF-0.0002))),which(abs(ho$MAF-0.0004)==min(abs(ho$MAF-0.0004))),which(abs(ho$MAF-0.0006)==min(abs(ho$MAF-0.0006))),which(abs(ho$MAF-0.001)==min(abs(ho$MAF-0.001))),which(abs(ho$MAF-0.0032)==min(abs(ho$MAF-0.0032))),which(abs(ho$MAF-0.01)==min(abs(ho$MAF-0.01))),which(abs(ho$MAF-0.0316)==min(abs(ho$MAF-0.0316))),which(abs(ho$MAF-0.1)==min(abs(ho$MAF-0.1))),which(abs(ho$MAF-0.5)==min(abs(ho$MAF-0.5)))),]
ho=melt(ho, id="MAF")
ho$`Reference Panel`=rep(0, nrow(ho))
ho$`Reference Panel`[which(ho$variable=="HO_1000G")]="1000G"
ho$`Reference Panel`[which(ho$variable=="HO_HRC")]="HRC"
ho$`Reference Panel`[which(ho$variable=="HO_TOPMed")]="TOPMed"
ho$`Reference Panel`=as.factor(ho$`Reference Panel`)
ho$Array=rep("Real array", nrow(ho))

metsim=fread("../analyses/percentileaggr2/metsim.bi.snv.mean.tab")
metsim=metsim[,c("MAF","OmniExp_1000G","OmniExp_HRC","OmniExp_TOPMed")]
metsim=metsim[c(which(abs(metsim$MAF-0.0002)==min(abs(metsim$MAF-0.0002))),which(abs(metsim$MAF-0.0004)==min(abs(metsim$MAF-0.0004))),which(abs(metsim$MAF-0.0006)==min(abs(metsim$MAF-0.0006))),which(abs(metsim$MAF-0.001)==min(abs(metsim$MAF-0.001))),which(abs(metsim$MAF-0.0032)==min(abs(metsim$MAF-0.0032))),which(abs(metsim$MAF-0.01)==min(abs(metsim$MAF-0.01))),which(abs(metsim$MAF-0.0316)==min(abs(metsim$MAF-0.0316))),which(abs(metsim$MAF-0.1)==min(abs(metsim$MAF-0.1))),which(abs(metsim$MAF-0.5)==min(abs(metsim$MAF-0.5)))),]
metsim=melt(metsim, id="MAF")
metsim$`Reference Panel`=rep(0, nrow(metsim))
metsim$`Reference Panel`[which(metsim$variable=="Core_1000G" | metsim$variable=="OmniExp_1000G" | metsim$variable=="Omni25_1000G" | metsim$variable=="MEGA_1000G")]="1000G"
metsim$`Reference Panel`[which(metsim$variable=="Core_HRC" | metsim$variable=="OmniExp_HRC" | metsim$variable=="Omni25_HRC"| metsim$variable=="MEGA_HRC")]="HRC"
metsim$`Reference Panel`[which(metsim$variable=="Core_TOPMed" | metsim$variable=="OmniExp_TOPMed" | metsim$variable=="Omni25_TOPMed" | metsim$variable=="MEGA_TOPMed")]="TOPMed"
metsim$`Reference Panel`=as.factor(metsim$`Reference Panel`)
metsim$Array=rep("WGS subset", nrow(metsim))

mean=rbind(ho, metsim)
mean$Measure=rep("mean", nrow(mean))

##Plotting
all=rbind(prop, mean)
all$`Reference Panel`=factor(all$`Reference Panel`, levels=c("TOPMed","HRC","1000G"))
all$Array=factor(all$Array, levels=c("Real array","WGS subset"))
ybrks<-c(0,0.2,0.4,0.6,0.8,1)
brks<-c(0.5, 0.1, 0.01,0.001,0.0001)
labs=c(expression(0.5),expression(10^-1), expression(10^-2), expression(10^-3),expression(10^-4))

a=ggplot(all[which(all$Measure=="proportion"),], aes(x=MAF,y=value, color=Array))+scale_x_log10(breaks=brks, labels=labs, limits=c(1e-4,0.5),expand=c(0,0))+scale_y_continuous(breaks=ybrks,limits=c(0,1), expand=c(0,0))+theme_bw()+facet_grid(`Reference Panel`~.)+theme(legend.position="bottom",legend.text=element_text(size=7),panel.grid.major=element_line(colour="#e7e7e7"),panel.spacing.y=unit(0.5, "lines"), strip.text=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=7),panel.border=element_rect(color="black"),strip.background=element_rect(fill="white",color=NA), legend.title=element_text(size=8), legend.key.width=unit(3,"lines"))+geom_line(aes(group=Array), size=0.8)+xlab("WGS Minor Allele Frequency")+ylab(expression(Proportion~of~well-imputed~SNVs))+guides(color=guide_legend(nrow=2))

b=ggplot(all[which(all$Measure=="mean"),],aes(x=MAF,y=value, color=Array))+scale_x_log10(breaks=brks, labels=labs, limits=c(1e-4,0.5),expand=c(0,0))+scale_y_continuous(breaks=ybrks,limits=c(0,1), expand=c(0,0))+theme_bw()+facet_grid(`Reference Panel`~.)+theme(legend.position="bottom",legend.text=element_text(size=7),panel.grid.major=element_line(colour="#e7e7e7"),panel.spacing.y=unit(0.5, "lines"), strip.text=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=7),panel.border=element_rect(color="black"),strip.background=element_rect(fill="white",color=NA), legend.title=element_text(size=8), legend.key.width=unit(3,"lines"))+geom_line(aes(group=Array), size=0.8)+xlab("WGS Minor Allele Frequency")+ylab(expression(Mean~observed~imputation~r^2))+guides(color=guide_legend(nrow=2))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend=get_legend(a)
a=a+theme(legend.position="none")
b=b+theme(legend.position="none")

tiff(file="~/FigS3.tiff", height=174, width=174, units="mm", res=300)
grid.arrange(ggplotGrob(a),ggplotGrob(b), legend, nrow=2, widths=c(1,1), heights=c(8,1), layout_matrix=rbind(c(1,2),c(3)))
dev.off()



