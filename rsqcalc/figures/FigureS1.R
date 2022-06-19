library(data.table)
library(ggplot2)
library(reshape2)
library(gridExtra)

metsim=fread("../analyses/downsample/metsim.all_mean.tab")
metsim=metsim[c(which(abs(metsim$MAF-0.0002)==min(abs(metsim$MAF-0.0002))),which(abs(metsim$MAF-0.0004)==min(abs(metsim$MAF-0.0004))),which(abs(metsim$MAF-0.0006)==min(abs(metsim$MAF-0.0006))),which(abs(metsim$MAF-0.001)==min(abs(metsim$MAF-0.001))),which(abs(metsim$MAF-0.0032)==min(abs(metsim$MAF-0.0032))),which(abs(metsim$MAF-0.01)==min(abs(metsim$MAF-0.01))),which(abs(metsim$MAF-0.0316)==min(abs(metsim$MAF-0.0316))),which(abs(metsim$MAF-0.1)==min(abs(metsim$MAF-0.1))),which(abs(metsim$MAF-0.316)==min(abs(metsim$MAF-0.316))),which(abs(metsim$MAF-0.5)==min(abs(metsim$MAF-0.5)))),]
metsim=melt(metsim, id="MAF")
metsim$Ancestry=rep("Finnish", nrow(metsim))

biome=fread("../analyses/downsample/biome.all_mean.tab")
biome=biome[c(which(abs(biome$MAF-0.0002)==min(abs(biome$MAF-0.0002))),which(abs(biome$MAF-0.0004)==min(abs(biome$MAF-0.0004))),which(abs(biome$MAF-0.0006)==min(abs(biome$MAF-0.0006))),which(abs(biome$MAF-0.001)==min(abs(biome$MAF-0.001))),which(abs(biome$MAF-0.0032)==min(abs(biome$MAF-0.0032))),which(abs(biome$MAF-0.01)==min(abs(biome$MAF-0.01))),which(abs(biome$MAF-0.0316)==min(abs(biome$MAF-0.0316))),which(abs(biome$MAF-0.1)==min(abs(biome$MAF-0.1))),which(abs(biome$MAF-0.316)==min(abs(biome$MAF-0.316))),which(abs(biome$MAF-0.5)==min(abs(biome$MAF-0.5)))),]
biome=melt(biome, id="MAF")
biome$Ancestry=rep("Hispanic/Latino", nrow(biome))

inpsyght=fread("../analyses/downsample/inpsyght.all_mean.tab")
inpsyght=inpsyght[c(which(abs(inpsyght$MAF-0.0002)==min(abs(inpsyght$MAF-0.0002))),which(abs(inpsyght$MAF-0.0004)==min(abs(inpsyght$MAF-0.0004))),which(abs(inpsyght$MAF-0.0006)==min(abs(inpsyght$MAF-0.0006))),which(abs(inpsyght$MAF-0.001)==min(abs(inpsyght$MAF-0.001))),which(abs(inpsyght$MAF-0.0032)==min(abs(inpsyght$MAF-0.0032))),which(abs(inpsyght$MAF-0.01)==min(abs(inpsyght$MAF-0.01))),which(abs(inpsyght$MAF-0.0316)==min(abs(inpsyght$MAF-0.0316))),which(abs(inpsyght$MAF-0.1)==min(abs(inpsyght$MAF-0.1))),which(abs(inpsyght$MAF-0.316)==min(abs(inpsyght$MAF-0.316))),which(abs(inpsyght$MAF-0.5)==min(abs(inpsyght$MAF-0.5)))),]
inpsyght=melt(inpsyght, id="MAF")
inpsyght$Ancestry=rep("African", nrow(inpsyght))

mlof=fread("../analyses/downsample/mlof.all_mean.tab")
mlof=mlof[c(which(abs(mlof$MAF-0.0002)==min(abs(mlof$MAF-0.0002))),which(abs(mlof$MAF-0.0004)==min(abs(mlof$MAF-0.0004))),which(abs(mlof$MAF-0.0006)==min(abs(mlof$MAF-0.0006))),which(abs(mlof$MAF-0.001)==min(abs(mlof$MAF-0.001))),which(abs(mlof$MAF-0.0032)==min(abs(mlof$MAF-0.0032))),which(abs(mlof$MAF-0.01)==min(abs(mlof$MAF-0.01))),which(abs(mlof$MAF-0.0316)==min(abs(mlof$MAF-0.0316))),which(abs(mlof$MAF-0.1)==min(abs(mlof$MAF-0.1))),which(abs(mlof$MAF-0.316)==min(abs(mlof$MAF-0.316))),which(abs(mlof$MAF-0.5)==min(abs(mlof$MAF-0.5)))),]
mlof=melt(mlof, id="MAF")
mlof$Ancestry=rep("European", nrow(mlof))

all=rbind(biome, mlof, inpsyght,metsim)
all$`Sample Size`=as.factor(as.numeric(substr(as.character(all$variable),2,5)))
all=all[which(!is.element(all$`Sample Size`, c(4000,5000,6000))),]
ybrks<-c(0,0.2,0.4,0.6,0.8,1)
brks<-c(0.5, 0.1, 0.01,0.001,0.0001)
labs=c(expression(0.5),expression(10^-1), expression(10^-2), expression(10^-3),expression(10^-4))
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))

#Panel B
b=ggplot(all, aes(x=MAF, y=value, color=`Sample Size`))+scale_x_log10(breaks=brks, labels=labs, limits=c(1e-4,0.5),expand=c(0,0))+scale_y_continuous(breaks=ybrks,limits=c(0,1), expand=c(0,0))+theme_bw()+facet_grid(.~Ancestry)+theme(panel.grid.major=element_line(colour="#e7e7e7"),panel.grid.minor=element_blank(),strip.text=element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=7), panel.spacing.y=unit(0.7, "lines"), panel.spacing.x=unit(0.8,"lines"),panel.border=element_rect(color="black"),strip.background=element_rect(fill="white",color=NA),plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"), legend.position="none")+geom_line(aes(group=`Sample Size`), size=0.6)+xlab("WGS Minor Allele Frequency")+ylab(expression(atop(Mean~observed,paste(imputation~r^2))))+scale_color_manual(values=c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#abdda4","#66c2a5","#3288bd","#5e4fa2"))

#####
metsim=fread("../analyses/downsample/metsim.all_perc8.tab")
metsim=metsim[c(which(abs(metsim$MAF-0.0002)==min(abs(metsim$MAF-0.0002))),which(abs(metsim$MAF-0.0004)==min(abs(metsim$MAF-0.0004))),which(abs(metsim$MAF-0.0006)==min(abs(metsim$MAF-0.0006))),which(abs(metsim$MAF-0.001)==min(abs(metsim$MAF-0.001))),which(abs(metsim$MAF-0.0032)==min(abs(metsim$MAF-0.0032))),which(abs(metsim$MAF-0.01)==min(abs(metsim$MAF-0.01))),which(abs(metsim$MAF-0.0316)==min(abs(metsim$MAF-0.0316))),which(abs(metsim$MAF-0.1)==min(abs(metsim$MAF-0.1))),which(abs(metsim$MAF-0.316)==min(abs(metsim$MAF-0.316))),which(abs(metsim$MAF-0.5)==min(abs(metsim$MAF-0.5)))),]
metsim=melt(metsim, id="MAF")
metsim$Ancestry=rep("Finnish", nrow(metsim))

biome=fread("../analyses/downsample/biome.all_perc8.tab")
biome=biome[c(which(abs(biome$MAF-0.0002)==min(abs(biome$MAF-0.0002))),which(abs(biome$MAF-0.0004)==min(abs(biome$MAF-0.0004))),which(abs(biome$MAF-0.0006)==min(abs(biome$MAF-0.0006))),which(abs(biome$MAF-0.001)==min(abs(biome$MAF-0.001))),which(abs(biome$MAF-0.0032)==min(abs(biome$MAF-0.0032))),which(abs(biome$MAF-0.01)==min(abs(biome$MAF-0.01))),which(abs(biome$MAF-0.0316)==min(abs(biome$MAF-0.0316))),which(abs(biome$MAF-0.1)==min(abs(biome$MAF-0.1))),which(abs(biome$MAF-0.316)==min(abs(biome$MAF-0.316))),which(abs(biome$MAF-0.5)==min(abs(biome$MAF-0.5)))),]
biome=melt(biome, id="MAF")
biome$Ancestry=rep("Hispanic/Latino", nrow(biome))

inpsyght=fread("../analyses/downsample/inpsyght.all_perc8.tab")
inpsyght=inpsyght[c(which(abs(inpsyght$MAF-0.0002)==min(abs(inpsyght$MAF-0.0002))),which(abs(inpsyght$MAF-0.0004)==min(abs(inpsyght$MAF-0.0004))),which(abs(inpsyght$MAF-0.0006)==min(abs(inpsyght$MAF-0.0006))),which(abs(inpsyght$MAF-0.001)==min(abs(inpsyght$MAF-0.001))),which(abs(inpsyght$MAF-0.0032)==min(abs(inpsyght$MAF-0.0032))),which(abs(inpsyght$MAF-0.01)==min(abs(inpsyght$MAF-0.01))),which(abs(inpsyght$MAF-0.0316)==min(abs(inpsyght$MAF-0.0316))),which(abs(inpsyght$MAF-0.1)==min(abs(inpsyght$MAF-0.1))),which(abs(inpsyght$MAF-0.316)==min(abs(inpsyght$MAF-0.316))),which(abs(inpsyght$MAF-0.5)==min(abs(inpsyght$MAF-0.5)))),]
inpsyght=melt(inpsyght, id="MAF")
inpsyght$Ancestry=rep("African", nrow(inpsyght))

mlof=fread("../analyses/downsample/mlof.all_perc8.tab")
mlof=mlof[c(which(abs(mlof$MAF-0.0002)==min(abs(mlof$MAF-0.0002))),which(abs(mlof$MAF-0.0004)==min(abs(mlof$MAF-0.0004))),which(abs(mlof$MAF-0.0006)==min(abs(mlof$MAF-0.0006))),which(abs(mlof$MAF-0.001)==min(abs(mlof$MAF-0.001))),which(abs(mlof$MAF-0.0032)==min(abs(mlof$MAF-0.0032))),which(abs(mlof$MAF-0.01)==min(abs(mlof$MAF-0.01))),which(abs(mlof$MAF-0.0316)==min(abs(mlof$MAF-0.0316))),which(abs(mlof$MAF-0.1)==min(abs(mlof$MAF-0.1))),which(abs(mlof$MAF-0.316)==min(abs(mlof$MAF-0.316))),which(abs(mlof$MAF-0.5)==min(abs(mlof$MAF-0.5)))),]
mlof=melt(mlof, id="MAF")
mlof$Ancestry=rep("European", nrow(mlof))

all=rbind(biome, mlof, inpsyght,metsim)
all$`Sample Size`=as.factor(as.numeric(substr(as.character(all$variable),2,5)))
all=all[which(!is.element(all$`Sample Size`, c(4000,5000,6000))),]
ybrks<-c(0,0.2,0.4,0.6,0.8,1)
brks<-c(0.5, 0.1, 0.01,0.001,0.0001)
labs=c(expression(0.5),expression(10^-1), expression(10^-2), expression(10^-3),expression(10^-4))
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))

#Panel A
a=ggplot(all, aes(x=MAF, y=value, color=`Sample Size`))+scale_x_log10(breaks=brks, labels=labs, limits=c(1e-4,0.5),expand=c(0,0))+scale_y_continuous(breaks=ybrks,limits=c(0,1), expand=c(0,0))+theme_bw()+facet_grid(.~Ancestry)+theme(panel.grid.major=element_line(colour="#e7e7e7"),panel.grid.minor=element_blank(),panel.spacing.y=unit(0.7, "lines"), panel.spacing.x=unit(0.8,"lines"),strip.text=element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=7), panel.border=element_rect(color="black"),strip.background=element_rect(fill="white",color=NA), plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"), legend.margin=unit(c(0.5,0.5,0.5,0.5),"lines"), legend.key.width=unit(2, "lines"),legend.text=element_text(size=8), legend.title=element_text(size=8))+geom_line(aes(group=`Sample Size`), size=0.6)+xlab("WGS Minor Allele Frequency")+ylab(expression(atop(Proportion~of,well-imputed~SNVs)))+scale_color_manual(values=c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#abdda4","#66c2a5","#3288bd","#5e4fa2"))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend=get_legend(a)
a=a+theme(legend.position="none")
tiff(file="~/FigS1.tiff", height=100, width=174, units="mm", res=300)
grid.arrange(ggplotGrob(a),ggplotGrob(b), legend, layout_matrix=rbind(c(1,1,3),c(2,2,3)), widths=c(1,1,0.3))
dev.off()



