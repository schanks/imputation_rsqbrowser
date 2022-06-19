library(data.table)
library(ggplot2)
library(gridExtra)

biome=fread("../analyses/VEP/biome.vep.perc8.tab")
colnames(biome)[1]="IMPACT"
biome=biome[,c("MAF","IMPACT","OmniExp_TOPMed")]
high=biome[which(biome$IMPACT=="IMPACT=HIGH"),]
moderate=biome[which(biome$IMPACT=="IMPACT=MODERATE"),]
low=biome[which(biome$IMPACT=="IMPACT=LOW"),]
modifier=biome[which(biome$IMPACT=="IMPACT=MODIFIER"),]
high=high[c(which(abs(high$MAF-0.0002)==min(abs(high$MAF-0.0002))),which(abs(high$MAF-0.0004)==min(abs(high$MAF-0.0004))),which(abs(high$MAF-0.0006)==min(abs(high$MAF-0.0006))),which(abs(high$MAF-0.001)==min(abs(high$MAF-0.001))),which(abs(high$MAF-0.0032)==min(abs(high$MAF-0.0032))),which(abs(high$MAF-0.01)==min(abs(high$MAF-0.01))),which(abs(high$MAF-0.0316)==min(abs(high$MAF-0.0316))),which(abs(high$MAF-0.1)==min(abs(high$MAF-0.1))),which(abs(high$MAF-0.316)==min(abs(high$MAF-0.316))),which(abs(high$MAF-0.5)==min(abs(high$MAF-0.5)))),]
moderate=moderate[c(which(abs(moderate$MAF-0.0002)==min(abs(moderate$MAF-0.0002))),which(abs(moderate$MAF-0.0004)==min(abs(moderate$MAF-0.0004))),which(abs(moderate$MAF-0.0006)==min(abs(moderate$MAF-0.0006))),which(abs(moderate$MAF-0.001)==min(abs(moderate$MAF-0.001))),which(abs(moderate$MAF-0.0032)==min(abs(moderate$MAF-0.0032))),which(abs(moderate$MAF-0.01)==min(abs(moderate$MAF-0.01))),which(abs(moderate$MAF-0.0316)==min(abs(moderate$MAF-0.0316))),which(abs(moderate$MAF-0.1)==min(abs(moderate$MAF-0.1))),which(abs(moderate$MAF-0.316)==min(abs(moderate$MAF-0.316))),which(abs(moderate$MAF-0.5)==min(abs(moderate$MAF-0.5)))),]
low=low[c(which(abs(low$MAF-0.0002)==min(abs(low$MAF-0.0002))),which(abs(low$MAF-0.0004)==min(abs(low$MAF-0.0004))),which(abs(low$MAF-0.0006)==min(abs(low$MAF-0.0006))),which(abs(low$MAF-0.001)==min(abs(low$MAF-0.001))),which(abs(low$MAF-0.0032)==min(abs(low$MAF-0.0032))),which(abs(low$MAF-0.01)==min(abs(low$MAF-0.01))),which(abs(low$MAF-0.0316)==min(abs(low$MAF-0.0316))),which(abs(low$MAF-0.1)==min(abs(low$MAF-0.1))),which(abs(low$MAF-0.316)==min(abs(low$MAF-0.316))),which(abs(low$MAF-0.5)==min(abs(low$MAF-0.5)))),]
modifier=modifier[c(which(abs(modifier$MAF-0.0002)==min(abs(modifier$MAF-0.0002))),which(abs(modifier$MAF-0.0004)==min(abs(modifier$MAF-0.0004))),which(abs(modifier$MAF-0.0006)==min(abs(modifier$MAF-0.0006))),which(abs(modifier$MAF-0.001)==min(abs(modifier$MAF-0.001))),which(abs(modifier$MAF-0.0032)==min(abs(modifier$MAF-0.0032))),which(abs(modifier$MAF-0.01)==min(abs(modifier$MAF-0.01))),which(abs(modifier$MAF-0.0316)==min(abs(modifier$MAF-0.0316))),which(abs(modifier$MAF-0.1)==min(abs(modifier$MAF-0.1))),which(abs(modifier$MAF-0.316)==min(abs(modifier$MAF-0.316))),which(abs(modifier$MAF-0.5)==min(abs(modifier$MAF-0.5)))),]
biome=rbind(high, moderate, low, modifier)
biome$Ancestry=rep("Hispanic/Latino", nrow(biome))

inpsyght=fread("../analyses/VEP/inpsyght.vep.perc8.tab")
colnames(inpsyght)[1]="IMPACT"
inpsyght=inpsyght[,c("MAF","IMPACT","OmniExp_TOPMed")]
high=inpsyght[which(inpsyght$IMPACT=="IMPACT=HIGH"),]
moderate=inpsyght[which(inpsyght$IMPACT=="IMPACT=MODERATE"),]
low=inpsyght[which(inpsyght$IMPACT=="IMPACT=LOW"),]
modifier=inpsyght[which(inpsyght$IMPACT=="IMPACT=MODIFIER"),]
high=high[c(which(abs(high$MAF-0.0002)==min(abs(high$MAF-0.0002))),which(abs(high$MAF-0.0004)==min(abs(high$MAF-0.0004))),which(abs(high$MAF-0.0006)==min(abs(high$MAF-0.0006))),which(abs(high$MAF-0.001)==min(abs(high$MAF-0.001))),which(abs(high$MAF-0.0032)==min(abs(high$MAF-0.0032))),which(abs(high$MAF-0.01)==min(abs(high$MAF-0.01))),which(abs(high$MAF-0.0316)==min(abs(high$MAF-0.0316))),which(abs(high$MAF-0.1)==min(abs(high$MAF-0.1))),which(abs(high$MAF-0.316)==min(abs(high$MAF-0.316))),which(abs(high$MAF-0.5)==min(abs(high$MAF-0.5)))),]
moderate=moderate[c(which(abs(moderate$MAF-0.0002)==min(abs(moderate$MAF-0.0002))),which(abs(moderate$MAF-0.0004)==min(abs(moderate$MAF-0.0004))),which(abs(moderate$MAF-0.0006)==min(abs(moderate$MAF-0.0006))),which(abs(moderate$MAF-0.001)==min(abs(moderate$MAF-0.001))),which(abs(moderate$MAF-0.0032)==min(abs(moderate$MAF-0.0032))),which(abs(moderate$MAF-0.01)==min(abs(moderate$MAF-0.01))),which(abs(moderate$MAF-0.0316)==min(abs(moderate$MAF-0.0316))),which(abs(moderate$MAF-0.1)==min(abs(moderate$MAF-0.1))),which(abs(moderate$MAF-0.316)==min(abs(moderate$MAF-0.316))),which(abs(moderate$MAF-0.5)==min(abs(moderate$MAF-0.5)))),]
low=low[c(which(abs(low$MAF-0.0002)==min(abs(low$MAF-0.0002))),which(abs(low$MAF-0.0004)==min(abs(low$MAF-0.0004))),which(abs(low$MAF-0.0006)==min(abs(low$MAF-0.0006))),which(abs(low$MAF-0.001)==min(abs(low$MAF-0.001))),which(abs(low$MAF-0.0032)==min(abs(low$MAF-0.0032))),which(abs(low$MAF-0.01)==min(abs(low$MAF-0.01))),which(abs(low$MAF-0.0316)==min(abs(low$MAF-0.0316))),which(abs(low$MAF-0.1)==min(abs(low$MAF-0.1))),which(abs(low$MAF-0.316)==min(abs(low$MAF-0.316))),which(abs(low$MAF-0.5)==min(abs(low$MAF-0.5)))),]
modifier=modifier[c(which(abs(modifier$MAF-0.0002)==min(abs(modifier$MAF-0.0002))),which(abs(modifier$MAF-0.0004)==min(abs(modifier$MAF-0.0004))),which(abs(modifier$MAF-0.0006)==min(abs(modifier$MAF-0.0006))),which(abs(modifier$MAF-0.001)==min(abs(modifier$MAF-0.001))),which(abs(modifier$MAF-0.0032)==min(abs(modifier$MAF-0.0032))),which(abs(modifier$MAF-0.01)==min(abs(modifier$MAF-0.01))),which(abs(modifier$MAF-0.0316)==min(abs(modifier$MAF-0.0316))),which(abs(modifier$MAF-0.1)==min(abs(modifier$MAF-0.1))),which(abs(modifier$MAF-0.316)==min(abs(modifier$MAF-0.316))),which(abs(modifier$MAF-0.5)==min(abs(modifier$MAF-0.5)))),]
inpsyght=rbind(high, moderate, low, modifier)
inpsyght$Ancestry=rep("African", nrow(inpsyght))

mlof=fread("../analyses/VEP/mlof.vep.perc8.tab")
colnames(mlof)[1]="IMPACT"
mlof=mlof[,c("MAF","IMPACT","OmniExp_TOPMed")]
high=mlof[which(mlof$IMPACT=="IMPACT=HIGH"),]
moderate=mlof[which(mlof$IMPACT=="IMPACT=MODERATE"),]
low=mlof[which(mlof$IMPACT=="IMPACT=LOW"),]
modifier=mlof[which(mlof$IMPACT=="IMPACT=MODIFIER"),]
high=high[c(which(abs(high$MAF-0.0002)==min(abs(high$MAF-0.0002))),which(abs(high$MAF-0.0004)==min(abs(high$MAF-0.0004))),which(abs(high$MAF-0.0006)==min(abs(high$MAF-0.0006))),which(abs(high$MAF-0.001)==min(abs(high$MAF-0.001))),which(abs(high$MAF-0.0032)==min(abs(high$MAF-0.0032))),which(abs(high$MAF-0.01)==min(abs(high$MAF-0.01))),which(abs(high$MAF-0.0316)==min(abs(high$MAF-0.0316))),which(abs(high$MAF-0.1)==min(abs(high$MAF-0.1))),which(abs(high$MAF-0.316)==min(abs(high$MAF-0.316))),which(abs(high$MAF-0.5)==min(abs(high$MAF-0.5)))),]
moderate=moderate[c(which(abs(moderate$MAF-0.0002)==min(abs(moderate$MAF-0.0002))),which(abs(moderate$MAF-0.0004)==min(abs(moderate$MAF-0.0004))),which(abs(moderate$MAF-0.0006)==min(abs(moderate$MAF-0.0006))),which(abs(moderate$MAF-0.001)==min(abs(moderate$MAF-0.001))),which(abs(moderate$MAF-0.0032)==min(abs(moderate$MAF-0.0032))),which(abs(moderate$MAF-0.01)==min(abs(moderate$MAF-0.01))),which(abs(moderate$MAF-0.0316)==min(abs(moderate$MAF-0.0316))),which(abs(moderate$MAF-0.1)==min(abs(moderate$MAF-0.1))),which(abs(moderate$MAF-0.316)==min(abs(moderate$MAF-0.316))),which(abs(moderate$MAF-0.5)==min(abs(moderate$MAF-0.5)))),]
low=low[c(which(abs(low$MAF-0.0002)==min(abs(low$MAF-0.0002))),which(abs(low$MAF-0.0004)==min(abs(low$MAF-0.0004))),which(abs(low$MAF-0.0006)==min(abs(low$MAF-0.0006))),which(abs(low$MAF-0.001)==min(abs(low$MAF-0.001))),which(abs(low$MAF-0.0032)==min(abs(low$MAF-0.0032))),which(abs(low$MAF-0.01)==min(abs(low$MAF-0.01))),which(abs(low$MAF-0.0316)==min(abs(low$MAF-0.0316))),which(abs(low$MAF-0.1)==min(abs(low$MAF-0.1))),which(abs(low$MAF-0.316)==min(abs(low$MAF-0.316))),which(abs(low$MAF-0.5)==min(abs(low$MAF-0.5)))),]
modifier=modifier[c(which(abs(modifier$MAF-0.0002)==min(abs(modifier$MAF-0.0002))),which(abs(modifier$MAF-0.0004)==min(abs(modifier$MAF-0.0004))),which(abs(modifier$MAF-0.0006)==min(abs(modifier$MAF-0.0006))),which(abs(modifier$MAF-0.001)==min(abs(modifier$MAF-0.001))),which(abs(modifier$MAF-0.0032)==min(abs(modifier$MAF-0.0032))),which(abs(modifier$MAF-0.01)==min(abs(modifier$MAF-0.01))),which(abs(modifier$MAF-0.0316)==min(abs(modifier$MAF-0.0316))),which(abs(modifier$MAF-0.1)==min(abs(modifier$MAF-0.1))),which(abs(modifier$MAF-0.316)==min(abs(modifier$MAF-0.316))),which(abs(modifier$MAF-0.5)==min(abs(modifier$MAF-0.5)))),]
mlof=rbind(high, moderate, low, modifier)
mlof$Ancestry=rep("European", nrow(mlof))

metsim=fread("../analyses/VEP/metsim.vep.perc8.tab")
colnames(metsim)[1]="IMPACT"
metsim=metsim[,c("MAF","IMPACT","OmniExp_TOPMed")]
high=metsim[which(metsim$IMPACT=="IMPACT=HIGH"),]
moderate=metsim[which(metsim$IMPACT=="IMPACT=MODERATE"),]
low=metsim[which(metsim$IMPACT=="IMPACT=LOW"),]
modifier=metsim[which(metsim$IMPACT=="IMPACT=MODIFIER"),]
high=high[c(which(abs(high$MAF-0.0002)==min(abs(high$MAF-0.0002))),which(abs(high$MAF-0.0004)==min(abs(high$MAF-0.0004))),which(abs(high$MAF-0.0006)==min(abs(high$MAF-0.0006))),which(abs(high$MAF-0.001)==min(abs(high$MAF-0.001))),which(abs(high$MAF-0.0032)==min(abs(high$MAF-0.0032))),which(abs(high$MAF-0.01)==min(abs(high$MAF-0.01))),which(abs(high$MAF-0.0316)==min(abs(high$MAF-0.0316))),which(abs(high$MAF-0.1)==min(abs(high$MAF-0.1))),which(abs(high$MAF-0.316)==min(abs(high$MAF-0.316))),which(abs(high$MAF-0.5)==min(abs(high$MAF-0.5)))),]
moderate=moderate[c(which(abs(moderate$MAF-0.0002)==min(abs(moderate$MAF-0.0002))),which(abs(moderate$MAF-0.0004)==min(abs(moderate$MAF-0.0004))),which(abs(moderate$MAF-0.0006)==min(abs(moderate$MAF-0.0006))),which(abs(moderate$MAF-0.001)==min(abs(moderate$MAF-0.001))),which(abs(moderate$MAF-0.0032)==min(abs(moderate$MAF-0.0032))),which(abs(moderate$MAF-0.01)==min(abs(moderate$MAF-0.01))),which(abs(moderate$MAF-0.0316)==min(abs(moderate$MAF-0.0316))),which(abs(moderate$MAF-0.1)==min(abs(moderate$MAF-0.1))),which(abs(moderate$MAF-0.316)==min(abs(moderate$MAF-0.316))),which(abs(moderate$MAF-0.5)==min(abs(moderate$MAF-0.5)))),]
low=low[c(which(abs(low$MAF-0.0002)==min(abs(low$MAF-0.0002))),which(abs(low$MAF-0.0004)==min(abs(low$MAF-0.0004))),which(abs(low$MAF-0.0006)==min(abs(low$MAF-0.0006))),which(abs(low$MAF-0.001)==min(abs(low$MAF-0.001))),which(abs(low$MAF-0.0032)==min(abs(low$MAF-0.0032))),which(abs(low$MAF-0.01)==min(abs(low$MAF-0.01))),which(abs(low$MAF-0.0316)==min(abs(low$MAF-0.0316))),which(abs(low$MAF-0.1)==min(abs(low$MAF-0.1))),which(abs(low$MAF-0.316)==min(abs(low$MAF-0.316))),which(abs(low$MAF-0.5)==min(abs(low$MAF-0.5)))),]
modifier=modifier[c(which(abs(modifier$MAF-0.0002)==min(abs(modifier$MAF-0.0002))),which(abs(modifier$MAF-0.0004)==min(abs(modifier$MAF-0.0004))),which(abs(modifier$MAF-0.0006)==min(abs(modifier$MAF-0.0006))),which(abs(modifier$MAF-0.001)==min(abs(modifier$MAF-0.001))),which(abs(modifier$MAF-0.0032)==min(abs(modifier$MAF-0.0032))),which(abs(modifier$MAF-0.01)==min(abs(modifier$MAF-0.01))),which(abs(modifier$MAF-0.0316)==min(abs(modifier$MAF-0.0316))),which(abs(modifier$MAF-0.1)==min(abs(modifier$MAF-0.1))),which(abs(modifier$MAF-0.316)==min(abs(modifier$MAF-0.316))),which(abs(modifier$MAF-0.5)==min(abs(modifier$MAF-0.5)))),]
metsim=rbind(high, moderate, low, modifier)
metsim$Ancestry=rep("Finnish", nrow(metsim))

all=rbind(biome, inpsyght, metsim, mlof)
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))
all$IMPACT[which(all$IMPACT=="IMPACT=HIGH")]="HIGH"
all$IMPACT[which(all$IMPACT=="IMPACT=MODERATE")]="MODERATE"
all$IMPACT[which(all$IMPACT=="IMPACT=LOW")]="LOW"
all$IMPACT[which(all$IMPACT=="IMPACT=MODIFIER")]="MODIFIER"
ybrks<-c(0,0.2,0.4,0.6,0.8,1)
brks<-c(0.5, 0.1, 0.01,0.001,0.0001)
labs=c(expression(0.5),expression(10^-1), expression(10^-2), expression(10^-3),expression(10^-4))

tiff(file="~/FigS11.tiff", height=75, width=174, units="mm", res=300)
ggplot(all, aes(x=MAF, y=OmniExp_TOPMed, color=IMPACT))+scale_x_log10(breaks=brks, labels=labs, limits=c(1e-4,0.5),expand=c(0,0))+scale_y_continuous(breaks=ybrks,limits=c(0,1), expand=c(0,0))+theme_bw()+theme(legend.position="bottom",panel.grid.major=element_line(colour="#e7e7e7"),panel.spacing=unit(0.8, "lines"),axis.title=element_text(size=8),axis.text=element_text(size=7),panel.border=element_rect(color="black"), legend.key.width=unit(3, "lines"),legend.text=element_text(size=7),legend.title=element_text(size=8), strip.background=element_rect(fill="white", color=NA), strip.text=element_text(size=8))+geom_line(aes(group=IMPACT), size=0.6)+facet_grid(.~Ancestry)+scale_color_manual(values=c("#1b9e77", "#e6ab02","#7570b3","#e7298a"))+xlab("WGS Minor Allele Frequency")+ylab(expression(atop(Proportion~of~TOPMed~SNVs,well-imputed~(r^2>0.8))))+guides(color=guide_legend(nrow=2))
dev.off()


