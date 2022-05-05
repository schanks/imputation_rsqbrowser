library(ggplot2)
library(gridExtra)
library(data.table)
library(scales)

reverselog_trans <- function(base = exp(1)) {
	trans <- function(x) -log(x, base)
	inv <- function(x) base^(-x)
	trans_new(paste0("reverselog-", format(base)), trans, inv,
		  log_breaks(base = base),
		  domain = c(1e-100, Inf))
}

#Panel A
irare=fread("inpsyght.tot.common.ind")
irare$Ancestry=rep("African", nrow(irare))
irare$MAF=rep("Rare", nrow(irare))
mrare=fread("mlof.tot.common.ind")
mrare$Ancestry=rep("European", nrow(mrare))
mrare$MAF=rep("Rare", nrow(mrare))
mtrare=fread("metsim.tot.common.ind")
mtrare$Ancestry=rep("Finnish",nrow(mtrare))
mtrare$MAF=rep("Rare", nrow(mtrare))
brare=fread("biome.tot.common.ind")
brare$Ancestry=rep("Hispanic/Latino", nrow(brare))
brare$MAF=rep("Rare",nrow(brare))
mtrare$V1=as.character(mtrare$V1)
rare=rbind(irare, mrare, mtrare, brare)
rare=rare[,c(1,11,12,13,18,19)]
colnames(rare)=c("ID","Ref","Het","Alt","Ancestry","MAF")
rare$Tot=rare$Ref+rare$Het+rare$Alt
rare$Discordance=1-rare$Het/(rare$Tot)
rare$Discordance_Quantile=rep("Q1", nrow(rare))
rare$Discordance_Quantile[which(rare$Discordance>quantile(rare$Discordance, 0.2,na.rm=TRUE) & rare$Discordance<quantile(rare$Discordance, 0.4, na.rm=TRUE))]="Q2"
rare$Discordance_Quantile[which(rare$Discordance>quantile(rare$Discordance, 0.4, na.rm=TRUE) & rare$Discordance<quantile(rare$Discordance, 0.6,na.rm=TRUE))]="Q3"
rare$Discordance_Quantile[which(rare$Discordance>quantile(rare$Discordance, 0.6, na.rm=TRUE) & rare$Discordance<quantile(rare$Discordance,0.8, na.rm=TRUE))]="Q4"
rare$Discordance_Quantile[which(rare$Discordance>quantile(rare$Discordance,0.8, na.rm=TRUE))]="Q5"
rare$Discordance_Quantile=factor(rare$Discordance_Quantile, levels=c("Q1","Q2","Q3","Q4","Q5"))
ipcs=fread("inpsyght.study_pc.txt")
colnames(ipcs)[1]="ID"
ipcs$Ancestry=rep("African", nrow(ipcs))
mpcs=fread("mlof.study_pc.txt")
colnames(mpcs)[1]="ID"
mpcs$Ancestry=rep("European", nrow(mpcs))
mtpcs=fread("metsim.study_pc.txt")
colnames(mtpcs)[1]="ID"
mtpcs$ID=as.character(mtpcs$ID)
mtpcs$indivID=as.character(mtpcs$indivID)
mtpcs$Ancestry=rep("Finnish", nrow(mtpcs))
bpcs=fread("biome.study_pc.txt")
colnames(bpcs)[1]="ID"
bpcs$Ancestry=rep("Hispanic/Latino", nrow(bpcs))
pcs=rbind(ipcs, mpcs, mtpcs, bpcs)
rare=merge(rare, pcs, by=c("ID","Ancestry"))
all=rbind(rare)
set.seed(1)
all=all[sample(nrow(all)),]
all$Concordance=1-all$Discordance
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))
all$`Concordance Quintile`=as.character(all$Discordance_Quantile)
all$`Concordance Quintile`[which(all$Discordance_Quantile=="Q1")]="Q5 (High)"
all$`Concordance Quintile`[which(all$Discordance_Quantile=="Q2")]="Q4"
all$`Concordance Quintile`[which(all$Discordance_Quantile=="Q4")]="Q2"
all$`Concordance Quintile`[which(all$Discordance_Quantile=="Q5")]="Q1 (Low)"
all$Panel=rep(1, nrow(all))

insetbreaks=c(0.9,0.95,1.0)

a=ggplot(all, aes(y=Concordance, x=Ancestry, fill=Ancestry))+geom_violin(aes(color=Ancestry))+geom_boxplot(coef=0,lwd=0.2, outlier.shape=NA, alpha=0)+scale_y_continuous(limits=c(0.3,1))+scale_color_manual(values=c("#cb181d","#238b45","#2171b5","#6a51a3"))+scale_fill_manual(values=c("#cb181d","#238b45","#2171b5","#6a51a3"))+ylab("Concordance")+theme_bw()+ggtitle("rates")+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"),plot.title=element_text(color="white",vjust=-55),legend.position="none",axis.title.x=element_text(color="white", size=8),axis.title.y=element_text(size=8), axis.text=element_text(size=6))+guides(color=guide_legend(nrow=2))

ainset=ggplot(all, aes(y=Concordance, x=Ancestry, fill=Ancestry))+geom_violin(aes(color=Ancestry))+geom_boxplot(coef=0,lwd=0.2, outlier.shape=NA, alpha=0)+scale_y_continuous(limits=c(0.9,1), breaks=insetbreaks)+scale_color_manual(values=c("#cb181d","#238b45","#2171b5","#6a51a3"))+scale_fill_manual(values=c("#cb181d","#238b45","#2171b5","#6a51a3"))+ylab("Concordance")+theme_bw()+theme(legend.position="none", axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(), axis.text.x=element_blank(),axis.text.y=element_text(size=6))+guides(color=guide_legend(nrow=2))

ainsetted=a+annotation_custom(ggplotGrob(ainset), ymin=0.3,ymax=0.8, xmin=2, xmax=4.5)

#Panel B
global=fread("inpsyght.percent_african.txt")
inpsyght=all[which(all$Ancestry=="African"),]
colnames(global)[1]="ID"
glob=merge(global, inpsyght, by="ID")
glob$afr=round(glob$african_percent, digits=1)
glob$`Proportion African ancestry`=as.factor(glob$afr)

b=ggplot(glob, aes(x=`Proportion African ancestry`, y=Concordance))+scale_y_continuous(limits=c(0.3,1))+ggtitle("African ancestry")+geom_violin(fill="#cb181d", color="#cb181d")+geom_boxplot(coef=0, lwd=0.2, outlier.shape=NA, alpha=0)+theme_bw()+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"),plot.title=element_text(vjust=-67,size=7),axis.text=element_text(size=6), axis.title.y=element_blank(), axis.title.x=element_text(size=8))+ylab("Concordance")

binset=ggplot(glob, aes(x=`Proportion African ancestry`, y=Concordance))+scale_y_continuous(limits=c(0.9,1), breaks=insetbreaks)+geom_violin(fill="#cb181d", color="#cb181d")+geom_boxplot(coef=0, lwd=0.2, outlier.shape=NA, alpha=0)+theme_bw()+theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=6), axis.title=element_blank())

binsetted=b+annotation_custom(ggplotGrob(binset), ymin=0.3, ymax=0.8, xmin=2.5, xmax=8.5)

#Panel C
pop=fread("biome.ancestry.txt")
map=fread("biome.map.txt")
pop=pop[,c("SUBJECT_ID","IBD_Community_Label","SUBCONTINENT")]
pop=merge(pop, map, by=c())
biome=all[which(all$Ancestry=="Hispanic/Latino"),]
colnames(pop)[4]="ID"
biome=merge(biome, pop, by="ID")
biome$Population=rep("Non-Caribbean", nrow(biome))
biome$Population[which(is.na(biome$SUBCONTINENT))]=NA
biome$Population[which(biome$SUBCONTINENT=="Caribbean")]="Caribbean"
biome$Population[which(biome$IBD_Community_Label=="Dominican" | biome$IBD_Community_Label=="PuertoRican")]="Caribbean"
biome$Population=factor(biome$Population, levels=c("Non-Caribbean","Caribbean"))

c=ggplot(biome[which(!is.na(biome$Population)),], aes(x=Population, y=Concordance))+scale_y_continuous(limits=c(0.3,1))+geom_violin(fill="#238b45", color="#238b45")+geom_boxplot(coef=0, lwd=0.2, outlier.shape=NA, alpha=0)+theme_bw()+ggtitle("Hispanic/Latino ancestry")+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"),axis.text=element_text(size=6), axis.title.y=element_blank(), plot.title=element_text(vjust=-67,size=7),axis.title.x=element_text(size=8))

cinset=ggplot(biome[which(!is.na(biome$Population)),], aes(x=Population, y=Concordance))+scale_y_continuous(limits=c(0.9,1), breaks=insetbreaks)+geom_violin(fill="#238b45", color="#238b45")+geom_boxplot(coef=0, lwd=0.2, outlier.shape=NA, alpha=0)+theme_bw()+theme(axis.ticks.x=element_blank(),axis.text.y=element_text(size=6), axis.text.x=element_blank(),axis.title=element_blank())

cinsetted=c+annotation_custom(ggplotGrob(cinset),ymin=0.3,ymax=0.8, xmin=1, xmax=2.5)

abc=grid.arrange(ainsetted,binsetted,cinsetted, nrow=1, widths=c(2.2,2.4,1.4))

#Panel D
d=ggplot(all, aes(x=PC1, y=PC2))+geom_point(aes(col=Ancestry),size=0.5)+theme_bw()+facet_grid(`Concordance Quintile`~Ancestry)+theme(strip.background=element_rect(fill="white", color=NA),strip.text=element_text(size=8),legend.position="none",axis.text=element_text(size=7), axis.title=element_text(size=8))+scale_color_manual(values=c("#cb181d","#238b45","#2171b5","#6a51a3"))

tiff(file="~/FigS6.tiff", height=174, width=174, units="mm", res=300)
grid.arrange(abc, d, nrow=2, heights=c(1,1.5))
dev.off()



