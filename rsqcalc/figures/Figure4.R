library(data.table)
library(ggplot2)
library(gridExtra)

biome=fread("../analyses/genomicfeatures/biome.logit.Omni25_TOPMed.coefs")
biome$Ancestry=rep("Hispanic/Latino", nrow(biome))
inpsyght=fread("../analyses/genomicfeatures/inpsyght.logit.Omni25_TOPMed.coefs")
inpsyght$Ancestry=rep("African", nrow(inpsyght))
metsim=fread("../analyses/genomicfeatures/metsim.logit.Omni25_TOPMed.coefs")
metsim$Ancestry=rep("Finnish", nrow(metsim))
mlof=fread("../analyses/genomicfeatures/mlof.logit.Omni25_TOPMed.coefs")
mlof$Ancestry=rep("European", nrow(mlof))

all=rbind(biome, mlof, inpsyght, metsim)
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))

ybrks=c(1,2,3,4,5,6)
ylabs=c("Recombination\nrate","GC content","Distance to\narray marker","Structural\nvariants","Segmental\nduplication","Repeats")

#% variance explained calculated in ../analyses/genomicfeatures/logit_feat.R
rsq=as.data.frame(matrix(c("MAF only",0.2832549,"African","MAF and all features",0.28379140,"African","MAF only",0.3029321,"Hispanic/Latino","MAF and all features",0.3034948,"Hispanic/Latino","MAF only",0.1675536,"European","MAF and all features",0.1686188,"European","MAF only",0.3745694,"Finnish","MAF and all features",0.3757176,"Finnish"), byrow=TRUE, ncol=3))
colnames(rsq)=c("Model","value","Ancestry")
rsq$value=as.numeric(rsq$value)
rsq$Model=factor(rsq$Model, levels=c("MAF only","MAF and all features"))
rsq$Ancestry=factor(rsq$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))

b=ggplot(rsq, aes(x=Ancestry, alpha=Model,y=value, fill=Ancestry))+ylab("Proportion of variance explained")+theme_bw()+scale_fill_manual(values=c("#e41a1c","#4daf4a","#377eb8","#984ea3"))+geom_bar(stat="identity", position="dodge")+scale_alpha_manual(values=c(0.7,1))+theme(axis.text.y=element_text(size=7), axis.title=element_text(size=8), legend.title=element_text(size=8), legend.text=element_text(size=6),axis.text.x=element_blank(), axis.ticks.x=element_blank())

all$Feature=factor(all$Feature, levels=c("RECOMB","GC","DIST","SV","SEG","REP"))

all$yval=as.numeric(as.factor(all$Feature))
all$yval[which(all$Ancestry=="Hispanic/Latino")]=all$yval[which(all$Ancestry=="Hispanic/Latino")]+0.1
all$yval[which(all$Ancestry=="African")]=all$yval[which(all$Ancestry=="African")]+0.3
all$yval[which(all$Ancestry=="European")]=all$yval[which(all$Ancestry=="European")]-0.1
all$yval[which(all$Ancestry=="Finnish")]=all$yval[which(all$Ancestry=="Finnish")]-0.3

a=ggplot(all, aes(y=yval))+geom_vline(xintercept=1)+geom_point(aes(x=OR, col=Ancestry))+geom_segment(aes(y=yval, yend=yval, x=LB, xend=UB, col=Ancestry))+scale_x_log10()+scale_y_continuous(breaks=ybrks, labels=ylabs)+theme_bw()+xlab("Odds Ratio")+ylab("Feature")+theme(panel.grid.minor=element_blank(),legend.position="none",strip.background=element_rect(fill="white",color=NA), axis.title=element_text(size=8), axis.text.y=element_text(size=7), axis.text.x=element_text(size=6))+scale_color_manual(values=c("#e41a1c","#4daf4a","#377eb8","#984ea3"))

tiff(file="~/Fig4.tiff", height=90, width=174, units="mm", res=300)
grid.arrange(a,b, nrow=1, widths=c(1,1.5))
dev.off()


