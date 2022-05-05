library(data.table)
library(ggplot2)
library(gridExtra)

biome=fread("../analyses/genomicfeatures/biome.represults.Omni25_TOPMed.coefs")
biome$Ancestry=rep("Hispanic/Latino", nrow(biome))
inpsyght=fread("../analyses/genomicfeatures/inpsyght.represults.Omni25_TOPMed.coefs")
inpsyght$Ancestry=rep("African", nrow(inpsyght))
metsim=fread("../analyses/genomicfeatures/metsim.represults.Omni25_TOPMed.coefs")
metsim$Ancestry=rep("Finnish", nrow(metsim))
mlof=fread("../analyses/genomicfeatures/mlof.represults.Omni25_TOPMed.coefs")
mlof$Ancestry=rep("European", nrow(mlof))

all=rbind(biome, mlof, inpsyght, metsim)
all$yval=as.numeric(as.factor(all$Feature))
all$yval[which(all$Ancestry=="Hispanic/Latino")]=all$yval[which(all$Ancestry=="Hispanic/Latino")]+0.1
all$yval[which(all$Ancestry=="African")]=all$yval[which(all$Ancestry=="African")]+0.3
all$yval[which(all$Ancestry=="European")]=all$yval[which(all$Ancestry=="European")]-0.1
all$yval[which(all$Ancestry=="Finnish")]=all$yval[which(all$Ancestry=="Finnish")]-0.3
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))

ybrks=seq(1,10)
ylabs=c("DNA","LINE","LowComplex","LTR","RC","RNA","Satellite","Simple","SINE","Unknown")

tiff(file="~/FigS9.tiff", height=114, width=114, units="mm", res=300)
ggplot(all, aes(y=yval))+geom_vline(xintercept=1)+geom_point(aes(x=OR, col=Ancestry))+geom_segment(aes(y=yval, yend=yval, x=LB, xend=UB, col=Ancestry))+scale_y_continuous(breaks=ybrks, labels=ylabs)+theme_bw()+xlab("Odds Ratio")+ylab("Repeat")+theme(strip.text=element_text(size=8), legend.title=element_text(size=8), legend.text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=8),panel.grid.minor=element_blank(),strip.background=element_rect(fill="white",color=NA))+scale_color_manual(values=c("#e41a1c","#4daf4a","#377eb8","#984ea3"))
dev.off()

