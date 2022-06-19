library(data.table)
library(ggplot2)

biome=fread("../analyses/genomicfeatures/biome.zoibresults.Omni25_TOPMed.coefs")
biome$Ancestry=rep("Hispanic/Latino", nrow(biome))
inpsyght=fread("../analyses/genomicfeatures/inpsyght.zoibresults.Omni25_TOPMed.coefs")
inpsyght$Ancestry=rep("African", nrow(inpsyght))
metsim=fread("../analyses/genomicfeatures/metsim.zoibresults.Omni25_TOPMed.coefs")
metsim$Ancestry=rep("Finnish", nrow(metsim))
mlof=fread("../analyses/genomicfeatures/mlof.zoibresults.Omni25_TOPMed.coefs")
mlof$Ancestry=rep("European", nrow(mlof))

all=rbind(biome, mlof, inpsyght, metsim)
all$yval=as.numeric(as.factor(all$Feature))
all$yval[which(all$Ancestry=="Hispanic/Latino")]=all$yval[which(all$Ancestry=="Hispanic/Latino")]+0.1
all$yval[which(all$Ancestry=="African")]=all$yval[which(all$Ancestry=="African")]+0.3
all$yval[which(all$Ancestry=="European")]=all$yval[which(all$Ancestry=="European")]-0.1
all$yval[which(all$Ancestry=="Finnish")]=all$yval[which(all$Ancestry=="Finnish")]-0.3
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino","European","Finnish"))

ybrks=c(1,2,3,4,5,6)
ylabs=c("Distance to\narray marker","GC content","Recombination\nrate","Repeats","Segmental\nduplication","Structural\nvariants")

all$Parameter=tolower(all$Parameter)

tiff(file="~/FigS10.tiff", height=174, width=114, units="mm", res=300)
ggplot(all, aes(y=yval))+geom_vline(xintercept=0)+geom_point(aes(x=Beta, col=Ancestry))+geom_segment(aes(y=yval, yend=yval, x=Beta-1.96*SE, xend=Beta+1.96*SE, col=Ancestry))+facet_grid(Parameter~., scales="free", labeller=label_parsed)+scale_y_continuous(breaks=ybrks, labels=ylabs)+theme_bw()+scale_color_manual(values=c("#e41a1c","#4daf4a","#377eb8","#984ea3"))+theme(legend.text=element_text(size=7), legend.title=element_text(size=8), strip.text=element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=7),panel.grid.minor=element_blank(),strip.background=element_rect(fill="white",color=NA))+ylab("Predictor")
dev.off()




