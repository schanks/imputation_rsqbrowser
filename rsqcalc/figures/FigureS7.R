library(ggplot2)
library(data.table)
library(gridExtra)

irare=fread("inpsyght.tot.rare.ind")
irare$Ancestry=rep("African", nrow(irare))
mrare=fread("mlof.tot.rare.ind")
mrare$Ancestry=rep("European", nrow(mrare))
mtrare=fread("metsim.tot.rare.ind")
mtrare$Ancestry=rep("Finnish",nrow(mtrare))
brare=fread("biome.tot.rare.ind")
brare$Ancestry=rep("Hispanic/Latino", nrow(brare))
mtrare$V1=as.character(mtrare$V1)
rare=rbind(irare, mrare, mtrare, brare)
rare=rare[,c("V1","Ancestry")]
colnames(rare)[1]="ID"

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
rare=rare[,-c("L","K","t","Z")]
rare$Samples=rep("WGS", nrow(rare))

ref=fread("reference_pc.txt")
colnames(ref)[1]="ID"
ref$Samples=rep("HGDP: Other", nrow(ref))
ref$Samples[which(is.element(ref$ID, c("Adygei","Basque","French","Italian","Orcadian","Russian","Sardinian","Tuscan")))]="HGDP: Europe"
ref$Samples[which(is.element(ref$ID, c("BantuKenya","BantuSouthAfrica","BiakaPygmy","Mandenka","MbutiPygmy","Mozabite","San","Yoruba")))]="HGDP: Africa"
ref$Samples[which(is.element(ref$ID, c("Colombian","Karitiana","Maya","Pima","Surui")))]="HGDP: Native America"
ref=ref[which(ref$Samples!="HGDP: Other"),]
ref$Ancestry=rep("African", nrow(ref))

all=rbind(rare, ref)
ref$Ancestry=rep("Hispanic/Latino", nrow(ref))
all=rbind(all, ref)
ref$Ancestry=rep("European", nrow(ref))
all=rbind(all, ref)
ref$Ancestry=rep("Finnish", nrow(ref))
all=rbind(all, ref)
all$Ancestry=factor(all$Ancestry, levels=c("African","Hispanic/Latino", "European","Finnish"))
all$Samples=factor(all$Samples, levels=c("WGS","HGDP: Africa","HGDP: Europe", "HGDP: Native America"))

tiff(file="~/FigS7.tiff", height=45, width=174, units="mm", res=300)
ggplot(all, aes(x=PC1, y=PC2, color=Samples, shape=Samples))+geom_point(size=0.5)+theme_bw()+facet_grid(.~Ancestry)+scale_color_manual(values=c("gray","#ff7f00","deepskyblue2","forestgreen"))+scale_shape_manual(values=c(16,17,15,3))+theme(legend.key.size=unit(1, "lines"), strip.background=element_rect(fill="white",color=NA), strip.text=element_text(size=8), legend.title=element_text(size=8), legend.text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=8))+guides(color=guide_legend(override.aes=list(size=2)))
dev.off()

