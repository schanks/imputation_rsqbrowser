library(data.table)
library(gamlss)

results=fread("../../metsim.bi.snv.tab")
results=results[which(results$CHR=="chr20"),]
results$MAF=as.numeric(results$AF)
results$MAF[which(results$MAF>0.5)]=1-results$MAF[which(results$MAF>0.5)]
results=results[which(results$MAF>0),]
results$Omni25_TOPMed[which(is.na(results$Omni25_TOPMed))]=0

indep=fread("metsim.bi.snv.indep.sites")
results=merge(results, indep, by=c("CHR","POS","REF","ALT","AF"))

gc=fread("chr20.gc")
sv=fread("chr20.sv")
seg=fread("chr20.segdup")
rep=fread("chr20.rep")
rec=fread("chr20.recomb.mean")
dist=fread("chr20.omni25")

results=merge(results, rep[,c(1,2,3,4,9)], by=c("CHR","POS","REF","ALT"))
results=merge(results, gc[,c(1,2,3,4,9)], by=c("CHR","POS","REF","ALT"))
results=merge(results, sv[,c(1,2,3,4,9)], by=c("CHR","POS","REF","ALT"))
results=merge(results, seg[,c(1,2,3,4,9)], by=c("CHR","POS","REF","ALT"))
results=merge(results, rec[,c(1,2,3,4,10)], by=c("CHR","POS","REF","ALT"))
results=merge(results, dist[,c(1,2,3,4,7)], by=c("CHR","POS","REF","ALT"))
results=results[which(!is.na(results$RECOMB_REG)),]

results=results[,c("Omni25_TOPMed","MAF","REP","GC","SV","SEG","RECOMB_REG","DIST")]

results$REP=scale(results$REP)
results$GC=scale(results$GC)
results$SV=scale(results$SV)
results$RECOMB_REG=scale(results$RECOMB_REG)
results$DIST=scale(results$DIST)

breaks=c(0,0.00025,0.0005,0.00075, 0.001,0.0032,0.01,0.032,0.1,0.5)
results$bin=cut(results$MAF, breaks=breaks, labels=FALSE)
results$MAF1=as.numeric(results$bin==1)
results$MAF2=as.numeric(results$bin==2)
results$MAF3=as.numeric(results$bin==3)
results$MAF4=as.numeric(results$bin==4)
results$MAF5=as.numeric(results$bin==5)
results$MAF6=as.numeric(results$bin==6)
results$MAF7=as.numeric(results$bin==7)
results$MAF8=as.numeric(results$bin==8)
results$MAF9=as.numeric(results$bin==9)

zoibresults=as.data.frame(matrix(c(rep("REP",4),rep("GC",4),rep("SV",4),rep("SEG",4),rep("RECOMB",4),rep("DIST",4),rep(c("Mu","Sigma","Nu","Tau"),6),rep(0, 4*6*4)), byrow=FALSE, nrow=6*4))
colnames(zoibresults)=c("Feature","Parameter","Beta","SE","P","NR2")

for (i in c(3:8)){
	fit=gamlss(reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8",colnames(results)[i]), "Omni25_TOPMed"),sigma.formula=reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8",colnames(results)[i]), "Omni25_TOPMed"), nu.formula=reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8",colnames(results)[i]), "Omni25_TOPMed"),tau.formula=reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8",colnames(results)[i]), "Omni25_TOPMed"), family=BEINF, data=results)
	coefs=as.data.frame(summary(fit))
	rsq=Rsq(fit)
	zoibresults[(i-2)*4-3,3]=as.numeric(coefs$Estimate[10])
	zoibresults[(i-2)*4-3,4]=as.numeric(coefs$`Std. Error`[10])
	zoibresults[(i-2)*4-3,5]=as.numeric(coefs$`Pr(>|t|)`[10])
	zoibresults[(i-2)*4-3,6]=rsq
	zoibresults[(i-2)*4-2,3]=as.numeric(coefs$Estimate[20])
        zoibresults[(i-2)*4-2,4]=as.numeric(coefs$`Std. Error`[20])
        zoibresults[(i-2)*4-2,5]=as.numeric(coefs$`Pr(>|t|)`[20])
        zoibresults[(i-2)*4-2,6]=rsq
	zoibresults[(i-2)*4-1,3]=as.numeric(coefs$Estimate[30])
        zoibresults[(i-2)*4-1,4]=as.numeric(coefs$`Std. Error`[30])
        zoibresults[(i-2)*4-1,5]=as.numeric(coefs$`Pr(>|t|)`[30])
        zoibresults[(i-2)*4-1,6]=rsq
	zoibresults[(i-2)*4,3]=as.numeric(coefs$Estimate[40])
        zoibresults[(i-2)*4,4]=as.numeric(coefs$`Std. Error`[40])
        zoibresults[(i-2)*4,5]=as.numeric(coefs$`Pr(>|t|)`[40])
        zoibresults[(i-2)*4,6]=rsq
	print(i)
}

write.table(zoibresults, "metsim.zoibresults.Omni25_TOPMed.coefs", sep="\t",quote=FALSE, row.names=FALSE)

fit_maf=gamlss(reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8"), "Omni25_TOPMed"),sigma.formula=reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8"), "Omni25_TOPMed"), nu.formula=reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8"), "Omni25_TOPMed"),tau.formula=reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8"), "Omni25_TOPMed"), family=BEINF, data=results)
fit_full=gamlss(reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8","REP","GC","SV","SEG","RECOMB_REG","DIST"), "Omni25_TOPMed"),sigma.formula=reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8","REP","GC","SV","SEG","RECOMB_REG","DIST"), "Omni25_TOPMed"), nu.formula=reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8","REP","GC","SV","SEG","RECOMB_REG","DIST"), "Omni25_TOPMed"),tau.formula=reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8","REP","GC","SV","SEG","RECOMB_REG","DIST"), "Omni25_TOPMed"), family=BEINF, data=results)
Rsq(fit_maf)
Rsq(fit_full)


