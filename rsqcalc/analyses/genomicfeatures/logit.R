library(data.table)
library(fmsb)

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
class=fread("chr20.rep.class")

results=merge(results, rep[,c(1,2,3,4,9)], by=c("CHR","POS","REF","ALT"))
results=merge(results, gc[,c(1,2,3,4,9)], by=c("CHR","POS","REF","ALT"))
results=merge(results, sv[,c(1,2,3,4,9)], by=c("CHR","POS","REF","ALT"))
results=merge(results, seg[,c(1,2,3,4,9)], by=c("CHR","POS","REF","ALT"))
results=merge(results, rec[,c(1,2,3,4,10)], by=c("CHR","POS","REF","ALT"))
results=merge(results, dist[,c(1,2,3,4,7)], by=c("CHR","POS","REF","ALT"))
results=merge(results, class[,-c(5,6,7,8)], by=c("CHR","POS","REF","ALT"))
results=results[which(!is.na(results$RECOMB_REG)),]

results=results[,c("Omni25_TOPMed","MAF","REP","GC","SV","SEG","RECOMB_REG","DIST","DNA","SINE","LINE","LTR","Simple","LowComplex","Satellite","RNA","RC","Unknown")]

results$REP=scale(results$REP)
results$GC=scale(results$GC)
results$SV=scale(results$SV)
results$RECOMB_REG=scale(results$RECOMB_REG)
results$DIST=scale(results$DIST)

results$DNA=scale(results$DNA)
results$SINE=scale(results$SINE)
results$LINE=scale(results$LINE)
results$LTR=scale(results$LTR)
results$Simple=scale(results$Simple)
results$LowComplex=scale(results$LowComplex)
results$Satellite=scale(results$Satellite)
results$RNA=scale(results$RNA)
results$RC=scale(results$RC)
results$Unknown=scale(results$Unknown)

results$wellimp=as.numeric(results$Omni25_TOPMed>0.8)
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

logitresults=as.data.frame(matrix(c("REP","GC","SV","SEG","RECOMB","DIST",rep(0,7*6)), byrow=FALSE, nrow=6))
represults=as.data.frame(matrix(c("DNA","SINE","LINE","LTR","Simple","LowComplex","Satellite","RNA","RC","Unknown",rep(0,7*10)), byrow=FALSE, nrow=10))
colnames(logitresults)=c("Feature","Beta","OR","SE","LB","UB","P","NR2")
colnames(represults)=c("Feature","Beta","OR","SE","LB","UB","P","NR2")

for (i in c(3:8)){
	g=glm(reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8",colnames(results)[i]),"wellimp"), data=results, family="binomial")
	gsum=summary(g)
	logitresults[i-2,2]=as.numeric(gsum$coef[10,1])
	logitresults[i-2,3]=exp(as.numeric(gsum$coef[10,1]))
	logitresults[i-2,4]=as.numeric(gsum$coef[10,2])
	logitresults[i-2,5]=exp(as.numeric(gsum$coef[10,1])-1.96*as.numeric(gsum$coef[10,2]))
	logitresults[i-2,6]=exp(as.numeric(gsum$coef[10,1])+1.96*as.numeric(gsum$coef[10,2]))
	logitresults[i-2,7]=as.numeric(gsum$coef[10,4])
	logitresults[i-2,8]=as.numeric(NagelkerkeR2(g)$R2)
	print(i)
}

for (i in c(9:18)){
        g=glm(reformulate(c("MAF1","MAF2","MAF3","MAF4","MAF5","MAF6","MAF7","MAF8",colnames(results)[i]),"wellimp"), data=results, family="binomial")
        gsum=summary(g)
        represults[i-8,2]=as.numeric(gsum$coef[10,1])
        represults[i-8,3]=exp(as.numeric(gsum$coef[10,1]))
        represults[i-8,4]=as.numeric(gsum$coef[10,2])
        represults[i-8,5]=exp(as.numeric(gsum$coef[10,1])-1.96*as.numeric(gsum$coef[10,2]))
        represults[i-8,6]=exp(as.numeric(gsum$coef[10,1])+1.96*as.numeric(gsum$coef[10,2]))
        represults[i-8,7]=as.numeric(gsum$coef[10,4])
        represults[i-8,8]=as.numeric(NagelkerkeR2(g)$R2)
        print(i)
}

write.table(logitresults, "metsim.logit.Omni25_TOPMed.coefs", quote=FALSE, row.names=FALSE, sep="\t")

write.table(represults, "metsim.represults.Omni25_TOPMed.coefs", quote=FALSE, row.names=FALSE, sep="\t")

#Calculate proportion of variance explained for Figure 4B
g=glm(wellimp~MAF1+MAF2+MAF3+MAF4+MAF5+MAF6+MAF7+MAF8, data=results, family="binomial")
NagelkerkeR2(g)

g2=glm(wellimp~MAF1+MAF2+MAF3+MAF4+MAF5+MAF6+MAF7+MAF8+REP+GC+RECOMB_REG+SV+SEG+DIST, data=results, family="binomial")
NagelkerkeR2(g2)

