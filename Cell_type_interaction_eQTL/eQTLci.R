options(stringsAsFactors = FALSE);
library(lme4)
library(pbkrtest)
load("$QTLDATA.RData");
load("$COVARIANTS.RData");
as.data.frame(t(cov2))->cov3
paste0("PEER",c(1:35))->subpeer
cov3[,subpeer]->peer
as.matrix(peer)->peer2
covs<-$PHENOTYPE_MATRIX[,c("Pairend","Batch","Age","Gender","Site","RIN","PC1","PC2","PC3","PC4","PC5","Fibrosis","DCT","CNT","ALOH","BIC","DC","Granul","Treg","NK","Th17","CD8T","CD4T","Macro","PT")]
allcovs<-matrix(as.numeric(unlist(cov2)),nrow=356)
read.table("$DOSAGE_FILE.txt",header=T,row.name=1)->snp
read.table("$INT_TRANSFORMED_TPM_FILE.txt",header=T,row.name=1)->gene
load("$OUTPUT_FROM_MATRIXEQTL.RData")
ordered_eqtls$gene<-as.character(ordered_eqtls$gene)
ordered_eqtls$snps<-as.character(ordered_eqtls$snps)
pvalues<-NULL
diff<-NULL
title<-NULL
for(j in 1:dim(gene)[1]){
rownames(gene[j,])->subgene
tpmexp<-t(gene[subgene,])
ordered_eqtls[ordered_eqtls$gene %in% subgene,]->tpmeqtl
dim(tpmeqtl)
for(m in 1:dim(tpmeqtl)[1]){
tpmeqtl[m,]$snps->subsnp
t(snp[subsnp,])->tpmdos
tpmmodel1<-lmer(tpmexp~tpmdos+(1|Batch)+Pairend+Age+(1|Gender)+(1|Site)+RIN+PC1+PC2+PC3+PC4+PC5+Fibrosis+DCT+CNT+ALOH+BIC+DC+Granul+Treg+NK+Th17+CD8T+CD4T+Macro+PT+peer2,data=covs,REML=FALSE)
summary(tpmmodel1)
tpmmodellarge<-lmer(tpmexp~PT:tpmdos+tpmdos+(1|Batch)+Pairend+Age+(1|Gender)+(1|Site)+RIN+PC1+PC2+PC3+PC4+PC5+Fibrosis+DCT+CNT+ALOH+BIC+DC+Granul+Treg+NK+Th17+CD8T+CD4T+Macro+PT+peer2,data=covs,REML=FALSE)
summary(tpmmodellarge)
anova(tpmmodellarge,tpmmodel1)->comp1
comp1$AIC[2]-comp1$AIC[1]->sdiff
compare<-PBmodcomp(tpmmodellarge,tpmmodel1,nsim=1000)
summary(compare)$test$p.value[3]->bpval
pvalues<-rbind(pvalues,bpval)
diff<-rbind(diff,sdiff)
stitle<-paste0(subgene,",",subsnp)
title<-rbind(title,stitle)
}
}
cbind(diff,pvalues)->res
rownames(res)<-title
colnames(res)<-c("AICDiff","chiSquarePvalue")
