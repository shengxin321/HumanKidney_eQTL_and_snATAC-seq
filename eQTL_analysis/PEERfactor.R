library(peer)
#load phenotype data and TPM data matrix (INT_Transformed)
load("$QTLDATAFILENAME.RData")
$INT_TRANSFORMED_TPM_DAT_AMATRIX->dat.int
exp <- dat.int
covs<-$PHENOTYPE_DATA_MATRIX[,c("Pairend","Batch","Age","Gender","Site","RIN","PC1","PC2","PC3","PC4","PC5","Fibrosis","DCT","CNT","ALOH","BIC","DC","Granul","Treg","NK","Th17","CD8T","CD4T","Macro","PT")]
expr <- t(exp)
model = PEER()
PEER_setPhenoMean(model, as.matrix(expr))
PEER_setNk(model,50)
del<-c("Pairend","Gender","Site","Batch")
covs[,!(colnames(covs)%in% del)]->covleft
model.matrix(~Pairend+Gender+Site+Batch,data=covs)[,-1]->add2
summary($PHENOTYPE_DATA_MATRIX$Batch)<5->yes
c(1:25)->bat
paste0("Batch",bat[yes])->delbatch
summary($PHENOTYPE_DATA_MATRIX$Site)==0->yes2
c(1:7)->sit
paste0("Site",bat[yes2])->delsite
add2[,!(colnames(add2) %in% c(delsite,delbatch))]->add
cbind(covleft,add)->allcov
as.matrix(allcov)->finalcov
PEER_setCovariates(model, finalcov)
PEER_setAdd_mean(model, TRUE)
PEER_update(model)
factors = PEER_getX(model)
residuals = PEER_getResiduals(model)
dat <- t(residuals)
dat.int <- matrix(,nrow(dat),ncol(dat))
for (i in 1:nrow(dat)){
  dat.int[i,] <- qqnorm(dat[i,],plot.it=F)$x
}
dat.int<-data.frame(dat.int)
rownames(dat.int)<-rownames(exp)
colnames(dat.int)<-colnames(exp)
dat.int->$TPMEXP
