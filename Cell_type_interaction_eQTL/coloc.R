library(coloc)
read.table("$GWAS_SNP_POS_INFO",header=F)->gwassnps
read.table("$eQTL_INFO.txt",header=T)->eqtl
read.table("$GWAS_INFO.txt",header=T)->gwas
gwassnps$snp<-as.character(gwassnps$snp)
gwassnps$gwas<-as.character(gwassnps$gwas)
eqtl$SNP<-as.character(eqtl$SNP)
eqtl$MARKLOC<-as.character(eqtl$MARKLOC)
eqtl$eGene<-as.character(eqtl$eGene)
gwas$MARKLOC<-as.character(gwas$MARKLOC)
pp4<-NULL
num<-NULL
snpPmax<-NULL
bestSNP<-NULL
for(j in 1:dim(allconnect)[1]){
tpmgwas<-allconnect[j,]$gwas
tpmgene<-allconnect[j,]$egene
eqtl[eqtl$eGene %in% tpmgene,]->qtlInfo
gwasmark<-qtlInfo[qtlInfo$SNP %in% gwasNeed,]$MARKLOC
rownames(gwas1[grep("\\.",gwas1$SNP),])->del
gwas1[!(rownames(gwas1) %in% del),]->gwasF
gwasF[xu,]->gwasF
eqtlF[xu,]->eqtlF
all(rownames(eqtlF)==rownames(gwasF))
gwasF$SNP<-gwasF$MARKLOC
eqtlF$SNP<-eqtlF$MARKLOC
gwasF$Var<-gwasF$STD^2
eqtlF$Var<-eqtlF$STD^2
my.res<-coloc.abf(dataset1=list(beta=gwasF$BETA,varbeta=gwasF$Var,N=gwasF$N,type="quant"),dataset2=list(beta=eqtlF$BETA,varbeta=eqtlF$Var,N=356,type="quant"),MAF=gwasF$MAF)#coloc 
my.res$summary[6]->colocValue
SNPNum<-dim(eqtlF)[1]
t(as.matrix(my.res[[1]]))->abf
my.res[[2]]->snpsPP
max(snpsPP$SNP.PP.H4)->snpPPmax
snpsPP[snpsPP$SNP.PP.H4==snpPPmax,]$snp->snpid
strsplit(snpid,".",fixed=TRUE)[[1]][2]->ci
ci<-as.numeric(as.character(ci))
eqtlF$SNP[ci]->topSNP
}
allconnect$mark<-paste0(allconnect$gwas,",",allconnect$egene)
rownames(pp4)<-allconnect$mark
rownames(num)<-allconnect$mark
cbind(num,pp4,snpPmax,bestSNP)->res
colnames(res)<-c("NumSNPs","PP4","bestcolocSNP_PP4","bestcolocSNP")
