library(susieR)
options(stringsAsFactors = FALSE);
load("$DOSAGE_FILE.RData")
load("$SNP.RData")
t(exp)->Y
read.table("$SIGPAIRS.txt",header=T,sep="\t",row.name=1)->subinfo
rownames(exp)->genes
info<-NULL
for(j in 1:length(genes))
{subinfo[subinfo$gene %in% genes[j],]->left
unique(left$snps)->subsnp
t(dos)->X
all(rownames(X)==rownames(Y))
colnames(X)->snps
fitted<-susie(X,Y[,j],L=10,estimate_residual_variance=TRUE,estimate_prior_variance=FALSE,scaled_prior_variance=0.1,verbose=TRUE)
print(j)
print(fitted$sets)
if(length(fitted$sets$cs)!=0){
sets <- susie_get_cs(fitted,
                     X = X,
             coverage = 0.9,
                     min_abs_corr = 0.1)
ssnp<-NULL
spip<-NULL
for(i in 1:allres){
sets$cs[[i]]->pos
snps[pos]->tmpsnp
print(tmpsnp)
ssnp<-c(ssnp,tmpsnp)
}
cbind(ssnp,spip)->tmpinfo
as.data.frame(tmpinfo)->tmpinfo
if(dim(tmpinfo)[1]>0){tmpinfo$Gene<-genes[j]}
colnames(info)<-c("SNP","Pip","Gene")
