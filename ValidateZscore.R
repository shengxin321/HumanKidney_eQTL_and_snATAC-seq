options(stringsAsFactors = FALSE);
Args <- commandArgs()
ct<-Args[6]
read.table(paste0(ct,"/GREGOR/StatisticSummaryFile.txt"),header=T)->Stat
read.table(paste0(ct,"/GREGOR/index_SNP/annotated.index.snp.txt"),header=T)->snps
dim(snps)[1]->Stat$Total
pvalues<-NULL
std<-NULL
zscore<-NULL
expP<-NULL
for(i in 1:dim(Stat)[1]){
pval=Stat$ExpectNum_of_InBed_SNP[i]/Stat$Total[i]
obsP<-Stat$InBed_Index_SNP[i]/Stat$Total[i]
sz<-(sqrt(Stat$Total[i])*(obsP-pval))/sqrt(obsP*(1-obsP))
sp<-binom.test(x=Stat$InBed_Index_SNP[i],n=Stat$Total[i],p=pval)$p.val
sstd<-sqrt(Stat$Total[i]*pval*(1-pval))
expP<-rbind(expP,pval)
pvalues<-rbind(pvalues,sp)
std<-rbind(std,sstd)
zscore<-rbind(zscore,sz)
}
