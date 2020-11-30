options(stringsAsFactors = FALSE);
library(mashr)
read.table("$SIG_BETA.txt",row.name=1,header=T,sep="\t")->beta
read.table("$SIG_SE.txt",row.name=1,header=T,sep="\t")->std
sig=mash_set_data(beta,std)
mSig = mash(sig, g=get_fitted_g(m), fixg=TRUE)
