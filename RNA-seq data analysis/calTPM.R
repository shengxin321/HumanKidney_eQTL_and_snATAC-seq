# Here we used read counts (output of RSEM) and effective gene length (output of RSEM) to calculate the TPM

#load RNA-seq and their matching phenotype data
load("$FILENAME.RData")
library(Biobase)
library(GenomicFeatures)
library(edgeR)
#Check whether the rownames of phenotype matrix are consistent with the colnames of the read counts matrix
all(rownames($PHENOTYPE_MATRIX)==colnames($READ_COUNTS_MATRIX)))

#Round the read counts estimated by RSEM
$READ_COUNTS_MATRIX <- data.frame(apply($READ_COUNTS_MATRIX, 2, function(x) as.numeric(as.character(x))))
round($READ_COUNTS_MATRIX)->$READ_COUNTS_MATRIX
rownames($READ_COUNTS_MATRIX)<-rownames($READ_COUNTS_MATRIX)

#Separated the read count matrix and phenotype matrix of the pair-end and single-end RNA-seq data
$PHENOTYPE_MATRIX[$PHENOTYPE_MATRIX$Pairend==1,]->$PHENOTYPE_MATRIX_PE
$PHENOTYPE_MATRIX[$PHENOTYPE_MATRIX$Pairend==0,]->$PHENOTYPE_MATRIX_SE

$READ_COUNTS_MATRIX[,rownames($PHENOTYPE_MATRIX_PE)]->$READ_COUNTS_MATRIX_PE
$READ_COUNTS_MATRIX[,rownames($PHENOTYPE_MATRIX_SE)]->$READ_COUNTS_MATRIX_SE




