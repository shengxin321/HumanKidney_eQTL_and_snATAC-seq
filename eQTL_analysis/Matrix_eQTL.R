library(MatrixEQTL);
load("$FILENAME.RData");#load covariant file
useModel = modelLINEAR;
SNP_file_name = paste("$DOSAGE_FILE.txt", sep="");
snps_location_file_name = paste("$SNP_POSITION_FILE.txt", sep="");
expression_file_name = paste("$INT_TRANSFORMED_TPM_MATRIX.txt", sep="");
gene_location_file_name = paste("$GENE_POSITION_FILE.txt", sep="");
pvOutputThreshold_cis = 5e-2;
pvOutputThreshold_tra = 0;
errorCovariance = numeric();
cisDist = 1e6;
snps = SlicedData$new();
snps$fileDelimiter = "\t";
snps$fileOmitCharacters="NA";
snps$fileSkipRows=1;
snps$fileSkipColumns=1;
snps$fileSliceSize = 2000;
snps$LoadFile(SNP_file_name);
gene = SlicedData$new();
gene$fileDelimiter = "\t";      
gene$fileOmitCharacters = "NA"; 
gene$fileSkipRows = 1;          
gene$fileSkipColumns = 1;       
gene$fileSliceSize = 2000;      
gene$LoadFile(expression_file_name);
#save(gene,file="MvalueTrans.RData")
cvrt = SlicedData$new();
cvrt=covariate35;
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE,sep="\t");
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE,sep="\t");
me = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE, 
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = TRUE,
noFDRsaveMemory = FALSE);
