library(cicero)
options(stringsAsFactors = FALSE);
library(GenomicRanges);
library(chromVAR);
fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
input_cds <-  suppressWarnings(newCellDataSet(indata3,
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily=VGAM::binomialff(),
                            lowerDetectionLimit=0))
input_cds@expressionFamily@vfamily <- "binomialff"
input_cds <- monocle::detectGenes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 
paste0(x.after.sp@metaData$barcode,x.after.sp@metaData$sample)->barSam
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = myumap)
