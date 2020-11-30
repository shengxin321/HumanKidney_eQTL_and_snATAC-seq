library("BSgenome.Hsapiens.scRef.hg19")
library(chromVAR);
library(motifmatchr);
library(SummarizedExperiment);
library(SnapATAC)
library(ggplot2)
library(GenomicRanges);
library(BSgenome)
options(stringsAsFactors = FALSE);
library(grid)
library(Seurat)
library(Matrix)
library(GenomicFiles)
library(GenomicAlignments)
library(plot3D)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(Rtsne)
library(GenomicRanges)
library(org.Hs.eg.db)
library(gridExtra)
library(pheatmap)
library(dplyr)
library(Seurat)
library(patchwork)
options(stringsAsFactors = FALSE);
rse<-SummarizedExperiment(assays=list(counts=t(data.use)),rowRanges=clusterDARpeaks,colData=DataFrame(Cell_Type=1:nrow(data.use),depth=Matrix::rowSums(data.use)))
rse<-addGCBias(rse,genome=BSgenome.Hsapiens.scRef.hg19)
library("TFBSTools")
library(JASPAR2018)
opts <- list()
opts[["species"]] <- "Homo sapiens"
opts[["all"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
PFMatrixList->motifs
species="Homo sapiens"
#motifs<-getJasparMotifs(collection="CORE",species=species)
motif_mm<-matchMotifs(motifs,rse,genome=BSgenome.Hsapiens.scRef.hg19)
save(motif_mm,file="motif2018DAR3.RData")
motif_pp<-matchMotifs(motifs,rse,genome=BSgenome.Hsapiens.scRef.hg19,out="positions")
motif_ss<-matchMotifs(motifs,rse,genome=BSgenome.Hsapiens.scRef.hg19,out="scores")
rownames(ordered[ordered$name %in% xuan,])->search
load("zscoreAllDAS3.RData")
zscore[search,]->info
colnames(zscore)->oldorder
all(oldorder==meta$barcode)
meta[,c("barcode","predicted.id")]->clusterinfo
clusterinfo[shun,]->new
info[,new[new$order %in% include,]$barcode]->gg
#info[,new$barcode]->gg
gg->info
as.character(ordered[rownames(info),]$name)->rownames(info)
rownames(info)
rownames(info)[duplicated(rownames(info))]
rownames(info)[duplicated(rownames(info))]->del
info[!duplicated(rownames(info)),]->info2
info2[xuan,]->info
