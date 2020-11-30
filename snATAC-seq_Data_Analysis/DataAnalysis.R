#Quality Control
library(SnapATAC)
library(ggplot2)
options(stringsAsFactors = FALSE);
snap.files=c(
#load snap files here
);
sample.names=c(
#add sample names here
);
barcode.files=c(
#barcode.files (output of cellranger ATAC)
);
x.sp.ls = lapply(seq(snap.files), function(i){
    createSnap(
        file=snap.files[i],
        sample=sample.names[i]
    );
  })
names(x.sp.ls)=sample.names;
barcode.ls = lapply(seq(snap.files), function(i){
    barcodes = read.csv(
        barcode.files[i], 
        head=TRUE
    );
    barcodes = barcodes[2:nrow(barcodes),];
    barcodes$logUMI = log10(barcodes$passed_filters + 1);
    barcodes$promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
    barcodes
  })
plots = lapply(seq(snap.files), function(i){
    p1 = ggplot(
        barcode.ls[[i]], 
        aes(x=logUMI, y=promoter_ratio)) + 
        geom_point(size=0.3, col="grey") +
        theme_classic()	+
        ggtitle(sample.names[[i]]) +
        ylim(0, 1) + xlim(0, 6) + 
        labs(x = "log10(UMI)", y="promoter ratio")
        p1+geom_hline(yintercept=0.25, linetype="dashed", color = "red") +
    geom_hline(yintercept=0.6, linetype="dashed", color = "red")+
    geom_vline(xintercept=3, linetype="dashed", color = "red")+
    geom_vline(xintercept=5, linetype="dashed", color = "red")

    })
pdf("$OUTPUTNAME.pdf")
plots
dev.off()
# for both datasets, we identify usable barcodes using [3-5] for log10(UMI) and [0.2-0.8] for promoter ratio as cutoff.
cutoff.logUMI.low = c(3, 3);
cutoff.logUMI.high = c(5, 5);
cutoff.FRIP.low = c(0.25, 0.25);
cutoff.FRIP.high = c(0.6, 0.6);
barcode.ls = lapply(seq(snap.files), function(i){
    barcodes = barcode.ls[[i]];
    idx = which(
        barcodes$logUMI >= cutoff.logUMI.low[i] & 
        barcodes$logUMI <= cutoff.logUMI.high[i] & 
        barcodes$promoter_ratio >= cutoff.FRIP.low[i] &
        barcodes$promoter_ratio <= cutoff.FRIP.high[i]
    );
    barcodes[idx,]
  });
x.sp.ls = lapply(seq(snap.files), function(i){
    barcodes = barcode.ls[[i]];
    x.sp = x.sp.ls[[i]];
    barcode.shared = intersect(x.sp@barcode, barcodes$barcode);
    x.sp = x.sp[match(barcode.shared, x.sp@barcode),];
    barcodes = barcodes[match(barcode.shared, barcodes$barcode),];
    x.sp@metaData = barcodes;
    x.sp
  })
x.sp.ls
barcode.ls = lapply(seq(snap.files), function(i){
    barcodes = barcode.ls[[i]];
      barcodes$mito_ratio=(barcodes$mitochondrial)/(barcodes$total);
      mito.cutoff = quantile(barcodes$mito_ratio[barcodes$mito_ratio>0], 0.95);
      print(mito.cutoff)
    idx = which(barcodes$mito_ratio <= mito.cutoff);
    barcodes[idx,]
  });
x.sp.ls = lapply(seq(snap.files), function(i){
    barcodes = barcode.ls[[i]];
    x.sp = x.sp.ls[[i]];
    barcode.shared = intersect(x.sp@barcode, barcodes$barcode);
    x.sp = x.sp[match(barcode.shared, x.sp@barcode),];
    barcodes = barcodes[match(barcode.shared, barcodes$barcode),];
    x.sp@metaData = barcodes;
    x.sp
  })
x.sp = Reduce(snapRbind, x.sp.ls);
x.sp@metaData["sample"] = x.sp@sample;
x.sp
table(x.sp@sample);
save(x.sp,file="filterMitoPromoter.RData")
x.sp = addBmatToSnap(x.sp, bin.size=5000);
x.sp = makeBinary(x.sp, mat="bmat");
save(x.sp,file="addBin.RData")
library(GenomicRanges);
black_list = read.table("wgEncodeHg19ConsensusSignalArtifactRegions.bed");
black_list$V1<-as.character(black_list$V1)
sub('^...','',black_list$V1)->black_list$V1
#black_list = read.table("hg19.blacklist.bed.gz");
black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
  );
idy = queryHits(
    findOverlaps(x.sp@feature, black_list.gr)
  );
if(length(idy) > 0){
    x.sp = x.sp[,-idy, mat="bmat"];
  };
x.sp

chr.exclude = seqlevels(x.sp@feature)[grep("hs37d5|GL|NC|MT|Y|X", seqlevels(x.sp@feature))];
c("M",chr.exclude)->chr.exclude
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){
    x.sp = x.sp[,-idy, mat="bmat"]
  };
x.sp

bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
hist(
    bin.cov[bin.cov > 0], 
    xlab="log10(bin cov)", 
    main="log10(Bin Cov)", 
    col="lightblue", 
    xlim=c(0, 5)
  );
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
idx = which(Matrix::rowSums(x.sp@bmat) > 1000);
x.sp = x.sp[idx,];
x.sp = runDiffusionMaps(
    obj=x.sp,
    input.mat="bmat", 
    num.eigs=50
  );

pdf("$OUTPUTNAME.pdf")
plotDimReductPW(
    obj=x.sp, 
    eigs.dims=1:50,
    point.size=0.3,
    point.color="grey",
    point.shape=19,
    point.alpha=0.6,
    down.sample=5000,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7
  );
dev.off()
#Batch effect correction
library(harmony)
x.after.sp = runHarmony(
    obj=x.sp, 
    eigs.dim=1:22, 
    meta_data=x.sp@sample # sample index
  );

x.after.sp = runKNN(
    obj=x.after.sp,
    eigs.dims=1:22,
    k=18
  );
library(leiden);
x.after.sp=runCluster(
    obj=x.after.sp,
    tmp.folder=tempdir(),
    louvain.lib="leiden",
    seed.use=10,
    resolution=0.5
  );
x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:22,
    k=18
  );
library(leiden);
x.sp=runCluster(
    obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="leiden",
    seed.use=10,
    resolution=0.5
  );
library(umap);
x.after.sp = runViz(
    obj=x.after.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:22, 
    method="umap",
    seed.use=10
  );
x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:22, 
    method="umap",
    seed.use=10
  );
pdf("$OUTPUTNAME.pdf",width=12)
par(mfrow = c(1, 2));
plotViz(
    obj= x.sp,
    method="umap", 
    main="Before Harmony",
    point.size=0.2, 
    point.shape=19,
    point.alpha=0.2, 
    point.color=x.sp@sample, 
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=10000,
    legend.add=TRUE
    );
plotViz(
    obj= x.after.sp,
    method="umap", 
    main="After Harmony",
    point.size=0.2, 
    point.shape=19,
    point.alpha=0.2, 
    point.color=x.after.sp@sample, 
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=10000,
    legend.add=TRUE
    );
dev.off()
pdf("$OUTPUTNAME.pdf")
plotViz(
    obj= x.after.sp,
    method="umap", 
    main="Cluster",
    point.color=x.after.sp@cluster, 
    point.size=0.2, 
    point.shape=19, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  );
dev.off()
library(ggplot2)
pdf("$OUTPUTNAME.pdf")
ggplot(data=samples.after,aes(x=Cluster,y=Counts,fill=as.character(Sample),group=Sample))+geom_bar(stat="identity", position=position_dodge())+theme_bw()+ggtitle("After Harmony")
dev.off()
rbind(sample1.before,sample2.before)->samples.before
as.data.frame(samples.before)->samples.before
colnames(samples.before)<-c("Counts","Sample","Cluster")
as.numeric(as.character(samples.before$Counts))->samples.before$Counts
as.factor(as.numeric(as.character(samples.before$Cluster)))->samples.before$Cluster
library(ggplot2)
pdf("$OUTPUTNAME.pdf")
ggplot(data=samples.before,aes(x=Cluster,y=Counts,fill=as.character(Sample),group=Sample))+geom_bar(stat="identity", position=position_dodge())+theme_bw()+ggtitle("Before Harmony")
dev.off()
pdf("$OUTPUTNAME.pdf")
plotViz(
    obj= x.after.sp,
    method="umap", 
    main="Cluster",
    point.color=x.after.sp@cluster, 
    point.size=0.2, 
    point.shape=19, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  );
dev.off()
pdf("$OUTPUTNAME.pdf")
plotFeatureSingle(
    obj=x.sp,
    feature.value=x.sp@metaData[,"logUMI"],
    method="umap", 
    main="Read Depth",
    point.size=0.2, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0.01, 0.99)
  );
dev.off()
pdf("$OUTPUTNAME.pdf")
par(mfrow = c(2, 2));
plotFeatureSingle(
    obj=x.after.sp,
    feature.value=x.after.sp@metaData[,"logUMI"],
    method="umap", 
    main="Read Depth",
    point.size=0.2, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0.01, 0.99)
  );
plotFeatureSingle(
    obj=x.after.sp,
    feature.value=x.after.sp@metaData$peak_region_fragments / x.after.sp@metaData$passed_filters,
    method="umap", 
    main="Kidney FRiP",
    point.size=0.2, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0.01, 0.99) # remove outliers
  );
  plotFeatureSingle(
    obj=x.after.sp,
    feature.value=x.after.sp@metaData$mitochondrial / x.after.sp@metaData$total,
    method="umap", 
    main="Kidney Mito%",
    point.size=0.2, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0.01, 0.99) # remove outliers
  );
plotFeatureSingle(
    obj=x.after.sp,
    feature.value=x.after.sp@metaData$duplicate / x.after.sp@metaData$total,
    method="umap", 
    main="Kidney Duplicate",
    point.size=0.2, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0.01, 0.99) # remove outliers
  );
dev.off()
ensemble.ls = lapply(split(seq(length(x.after.sp@cluster)), x.after.sp@cluster), function(x){
	SnapATAC::colMeans(x.after.sp[x,], mat="bmat");
	})
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
pdf("$OUTPUTNAME.pdf")
plot(hc, hang=-1, xlab="");
dev.off()
library(GenomicRanges);
genes = read.table("~/hg19sc.gene.bed");
genes.gr = GRanges(genes[,1], 
    IRanges(genes[,2], genes[,3]), name=genes[,4]
  );
#marker genes
marker.genes = c(
    "NRP1", "KDR", "NPHS1","NPHS2", "SLC27A2", "LRP2","SLC12A1","UMOD","C1QA","C1QB","S100A8", "S100A9","CD79A","CD79B","LTB","CXCR6","GZMA","NKG7","STMN1",
    "SLC12A3", "PVALB", "AQP2","HSD11B2","ATP6V1G3","ATP6V0D2","INSRR","RHBG","MKI67","CDCA3","PLAC8","S100A4","TRPV5"
  );
genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)];
x.after.sp = addBmatToSnap(x.after.sp);
x.after.sp = createGmatFromMat(
    obj=x.after.sp, 
    input.mat="bmat",
    genes=genes.sel.gr,
    do.par=TRUE,
    num.cores=10
  );
x.after.sp = scaleCountMatrix(
    obj=x.after.sp, 
    cov=x.after.sp@metaData$passed_filters + 1,
    mat="gmat",
    method = "RPM"
  );
x.after.sp = runMagic(
    obj=x.after.sp,
    input.mat="gmat",
    step.size=3
  );
pdf("$OUTPUTNAME.pdf",width=10)
par(mfrow = c(2, 3));
for(i in 1:length(marker.genes)){
    plotFeatureSingle(
        obj=x.after.sp,
        feature.value=x.after.sp@gmat[, marker.genes[i]],
        method="umap", 
        main=marker.genes[i],
        point.size=0.1, 
        point.shape=19, 
        down.sample=10000,
        quantiles=c(0.01, 0.99)
  )};
dev.off()
clusters.sel = names(table(x.after.sp@cluster))[which(table(x.after.sp@cluster) > 100)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[i]);
    runMACS(
        obj=x.after.sp[which(x.after.sp@cluster==clusters.sel[i]),], 
        output.prefix=paste0("HumanKidney.", gsub(" ", "_", clusters.sel)[i]),
        path.to.snaptools="$PATH_TO_SNAPTOOLS",
        path.to.macs="$PATH_TO_MACS2",
        gsize="hs", # mm, hs, etc
        buffer.size=500, 
        num.cores=10,
        macs.options="--nomodel --keep-dup all --shift 100 --ext 200 --qval 5e-2 -B --SPMR --call-summits --nolambda",
        tmp.folder="$PATH_TO_SAVE_TPM_FILES"
   );
 }, mc.cores=10);
