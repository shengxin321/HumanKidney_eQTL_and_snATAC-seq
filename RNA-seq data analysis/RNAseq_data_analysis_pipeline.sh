# For our pair-end RNA-seq samples
# Quality Control
trim_galore --fastqc --retain_unpaired --paired $FILENAME_R1.fastq.gz $FILENAME_R2.fastq.gz -o .

# Alignment for our strand specific pair-end RNA-seq samples
STAR --runMode alignReads --runThreadN 1 --genomeDir $REFERENCE_TO_STAR_Index --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --readFilesIn $FILENAME_AFTER_QC_1.fq.gz $FILENAME_AFTER_QC_2.fq.gz --readFilesCommand zcat --outFileNamePrefix $OUTPUT_PREFIX --alignSoftClipAtReferenceEnds Yes --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted SortedByCoordinate --outSAMunmapped Within KeepPairs --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimOutType WithinBAM SoftClip --outSAMattributes NH HI AS nM NM MD jM jI XS --outSAMattrRGline ID:rg1 SM:sm1

# Quantification
rsem-calculate-expression --num-threads 1 --fragment-length-max 1000 --no-bam-output --paired-end --estimate-rspd --forward-prob 0 --bam $STAR_OUTPUT_FILENAME.toTranscriptome.out.bam $PATH_TO_RSEM_INDEX $OUTPUT_PREFIX


# For our single-end RNA-seq samples
# Quality Control
trim_galore --fastqc $FILENAME.fastq.gz -o $OUTPUT_FOLDER_NAME

# Alignment
STAR --runMode alignReads --runThreadN 1 --genomeDir $REFERENCE_TO_STAR_Index --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --readFilesIn $FILENAME_AFTER_QC.fq.gz --readFilesCommand zcat --outFileNamePrefix $OUTPUT_PREFIX --alignSoftClipAtReferenceEnds Yes --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted SortedByCoordinate --outSAMunmapped Within --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimOutType WithinBAM SoftClip --outSAMattributes NH HI AS nM NM MD jM jI XS --outSAMattrRGline ID:rg1 SM:sm1

# Quantification
rsem-calculate-expression --num-threads 1 --no-bam-output --estimate-rspd --bam $STAR_OUTPUT_FILENAME.toTranscriptome.out.bam $PATH_TO_RSEM_INDEX $OUTPUT_PREFIX
