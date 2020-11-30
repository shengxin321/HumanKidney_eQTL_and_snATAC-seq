for j in $(seq 1 100); do
   echo "fastQTL.static --vcf $DOSAGE_FILE_NAME.vcf.gz --bed $INT_TRANSFORMED_GENE_EXP_MATRIX.bed.gz --permute 10000 --out permutations.chunk$j.txt.gz --window 1e6 --cov $COVARIATES.txt.gz --chunk $j 100"
done
