magma --bfile $KGP_PLINK_FILE synonyms=dbsnp151.synonyms --gene-annot $ANNOTATION_FILE --pval $GWAS_PVALUE.txt 'use=rsid,pval' ncol=N --out $OUTPUT

magma --gene-results Genes.$GWAS_TRAIT.genes.raw --gene-covar $CELL_TYPE_SPECIFICITY_RANK --out $OUTPUT
