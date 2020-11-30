#snATAC-seq
python ldsc.py --h2 $sumstats_formatted_filename.gz --w-ld $KGP_LD_score_filename --ref-ld $Celltype_ATAC-seq_LDscore,$BASELINE_REF_LDSCORE --overlap-annot --frqfile $KGP.frq --out $OUTPUT --print-coefficients

#scRNA-seq
#cov_ldsc for multiethinic GWAS data, ldsc for european GWAS data
#get annoataion ldsc scores
python ldsc.py --bfile $KGP_PLINK_FILE --l2 --ld-wind-cm 20 --cov $KGPCOV_FILENAME --out $OUTPUT --annot $ANNOTATION.annot.gz

#ldsc
python cov-ldsc/ldsc.py --h2 $GWAS_TRAIT.sumstats.gz --w-ld $KGP_LD_Score_file --ref-ld $Celltype_scRNA-seq_Per_BIN_LDscore,$BASELINE_REF_LDSCORE --cov $KGPCOVFILE.txt --overlap-annot --frqfile $KGP.frq --out $OUTPUT --print-coefficients

