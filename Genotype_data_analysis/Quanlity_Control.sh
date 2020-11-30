#Remove samples with age <20
plink --bfile $GENOTYPE_DATA_NAME --remove $YoungSampleList --make-bed --out $OUTPUT1 --noweb

#Remove samples with transplant kidney
plink --bfile $OUTPUT1 --remove $TransplantSampleList.txt --make-bed --out $OUTPUT2 --noweb

#Remove samples with inconsistent gender  
plink --bfile $OUTPUT2 --remove $MISGENDERSampleList.txt --make-bed --out $OUTPUT3 --noweb

#Remove monosnps
plink --bfile $OUTPUT3 --exclude $MONOSNPS.txt --make-bed -out $OUTPUT4

#Remove raresnps
plink --bfile $OUTPUT4 --exclude $RARESNPS.txt --make-bed -out $OUTPUT5

#Remove SNPs with Hardy-Weinberg equilibrium P <1×10-6
plink --bfile $OUTPUT5 --exclude $HARDY_FILTERED_SNPS.txt --make-bed -out $OUTPUT6 --noweb

#Reomve low call rate SNPs
plink --bfile $OUTPUT6  --exclude $LOW_CALL_RATE_SNPS.txt --make-bed -out $OUTPUT7 --noweb

#Remove low call rate samples
plink --bfile $OUTPUT7 --remove $LOW_CALL_RATE_SAMPLES.txt --make-bed --out $OUTPUT8 --noweb

#Remove SNPs associated with chemistry plate ID 
plink --bfile $OUTPUT8 --loop-assoc $BATCHES.list --assoc --pfilter 1e-3 --allow-no-sex --noweb
plink --bfile $OUTPUT8 --exclude $SNPs_PLATE_LIST.txt --make-bed -out $OUTPUT9 --noweb

#Remove contaminated samples （inbreeding coefficient cutoff: heterozygosity rate ± 3-fold standard deviations from the mean)
plink --bfile $OUTPUT9 --remove $CONTAMINATED_SAMPLE_LIST.txt --make-bed --out $OUTPUT10 --noweb

#Remove cryptic samples
plink --bfile $OUTPUT10 --genome --out $OUTPUT11 --noweb
plink --bfile cleanSamples --remove $RELATIVE_SAMPLES.txt --make-bed --out $OUTPUT12 --noweb

#Sex chromosome
plink --bfile $OUTPUT12 --chr X --make-bed --out chrx
plink --bfile $OUTPUT12 --chr XY --make-bed --out chrxy
plink --bfile $OUTPUT12 --chr Y --make-bed --out chry

#Remove SNPs on sex chromosomes
plink --bfile $OUTPUT13 --exclude $SNPs_on_SEX_CHROM --make-bed --out $OUTPUT14 --noweb

#Combine with KGP reference 
#Remove SNPs in the long LD region of hg19
plink --bfile $OUTPUT14 --make-set longLDRhg19.txt --write-set --out $OUTPUT15
plink --bfile $OUTPUT15 --exclude hild.set --make-bed --out $OUTPUT16
plink --bfile $OUTPUT16 --indep-pairwise 50 5 0.2
plink --bfile $OUTPUT16 --extract plink.prune.in --make-bed --out $OUTPUT17

#Extract shared SNPs in KGP
plink --bfile $OUTPUT17 --extract $SHARED.txt --make-bed --out $OUTPUT18
#Combine
plink --bfile $OUTPUT18 --bmerge $KGP.bed $KGP.bim $KGP.fam --recode --out $OUTPUT19


#PCA by EIGENSTRAT
convertf -p $TRANSFORMAT.conf
smartpca -p $SMARTPCA.conf
















