#Remove samples with age <20
plink --bfile $GENOTYPE_DATA_NAME --remove $YoungSampleList --make-bed --out $OUTPUT1 --noweb

#Remove samples with transplant kidney
plink --bfile $OUTPUT1 --remove $TransplantSampleList.txt --make-bed --out $OUTPUT2 --noweb

#Remove samples with inconsistent gender  
plink --bfile $OUTPUT2 --remove $MISGENDERSampleList.txt --make-bed --out $OUTPUT3 --noweb

#Remove monosnps
plink --bfile $OUTPUT3 --exclude $MONOSNPS.txt --make-bed -out $OUTPUT4





