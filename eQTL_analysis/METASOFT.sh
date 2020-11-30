java -jar Metasoft.jar -input $FILENAME.txt -mvalue -mvalue_method mcmc -pvalue_table HanEskinPvalueTable.txt -output $OUTPUT

#P-M plot
ython pmplot.py $INPUT.txt $OUTPUT.txt title.txt studyorder.txt $SNP $GENE $OUTPUT.pdf 
