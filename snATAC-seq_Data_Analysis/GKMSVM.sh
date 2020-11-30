#Model training
gkmtrain -t 3 -c 10 -g 2 -x 10 $INPUT_POSITIVE.fa $INPUT_NEGATIVE.fa $CELLTYPE.gkmtrain -i 1

#gkm predict
gkmpredict 11mer.fa $CELLTYPE.gkmtrain.model.txt $CELLTYPE.gkmpredict.txt

#DeltaSVM
perl deltasvm.pl refsnp.fa altsnp.fa $CELLTYPE.gkmpredict.txt $CELLTYPE.Delta.output.txt

#Permutation
perl deltasvm.pl refsnp.fa altsnp.fa model$CELLTYPE.Permute.$PERMUTATION_COUNT.txt permout/model$CELLTYPE_per$PERMUTATION_COUNT_out.txt


