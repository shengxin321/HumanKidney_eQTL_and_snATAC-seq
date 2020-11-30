#Prephasing, $i represents for chromosome ID, $j represents for each chunk
shapeit -B chr$i.filtered -M genetic_map_chr$i_combined_b37.txt --exclude-snp chr$i.alignments.snp.strand.exclude -O chr$i_noref.phased --thread 1

#Impute
impute2 -use_prephased_g -known_haps_g chr$i_noref.phased.haps -h 1000GP_Phase3_chr$i.hap.gz -l  1000GP_Phase3_chr$i.legend.gz -m genetic_map_chr$i_combined_b37.txt -int 0 5000000 -Ne 20000 -o chr$i_chunk$j.imputed

