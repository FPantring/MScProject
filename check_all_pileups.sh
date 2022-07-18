#!/bin/bash

# This script lists the allele pileup at all True Negative Locations in both the "truth files" and the realigned bam file from BWA-MEM. 
# It also prints information (read name, true origin position according to the personalised coordinates and read sequence) for all reads incorrectly aligned to the locus
# and all reads originating from the locus that were incorrectly aligned to another location.
# The script needs to be submitted with paths and names of all required input files which are (in that order):
# "TN_and_FL.csv": a file produced by the script Result_Analysis.py which includes the names of all True Negative Locations
# two liftover files: one for each haplotype, which enable the liftover of the reference coordinates to the corresponding personalised coordinates used in the simulation
# "chr_list.txt": a list with line separated paths and names of "truth files" in a format similar to T_X1_chr1.bam, with X1 or X2 specifies the haplotype and chr1 the chromosome. 
# 		  There should be one file for each chromosome and each of the two haplotypes respectively
# the bam file
# "TN_FL_locs.txt": a file produced by the script Result_Analysis.py which lists the locus and alternative locus of all possible False Location variants, 
# 		    i.e. mapping quality variants that have a variant present at one of the alternative mapping locations. 
#		    These loci will be lifted-over to the personlised coordinates of both haplotypes at the end of the script.
# the name of the output file that will be produced (usually pileup_infos.txt)

# Thus, the script needs to be called with a command resembling: 
# sbatch check_pileups.sh TN_and_FL.csv liftover_X1_HG00110.condense.txt liftover_X2_HG00110.condense.txt chr_list.txt T_WGS.bam TN_FL_locs.txt pileup_infos.txt


# If a file with this name exists it has to be removed, since this file will later be appended with information
rm -f $7

cut -f1 $1 | tail -n +2 | while read loc; do
# liftover the locus of interest to the locus in the simulated genomes:
x1_loc=($(./2wayLiftover r $loc $2))
x2_loc=($(./2wayLiftover r $loc $3))
locnum=$( cut -d ":" -f2 <<< $loc ) #create a variable with the location without the chromosome ID before that
locchr=$( cut -d ":" -f1 <<< $loc )
echo $loc $locnum $locchr $(echo "X1_"$locchr"_") 

# get the IDs of all reads that were created at the locus of interest:
/data/fpantring/samtools mpileup --no-BAQ --output-QNAME $(grep $(echo "X1_"$locchr"_") $4) -r $(echo "${x1_loc[0]}"":""${x1_loc[1]}""-""${x1_loc[1]}") | cut -f7 | tr ',' '\012' | sort > gt.lst
/data/fpantring/samtools mpileup --no-BAQ --output-QNAME $(grep $(echo "X2_"$locchr"_") $4) -r $(echo "${x2_loc[0]}"":""${x2_loc[1]}""-""${x2_loc[1]}") | cut -f7 | tr ',' '\012' | sort >> gt.lst
# get the IDs of all reads that were allocated to the locus of interest:
/data/fpantring/samtools mpileup --no-BAQ --output-QNAME $5 -r $(echo "$loc""-""$locnum") | tr ',' '\012' | tail -n +2  | sort > called.lst

# check the allele pileups at the locus in the truth files and the bwa realigned file:
echo -e " \n"$loc" \nTrue locus allele pileup:" >> $7
/data/fpantring/samtools mpileup --no-BAQ $(grep $(echo "X1_"$locchr"_") $4) -r $(echo "${x1_loc[0]}"":""${x1_loc[1]}""-""${x1_loc[1]}") >> $7
/data/fpantring/samtools mpileup --no-BAQ $(grep $(echo "X2_"$locchr"_") $4) -r $(echo "${x2_loc[0]}"":""${x2_loc[1]}""-""${x2_loc[1]}") >> $7
echo -e " \n"$loc" \nLocus allele pileup according to bwa allocation:" >> $7
/data/fpantring/samtools mpileup --no-BAQ $5 -r $(echo "$loc""-""$locnum") >> $7 


# print the differences in read ID at the pileups:
echo -e " \n"$loc" \nIDs of reads differing in the pileups:" >> $7
echo $(diff called.lst gt.lst) >> $7


# check which reads were wrongly allocated and where
echo -e " \n"$loc" \nReads falsely allocated to the locus:" >> $7
echo $(diff called.lst gt.lst)

# because it is possible that reads from a different chromosome were falsely aligned to the locus, it is necessary to loop through all chromosomes from which read originated:
diff called.lst gt.lst | grep "<"| grep "X1"| sed -e 's/< //1' > diffs_left_X1.txt
[[ -f diffs_left_X1.txt ]] && sed -e 's/T_X1_\(chr[[:digit:]]\+\)-[[:digit:]]\+|*.*/\1/' diffs_left_X1.txt | sort --version-sort | uniq | while read chr; do
echo $chr
/data/fpantring/samtools view $(grep $(echo "X1_"$chr"_") $4) $(echo $chr) | grep -P "$(grep $(echo $chr"-") diffs_left_X1.txt | paste -s -d"|" -)" | cut -f1-4,10 >> $7
done

# do the same for all X2 reads:
diff called.lst gt.lst | grep "<"| grep "X2"| sed -e 's/< //1' > diffs_left_X2.txt
[[ -f diffs_left_X2.txt ]] && sed -e 's/T_X2_\(chr[[:digit:]]\+\)-[[:digit:]]\+|*.*/\1/' diffs_left_X2.txt | sort --version-sort | uniq | while read chr; do
/data/fpantring/samtools view $(grep $(echo "X2_"$chr"_") $4) $(echo $chr) | grep -P "$(grep $(echo $chr"-") diffs_left_X2.txt | paste -s -d"|" -)" | cut -f1-4,10 >> $7
done


# Repeat the same as above for all reads falsely allocated to alternative loci:
echo -e " \n"$loc" \nReads falsely allocated to alternative loci:" >> $7
diff called.lst gt.lst | grep ">"| grep "X1"| sed -e 's/> //1' > diffs_right_X1.txt
[[ -f diffs_right_X1.txt ]] && sed -e 's/T_X1_\(chr[[:digit:]]\+\)-[[:digit:]]\+|*.*/\1/' diffs_right_X1.txt | sort --version-sort | uniq | while read chr; do
/data/fpantring/samtools view $(grep $(echo "X1_"$chr"_") $4) $(echo $chr) | grep -P "$(grep $(echo $chr"-") diffs_right_X1.txt | paste -s -d"|" -)" | cut -f1-4,10 >> $7
done

diff called.lst gt.lst | grep ">"| grep "X2" | sed -e 's/> //1' > diffs_right_X2.txt
[[ -f diffs_right_X2.txt ]] && sed -e 's/T_X2_\(chr[[:digit:]]\+\)-[[:digit:]]\+|*.*/\1/' diffs_right_X2.txt | sort --version-sort | uniq | while read chr; do
/data/fpantring/samtools view $(grep $(echo "X2_"$chr"_") $4) $(echo $chr) | grep -P "$(grep $(echo $chr"-") diffs_right_X2.txt | paste -s -d"|" -)" | cut -f1-4,10 >> $7
done

done

rm tmp1.txt tmp2.txt diffs_left_X1.txt diffs_left_X2.txt diffs_right_X1.txt diffs_right_X2.txt
cut -f1 $6 | while read loc; do
./2wayLiftover r $loc $2 >> tmp1.txt 
./2wayLiftover r $loc $3 >> tmp2.txt
done
paste $6 tmp1.txt tmp2.txt > TN_FL_liftover.txt


rm gt.lst called.lst

