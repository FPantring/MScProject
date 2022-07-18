#!/bin/bash

# This script checks for all mapping quality variants and all alternative mapping locations whether germline variants are present.
# It needs to be called with the paths of three files: 
# "Alt_locs.csv": a file produced by Processing.py including all alternative mapping locations annotated in the reads of the locus pileups
# a .vcf file containing the germline variants
# a ground truth map file where the names of all mapping quality variants are listed
# sbatch germline_map.sh Alt_locs.csv HG00110.vcf WGS_gtMap.txt

# If a file with this name exists it has to be removed, since this file will later be appended with information
rm germline_alt_loc.txt germline_vars.txt

# check all alternative mapping locations for germline variants:
cut -f1,2 $1 | tail -n +2 | awk '{print $1"[[:space:]]"$2}' | while read loc; do
if grep "$loc[[:space:]]" $2 > /dev/null; then  
echo -e "$(grep -P "$loc[[:space:]]" $1)\t$(grep -P "$loc[[:space:]]" $2)" >> germline_alt_loc.txt 
fi
done

# check for germline variants at the initial locus:
cut -f1 $3 | sed -e 's/:/[[:space:]]/1' | while read loc; do
grep "$loc[[:space:]]" $2 >> germline_vars.txt
done
