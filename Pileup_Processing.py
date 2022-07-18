#!/usr/bin/env python
# coding: utf-8

# This script extracts all alternative mapping locations from the annotations of the reads in the pileup of all mapping quality variant loci.
# It needs to be called with the paths and names to three input files: the .bam file, the .fa reference file and 
# the ground truth map file which lists all mapping quality variants with including information about variants spiked in at that location.
# When running this script with real data instead of simulated, instead of the ground truth map a tap separated file should be uploaded,
# which features the mapping quality variant positions in the first column in the format chr1:12345678 and the single base substitution at that position
# in the second column in the format R>V, with R being the reference allele at this position and V the variant allele. Note that this file should include a header, 
# otherwise the first row of the file will be missed in the analysis.
# Example call of this script: python3 Pileup_Processing.py T_WGS.bam GRCh38.d1.vd1.fa WGS_gtMap.txt > output.txt

# This script creates a lot of output which should be saved in a log file (using the command shown above) and can be manually reviewed for troubleshooting purposes
# However, truly relevant for the Analysis Pipeline are only the files Alt_locs.csv and excluded.txt which will be produced at the end of the script

# Import packages:
import pysam
import sys
import regex as re
import pandas as pd
from Bio.Seq import Seq
from fuzzywuzzy import process
from fuzzywuzzy import fuzz
# Import pairwise2 module
from Bio import pairwise2
# Import format_alignment method
from Bio.pairwise2 import format_alignment
import csv


# Open connection to the alignment file
print(sys.argv, sys.argv[1])
samfile = pysam.AlignmentFile(sys.argv[1], "rb")

from pyfaidx import Faidx
# Open connection to reference file
fa = Faidx(sys.argv[2])


# create empty dictionaries that will be filled in later loops
poslist = []
refdict = {}
alleledict = {}
appended_data = {}


# define a helper function to determine the position of variant bases. Arguments:
# l: the variant position
# startpos: the start position of the read
# cigar: the CIGAR string
# orig: boolean, select True if the function is used to adjust according to the CIGAR string of the original alignment, 
# 	select False if the function is applied to adjust according to the CIGAR string of an alternative alignment location.
def base_pos(l, startpos, cigar, orig):
    # access the position of the bases at the variant position of interest and the corresponding reference position
    # to get this position it is necessary to subtract the start position of the read alignment to the reference
    # from the position of the variant (the second entry of the list loc -1 (subtract one to correct to 0-based))
    # Get base position at query_alignment_sequence:
    distance = l-startpos
    # also correct the base position if indels are present:
    if re.findall('I|D|S', cigar) != []: # check if insertion or deletion present
        start = 0
        dl = []
        # separate the string at the number before a D or I to extract all deletions and insertions:
        for num1, i_or_d in re.findall('(\d+)([IDSM])', cigar):
            if i_or_d not in 'IDS':
                # add the number of matches and substitutions before the indel 
                # to calculate the indels positions within the cigar string:
                start += int(num1)
                continue
            dl.append([num1, i_or_d, start])
            if num1:
                start += int(num1)
        for j in dl:
            # if the indel occurred after the loci in question, no adjustment of the position is necessary)
            if distance < j[2]:
                continue
            # if a region in the beginning is sofclipped, this is likely a result of the variant in the read,
            # because due to the mismatch in the beginning the alignment software assumes softclipping.
            # Since the "coordinate shown in SAM is the position of the first aligned base" (Li, Heng et al. 
            # “The Sequence Alignment/Map format and SAMtools.” Bioinformatics (Oxford, England) 
            # vol. 25,16 (2009): 2078-9. doi:10.1093/bioinformatics/btp352), this coordinate has to be adjusted 
            # by the number of soft clipped bases: 
            elif j[2] == 0 and j[1] == "S":
                if orig == False:
                    distance = distance - int(j[0])
                else:
                    distance = distance + int(j[0])
            # check whether the position of interest is within a deleted or softclipped region:
            elif distance in range(int(j[2]), int(j[2])+ int(j[0])) and (j[1] == "D" or j[1] == "S"):
                # if there is a deletion present at the position of interest, skip this read later:
                # note that positions within softclipped regions are treated as deletions 
                # because they usually diverge clearly from the original sequence. Thus it is unlikely
                # that reads from this alternative softclipped regions would be allocated to the initial position,
                # while the probability of identifying the variant allele as the reference allele at this position
                # can be suspected to be close to 1/4, thus leading to falsely excluded variants 
                print("\n", "DELETION!!!!!!!", cigar, distance, "\n")
                return "DELETION"
                    
            # if the indel occurred before the loci in question (or at the loci in the case of an insertion), 
            # adjust the position number:
            else:
                # in case of an insertion, subtract the number of inserted bases from the current position 
                # when adjusting for an alternative location. 
                if j[1] == "I":
                    if orig == False:
                        distance = distance - int(j[0])
                    else: # When adjusting for the original CIGAR string, add the number of inserted bases
                        distance = distance + int(j[0])
                # in case of a deletion, add the number of deleted bases to the current position when adjusting for an alternative location
                elif j[1] == "D":
                    if orig == False:
                        distance = distance + int(j[0])
                    else: # When adjusting for the original CIGAR string, subtract the number of deleted bases
                        distance = distance - int(j[0])
    if cigar != "76M":
       # if cigar not in appended_data.keys():
        appended_data[cigar] = [l-startpos, distance]
    return distance


with open(sys.argv[3]) as gtMap:
    next(gtMap) # skip the header line of the txt file
    for line in csv.reader(gtMap, dialect="excel-tab"):
        # assign the suspected reference and variant allele to dictionaries:
        VafLoc, alleles = line[0], line[1].split(">")
        refdict[VafLoc], alleledict[VafLoc] = alleles[0], alleles[1]
        # make the string location accessible:
        loc = VafLoc.split(":")

        # access the read pileup at a specific position:
        # REMEMBER: take (position-1, position) instead of (position, position+1), because coordinates are 0-based
        # access the alternative hits with .get_tag("XA") and concatenate all alternative hits to one list
        # each list element then includes: chromosome and position of the read alignment start position at the alternative location, 
        # CIGAR string of the alignment at the alternative location, edit distance (number of nucleotide changes), 
        # the position of the variant within the read, the reverse position within the read (position when counted from the end of the read),
        # the initial location of the read in the format chr:positionnumber, the read sequence and the CIGAR string of the original alignment
        element=[re.sub("[+-]", "", s) + "," + str(base_pos(l = (int(loc[1]) - 1), startpos = pileupread.alignment.reference_start,
                                                            cigar = pileupread.alignment.cigarstring, orig = True)) + "," +
                 str(len(pileupread.alignment.query_sequence)-base_pos(l = (int(loc[1]) - 1), startpos = pileupread.alignment.reference_start,
                                                                       cigar = pileupread.alignment.cigarstring, orig = True) - 1) +
                 "," + loc[0] + ":" + loc[1] + "," + str(pileupread.alignment.query_sequence) + "," + str(pileupread.alignment.cigarstring)
                 for pileupcolumn in samfile.pileup(loc[0], (int(loc[1]) - 1), int(loc[1])) for pileupread in pileupcolumn.pileups
                 # check whether the position in question has no alternative alignments in the genome
                 # and whether the position is deleted. If so, this read is not relevant:
                 if pileupread.alignment.has_tag("XA") and base_pos(l = (int(loc[1]) - 1), startpos = pileupread.alignment.reference_start,
                                                                    cigar = pileupread.alignment.cigarstring, orig = True) != "DELETION" and \
                 base_pos(l = (int(loc[1]) - 1), startpos = pileupread.alignment.reference_start,
                          cigar = pileupread.alignment.cigarstring, orig = True) != "DELETION"
                 for s in list(filter(None, re.split(";", pileupread.alignment.get_tag("XA"))))]
        poslist = poslist + element

                     
# Convert the list into a data frame where the information is separated into columns
# dis is short for edit distance, org_cig stands for the CIGAR string of the original alignment
posframe = pd.DataFrame(columns = ["chromosome", "position", "CIGAR", "dis", "pos", "reverse", "initial pos", "seq", "org_cig"], 
                        data=[col.split(",") for col in poslist[0:]])
posframe = posframe.drop_duplicates(keep = "first")

# transform the columns pos and reverse to adjust the variant position according to the CIGAR string at the alternative location:
posframe["pos"] = posframe.apply(lambda x: base_pos(l = int(x["pos"]), startpos = 0, cigar = x["CIGAR"], orig = False), axis = 1)
posframe["reverse"] = posframe.apply(lambda x: base_pos(l = int(x["reverse"]), startpos = 0, cigar = x["CIGAR"], orig = False), axis = 1)

# drop all rows whose variant position is within a deleted or softclipped region:
posframe = posframe.drop(posframe[(posframe["pos"] == "DELETION") | (posframe["reverse"] == "DELETION")].index)


# check for all alternative positions whether the reference corresponds to a reads forward, reverse, complement or reverse complement strand:
def revcheck(i):
    # it is important to adjust the position according to the CIGAR string, because if some bases are annotated as soft clipped
    # the alignment of the read will be shifted by the amount of softclipped bases in the beginning
    if (int(i["position"]) + base_pos(0, 0, i["CIGAR"], orig = False)) > 0 and (int(i["position"]) + base_pos(0, 0, i["CIGAR"], orig = False)) > 0:
        choices = {str(fa.fetch(i["chromosome"], (int(i["position"]) + base_pos(0, 0, i["CIGAR"], orig = False)), 
                                int(i["position"]) + base_pos(0, 0, i["CIGAR"], orig = False) + len(i["seq"]) - 1)): "forward",
                   Seq(str(fa.fetch(i["chromosome"], (int(i["position"]) + base_pos(0, 0, i["CIGAR"], orig = False)),
                                    int(i["position"]) + base_pos(0, 0, i["CIGAR"], orig = False) + len(i["seq"]) - 1))).reverse_complement(): "rev_comp"}
        return choices[process.extractOne(i["seq"], choices.keys())[0]]
    else:
        return "DELETION"

posframe['orientation'] = posframe.apply(revcheck, axis = 1)
posframe = posframe[posframe['orientation'] != "DELETION"]


# Using the orientation of the read, it is finally possible to calculate the exact variant location at the alternative position
# by adding the position within the read to the read alignment start position at the the alternative location. 
# Define the helper function for that purpose:   
def poscalc(i):
    if i["orientation"] == "forward":
        return pd.to_numeric(i["position"]) + pd.to_numeric(i["pos"])
    elif i["orientation"] == "rev_comp":
        return pd.to_numeric(i["position"]) + pd.to_numeric(i["reverse"])
    
print(posframe.apply(revcheck, axis = 1))

# create an additional column in the data frame, which includes the actual variant position:
posframe['actual pos'] = posframe.apply(poscalc, axis = 1)

# drop duplicated rows to minimise computational time:
posframe = posframe.drop_duplicates(subset = ["chromosome", "actual pos", "orientation", "initial pos"], keep = "first")



# go through all alternative loci and check the base at the variant position. 
discards = []
def varbase(i): 
    if i["orientation"] == "rev_comp":
        allelepos = Seq(str(fa.fetch(i["chromosome"], (int(i["actual pos"])), int(i["actual pos"])))).complement()
        if int(i["actual pos"]) > 11:
            print(str(fa.fetch(i["initial pos"].split(":")[0], (int(i["initial pos"].split(":")[1])-10), int(i["initial pos"].split(":")[1]) + 10)), i["initial pos"], "initial")
            print(Seq(str(fa.fetch(i["chromosome"], (int(i["actual pos"])-10), int(i["actual pos"]) + 10))).reverse_complement(), i["CIGAR"], "comp", i["initial pos"])
        if str(allelepos) in alleledict[i["initial pos"]]:
            print("\n", "REFERENCE BASE!!!!", "\n", Seq(str(fa.fetch(i["chromosome"], (int(i["actual pos"]) - 10), int(i["actual pos"]) + 10))).reverse_complement(), 
                  i["initial pos"], i["actual pos"], i["orientation"])
    else:
        allelepos = fa.fetch(i["chromosome"], (int(i["actual pos"])), int(i["actual pos"]))
        if int(i["actual pos"]) > 11:
            print(fa.fetch(i["initial pos"].split(":")[0], (int(i["initial pos"].split(":")[1])-10),
                           int(i["initial pos"].split(":")[1]) + 10), i["initial pos"], "initial")
            print(fa.fetch(i["chromosome"], (int(i["actual pos"]) -10), int(i["actual pos"]) + 10), i["CIGAR"], i["initial pos"])
        if str(allelepos) in alleledict[i["initial pos"]]:
            print("\n", "REFERENCE BASE!!!!", "\n", fa.fetch(i["chromosome"], (int(i["actual pos"]) - 10), int(i["actual pos"]) + 10), i["initial pos"], i["actual pos"], i["orientation"])
    if str(allelepos) in alleledict[i["initial pos"]]:
        print("\n", "REFERENCE BASE!!!!", "\n")
        discards.append(i["initial pos"])


            
# apply the function to the data frame of unique rows       
posframe.apply(varbase, axis = 1)  

# create a list including all loci for which no alternate loci with the variant allele exist
ls = list(refdict.keys())
[ls.remove(d) for d in set(discards)] 


# Extract the relevant reference positions by dropping duplicate positions
uniqs = posframe[["chromosome", "actual pos", "initial pos", "orientation"]].drop_duplicates()

# This list includes all variants whose alternate alleles were not present in the alternative loci:
print(ls)
print(len(ls), "variants had their variant alleles not present in the alternative loci. \nThe following",
     len(set(discards)), "variant(s) had their variant allele present at an alternate position:", sorted(set(discards), key=lambda x: int(x.partition('chr')[2].partition(':')[0])))
print(posframe)

print(appended_data)

uniqs.to_csv('Alt_locs.csv',index=False, sep='\t')
with open('excluded.txt', 'w') as f:
    f.write('\n'.join(sorted(set(discards), key=lambda x: int(x.partition('chr')[2].partition(':')[0]))))

samfile.close()
fa.close()

