#!/usr/bin/env python
# coding: utf-8

# This script checks for each mapping quality variant whether potential germline variants or spiked in somatic variants 
# at that locus have the same allele base change as the variant.
# Output: "TN_and_FL.csv": tab separated file with data frame entries to all True Negative (TN) mapping quality variants
#	  "FN.csv": tab separated file with the data frame entries to all False Negative mapping quality variants
#	  "TN_unexplained.csv": tab separated file with the data frame entries to all TN mapping quality variants that 
#	  are not recognised as False Location (FL) or as having the variant allele as reference base at on of the  
#	  alternative locations. These variants need to be manually examined to determine the origin.
#	  "TN_FL_locs.txt": the names of the potential FL mapping quality variants together with the name of the 
#	  alternative mapping location which holds the variant

import pandas as pd
from Bio.Seq import Seq
import matplotlib.pyplot as plt

alt_locs = []
# define a helper function to check for variants noted in the ground truth file whether their 
# allele base change corresponds to the allele base change marked as possible variant by MuTect2. 
# Arguments:
# x: the data frame including the columns SBS with the allele base change according to MuTect2 and 
#    Hap1SBS and Hap2SBS for spiked in allele base changes in the Haplotypes 1 and 2, if present
# reverse: either "forward" if the strand orientation is the same at the alternative position or 
#          "rev_comp" if the strand orientation is the reverse complement of the initial position.
#          Always select "forward" when applying the function to determine whether allele changes
#          are present at the initial position itself
# alt: boolean, whether or not the the function is applied to check an alternative locus
def varcheck(x, reverse, alt):
    for h in ["Hap1", "Hap2"]:
        # check whether a truth variant is present for the position:
        if x["".join([h, "Filter"])] in ["PASS", "MASKED"]:
            if reverse == "forward":
                # Extract the variant allele and check if it matches (one of) the spiked in variant allele(s)
                if x["SBS"].split(">")[1] in x["".join([h, "SBS"])].split(">")[1]:
                    if alt == True: #if checking for the variant position at an alternative mapping location
                        # append the list with the alternative location and the variant location:
                        alt_locs.append(x[["AltLocus", "VafLocus"]])
                    return True
            elif reverse == "rev_comp":
                # when dealing with a reverse complement orientation, check if the variant allele
                # matches the complement of the spiked in variant:
                comp = [str(Seq(i).complement()) for i in x["SBS"].split(">")]
                if comp[1] in x["".join([h, "SBS"])].split(">")[1]:
                    if alt == True:
                        alt_locs.append(x[["AltLocus", "VafLocus"]])
                    return True
            else:
                return "Wrong Input for Argument reverse"
    return False

# also create a slightly modified version of the varcheck function 
# to handle the format of the germline variant files. Arguments:
# x, reverse and alt as in varcheck, 
# y: the germline variant data set
# column: name of the column the values should be added to
def germlinecheck(x, y, reverse, column, alt):
    line = x.loc[x["VafLocus"] == y["VafLocus"]]
    if reverse == "forward":
        if line["SBS"].to_string(index=False).split(">")[1] in y["germSBS"].split(">")[1]:
            if alt == True:
                alt_locs.append(y[["AltLocus", "VafLocus"]])
            x.loc[line.index, column] = True
    elif reverse == "rev_comp":
        comp = [str(Seq(i).complement()) for i in line["SBS"].to_string(index=False).split(">")]
        if comp[1] in y["germSBS"].split(">")[1]:
            if alt == True:
                alt_locs.append(y[["AltLocus", "VafLocus"]])
            x.loc[line.index, column] = True
    else:
        return "Wrong Input for Argument reverse"

# Read in the ground truth Map:
gtMap_df = pd.read_csv('WGS_gtMap.txt', sep="\t")
print(gtMap_df)

# Create a data frame that will contain the results without unnecessary additional information:
gt = gtMap_df[["VafLocus", "SBS", "VafFreq"]]

# Add a boolean column determining whether a somatic variant was spiked in at the locus 
# and the variant is therefore false negative (FN):
gt.insert(3, "FN", gtMap_df.apply(lambda y: varcheck(y, reverse = "forward", alt =  False), axis = 1))
print(gt[["FN"]].apply(sum)) # This is the number of SOMATIC FN mutations

# Read in the germline variant files and transform it to the required format:
gv = pd.read_csv('germline_vars.txt', sep="\t", header = None)
gv["germSBS"] = gv[3] + ">" + gv[4]
gv["VafLocus"] = gv[0] + ":" + gv[1].astype(str)
gv = gv[["VafLocus", "germSBS"]]
gv.apply(lambda l: germlinecheck(gt, l, reverse = "forward", column = "FN", alt =  False), axis = 1)

# Read in the germline variant at alternative locations files and transform it to the required format:
gv_alt = pd.read_csv('germline_alt_loc.txt', sep="\t", header = None)
gv_alt["germSBS"] = gv_alt[7] + ">" + gv_alt[8]
gv_alt["AltLocus"] = gv_alt[0] + ":" + gv_alt[1].astype(str)
gv_alt = gv_alt[["AltLocus", 2, 3, "germSBS"]]
gv_alt.columns = ["AltLocus", "VafLocus", "Orientation", "germSBS"]

# Read in the information on alternative loci
altLoc_df = pd.merge(pd.read_csv('alt_loc_gt.txt', sep="\t").rename(columns = {"VafLocus":"AltLocus"}), 
                     gt[["VafLocus", "SBS", "VafFreq"]], left_on='OriginalLocus', right_on='VafLocus')

# Create a boolean column marking possible false locations (FL). 
# Whether the reads of this alternative locus were actually allocated to the initial locus has to be checked
gt.insert(4, "FL", altLoc_df.apply(lambda l: varcheck(l, reverse = l["Orientation"], alt =  True), axis = 1))
gv_alt.apply(lambda l: germlinecheck(gt, l, reverse = l["Orientation"], column = "FL", alt =  True), axis = 1)

# read in the file excluded.txt produced by the Pileup_Processing.py script:
with open("excluded.txt", "r") as f:
    ex = list(filter(None, f.read().split("\n") ))

# Insert the information from these files as columns into the data frame:
gt.insert(5, "excluded", [j in ex for j in gt["VafLocus"]])                 

# Create a list with all TN and potential FL variants:
TN_FL = []
for i in range(0, len(alt_locs)):
    if alt_locs[i][1] in set(gt.loc[(gt["FN"]==False) | (gt["FL"]==True)]["VafLocus"]):
        print(alt_locs[i][0], alt_locs[i][1], i)
        TN_FL.append([alt_locs[i][0], alt_locs[i][1]])
        

#  Write data frames and lists to files:

# Subset the data frame to only the TNs:
gt.loc[(gt["FN"]==False) | (gt["FL"]==True)].to_csv('TN_and_FL.csv',index=False, sep='\t') 

# Subset to only FNs:
gt.loc[(gt["FN"]==True) & (gt["FL"]==False)].to_csv('FN.csv',index=False, sep='\t')

# Subset to all TNs that are not FL or excluded when checking the reference alleles at alternative locations:
gt.loc[(gt["FN"]==False) & (gt["FL"]==False) & (gt["excluded"]==False)].to_csv('TN_unexplained.csv', index=False, sep='\t')

with open("TN_FL_locs.txt", "w") as f:
    for item in TN_FL:
        for i in item:
            f.write("%s\t" % i)
        f.write("\n")

