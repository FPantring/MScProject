#!/usr/bin/env python
# coding: utf-8

# This script checks for all suspected FL variants whether reads from the alternative mapping location
# which holds a variant have truly been incorrectly aligned at the locus

import re
import pandas as pd
from Result_Analysis import gt
from Result_Analysis import varcheck

# Create an empty dictionary which will be filled up with all reads falsely aligned to the locus of all suspected FL and TN: 
pileup_dict = {}

# Open the file created by check_all_pileups.sh which (among others) compiled all reads incorrectly aligned to the locus:
with open("pileup_infos.txt") as f:
    # Extract all falsely allocated reads and information on their origin into a dataframe with the name of the locus as dictonary key: 
    for result, chr_name in re.findall("Reads falsely allocated to the locus:\n(.*?)\n(chr\d+:\d+) \nReads falsely allocated to alternative loci:", f.read(), re.S):
        if result != " ":
            temp_df = pd.DataFrame([x.split('\t') for x in result.split('\n')], columns = ["Read_name", "Nr", "Pos_chr", "Position", "Sequence"])
            temp_df.drop(temp_df.tail(1).index,inplace=True)
            pileup_dict[chr_name] = temp_df
        else: # if no reads were incorrectly aligned to the locus, asign an empty dataframe to the dictonary
            pileup_dict[chr_name] = pd.DataFrame()

# Read in the file TN_FL_liftover.txt which includes the following information in the columns:
# VafLoc: the position of the FL or TN
# AltLoc: the reference position of the alternative mapping location of interest 
#         (reminder: only relevant alternative locations, i.e. locations where a variant or sequencing error
#         were spiked in are included in this list)
# Alt_X1: the lifted over position of that alternative mapping location in the simulated haplotype 1
# Alt_X2: the lifted over position of that alternative mapping location in the simulated haplotype 2
# the columns Alt_chr1 and Alt_chr2 hold the chromosomes corresponding to the positions Alt_X1 and Alt_X2
liftovers = pd.read_csv("TN_FL_liftover.txt", sep = "\t", header = None)
liftovers = liftovers[[0, 1, 3, 4, 5, 6]]
liftovers.columns = ["AltLoc", "VafLoc", "Alt_chr1", "Alt_X1", "Alt_chr2", "Alt_X2"]

# Check if for suspected FL or TN some of the reads of the alternative position with the variant or seq error
# were actually aligned to the MuTect2 variant position:
TN_FL_dict = {}
for i in range(0, len(liftovers["VafLoc"])):
    for index, row in pileup_dict[liftovers["VafLoc"][i]].iterrows():
        if "_X1_" in row["Read_name"] and row["Pos_chr"] == liftovers["Alt_chr1"][i] and int(row["Position"]) in range(int(liftovers["Alt_X1"][i])- 76, int(liftovers["Alt_X1"][i])+ 76):
            TN_FL_dict[liftovers["AltLoc"][i]] = liftovers["VafLoc"][i]
        elif "_X2_" in row["Read_name"] and row["Pos_chr"] == liftovers["Alt_chr2"][i] and int(row["Position"]) in range(int(liftovers["Alt_X2"][i])- 76, int(liftovers["Alt_X2"][i])+ 76):
            TN_FL_dict[liftovers["AltLoc"][i]] = liftovers["VafLoc"][i]

# Create a new data frame (tnfl_df) which is a subset of the data frame gt from the Result_Analysis.py script
# including only the TN and all suspected FL variants
tnfl_df = gt.loc[(gt["FN"]==False) | (gt["FL"]==True)]
print(tnfl_df) 

# Set all values in the columns FL to False and reassign the value True for all actual FL variants:
tnfl_df["FL"] = False
gt["FL"] = False
for v in TN_FL_dict.keys():
    tnfl_df.loc[tnfl_df["VafLocus"] == TN_FL_dict[v], "FL"] = True
    gt.loc[gt["VafLocus"] == TN_FL_dict[v], "FL"] = True



# Draw bar plot of the filter categories:
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt

# read in a ground truth file of all SBSs filtered by FilterMutectCalls
# this data frame includes a column called "Filter" which holds the category used to filter the SBS  
filters = pd.read_csv("gtMap_Filters.txt", sep = "\t")

# Insert a column determining whether the variant is a FN using the function varcheck from Results_Analysis.py: 
filters.insert(14, "FN", filters.apply(lambda y: varcheck(y, reverse = "forward", alt =  False), axis = 1))

# Change some of the category labels as needed:
for i, row in filters.iterrows():
    # whenever multiple filter categories are listed (separated with a semicolon), 
    # replace with the category with the label "multiple filters":
    if ";" in str(row["Filter"]):
        filters.loc[i,"Filter"] = "multiple filters"  
    # if a filter occurs only once, replace the category with the label "other filters":
    elif sum(filters["Filter"] == row["Filter"]) == 1:
        filters.loc[i,"Filter"] = "other filters"
df = filters.groupby(["FN", "Filter"]).size()


# Set the variables needed to create the bar plot:
index = np.arange(len(set(filters["Filter"])))
bar_width = 0.42
opacity = 0.8

# insert first set of bars, displaying the number of FNs in each category
rects1 = plt.bar(index, (15, 2, 215, 135, 694, 0, 815, 275), bar_width,
alpha=opacity,
color='red',
label='False Negative')

# insert second set of bars, displaying the number of TNs in each category
rects2 = plt.bar(index + bar_width, (0, 1, 28, 1908, 116, 3, 22, 150), bar_width,
alpha=opacity,
color='green',
label='True Negative')

# Add the observed numbers as text above the bars 
for i, n in enumerate((15, 2, 215, 135, 694, 0, 815, 275)):
    plt.text(i, n+10, n, horizontalalignment='center', color='red', fontsize=9)
    
for i, n in enumerate((0, 1, 28, 1908, 116, 3, 22, 150)):
    plt.text(i + bar_width, n+10, n, horizontalalignment='center', color='green', fontsize=9)

# Set axis labels:
plt.xlabel('Filter category', fontsize=12)
plt.ylabel('Number of SBSs filtered', fontsize=12)
plt.xticks(index + (bar_width/2), sorted(set(filters["Filter"])), rotation=60, horizontalalignment='right', fontsize=12)
plt.legend()

# Create and save the graphic:
plt.gcf().set_size_inches(6, 6)
plt.tight_layout()
plt.savefig('Filter_cats.pdf')



# Final statistics on the occurrence of FN and FL variants
print(gt[["FN", "FL"]].apply(sum))
len(gt["VafLocus"])

