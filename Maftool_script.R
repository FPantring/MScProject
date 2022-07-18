# This script was used to create the Maftool plots and to extract the gene lengths
# of genes with mapping quality variants to calculate the proportion of mutations ´

library(maftools)                 

# read in a file containing all SNVs excluded due to mapping quality from the BRCA data set
somvars <- read.maf("brca.map_qual_SNV.mav")
plotmafSummary(somvars) # overview of the data

# Create the Oncoplot and export it as a PDF:
pdf("Oncoplot.pdf", 6, 4)
par(oma = c(0, 1, 2, 0))
oncoplot(maf = somvars, top = 10)
dev.off()

# Create Lollipop plots for all genes with many non-synonymous mutations
# and export the plots as PDFs:


lollipopPlot(maf = somvars, gene = "HAUS8", showDomainLabel = FALSE) # fails with error:
# Error in lollipopPlot(maf = somvars, gene = "HAUS8", showDomainLabel = FALSE) : 
# Structure for protein HAUS8 not found.

pdf("TTN.pdf", 6, 4)
lollipopPlot(maf = somvars, gene = "TTN", showDomainLabel = FALSE) 
dev.off() 

pdf("HAS2.pdf", 6, 4)
lollipopPlot(maf = somvars, gene = "HAS2", showDomainLabel = FALSE) 
dev.off() 

pdf("KIF16B.pdf", 6, 4)
lollipopPlot(maf = somvars, gene = "KIF16B", showDomainLabel = FALSE) 
dev.off() 

pdf("TBK1.pdf", 6, 4)
lollipopPlot(maf = somvars, gene = "TBK1", showDomainLabel = FALSE) 
dev.off() 



# Extract the gene length of recurringly mutated genes:
library(biomaRt)
library(xlsx)
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Read in a list wiht the number of individuals for which the gene were mutated
# in the first column and the Gene name (hgnc symbol) in the second column:
gene_list <- read.table("genelist.txt")

# subset all genes to genes that were mutated in more than 3 individuals:
gene_list <- gene_list[which(gene_list[,1] > 3),]

# import Supplementary Table 1 from Lopes et al. (2021):
length_list <- readxl::read_xlsx("Table_1_Gene Size Matters An Analysis of Gene Length in the Human Genome.xlsx", skip = 1)

# get the Ensembl gene IDs:
attributes <- c("ensembl_gene_id","hgnc_symbol")
gene_coords <- getBM(attributes=attributes, filter="hgnc_symbol", values=gene_list[,2], mart=ensembl)

# Allocate the gene size to the genes:
gene_coords$size <- sapply(gene_coords$ensembl_gene_id, 
                           function(x) length_list$`Gene length (bp)`[which(length_list$Gene == x)])
# Calculate the rate of mutations proportional to the length of the gene: 
gene_coords$proportion <- apply(gene_coords, 1, 
                                function(x) gene_list[which(gene_list[,2] == x$hgnc_symbol),1]/as.numeric(x$size))

# Extract the genes with the highest proportion of mutations:
high_ratios <- gene_coords[which(gene_coords$proportion > 0.001),]
# format the data frame:
high_ratios$size <- unlist(high_ratios$size)
high_ratios$proportion <- unlist(high_ratios$proportion)

# Display the genes with the highest proportion sorted by proportion
# and export it to a csv file
high_ratios[order(high_ratios$proportion, decreasing = TRUE),]
write.csv(high_ratios[order(high_ratios$proportion, decreasing = TRUE),], 
           "gene_proportions.csv", row.names = FALSE)

