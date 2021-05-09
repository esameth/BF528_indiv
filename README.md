# Project Description
This project focuses on reproducing the analytical and biological aspects from O’Meara et al. The objectives of the study were to determine if myocytes revert the transcriptional phenotype to a less differentiated state during regeneration and to systematically interrogate the transcriptional data to identify and validate potential regulators of this process.

In reproducing the study differential expression values between P0_1 and Ad was interpreted using DAVID. By doing so, we were able to better understand how neonatal mice are able to regenerate their heart tissue but lose this ability later in development.

O’Meara et al. Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration. Circ Res. Feb 2015. PMID: 25477501

# Repository Contents
## histogram.R
* Input: **gene_exp.diff** from `Cuffdiff`
* Output: `up_reg.txt`, `down_reg.txt`
Identifies differentially expressed genes in the P0 and Ad samples. Genes are filtered using test status = OK and significant = yes. The distribution of the log2FC of these genes were then plotted and up-regulated (log2FC > 0) and down-regulated (log2FC < 0) genes were separated into their perspective files.

## GO_terms.R
Creates a bar plot of the top up- and down-regulated cluster terms from DAVID and their enrichment scores

## maturation_and_clustering.R
* Input: **[replicate].fpkm_tracking** for each replicate outputted by `Cufflinks`
Creates line graphs showing the change in FPKM during maturation (P0 to Ad) for genes associated with Sarcomere, Mitochondria, and Cell Cycle. Hierarchical clustering of the replicates was then performed.

