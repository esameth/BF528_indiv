#####################################################################
# Biologist section                                                 #
# Create line graphs of the FPKM values for each sample             #
# and perform hierarchical clustering of samples                    #
#####################################################################

library(gplots)
library(gridExtra)
library(reshape2)
library(tidyverse)
library(data.table)

# Genes of interest
sarcomere <- c('Pdlim5', 'Pygm', 'Myoz2', 'Des', 'Csrp3', 'Tcap', 'Cryab')
mitochondria <- c('Mpc1', 'Prdx3', 'Acat1', 'Echs1', 'Slc25a11', 'Phyh')
cell_cycle <- c('Cdc7','E2f8','Cdk7','Cdc26','Cdc6','E2f1','Cdc27','Bora','Cdc45','Rad51','Aurkb','Cdc23')

# Read in FPKM tracking tables for each of the replicates of the samples 
P0_1 <- data.frame(fread('samples/P0_1.fpkm_tracking', select=c('gene_short_name', 'FPKM'), header=TRUE))
P0_2 <- data.frame(fread('samples/P0_2.fpkm_tracking', select=c('gene_short_name', 'FPKM'), header=TRUE))
P4_1 <- data.frame(fread('samples/P4_1.fpkm_tracking', select=c('gene_short_name', 'FPKM'), header=TRUE))
P4_2 <- data.frame(fread('samples/P4_2.fpkm_tracking', select=c('gene_short_name', 'FPKM'), header=TRUE))
P7_1 <- data.frame(fread('samples/P7_1.fpkm_tracking', select=c('gene_short_name', 'FPKM'), header=TRUE))
P7_2 <- data.frame(fread('samples/P7_2.fpkm_tracking', select=c('gene_short_name', 'FPKM'), header=TRUE))
Ad_1 <- data.frame(fread('samples/Ad_1.fpkm_tracking', select=c('gene_short_name', 'FPKM'), header=TRUE))
Ad_2 <- data.frame(fread('samples/Ad_2.fpkm_tracking', select=c('gene_short_name', 'FPKM'), header=TRUE))

########### Line graph
# Get the average FPKM for each sample
average <- function(geneList, rep_1, rep_2) {
  df <- merge(subset(rep_1, gene_short_name %in% geneList), 
              subset(rep_2, gene_short_name %in% geneList), 
              by='gene_short_name')
  rownames(df) <- df$gene_short_name
  # Drop first column
  df <- subset(df, select=-c(1))
  df$avg <- rowMeans(df)
  return(subset(df, select=-c(1,2)))
}

# Plot the FPKM values
plot_fpkm <- function(df, title) {
  colnames(df) <- c('P0', 'P4', 'P7', 'Ad')
  
  # Add missing genes
  if (title == 'C) Cell Cycle') {
    df <-rbind(df, 'Bora'=c(0,0,0,0))
  }
  if (title == 'B) Mitochondria') {
    df <-rbind(df, 'Mpc1'=c(0,0,0,0))
  }
  
  # Transpose and melt
  fpkm <- melt(t(df), id.vars=colnames(df)) 
  ggplot(fpkm, aes(x=Var1, y=value, group=Var2)) + 
    geom_line(aes(color=Var2)) + 
    geom_point() + 
    labs(title = title, x='In vivo Maturation', y='FPKM', colour='Gene')
}


sarc <- data.frame(P0=average(sarcomere, P0_1, P0_2))
sarc$P4 <- average(sarcomere, P4_1, P4_2)
sarc$P7 <- average(sarcomere, P7_1, P7_2)
sarc$Ad <- average(sarcomere, Ad_1, Ad_2)
plot_fpkm(sarc, 'A) Sarcomere')

mito <- data.frame(P0=average(mitochondria, P0_1, P0_2))
mito$P4 <- average(mitochondria, P4_1, P4_2)
mito$P7 <- average(mitochondria, P7_1, P7_2)
mito$Ad <- average(mitochondria, Ad_1, Ad_2)
plot_fpkm(mito, 'B) Mitochondria')

cc <- data.frame(P0=average(cell_cycle, P0_1, P0_2))
cc$P4 <- average(cell_cycle, P4_1, P4_2)
cc$P7 <- average(cell_cycle, P7_1, P7_2)
cc$Ad <- average(cell_cycle, Ad_1, Ad_2)
plot_fpkm(cc, 'C) Cell Cycle')



########### Hierarchical Clustering
# Read in Cuffdiff output
cuffdiff <- read.table('cuffdiff_out/gene_exp.diff', header=TRUE)
cuffdiff <- subset(cuffdiff, status == "OK" & significant == 'yes')

# Determine duplicated genes in P0_1
duplicates <- P0_1[duplicated(P0_1$gene_short_name), 'gene_short_name']

# Get top 1000 DE genes in P0_1 and Ad that does not include the duplicates from P0_1
cuffdiff_no_dup <- subset(cuffdiff, !gene %in% duplicates)
genes <- cuffdiff[order(cuffdiff_no_dup$q_value), ]$gene[1:1000]

FPKM_matrix <- merge(subset(P0_1, gene_short_name %in% genes), 
                     subset(P0_2, gene_short_name %in% genes), 
                     by='gene_short_name')

FPKM_matrix <- merge(FPKM_matrix, P4_1, by='gene_short_name', all=TRUE)
FPKM_matrix <- merge(FPKM_matrix, P4_2, by='gene_short_name', all=TRUE)
FPKM_matrix <- merge(FPKM_matrix, P7_1, by='gene_short_name', all=TRUE)
FPKM_matrix <- merge(FPKM_matrix, P7_2, by='gene_short_name', all=TRUE)
FPKM_matrix <- merge(FPKM_matrix, Ad_1, by='gene_short_name', all=TRUE)
FPKM_matrix <- merge(FPKM_matrix, Ad_1, by='gene_short_name', all=TRUE)



