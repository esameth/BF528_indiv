#####################################################################
# Analyst section                                                   #
# Identifies differentially expressed genes in P0 and adult samples #
#####################################################################

# Read in Cuffdiff output
cuffdiff <- read.table('cuffdiff_out/gene_exp.diff', header=TRUE)

###### Data filtering
# Get summary of testing status
summary(as.factor(cuffdiff$status))
# Filter genes with test status successful
cuffdiff_OK <- subset(cuffdiff, status == "OK")

###### Top 10 differentially expressed genes according to q-value
cuffdiff_OK <- cuffdiff[order(cuffdiff_OK$q_value), ]

###### Histograms of log2FC
# Subset of each with significant = yes
cuffdiff_OK_sig <- subset(cuffdiff_OK, significant == "yes")

par(mfrow=c(1, 2), oma = c(0, 0, 2, 0))
hist(cuffdiff_OK$log2.fold_change, breaks=30, 
     main="A) Log2FC for All Genes", xlab='log2FC')
hist(cuffdiff_OK_sig$log2.fold_change, breaks=30, 
     main="B) Log2FC for Significant Genes", xlab='log2FC')
mtext("Genes with Test Status OK", outer = TRUE, cex = 1.5)

###### Separate up- and down-regulated gene
# Summary of number of up- and down-regulated (all)
table(sign(cuffdiff$log2.fold_change.))
# Summary of number of up- and down-regulated (OK)
table(sign(cuffdiff_OK$log2.fold_change.))
# Summary of number of up- and down-regulated (OK, significant)
table(sign(cuffdiff_OK_sig$log2.fold_change.))
# Summary of number of up- and down-regulated (OK, significant p < 0.01)
table(sign(subset(cuffdiff_OK, p_value < 0.01)$log2.fold_change.))

###### Write to files
up_reg <- subset(cuffdiff_OK_sig, log2.fold_change. > 0)
down_reg <- subset(cuffdiff_OK_sig, log2.fold_change. < 0)

write.table(up_reg$gene, 'up_reg.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(down_reg$gene, 'down_reg.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)

