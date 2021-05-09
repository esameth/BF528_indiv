#####################################################################
# Analyst section                                                   #
# Plot top GO terms                                                 #
#####################################################################

up_reg <- data.frame(Term='Extracellular Organelle', ES=5.98)
up_reg <- rbind(up_reg, data.frame(Term='Respiratory Chain', ES=6.24))
up_reg <- rbind(up_reg, data.frame(Term='Cellular Respiration', ES=6.7))
up_reg <- rbind(up_reg, data.frame(Term='Metabolic Process', ES=7.1))
up_reg <- rbind(up_reg, data.frame(Term='Mitochondria', ES=13.48))

down_reg <- data.frame(Term='Extracellular Matrix', ES=4.6)
down_reg <- rbind(down_reg, data.frame(Term='Neurogenesis', ES=4.8))
down_reg <- rbind(down_reg, data.frame(Term='Embryonic Morphogenesis', ES=5.23))
down_reg <- rbind(down_reg, data.frame(Term='Cell Proliferation', ES=6.64))
down_reg <- rbind(down_reg, data.frame(Term='Cardiovascular System', ES=7.87))

par(mar=c(5,16,4,2))
barplot(up_reg$ES, main = "Up-Regulated GO Terms", xlab = "Enrichment Score", 
        names.arg = up_reg$Term, 
        cex.names = 1.5, cex.main = 1.5, cex.lab = 1.5, horiz = TRUE, las=1)

barplot(down_reg$ES, 
        main = "Down-Regulated GO Terms", 
        xlab = "Enrichment Score", 
        names.arg = down_reg$Term,
        cex.names = 1.5,
        cex.main = 1.5,
        cex.lab = 1.5,
        horiz = TRUE, las=1)
