#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# set wd to source file
source("p27048_usefulFunctions.R") # read in function

# globals
fgczProject <- "p27048"
descri <- "LookAtDeamidation"

source("p27048_somePreprocessing.R")



# what do we have in datOK
colnames(datOK)
datOK$pep_var_mod

table(datOK$pep_var_mod)

# var mod pos string
# LTNNRYSDIDSIALIDK Y 2 Deamidated (NQ) 0.00110000000000000.0
# Deamidation == 1

datOK$DeamidationBool <- str_count(datOK$pep_var_mod_pos, pattern = "1") > 0
datOK$DeamidationBool[datOK$DeamidationBool == NA] <- FALSE
# text table
table(datOK$DeamidationBool, datOK$rawFile)

# plot
# q: how can I increase margin to read the labels
pdf("Deamidation_yes_no.pdf", 19,19)
par(mar = c(4.1, 14.4, 4.1, 1.9))
barplot(table(datOK$DeamidationBool, datOK$rawFile), horiz = TRUE, las=2, margin = 1, col = c("red", "blue"), main = "Deamidation", xlab = "Count", ylab = "RawFile")
dev.off()


pdf("Deamidation_yes_no_nonstacked.pdf", 19,28)
par(mar = c(4.1, 14.4, 4.1, 1.9))
barplot(table(datOK$DeamidationBool, datOK$rawFile), horiz = TRUE, las=2, margin = 1, beside = TRUE, col = c("red", "blue"), main = "Deamidation", xlab = "Count", ylab = "RawFile")
dev.off()
