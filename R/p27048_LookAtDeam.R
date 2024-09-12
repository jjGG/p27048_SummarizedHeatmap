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

table(is.na(datOK$pep_var_mod_pos))
datOK$DeamidationBool <- str_count(datOK$pep_var_mod_pos, pattern = "1") > 0
datOK$DeamidationBool[is.na(datOK$DeamidationBool)] <- FALSE
table(datOK$DeamidationBool)
#
datOK$pep_var_mod[datOK$DeamidationBool]

table(datOK$pep_var_mod[datOK$DeamidationBool])

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


# hunt the methyl
datOK$MethylBool <- str_count(datOK$pep_var_mod, pattern = "ethyl") > 0
datOK$MethylBool[is.na(datOK$MethylBool)] <- FALSE
table(datOK$MethylBool)
#
datOK$pep_var_mod[datOK$MethylBool]
# table(datOK$rawFile[datOK$MethylBool])
# > table(datOK$rawFile[datOK$MethylBool])
# 
# 20221202_010_S433299_75773 
# 213

str(datOK)
