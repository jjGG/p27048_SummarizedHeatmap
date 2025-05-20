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
colnames(datOK)
datOK[91, c("pep_seq", "pep_var_mod_pos")]
# q: I want to parse the sequence and pull the residue that is deamidated (indicated by 1 in pep_var_mod_pos)



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

# Create the function to extract residues at positions of '1'
extract_residues_at_1 <- function(sequence, value) {
  # Handle the case where the value is NA
  if (is.na(value)) {
    return(NA)  # Return NA if value is NA
  }
  
  # Extract the decimal part of the numeric value
  decimal_part <- sub("^[^.]*\\.", "", value)  # Get the part after the decimal point
  
  # Find all positions of '1' in the decimal part
  positions_of_1 <- gregexpr("1", decimal_part)[[1]]  # Returns positions of all '1's
  
  # If '1' is found at any position
  if (length(positions_of_1) > 0) {
    # Extract residues at the corresponding positions
    residues <- sapply(positions_of_1, function(pos) {
      if (pos <= nchar(sequence)) {
        return(substr(sequence, pos, pos))
      } else {
        return(NA)  # Return NA if position exceeds sequence length
      }
    })
    return(paste(residues, collapse = ", "))  # Return residues as a string (separated by commas)
  } else {
    return(NA)  # Return NA if no '1' found
  }
}

# # Example dataframe
# df <- data.frame(
#   Sequence = c("LTNNRYSDIDSIALIDK", "AGCTGGA", "TTTTAA"),
#   Value = c("0.00110000000000000.0", "0.020000000000000.0", "0.0010000000000000.0")
# )
# 
# # Apply the function to each row of the dataframe
# df$Residues <- mapply(extract_residues_at_1, df$Sequence, df$Value)
# 
# # View the updated dataframe
# print(df)

datOK$ModResidues <- mapply(extract_residues_at_1, datOK$pep_seq, datOK$pep_var_mod_pos)




