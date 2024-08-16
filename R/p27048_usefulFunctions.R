#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# libraries
library(dplyr)
#library(prolfqua)
library(stringr)
library(readr)
library(ggplot2)
library(gplots)
#library(heatmap3)
#library(pheatmap)
library(heatmap.plus)


#functions
dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1
  return(edist)
}

get_protFDR <-  function(tableMat, revString = "REV_") {
  protFDR <- round(sum((str_count(rownames(tableMat), pattern = revString) > 0))/nrow(tableMat)*100,3)
  return(protFDR)
}


