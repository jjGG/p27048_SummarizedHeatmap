#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#set working directory
#setwd("~/Dropbox/GITHUB/OES beads/p27048_SummarizedHeatmap/R")
# read in csv
# see how many csv files are in folder, read in all and concatenate! # the MascotExport directory should only contain files of interest
#allResultFiles <- system("ls ../2024-01-30_withMascotExports/*.csv")
allFiles <- dir(path = "../2024-01-30_withMascotExports/", full.names = TRUE)
myCSVs <- allFiles[grepl(x = allFiles, pattern = "All.*.csv")]
myCSVs

# get in the first
dat <- read_csv(file = myCSVs[1])
for (i in 2:length(myCSVs)) {
  myDat <- read_csv(file = myCSVs[i])
  dat <- rbind(myDat, dat)
}

# do some filtering -> most BUT NOT ALL are having evalues < 0.05.. -> still do some filtering??
# rank = 1
# bold_red = 1

# maybe further pep_expect or pep_score, maybe decoy hits (REV_) -> no we keep it at this point to calculate protFDR afterwards!
datOK <- dat |> filter(pep_rank == 1 &  pep_isbold == 1)
datOK$rawFile <- sapply(strsplit(sapply(strsplit((datOK$pep_scan_title), split = "\\\\"), function(x)x[6]), split = "\\."), function(x)x[1])

# clip off the first 21 characters for better sample names! #str_sub
datOK$nameTag <- str_sub(string = datOK$rawFile, start = 22, end = nchar(datOK$rawFile))

# every scan only once?
datOK |> dim()
datOK |> distinct() |> dim()
# 2024-08-22 -> all fine at the moment.. (only the 6 header lines are doublicated since we read in 6 same files)

# make all lines unique
datOK <- datOK |> distinct()


