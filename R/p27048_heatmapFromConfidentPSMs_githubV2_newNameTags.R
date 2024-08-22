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
descri <- "newNameTagFile"

# read in data and do some preprocessing
source("p27048_somePreprocessing.R")


# how does it look 
# here nameTag from rawFileName
table(datOK$nameTag)

#filtering for and evalue
# filtering for unique lines is done in the preprocessing (only header lines are ok to be doublicates)
pExpectThreshold <- 0.05
datOK |> dim()
datOK <- datOK |> filter(pep_expect < pExpectThreshold) # at this step we looks actually quite a bit of psms
datOK |> dim()
# before eV filtering: 181954
# after filtering (0.05): 133880



# joining back better sample names and col annotation from Shevan!
# sampleAnno <- read_tsv(file = "p27048_sampleNamesNrawFile-V6.txt")
sampleAnno <- read_tsv(file = "newSampleNamesNrawFile_annotation_Aug2024-V9-newNameTags.txt")
colnames(sampleAnno) <- c("rawFile", "nameTag", "Time", "Location")

# clean nameTags
sampleAnno$nameTag <- gsub(x = sampleAnno$nameTag, pattern = " ", replacement = "_")

# protein annotation for better classification w/ Shevan
protAnno <- read.csv(file = "Human_proteins_classifications-V2_2024-05-27_noMoreDoubs.csv")
colnames(protAnno) <- c("prot_acc","TrivialName", "Species", "Primary", "Secondary")
table(protAnno$Primary)
table(protAnno$Secondary)

# decoy tag lost for trivial, bring it back
decoyString <- rep("",nrow(protAnno))
decoyBool <- str_count(string = protAnno$prot_acc,pattern = "REV_") >  0
table(decoyBool)
decoyString[decoyBool] <- "REV_"
# in trivialName 2 we can still have the decoy-tag in
protAnno$TrivialName2 <- paste(decoyString,  protAnno$TrivialName, sep = "")

# cosmetics
protAnno$TrivialName2 <- gsub(x = protAnno$TrivialName2, pattern = " ",  replacement = "_")

# first join
datJoined_withSampleAnno <- left_join(x = datOK, y= sampleAnno, by="rawFile")
# now not really matched
sum(datJoined_withSampleAnno$nameTag.x == datJoined_withSampleAnno$nameTag.y, na.rm = TRUE)/nrow(datJoined_withSampleAnno)
#datJoined_withSampleAnno$nameTag <- datJoined_withSampleAnno$nameTag.x
unique(datJoined_withSampleAnno$nameTag.y)
datJoined_withSampleAnno$nameTag <- datJoined_withSampleAnno$nameTag.y

# which ones are missing now from the annotation?
newSampleNamesNrawFile_annotation_June2024 <-  datJoined_withSampleAnno |> select(rawFile, nameTag, Location, Time) |> distinct()
#write_tsv(newSampleNamesNrawFile_annotation_June2024, "newSampleNamesNrawFile_annotation_June2024.tsv")

# these will be taken out later (faulty injections or blanks) at line round 286 (when looking for NA~NA) 
newSampleNamesNrawFile_annotation_June2024[is.na(newSampleNamesNrawFile_annotation_June2024$Time),]

# second join with protein annotation (curated manually and categorized)
# in some we find a space at the end before newline.. -> difficult to spot
protAnno$TrivialName <- gsub(x = protAnno$TrivialName, pattern = " $", replacement = "")

datJoined_withSampleAnno <- as.data.frame(datJoined_withSampleAnno)
protAnno <- as.data.frame(protAnno)


# all left joined
datJoined <- left_join(x = datJoined_withSampleAnno, y = protAnno, by = c("prot_acc" = "prot_acc"))

# replace spaces with underscore in Species
datJoined$Species <- gsub(x = datJoined$Species, pattern = " ", replacement = "_")

# work more on protNames and sampleNames and concatenated with ~ tilde
datJoined$myProteins <- paste(datJoined$TrivialName2, datJoined$Species, sep = "~")
datJoined$mySampleName <- paste(datJoined$nameTag, datJoined$Location, datJoined$Time, sep="~")

# Now summarization to get heatmap
# summarize psms by table -> all non human proteins will be NA~NA
# we have to take species into protein name only like this makes sense
dim(relevantTableMat <- table(datJoined$myProteins, datJoined$mySampleName))
#[1] 524 104

# get sum column
head(relevantTableMat)
myPSMsum <- rowSums(relevantTableMat)

# protFDR before filtering
dim(relevantTableMat) # now 526 rows -> 524
get_protFDR(relevantTableMat)
# 1.711 -> 1.718

# filtering on rowSum with PSM threshold
PSMthreshold <- 5

bool_keep <- myPSMsum >= PSMthreshold

# filter mat according to total PSM sum threshold
filtMat <- as.matrix(relevantTableMat[bool_keep,])
dim(filtMat) # back to 366  without  filtering for eValue
# 346 104 with eValue filtering
get_protFDR(filtMat) # 1.093 -> 1.156%



# go for heatmap plotting w gplots
# df[df == 0] <- NA
filtMat[filtMat == 0] <- NA

# split for HUMAN
myKeyword <- "Homo" # maybe here also use HUMAN to not loose contaminants..

bool_homo <- str_count(string = row.names(filtMat), pattern = myKeyword) > 0
sum(bool_homo) #  361 human proteins left -> 342 (365 w/o evalue filter due to new search results)


filtMat[which(bool_homo == FALSE),] # this one NA~NA gathers quite some  psms.. -> these are NON Human proteins

# work on protein row colors!!
myHumanProteins <-data.frame(rownames(filtMat)[bool_homo])

colnames(myHumanProteins) <- "FullDesc"
#myHumanProteins$desc <- sapply(strsplit(myHumanProteins$FullDesc, split = "~"), function(x)x[1])
# check in original file for classes
protAnno$desc <-  gsub(x = protAnno$TrivialName, pattern = " ", replacement = "_")
dim(protAnno)

#we have to make protAnno unique first based on our annotation without sp-Acc
protAnno_unique <- protAnno |> select(TrivialName, Species, Primary, Secondary) |> distinct()
protAnno_unique$TrivialName <-  gsub(x = protAnno_unique$TrivialName, pattern = " ", replacement = "_")

#  here is an issue -> issue solved by removing doubs from input, Keratin_32 classified  as Hair!
colnames(protAnno_unique)[1] <- "desc"
colnames(protAnno_unique)
protAnno_unique$Species <- gsub(x = protAnno_unique$Species,pattern = " ", replacement = "_")
protAnno_unique$FullDesc <- paste(protAnno_unique$desc,protAnno_unique$Species, sep = "~" )
protAnno_unique$desc <- NULL
protAnno_unique$Species <- NULL

myHumanProteins <- left_join(x = myHumanProteins, y = protAnno_unique) #Joining with `by = join_by(desc)`
dim(myHumanProteins)

# Now we have 345 proteins (human) left w/ eValue < 0.05
#
# working on archeological sites for col_side_colors
#

# now with colSideColors
# check unique length before assigning colors
myLocations <- unique(sapply(strsplit(colnames(filtMat), split = "~"), function(x)x[2]))
myTimes <- unique(sapply(strsplit(colnames(filtMat), split = "~"), function(x)x[3]))

# V0ector with 4 levels we have some NAs which ones?
Location <- sapply(strsplit(colnames(filtMat), split = "~"), function(x)x[2])
Times <- sapply(strsplit(colnames(filtMat), split = "~"), function(x)x[3])
table(Location)
table(Times)

# which ones are NAs -> what to do with it?
colnames(filtMat)[Times == "NA"] # this filter can be used
colnames(filtMat)[Location == "NA"] # we want to keep "090_Beth~NA~ARCH_SouthAfrica"

# do we still have the kenya samples?
sum(grepl(x = colnames(filtMat), "enya")) # yes 7


# Define colors for each level
unique(Location)
length(unique(Location))

# Map colors to levels in the data vector
color_vector_loc <- c("darkred", "tomato4", "tomato3", "salmon3", "brown")
assigned_colors_locations <- as.matrix(color_vector_loc[as.numeric(factor(Location))], ncol = 1)
color_vector_times <- c("navajowhite2", "lemonchiffon", "tan", "antiquewhite")
assigned_colors_times <- as.matrix(color_vector_times[as.numeric(factor(Times))], ncol = 1)

assigned_colors <- cbind(assigned_colors_locations, assigned_colors_times)

# here how the colors are defined for the samples
cbind(colnames(filtMat), assigned_colors)


# Create a matrix with multiple colors for each row
row_side_colors_matrix <- matrix(NA, nrow = nrow(filtMat[bool_homo, ]), ncol = 2)

# REVorFW <- c("TRUE", "FALSE")
REVorFW <- str_count(string = myHumanProteins$FullDesc, pattern = "REV_") == 0
sum(REVorFW)
numRev <- str_count(string = myHumanProteins$FullDesc, pattern = "REV_") > 0
sum(numRev)

# Define colors for each level
color_vector <- c("red", "white") # red for decoy! white for fw  --> this should probably be ommited for the final matrix
# Map colors to levels in the data vector
fwOrrevColor <- as.matrix(color_vector[as.numeric(factor(REVorFW))], ncol = 1)

# Now we should use our classifications to come up with better protein colors
specialProteinsClassification_bodyParts <- matrix("other",nrow = nrow(filtMat[bool_homo, ]), ncol = 1 )
table(myHumanProteins$Primary)
specialProteinsClassification_bodyParts[str_count(string = myHumanProteins$Primary, pattern = "Several") > 0] <- "Several"
specialProteinsClassification_bodyParts[str_count(string = myHumanProteins$Primary, pattern = "Saliva|Blood|Digestive|Milk|Urine") > 0] <- "BodyFluid"
specialProteinsClassification_bodyParts[str_count(string = myHumanProteins$Primary, pattern = "Skin|Hair|Muscle") > 0] <- "BodyPart"
table(specialProteinsClassification_bodyParts)

# Define colors for each level.
# BodyPart should be tan4, BodyFluid should be red, and Several should be khaki2, Other should be black
length(unique(specialProteinsClassification_bodyParts))

# several = khaki2
# BodyFluid = red
# Bodypart = tan4
# other = black
color_vector <- c("red2", "tan4", "black", "khaki2")

# Map colors to levels in the data vector
specialProteinsClassification_bodyParts <- as.matrix(color_vector[as.numeric(factor(specialProteinsClassification_bodyParts))], ncol = 1)

# work on second column
table(myHumanProteins$Secondary)

# @Shevan how shall we handle Male/Female case?
specialProteinsClassification_sex <- matrix("other",nrow = nrow(filtMat[bool_homo, ]), ncol = 1 )
specialProteinsClassification_sex[str_count(string = myHumanProteins$Secondary, pattern = "Male") > 0] <- "male"
specialProteinsClassification_sex[str_count(string = myHumanProteins$Secondary, pattern = "Female") > 0] <- "female"
specialProteinsClassification_sex[str_count(string = myHumanProteins$Secondary, pattern = "Male/Female") > 0] <- "Male_Female"
table(specialProteinsClassification_sex)
table(myHumanProteins$Secondary)


# Define colors for each level
## I would like "male" to be darklue, "female" to be gold and the male_female to be forestgreen, and blanks to be white.
table(myHumanProteins$Secondary)
length(unique(specialProteinsClassification_sex))
(unique(specialProteinsClassification_sex))
color_vector <- c("gold", "darkblue", "forestgreen", "white")

# Map colors to levels in the data vector
specialProteinsClassification_sex <- as.matrix(color_vector[as.numeric(factor(specialProteinsClassification_sex))], ncol = 1)



# bind together
row_side_colors_matrix <- cbind(fwOrrevColor, specialProteinsClassification_bodyParts, specialProteinsClassification_sex)
dim(row_side_colors_matrix)
dim(myHumanProteins)

# here to check what proteins are assigned to what colors w/ its categories before..
tt <- cbind(myHumanProteins, row_side_colors_matrix)
View(tt)

# color scheme for proteins
colorNproteins <- data.frame(cbind(row_side_colors_matrix, myHumanProteins))
#write_tsv(colorNproteins, file = "ColorSchemeForProteins_new.tsv")

# plotting
(fN <- paste(fgczProject, "_", descri, "_bigPage", ".pdf",sep = ""))
pdf(fN, 35,35)
# Create the heatmap with row side colors list
heatmap.plus(
  filtMat[bool_homo, ],
  Rowv = FALSE,
  distfun = dist_no_na,
  margins = c(20, 25),
  RowSideColors = row_side_colors_matrix,
  ColSideColors = assigned_colors,
  main="Assigned PSM heatmap from all Human Proteins w multiple row and colSide Colors"
)
dev.off()


# identify blanks and take them out -> key is NA~NA for location and times!
colnames(filtMat[bool_homo, ])
idx_takeout <- which(grepl(x = colnames(filtMat[bool_homo, ]), "NA~NA"))
# to be taken out
colnames(filtMat[bool_homo, idx_takeout]) 

# to keep
colnames(filtMat[bool_homo,-which(grepl(x = colnames(filtMat[bool_homo, ]), "NA"))])


# plotting
(fN <- paste(fgczProject, "_", descri, "_min5psms", ".pdf",sep = ""))
pdf(fN, 35,35)
# Create the heatmap with row side colors list
heatmap.plus(
  filtMat[bool_homo, -idx_takeout],
  Rowv = FALSE,
  distfun = dist_no_na,
  margins = c(20, 25),
  RowSideColors = row_side_colors_matrix,
  ColSideColors = assigned_colors[-idx_takeout,],
  main="Assigned PSM heatmap from all Human Proteins w multiple row and colSide Colors"
)
dev.off()



# do more on filter w/ psms
filtMat[which(bool_homo == FALSE),]

# reduce color vectors
dim(row_side_colors_matrix)
nonHumanProt <- which(bool_homo == FALSE)
row_side_colors_matrix2 <- row_side_colors_matrix[-nonHumanProt, ]
assigned_colors2 <- assigned_colors[-idx_takeout,]
# 

# generate a filtered matrix again
dim(filtMat)
newFiltMat <- filtMat[bool_homo, -idx_takeout]
dim(newFiltMat)

# filtering on rowSum with PSM threshold
nrow(newFiltMat)
myPSMsum <- rowSums(newFiltMat, na.rm = TRUE)
PSMthreshold_10 <- 10
bool_keep_10 <- myPSMsum >= PSMthreshold_10
sum(bool_keep_10)

# filter mat according to total PSM sum threshold
filtMat_min10 <- newFiltMat[bool_keep_10,]

# 
rscol2 <- row_side_colors_matrix[bool_keep_10, ]
(fN <- paste(fgczProject, "_", descri, "_noBlanks_min10psms", ".pdf",sep = ""))
pdf(fN, 31,31)
# Create the heatmap with row side colors list
heatmap.plus(
  filtMat_min10,
  Rowv = FALSE,
  distfun = dist_no_na,
  margins = c(20, 25),
  RowSideColors = rscol2,
  ColSideColors = assigned_colors2,
  main="Assigned PSM heatmap from all Human Proteins w multiple row and colSide Colors"
)
dev.off()


# Create the heatmap with row side colors list

# more filtering
PSMthreshold_20 <- 20
bool_keep_20 <- myPSMsum >= PSMthreshold_20
sum(bool_keep_20)

# filter mat according to total PSM sum threshold
filtMat_min20 <- newFiltMat[bool_keep_20,]


rscol3 <- row_side_colors_matrix[bool_keep_20, ]
(fN <- paste(fgczProject, "_", descri, "_noBlanks_min20psms", ".pdf",sep = ""))
pdf("p27048_heatmap_HUMANonly_betterColors_noBlanks_min20psms.pdf", 25,25)
heatmap.plus(
  filtMat_min20,
  Rowv = FALSE,
  distfun = dist_no_na,
  margins = c(20, 25),
  RowSideColors = rscol3,
  ColSideColors = assigned_colors2,
  main="Assigned PSM heatmap from all Human Proteins w multiple row and colSide Colors"
)
dev.off()











