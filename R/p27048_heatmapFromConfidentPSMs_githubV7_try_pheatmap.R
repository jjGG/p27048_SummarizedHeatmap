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
descri <- "SA_and_modern_reorderColNrows"

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
sampleAnno

# 2025-05-20: Take out some sites to only keep South Africa and Ethnographic
unique(sampleAnno$Location)
sampleAnno <- sampleAnno |> filter(Location == "ARCH_SouthAfrica" | Location == "Ethnographic")


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
# (fN <- paste(fgczProject, "_", descri, "_newSampleNamesNrawFile_annotation_", Sys.Date(), ".tsv", sep = ""))
#write_tsv(newSampleNamesNrawFile_annotation_June2024, "newSampleNamesNrawFile_annotation_June2024.tsv")
# these will be taken out later (faulty injections or blanks) at line round 286 (when looking for NA~NA) 
# 2025-05-20: Yes this is OK only Ethnographic and South Africa are still in the game
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
#[1] 524 70 (69 + NA)

# get sum column
head(relevantTableMat)
myPSMsum <- rowSums(relevantTableMat)

# protFDR before filtering
dim(relevantTableMat) # now 526 rows -> 524
get_protFDR(relevantTableMat)
# 1.711 -> 1.718 # 2025-05-20: still the same! (1.718)

# filtering on rowSum with PSM threshold
PSMthreshold <- 5
bool_keep <- myPSMsum >= PSMthreshold

# filter mat according to total PSM sum threshold
filtMat <- as.matrix(relevantTableMat[bool_keep,])

dim(filtMat) # back to 366  without  filtering for eValue
# 346 104 with eValue filtering
get_protFDR(filtMat) # 1.093 -> 1.156%
# 2025-05-20: 1.163 % now!



# now eliminate REVs
decoyBool <- grepl(x = rownames(filtMat), pattern = "REV_")
nondecoy_bool <- decoyBool == FALSE

dim(filtMat)
filtMat <- filtMat[nondecoy_bool, ]
dim(filtMat)


# go for heatmap plotting w gplots
# df[df == 0] <- NA
filtMat[filtMat == 0] <- NA

# split for HUMAN
myKeyword <- "Homo" # maybe here also use HUMAN to not loose contaminants..

bool_homo <- str_count(string = row.names(filtMat), pattern = myKeyword) > 0
sum(bool_homo) #  361 human proteins left -> 342 (365 w/o evalue filter due to new search results)
# 339 human proteins left

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

# Now we have 345 proteins (human) left w/ eValue < 0.05 # update: 2025-05-20: 339 proteins (human) left w/ eValue < 0.05 !
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
# fix one NA for one sample?? why is this NA?
Times[36] <- "not"
location_times <- cbind(Location, Times)


# which ones are NAs -> what to do with it?
colnames(filtMat)[Times == "NA"] # this filter can be used
colnames(filtMat)[Location == "NA"] # we want to keep "090_Beth~NA~ARCH_SouthAfrica"

# do we still have the kenya samples?
sum(grepl(x = colnames(filtMat), "enya")) # yes 7 --> 0 # No - not anymore


# Define colors for each level
unique(Location)
length(unique(Location))


#I would also like to make two copies of this heatmaps, one for the main text and one for the SI.
color_vector_loc <- c("darkred",  "tomato3", "brown")
assigned_colors_locations <- as.matrix(color_vector_loc[as.numeric(factor(Location))], ncol = 1)
color_vector_times <- c("navajowhite2", "lemonchiffon", "tan", "antiquewhite")
assigned_colors_times <- as.matrix(color_vector_times[as.numeric(factor(Times))], ncol = 1)

#assigned_colors <- cbind(assigned_colors_locations, assigned_colors_times)
assigned_colors <- cbind("white", assigned_colors_locations)

# here how the colors are defined for the samples
ss <- as.data.frame(cbind(colnames(filtMat), assigned_colors))
View(ss)

# to-be-done Last entry is NA~NA~NA .. take it out
# Create a matrix with multiple colors for each row
row_side_colors_matrix <- matrix(NA, nrow = nrow(filtMat[bool_homo, ]), ncol = 2)

# REVorFW <- c("TRUE", "FALSE")
REVorFW <- str_count(string = myHumanProteins$FullDesc, pattern = "REV_") == 0
sum(REVorFW) # 339
numRev <- str_count(string = myHumanProteins$FullDesc, pattern = "REV_") > 0
sum(numRev) # 0

# Define colors for each level
color_vector <- c("white", "red") # red for decoy! white for fw  --> this should probably be ommited for the final matrix
# Map colors to levels in the data vector
fwOrrevColor <- as.matrix(color_vector[as.numeric(factor(REVorFW))], ncol = 1)
table(fwOrrevColor)

# Now we should use our classifications to come up with better protein colors
specialProteinsClassification_bodyParts <- matrix("other",nrow = nrow(filtMat[bool_homo, ]), ncol = 1 )
table(myHumanProteins$Primary)
specialProteinsClassification_bodyParts[str_count(string = myHumanProteins$Primary, pattern = "Several") > 0] <- "Several"
specialProteinsClassification_bodyParts[str_count(string = myHumanProteins$Primary, pattern = "Saliva|Blood|Digestive|Milk|Urine") > 0] <- "BodyFluid"
specialProteinsClassification_bodyParts[str_count(string = myHumanProteins$Primary, pattern = "Skin|Hair|Muscle") > 0] <- "BodyPart"
table(specialProteinsClassification_bodyParts)



bodypart_order <- order(specialProteinsClassification_bodyParts)
# specialProteinsClassification_bodyParts[bodypart_order]
# Define colors for each level.
# BodyPart should be tan4, BodyFluid should be red, and Several should be khaki2, Other should be black

# several = khaki2
# BodyFluid = red
# Bodypart = tan4
# other = black
color_vector <- c("red2", "tan4", "black", "khaki2")

# Map colors to levels in the data vector
#specialProteinsClassification_bodyParts <- as.matrix(color_vector[as.numeric(factor(specialProteinsClassification_bodyParts))], ncol = 1)
#table(specialProteinsClassification_bodyParts)


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
#color_vector <- c("gold", "darkblue", "forestgreen", "white")

special_protein_class <- cbind(specialProteinsClassification_bodyParts, specialProteinsClassification_sex)

# Map colors to levels in the data vector
#specialProteinsClassification_sex <- as.matrix(color_vector[as.numeric(factor(specialProteinsClassification_sex))], ncol = 1)
table(specialProteinsClassification_sex)

head(location_times)
head(special_protein_class)


filtMat
df <- as.data.frame.matrix(filtMat)
View(df)

# get rid of NA~NA protein
NAprotein <- which(rownames(df) == "NA~NA")
df <- df[-NAprotein, ] # take out NA~NA protein

# get rid of NAsite
NAsite <- which(colnames(df) == "NA~NA~NA")
df <- df[, -NAsite] # take out NA~NA~NA site

# go for pheatmap
library(pheatmap)


# Customize NA color
pheatmap(df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         na_col = "gray")

# sort some
nums <- colSums(df, na.rm = TRUE)

dim(df)
# go for better sorting
col_order <- order(nums, decreasing = TRUE)
row_order <- order(specialProteinsClassification_bodyParts)

# Reorder your matrix
ordered_matrix <- df[row_order, col_order]

# Plot without clustering
pheatmap(ordered_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         na_col = "gray")


# simple color vector
# Just one annotation per axis
# fix minus in open-air
location_times[,2] <- gsub(x = location_times[,2], pattern = "-", replacement = "_")
head(location_times)
head(special_protein_class)
colnames(special_protein_class) <- c("Tissue", "Sexclassification")

row_colors <- data.frame(special_protein_class,
                         row.names = rownames(df))

col_colors <- data.frame(location_times[-NAsite,],
           row.names = colnames(df))

# col_colors <- data.frame(
#   Group = assigned_colors[-NAsite,2],
#   row.names = colnames(df))

pheatmap(df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = row_colors,
         annotation_col = col_colors,
         na_col = "gray")




# Define specific colors for each category
my_colors <- list(
  Tissue = c(other = "black", BodyPart = "khaki2", BodyFluid = "red2", Several = "tan4"),
  Location = c(Ethnographic = "tomato3", ARCH_SouthAfrica = "darkred"),
  Sexclassification = c(female="pink", male="lightblue", Male_Female="forestgreen", other="white"),
  Times = c(Curated = "darkgreen", not ="blue", Open_air = "lightblue", Rockshelter = "orange")
)

pheatmap(df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = row_colors,
         annotation_col = col_colors,
         annotation_colors = my_colors,
         na_col = "gray")

# save the heatmap
(fN <- paste(fgczProject, "_", descri, "_pheatmap_", Sys.Date(), ".pdf", sep = ""))
pdf(file = fN, width = 20, height =40)
pheatmap(df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = row_colors,
         annotation_col = col_colors,
         annotation_colors = my_colors,
         na_col = "gray")
dev.off()



# now do some reordering of columns and rows
# go for better sorting
col_order <- order(nums, decreasing = TRUE)
row_order <- order(specialProteinsClassification_bodyParts)

# Reorder your matrix
ordered_matrix <- df[row_order, col_order]

# Plot without clustering
(fN <- paste(fgczProject, "_", descri, "_pheatmap_reordered_", Sys.Date(), ".pdf", sep = ""))
pdf(file = fN, width = 20, height =40)
pheatmap(ordered_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = row_colors,
         annotation_col = col_colors,
         annotation_colors = my_colors,
         na_col = "gray")
dev.off()

# now apply some more psmThreshold
psmThreshold <- 10
# filtering on rowSum with PSM threshold
psmSumOrderedMat <- rowSums(ordered_matrix, na.rm = TRUE)
bool_keep <- psmSumOrderedMat >= psmThreshold
length(bool_keep)

# filter mat according to total PSM sum threshold
filtOrderedMat <- ordered_matrix[bool_keep,]

# Plot without clustering
(fN <- paste(fgczProject, "_", descri, "_pheatmap_reordered_psm",psmThreshold,"_", Sys.Date(), ".pdf", sep = ""))
pdf(file = fN, width = 20, height =30)
pheatmap(filtOrderedMat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = row_colors,
         annotation_col = col_colors,
         annotation_colors = my_colors,
         na_col = "gray")
dev.off()



# now apply some more psmThreshold
psmThreshold <- 25
# filtering on rowSum with PSM threshold
psmSumOrderedMat <- rowSums(ordered_matrix, na.rm = TRUE)
bool_keep <- psmSumOrderedMat >= psmThreshold
length(bool_keep)

# filter mat according to total PSM sum threshold
filtOrderedMat <- ordered_matrix[bool_keep,]

# Plot without clustering
(fN <- paste(fgczProject, "_", descri, "_pheatmap_reordered_psm",psmThreshold,"_", Sys.Date(), ".pdf", sep = ""))
pdf(file = fN, width = 20, height =25)
pheatmap(filtOrderedMat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = row_colors,
         annotation_col = col_colors,
         annotation_colors = my_colors,
         na_col = "gray")
dev.off()










