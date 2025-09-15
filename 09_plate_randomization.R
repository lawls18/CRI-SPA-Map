library(tidyverse)
library(readxl)
library(lme4)
library(qvalue)
library(lmerTest)
library(wpm)
library(ggplate)

#The purpose of this code is to randomize source plates for CRI-SPA-Map/Validation phenotyping.
#This code was made for a specific use case and is not adapted for general use. 
#Included in this code is the following functionality:
# -Reading original plate layouts, compiling them, inserting wild type strains, and then randomizing them with the WPM package
# -Checking for proper wild type distribution across paltes
# -Checking for proper distribution of isolates on edges
# -Renaming relevant well information
# -Creating PIXL rearray files to create the randomization plates (96 and 384 format)
# -Generating visual maps of randomization plates


#The path to a folder containing the original plate layouts of all plates being pinned from
origPath <- "/Users/samuelamidon/Desktop/Aim 1/Randomization/Randomization 3 From DWP/Orig_Layouts_TargetsOnly"



#Map Combiner - combines all original plates into a single data frame
origLayoutNames <- dir(origPath)
numFiles <- length(origLayoutNames)
origLayouts <- list()
for (i in 1:numFiles) {
  layoutPath <- str_c(origPath,'/',origLayoutNames[i])
  layout <- read_excel(layoutPath)
  layout$...5<-NULL
  layout$...6<-NULL
  layout <- layout %>% rename(Plate_ID="...7")
  origLayouts[[i]] <- layout
}
origLayoutsCombined <- origLayouts[[1]]
origLayoutsCombined$Isolate <- as.character(origLayoutsCombined$Isolate)
for (i in 2:numFiles) {
  layout <- origLayouts[[i]]
  layout$Isolate <- as.character(layout$Isolate)
  origLayoutsCombined <- bind_rows(origLayoutsCombined,layout)
}



#Sequence Matcher - filters out all blanks strains that won't be used
plateReaderMap <- read_excel("/Users/samuelamidon/Desktop/Aim 1/Randomization/Randomization 3 From DWP/Target_Genes.xlsx") #An Excel sheet containing the names of all strains to be randomized
plateReaderMap <- plateReaderMap %>% mutate(geneIsolate = str_c(Gene,"_",Isolate))
origSourceWells <- origLayoutsCombined %>% filter(!is.na(origLayoutsCombined$'Gene deletion'))
origSourceWells <- origSourceWells %>% mutate(geneIsolate = str_c(origSourceWells$'Gene deletion',"_",Isolate))
origSourceWells <- origSourceWells %>% filter(origSourceWells$geneIsolate %in% plateReaderMap$geneIsolate)



#Wild Type Inserter - disperses wild types among the original strains
wtMap <- read_excel("/Users/samuelamidon/Desktop/Aim 1/Randomization/Randomization 3 From DWP/wt_map.xlsx")
wtMap$Isolate <- as.character(wtMap$Isolate)
wtMap <- wtMap %>% mutate(geneIsolate = str_c(wtMap$'Gene deletion',"_",Isolate))
targetPlateCount <- 8 #Define the amount of wild types to insert
totalWells <- targetPlateCount*96
totalStrains <- nrow(origSourceWells)
availableWells <- totalWells - totalStrains
wtReps <- ceiling(availableWells/nrow(wtMap))
for (i in 1:wtReps) {
  origSourceWells <- bind_rows(origSourceWells,wtMap)
}
origSourceWells<-origSourceWells[1:totalWells,]



#Randomization Preparation - shuffles strains and creates source files for WPM randomization
targetFolder <-"/Users/samuelamidon/Desktop/Aim 1/Randomization/Randomization 3 From DWP/Randomization_Source_Files"
sampVec <- sample(1:totalWells)
origSourceWells <- origSourceWells[sampVec,] #Shuffle all plates before randomizing
wpmPrep <- data.frame(Identifier=str_c(origSourceWells$Wells,"_",origSourceWells$Isolate,"_",origSourceWells$Plate_ID))
wpmPrep$Gene <- origSourceWells$'Gene deletion'
for (i in 0:(targetPlateCount-1)) { #Separate shuffled strains into 96 strain groups for use in WPM package
  chunk <- wpmPrep[(i*96+1):((i+1)*96),]
  write.csv(chunk,str_c(targetFolder,"/","P",as.character(i+1),".csv"),row.names=FALSE)
}



#Randomizer - applies WPM randomization to source files
wpmPrepNames <- dir(targetFolder)
array <- c('A','B')
randArrays <- list()
arrayCount <- 1
for (i in 1:targetPlateCount) {
  csvName <- str_c(targetFolder,"/",wpmPrepNames[i])
  for (j in 1:2) {
    imported_csv <- wpm::convertCSV(csvName,
                                    gp_field = 'Gene',
                                    header = TRUE,
                                    sep = ',')
    wpm_result <- wpm::wrapperWPM(user_df = imported_csv$df_wpm,
                                  plate_dims = list(8,12),
                                  nb_plates = 1,
                                  spatial_constraint = 'NEWS',
                                  max_iteration = 50)
    wpm_result$Plate_ID <- str_c(as.character(i),"_",array[j])
    randArrays[[arrayCount]] <- wpm_result
    arrayCount<-arrayCount+1
  }
}



#Randomization Processer - processes output of WPM randomization
randArraysComp <- randArrays[[1]] #Initialize compilation of random arrays
for (i in 2:length(randArrays)) { #Iterate through all random arrays
  randArraysComp <- bind_rows(randArraysComp,randArrays[[i]]) #Add each random array to the compilation
}
randArraysComp <- randArraysComp %>% mutate(splitSample = strsplit(Sample,"_")) #Add split sample info column
range <- 1:nrow(randArraysComp) #Set range for sapply to iterate through
randArraysComp <- randArraysComp %>% mutate(Isolate = sapply(range,function(x) { #Iterate through range
  randArraysComp[x,]$splitSample[[1]][2] #Grab isolate number from split sample info
}))
randArraysComp <- randArraysComp %>% mutate(Source_Plate = sapply(range,function(x) { #Iterate through range
  randArraysComp[x,]$splitSample[[1]][3] #Grab plate ID from split sample info
}))




#Wild Type Fixer - redistribute wild types so they are well distributed across plates
wellSwapper <- function(df,tarSwap,sourSwap,tarPlate,sourPlate,swapPlates=FALSE) { #Function takes input data frame and two entire rows from that frame source and target and a plate ID. Swaps the source and target well information in df by plate. One swap per call.
  res <- df
  res$Well[res$Sample==tarSwap$Sample&res$Plate_ID==tarPlate&res$ID==tarSwap$ID] <- sourSwap$Well #Insert source well to target well
  res$Well[res$Sample==sourSwap$Sample&res$Plate_ID==sourPlate&res$ID==sourSwap$ID] <- tarSwap$Well #Insert target well to source well
  res$Row[res$Sample==tarSwap$Sample&res$Plate_ID==tarPlate&res$ID==tarSwap$ID] <- sourSwap$Row #Insert source row num to target row num
  res$Row[res$Sample==sourSwap$Sample&res$Plate_ID==sourPlate&res$ID==sourSwap$ID] <- tarSwap$Row #.......
  res$Column[res$Sample==tarSwap$Sample&res$Plate_ID==tarPlate&res$ID==tarSwap$ID] <- sourSwap$Column
  res$Column[res$Sample==sourSwap$Sample&res$Plate_ID==sourPlate&res$ID==sourSwap$ID] <- tarSwap$Column
  if (swapPlates == TRUE) {
    res$Plate_ID[res$Sample==tarSwap$Sample&res$Plate_ID==tarPlate&res$ID==tarSwap$ID] <- sourSwap$Plate_ID
    res$Plate_ID[res$Sample==sourSwap$Sample&res$Plate_ID==sourPlate&res$ID==sourSwap$ID] <- tarSwap$Plate_ID
  }
  res #Return a new data frame with swapped well info
}

numWT <- sum(randArraysComp$Group=="wt")/2 #Calculate total wild types
wtPerPlate <- floor(numWT/targetPlateCount) #Calculate possible wild types per plate

randArraysComp <- randArraysComp %>% mutate(Plate_Num = substr(Plate_ID,1,1)) #Create target plate number variable
randArraysWT <- randArraysComp[randArraysComp$Group=="wt",] #Extract all wild type rows
WTCounts <- as.data.frame(table(randArraysWT$Plate_Num)/2) #Get occurrences of wild types on each plate
lowWTCounts <- WTCounts[WTCounts$Freq<(wtPerPlate),] #Get list of wild type counts less than the possible distribution
highWTCounts <- WTCounts[WTCounts$Freq>(wtPerPlate),] #Get list of counts above the possible distribution

while (nrow(lowWTCounts)>=1) { #Repeat until no more low wild type counts exist
  randArraysComp <- randArraysComp %>% mutate(Plate_Num = substr(Plate_ID,1,1)) #Recalculate (can maybe only happen at end of loop)
  randArraysWT <- randArraysComp[randArraysComp$Group=="wt",] #Recalculate (can maybe only happen at end of loop)
  WTCounts <- as.data.frame(table(randArraysWT$Plate_Num)/2) #Recalculate (can maybe only happen at end of loop)
  lowWTCounts <- WTCounts[WTCounts$Freq<(wtPerPlate),] #Recalculate (can maybe only happen at end of loop)
  highWTCounts <- WTCounts[WTCounts$Freq>(wtPerPlate),] #Recalculate (can maybe only happen at end of loop)
  
  lowPlate <- lowWTCounts$Var1[1] #Get plate number of first low wild type count plate
  highPlate <- sample(highWTCounts$Var1,1) #Pick random plate from high count plates
  lowPlate_IDA <- str_c(lowPlate,"_A") #Concatenate A onto plate num
  highPlate_IDA <- str_c(highPlate,"_A") #see above
  
  lowWTSwapA <- slice_sample(randArraysComp[randArraysComp$Plate_ID==lowPlate_IDA & randArraysComp$Group!="wt",]) #Pick a random NON WT row from the plate that has low wt count to swap with that of high wt count
  highWTSwapA <- slice_sample(randArraysComp[randArraysComp$Plate_ID==highPlate_IDA & randArraysComp$Group=="wt",]) #Pick a random WT row from plate that has high wt count to swap with low wt count
  lowWTSwapB <- randArraysComp[randArraysComp$Sample==lowWTSwapA$Sample&randArraysComp$ID==lowWTSwapA$ID&randArraysComp$Plate_ID!=lowPlate_IDA,] #Repeat for plate B, use same gene that was picked for low plate A - do this so that isolates still exist on the same plate number, one on A one on B
  highWTSwapB <- randArraysComp[randArraysComp$Sample==highWTSwapA$Sample&randArraysComp$ID==highWTSwapA$ID&randArraysComp$Plate_ID!=highPlate_IDA,] #Repeat for plate B, use same gene that was picked for high plate A
  
  randArraysComp <- wellSwapper(randArraysComp,lowWTSwapA,highWTSwapA,str_c(lowPlate_IDA),str_c(highPlate_IDA),swapPlates=TRUE) #Swap plate A isolate well info and plate IDs
  randArraysComp <- wellSwapper(randArraysComp,lowWTSwapB,highWTSwapB,str_c(lowPlate,"_","B"),str_c(highPlate,"_","B"),swapPlates=TRUE) #Swap plate B isolate well info and plate IDs
  
  lowWTSwapA <- NULL #Reset all swap variables
  highWTSwapA <- NULL
  lowWTSwapB <- NULL
  highWTSwapB <- NULL
  
  randArraysComp <- randArraysComp %>% mutate(Plate_Num = substr(Plate_ID,1,1)) #Recalculate so while loop is run on proper lowWTCounts
  randArraysWT <- randArraysComp[randArraysComp$Group=="wt",]
  WTCounts <- as.data.frame(table(randArraysWT$Plate_Num)/2)
  lowWTCounts <- WTCounts[WTCounts$Freq<(wtPerPlate),]
}
randArraysComp <- randArraysComp %>% mutate(Plate_Num = substr(Plate_ID,1,1)) #Recalculate plate numbers to avoid issues
randArraysComp <- randArraysComp[order(randArraysComp$Plate_ID),] #Reorder by plate ID to put swapped rows back with the right plates



#Edge Duplicate Checker - finds cases where isolates are on an edge in both random arrays and redistributes them
#Setup
randArraysCompOrig <- randArraysComp #Save compilation before edge corrections
randEdges <- randArraysComp %>% filter((Row==1|Row==8|Column==1|Column==12)) #Get list of all edge positions from compilation
randCenters <- randArraysComp %>% filter(!(Row==1|Row==8|Column==1|Column==12)) #Get list of all center positions from compilation
duplicateEdges <- as.data.frame(table(randEdges$Sample)[table(randEdges$Sample)==2&!grepl("W",as.data.frame(table(randEdges$Sample))$Var1)]) #Get list of all isolates that are on an edge twice
#tables are used to count occurances, occurance of 2 means an isolate was on an edge twice
#grepl used to filter wild types with "plate" W
duplicateCenters <- as.data.frame(table(randCenters$Sample)[table(randCenters$Sample)==2&!grepl("W",as.data.frame(table(randCenters$Sample))$Var1)]) #Get list of all isolates that are on a center well twice
duplicateCenters <- randCenters[randCenters$Sample %in% duplicateCenters$Var1,] #Get additional info about center wells

for (i in 1:(nrow(duplicateEdges)-1)) { #Iterate through all but the last duplicate edge
  randEdges <- randArraysComp %>% filter((Row==1|Row==8|Column==1|Column==12)) #Get new list of edges
  randCenters <- randArraysComp %>% filter(!(Row==1|Row==8|Column==1|Column==12)) #Get new list of centers 
  duplicateEdges <- as.data.frame(table(randEdges$Sample)[table(randEdges$Sample)==2&!grepl("W",as.data.frame(table(randEdges$Sample))$Var1)]) #Get new list of duplicate edges
  duplicateCenters <- as.data.frame(table(randCenters$Sample)[table(randCenters$Sample)==2&!grepl("W",as.data.frame(table(randCenters$Sample))$Var1)]) #Get list of new duplicate centers
  duplicateCenters <- randCenters[randCenters$Sample %in% duplicateCenters$Var1,] #Get more info on centers
  
  dupSample <- as.character(duplicateEdges$Var1[1]) #Pick isolate thats on an edge twice to fix (always use first because edge list is updated every iteration)
  plate <- sample(randArraysComp$Plate_ID[randArraysComp$Sample==dupSample],1) #Randomly select occurance of dupSample to fix (always between A or B, which should both have dupSample on an edge...) Also get plate number, plate looks like 1_A
  targetInd <- sample(1:sum(duplicateCenters$Plate_ID==plate),1) #Pick a random index from the number of eligible swap positions (isolates on an edge twice with the correct plate)
  targetSwap <- duplicateCenters[duplicateCenters$Plate_ID==plate,][targetInd,] #Get all info about the chosen isolate to swap with the duplicate edge isolate, from plate and random index
  sourceSwap <- randEdges[randEdges$Sample==dupSample&randEdges$Plate_ID==plate,] #Get all info about the isolate that needs to be fixed
  
  randArraysComp <- wellSwapper(randArraysComp,targetSwap,sourceSwap,plate,plate) #Swap target and source isolate well information
  
  targetSwap <- NULL #Reset target isolate - issues have occurred where well info gets reset mid iteration. This may or may not have been the fix. Anyways, doesn't hurt
  sourceSwap <- NULL #Reset source isolate
}
#Finish fixing the last duplicate edge. Issue occurred where table function does not properly count occurrences if only one isolate is available to count. Because of this, last duplicate edge must be dealt with slightly differently.
dupSample <- as.character(duplicateEdges$Var1[2]) #Get duplicate edge sample info - duplicateEdges[1] has already been dealt with in the loop, now just need to fix duplicateEdges[2]

#Recalculate everything but duplicateEdges because that is the one that is depleted
randEdges <- randArraysComp %>% filter((Row==1|Row==8|Column==1|Column==12)) #Not really necessary to recalculate this
randCenters <- randArraysComp %>% filter(!(Row==1|Row==8|Column==1|Column==12)) #Need to recalculate this because randCenters changes during the last iteration above and is not recalculated in the loop
duplicateCenters <- as.data.frame(table(randCenters$Sample)[table(randCenters$Sample)==2&!grepl("W",as.data.frame(table(randCenters$Sample))$Var1)]) #Same here, a duplicate center exists in the old instance of this data frame that must be removed by recalculating
duplicateCenters <- randCenters[randCenters$Sample %in% duplicateCenters$Var1,] #Same here
plate <- sample(randArraysComp$Plate_ID[randArraysComp$Sample==dupSample],1) #All below is the same as in the loop
targetInd <- sample(1:sum(duplicateCenters$Plate_ID==plate),1)
targetSwap <- duplicateCenters[duplicateCenters$Plate_ID==plate,][targetInd,]
sourceSwap <- randEdges[randEdges$Sample==dupSample&randEdges$Plate_ID==plate,]
randArraysComp <- wellSwapper(randArraysComp,targetSwap,sourceSwap,plate,plate)



#Edge Frequency Checker - finds cases where isolates have higher frequency of edge positions and redistributes them
numIsolatesPerGene <- 7 #Number of isolates per gene
nIso <- numIsolatesPerGene #Set isolates oer gene to a shorter variable?
randEdges <- randArraysComp %>% filter((Row==1|Row==8|Column==1|Column==12)) #get new list of centers
randCenters <- randArraysComp %>% filter(!(Row==1|Row==8|Column==1|Column==12)) #Get new list of edges
maxEdgeCounts <- as.data.frame(table(randEdges$Group)[table(randEdges$Group)==nIso]) #Initial list of genes with max edge occurances (all isolates on an edge once)
lowEdgeCounts <- as.data.frame(table(randEdges$Group)[table(randEdges$Group)<=(nIso-4)]) #Initial list of genes with low edge occurances

while (nrow(maxEdgeCounts)>2) {
  randEdges <- randArraysComp %>% filter((Row==1|Row==8|Column==1|Column==12)) #Get new list of edges
  randCenters <- randArraysComp %>% filter(!(Row==1|Row==8|Column==1|Column==12)) #Get new list of centers
  
  if (sum(table(randEdges$Group)<=(nIso-4))>1) {
    lowEdgeCounts <- as.data.frame(table(randEdges$Group)[table(randEdges$Group)<=(nIso-4)])
  } else if (sum(table(randEdges$Group)<=(nIso-4))==1) {
    lowEdgeCounts <- as.data.frame(table(randEdges$Group)[table(randEdges$Group)<=(nIso-3)])
  }
  midEdgeCountsBackup <- as.data.frame(table(randEdges$Group)[table(randEdges$Group)==nIso-3]) #Get backup set of only mid edge counts of 4 in case no low edge counts of 1,2,3 match the chosen source swap's plate
  maxEdgeCounts <- as.data.frame(table(randEdges$Group)[table(randEdges$Group)==nIso]) #Get new list of genes with all isolates on an edge once
  
  duplicateCenters <- as.data.frame(table(randCenters$Sample)[table(randCenters$Sample)==2&!grepl("W",as.data.frame(table(randCenters$Sample))$Var1)]) #Get list of swappable positions (isolates that are in center positions twice)
  duplicateCenters <- randCenters[randCenters$Sample %in% duplicateCenters$Var1,] #Get isolate info for all eligible center positions
  
  sourceSwapGene <- maxEdgeCounts$Var1[1] #Get the first gene whose isolates are all on an edge once
  edgeCands <- randEdges[randEdges$Group==sourceSwapGene,] #List all instances of source gene that are on edges
  sourceSwap <- slice_sample(edgeCands,n=1) #Pick random gene instance from edge candidates
  plate <- sourceSwap$Plate_ID #Get plate ID of chosen gene instance
  
  n<-0 #Set iterating variable to 0. While loop proceeds until n is set to 1
  x<-1 #Set counter to prevent infinite iterations
  while (n==0 & x<=1000) { #Iterates until a suitable swap isolate has been chosen, or until too many iterations have passed.
    targetSwapGeneCand <- case_when( #Picking a swap isolate based on availability. Sample randomly from low or mid edge count data frames
      x>=100 ~ sample(midEdgeCountsBackup$Var1,1), #Only pick from edge counts of 4 if 2 and 3 have been depleted
      .default = sample(lowEdgeCounts$Var1,1) #Favored outcome. Want to replace high edge counts with low edge coutnts, 2 or 3.
    )
    if (x==100) {print("Exceeded 100 iterations without finding an eligible isolate, using mid edge counts")}
    # targetSwapGeneCand <- sample(lowEdgeCounts$Var1,1) #Pick a random gene candidate from all genes with low edge counts
    centerCands <- duplicateCenters[duplicateCenters$Group==targetSwapGeneCand,] #Get all isolates of candidate gene that are in center positions both instances
    centerCands <- centerCands[centerCands$Plate_ID==plate,] #Filter isolates that are not located on the same plate as sourceSwap isolate
    n <- case_when(nrow(centerCands)>=1 ~ 1, .default = 0) #If no isolates were found on the target plate, continue searching
    x <- x+1
    if (x==1000) {print("Exceeded 1000 iterations without finding an eligible isolate, even from backup group. Try restarting?")}
  }
  targetSwap <- slice_sample(centerCands,n=1) #Pick target isolate to swap with source from all eligible isolates
  
  randArraysComp <- wellSwapper(randArraysComp,targetSwap,sourceSwap,plate,plate) #Swap well info for source and target isolates
  
  targetSwap <- NULL #Reset target isolate
  sourceSwap <- NULL #Reset source isolate
}
sourceSwapGene <- maxEdgeCounts$Var1[2] #Repeat an iteration of the loop except with maxEdgeCounts[2] (maxEdgeCounts[1] already dealt with) - see notes with the first for loop.
#Need to recalculate everything except maxEdgeCounts
randEdges <- randArraysComp %>% filter((Row==1|Row==8|Column==1|Column==12))
randCenters <- randArraysComp %>% filter(!(Row==1|Row==8|Column==1|Column==12))
if (sum(table(randEdges$Group)<=(nIso-4))>1) {
  lowEdgeCounts <- as.data.frame(table(randEdges$Group)[table(randEdges$Group)<=(nIso-4)])
} else if (sum(table(randEdges$Group)<=(nIso-4))==1) {
  lowEdgeCounts <- as.data.frame(table(randEdges$Group)[table(randEdges$Group)<=(nIso-3)])
}
midEdgeCountsBackup <- as.data.frame(table(randEdges$Group)[table(randEdges$Group)==nIso-3])
edgeCands <- randEdges[randEdges$Group==sourceSwapGene,]
sourceSwap <- slice_sample(edgeCands,n=1)
plate <- sourceSwap$Plate_ID
n<-0
x<-1
while (n==0 & x<=1000) {
  targetSwapGeneCand <- case_when(
    x>=50 ~ sample(midEdgeCountsBackup$Var1,1),
    .default = sample(lowEdgeCounts$Var1,1)
  )
  if (x==50) {print("Exceeded 50 iterations without finding an eligible isolate, using mid edge counts")}
  centerCands <- duplicateCenters[duplicateCenters$Group==targetSwapGeneCand,]
  centerCands <- centerCands[centerCands$Plate_ID==plate,]
  n <- case_when(nrow(centerCands)>=1 ~ 1, .default = 0)
  x <- x+1
  if (x==1000) {print("Exceeded 1000 iterations without finding an eligible isolate, even from backup group. Try restarting?")}
}
targetSwap <- slice_sample(centerCands,n=1)
randArraysComp <- wellSwapper(randArraysComp,targetSwap,sourceSwap,plate,plate)



#Wild Type Metadata Renaming
randArraysComp<-randArraysComp %>% mutate(geneIso = str_c(Group,"_",Isolate))

range <- 1:nrow(randArraysComp)
randArraysComp <- randArraysComp %>% mutate(Source_Well = sapply(range, function(x) {
  randArraysComp[x,]$splitSample[[1]][1]
}))

randArraysComp[randArraysComp$geneIso=="wt_1",]$Source_Well <- "A2"
randArraysComp[randArraysComp$geneIso=="wt_1",]$Source_Plate <- "S3"
randArraysComp[randArraysComp$geneIso=="wt_2",]$Source_Well <- "B2"
randArraysComp[randArraysComp$geneIso=="wt_2",]$Source_Plate <- "S3"
randArraysComp[randArraysComp$geneIso=="wt_3",]$Source_Well <- "C2"
randArraysComp[randArraysComp$geneIso=="wt_3",]$Source_Plate <- "S3"

randArraysComp[randArraysComp$geneIso=="wt_4",]$Source_Well <- "A12"
randArraysComp[randArraysComp$geneIso=="wt_4",]$Source_Plate <- "S4"
randArraysComp[randArraysComp$geneIso=="wt_5",]$Source_Well <- "B12"
randArraysComp[randArraysComp$geneIso=="wt_5",]$Source_Plate <- "S4"
randArraysComp[randArraysComp$geneIso=="wt_6",]$Source_Well <- "C12"
randArraysComp[randArraysComp$geneIso=="wt_6",]$Source_Plate <- "S4"

randArraysComp[randArraysComp$geneIso=="wt_7",]$Source_Well <- "A7"
randArraysComp[randArraysComp$geneIso=="wt_7",]$Source_Plate <- "S6"
randArraysComp[randArraysComp$geneIso=="wt_8",]$Source_Well <- "B7"
randArraysComp[randArraysComp$geneIso=="wt_8",]$Source_Plate <- "S6"
randArraysComp[randArraysComp$geneIso=="wt_9",]$Source_Well <- "C7"
randArraysComp[randArraysComp$geneIso=="wt_9",]$Source_Plate <- "S6"
randArraysComp[randArraysComp$geneIso=="wt_10",]$Source_Well <- "D7"
randArraysComp[randArraysComp$geneIso=="wt_10",]$Source_Plate <- "S6"

randArraysComp[randArraysComp$geneIso=="wt_11",]$Source_Well <- "A9"
randArraysComp[randArraysComp$geneIso=="wt_11",]$Source_Plate <- "S8"
randArraysComp[randArraysComp$geneIso=="wt_12",]$Source_Well <- "B9"
randArraysComp[randArraysComp$geneIso=="wt_12",]$Source_Plate <- "S8"
randArraysComp[randArraysComp$geneIso=="wt_13",]$Source_Well <- "C9"
randArraysComp[randArraysComp$geneIso=="wt_13",]$Source_Plate <- "S8"

randArraysComp[randArraysComp$geneIso=="wt_14",]$Source_Well <- "A6"
randArraysComp[randArraysComp$geneIso=="wt_14",]$Source_Plate <- "S0"
randArraysComp[randArraysComp$geneIso=="wt_15",]$Source_Well <- "B6"
randArraysComp[randArraysComp$geneIso=="wt_15",]$Source_Plate <- "S0"
randArraysComp[randArraysComp$geneIso=="wt_16",]$Source_Well <- "C6"
randArraysComp[randArraysComp$geneIso=="wt_16",]$Source_Plate <- "S0"



#PIXL ReArray Prepper - creates PIXL rearray files from the randomized well positions
randArraysComp <- randArraysComp %>% rename(Target_Plate = "Plate_ID")

randArraysComp <- randArraysComp %>% mutate(Target_Row = substr(Well,1,1))
randArraysComp <- randArraysComp %>% mutate(Target_Col = substr(Well,2,3))
randArraysComp <- randArraysComp %>% mutate(Source_Row = substr(Source_Well,1,1))
randArraysComp <- randArraysComp %>% mutate(Source_Col = substr(Source_Well,2,3))

randArraysComp <- randArraysComp %>% mutate(Plate_Num = substr(Target_Plate,1,1))

targetFolder <-"/Users/samuelamidon/Desktop/Aim 1/Randomization/Randomization 3 From DWP/PIXL_Instruction_Files"
plateNums <- unique(randArraysComp$Plate_Num)
for (i in 1:length(plateNums)) {
  dat <- randArraysComp[randArraysComp$Plate_Num==plateNums[i],]
  dat <- dat[order(dat$Source_Plate),]
  datPrint <- dat %>% select(c(Source_Plate, Source_Row, Source_Col,Target_Plate,Target_Row,Target_Col))
  
  sourcePlates<-unique(dat$Source_Plate)
  targetPlates<-unique(dat$Target_Plate)
  header <- data.frame("Source_Plate"=sourcePlates[1],"Source_Row"="SBS","Source_Col"="96","Target_Plate"="Source","Target_Row"=NA,"Target_Col"=NA)
  for (j in 2:length(sourcePlates)) {
    newHead <- data.frame("Source_Plate"=sourcePlates[j],"Source_Row"="SBS","Source_Col"="96","Target_Plate"="Source","Target_Row"=NA,"Target_Col"=NA)
    header <- bind_rows(header,newHead)
  }
  n<-1
  for (j in (length(sourcePlates)+1):(length(sourcePlates)+length(targetPlates))) {
    newHead <- data.frame("Source_Plate"=targetPlates[n],"Source_Row"="SBS","Source_Col"="96","Target_Plate"="Target","Target_Row"=NA,"Target_Col"=NA)
    header <- bind_rows(header,newHead)
    n<-n+1
  }
  
  targetFile <- str_c(targetFolder,"/Plate",as.character(plateNums[i]),".csv")
  
  datWrite <- as.matrix(bind_rows(header,datPrint))
  colnames(datWrite)<-NULL
  write.table(datWrite, file=targetFile, quote=FALSE, sep=",", na="", col.names=FALSE, row.names=FALSE)
}



#Plate Layout Plotter - generates figures with randomized playe layouts
plates <- unique(randArraysComp$Target_Plate)
plots<-list(plate_plot(
  data = data_continuous_96,
  position = well,
  value = Value,
  plate_size = 96,
  plate_type = "round"
))
for (i in 1:length(plates)) {
  dat<-randArraysComp[randArraysComp$Target_Plate==plates[i],]
  plots[[i+1]]<-plate_plot(
    data=dat,
    position=Well,
    # value=Group,
    value=Source_Plate,
    # label=Source_Plate,
    # label=Sample,
    label=geneIso,
    plate_size=96,
    title=plates[i],
    plate_type="square",
    label_size=1.2,
    scale=1.5,
    legend_n_row=27
  )
}
pdf("/Users/samuelamidon/Desktop/Rand_Plate_Check.pdf")
for (i in 1:(length(plates)+1)) {
  print(plots[[i]])
}
dev.off()



#96->384Mapper - converts layouts in 96 density to 384 density layouts
plateNames <- unique(randArraysComp$Target_Plate)
rand384Map <- data.frame(Well="A1",Type="S",Strain="YNL001C",Isolate="1",Quartet="1",Plate_ID="1_A")
quartetVec <- c(1,2,3,4)

count384 <- 1
for (i in 1:nrow(randArraysComp)) {
  origRow <- which(randArraysComp[i,]$Target_Row==LETTERS)
  origCol <- randArraysComp[i,]$Target_Col
  rand384Map[count384,]$Well <- str_c(as.character(LETTERS[2*origRow-1]),as.character(2*origCol-1))
  rand384Map[count384+1,]$Well <- str_c(as.character(LETTERS[2*origRow-1]),as.character(2*origCol))
  rand384Map[count384+2,]$Well <- str_c(as.character(LETTERS[2*origRow]),as.character(2*origCol-1))
  rand384Map[count384+3,]$Well <- str_c(as.character(LETTERS[2*origRow]),as.character(2*origCol))
  rand384Map[count384:(count384+3),]$Type <- strsplit(randArraysComp[i,]$Source_Plate,split="")[[1]][1]
  rand384Map[count384:(count384+3),]$Strain <-strsplit(as.character(randArraysComp[i,]$Strain),split="_")[[1]][1]
  rand384Map[count384:(count384+3),]$Isolate <- strsplit(as.character(randArraysComp[i,]$Strain),split="_")[[1]][2]
  rand384Map[count384:(count384+3),]$Plate_ID <- randArraysComp[i,]$Target_Plate
  rand384Map[count384:(count384+3),]$Quartet <- quartetVec
  count384 <- count384+4
}



#384 Plate Layout Plotter - makes plate layout figures in 384 density
plates <- unique(rand384Map$Plate_ID)
plots<-list(plate_plot(
  data = data_continuous_96,
  position = well,
  value = Value,
  plate_size = 384,
  plate_type = "round"
))
for (i in 1:length(plates)) {
  dat<-rand384Map[rand384Map$Plate_ID==plates[i],]
  plots[[i+1]]<-plate_plot(
    data=dat,
    position=Well,
    # value=Group,
    value=Strain,
    # label=Source_Plate,
    # label=Sample,
    label=Quartet,
    plate_size=384,
    title=plates[i],
    plate_type="square"
  )
}
pdf("/Users/samuelamidon/Desktop/384Map.pdf")
for (i in 1:(length(plates)+1)) {
  print(plots[[i]])
}
dev.off()



#Saves plate maps separated into sheets by plate name
rand384Map$PP <- 1
rand384Map <- rand384Map %>% mutate(Row = substr(Well,1,1))
rand384Map <- rand384Map %>% mutate(Col = as.numeric(substr(Well,2,3)))
for (i in 1:nrow(rand384Map)) {
  rand384Map$PP[i] <- str_c("P",strsplit(rand384Map$Plate_ID[i],"_")[[1]][1],strsplit(rand384Map$Plate_ID[i],"_")[[1]][2],"")
}

wb <- createWorkbook("Plate_Map_384_SplitByPlate.xlsx")
plates <- unique(rand384Map$PP)
for (i in 1:length(plates)) {
  sheetName <- str_c(plates[i],"_map")
  dat <- rand384Map[rand384Map$PP == plates[i],]
  dat <- dat[order(dat$Col,dat$Row),]
  dat <- dat[,1:5]
  addWorksheet(wb, sheetName)
  writeData(wb, sheet = sheetName, dat)
}
saveWorkbook(wb, "/Users/samuelamidon/Desktop/Aim 1/Plate_Map_384_SplitByPlate.xlsx",overwrite=TRUE)
