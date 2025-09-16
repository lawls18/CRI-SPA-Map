library(dplyr)
library(openxlsx)
library(tidyverse)
library(readxl)

#The purpose of this code is to take folders with raw PIXL data and match up colonies with the correct strain information.
#Processed data is saved to excel file with sheets corresponding to conditions.

#ENTER FOLDER NAME CONTAINING RAW DATA
dataFolder <- "data/pixl_data"    

#ENTER CONDITION NAMES (Could also be just YPD)
conditions <- c('Caf', 'EtOH', 'LiAc', 'LowG', 'NaCl', 'SDS', 'Tun', 'YNB', 'YPD', 'YPD37C', 'YPG') #THESE NAMES MUST MATCH THE CONDITION IN THE PIXL FILE NAMES       

#ENTER FILE PATH FOR PLATE MAPS
mapName <- "data/CRI-SPA-Map_ValidationNo1_Map.xlsx"

#NEED THIS FILE PRESENT FOR CODE TO RUN
coordMap <- read_excel("data/coordMap384.xlsx")  

#ENTER FILE WHERE YOU WANT DATA SENT
targetDataFilePath <- "tables/processed_pixl_data.xlsx" 

fileNames <- dir(dataFolder) #Get all raw data folder names

sortedData <- list() #Initialize sorted data storage list
exDF <- read.table(str_c(dataFolder,"/",fileNames[1]),header=TRUE,sep="\t") #Get example data to read column names
for (condNum in 1:length(conditions)) { #Iterate through conditions
   sortedData[[condNum]] <- setNames(data.frame(matrix(ncol = (ncol(exDF)+10), nrow = 1)), c(colnames(exDF),"Row","Column","Well","Type","gene","Isolate","Replicate","Plate_ID","plate","Array")) #Initialize data frames to store sorted data in with correct column names and a row of NAs
   sortedData[[condNum]] <- sortedData[[condNum]] %>% rename(`Colony Area (mm^2)` = "Colony_Area_.mm2.") #Rename colony area column
}

for (fileNum in 1:length(fileNames)) { #Iterate through all data files
   dat <- read.table(str_c(dataFolder,"/",fileNames[fileNum]),header=TRUE,sep="\t") #Read data file
   curPlate <- str_split_i(fileNames[fileNum], '_AllCol', 1) #Extract plate identifying information
   
   #This section rotates coordinates for Lithium acetate plate 6B2, which was rotated during creation
   if (curPlate == "P6_B_2_LiAc") {
      x <- dat$Colony_X_.mm.
      y <- dat$Colony_Y_.mm.
      dat$Colony_X_.mm. <- x*cos(-0.02)-y*sin(-0.02) #Rotates every point by 0.02 radians
      dat$Colony_Y_.mm. <- y*cos(-0.02)+x*sin(-0.02)
   }
   
   #This section matches coordinates that should be the same but are not set equal to eachother by rounding
   dat$Colony_X_.mm. <- ceiling(dat$Colony_X_.mm.) #Round all x coords
   dat <- dat[order(dat$Colony_X_.mm.),] #Sort by x coords
   for (coordNum in 2:nrow(dat)) { #Iterate through colonies, starting at the second
      if ((dat$Colony_X_.mm.[coordNum] == (dat$Colony_X_.mm.[coordNum-1]+1)) | (dat$Colony_X_.mm.[coordNum] == (dat$Colony_X_.mm.[coordNum-1]-1))) #Check if n-1th x coord is close to nth coord (+/- 1 coord value)
         dat$Colony_X_.mm.[coordNum] <- dat$Colony_X_.mm.[coordNum-1] #If n-1 and nth coords are close, set them equal to eachother
   }
   #Repeat for y coords
   dat$Colony_Y_.mm. <- ceiling(dat$Colony_Y_.mm.) 
   dat <- dat[order(dat$Colony_Y_.mm.),]
   for (coordNum in 2:nrow(dat)) {
      if ((dat$Colony_Y_.mm.[coordNum] == (dat$Colony_Y_.mm.[coordNum-1]+1)) | (dat$Colony_Y_.mm.[coordNum] == (dat$Colony_Y_.mm.[coordNum-1]-1)))
         dat$Colony_Y_.mm.[coordNum] <- dat$Colony_Y_.mm.[coordNum-1]
   }
   
   dat <- dat[order(dat$Colony_Y_.mm.,dat$Colony_X_.mm.),] #Leveled sort by row then column to match coordinate map order
   
   #This section matches colony coords to coord map coords and extracts row letters and column numbers
   for (i in 1:nrow(dat)) { #Iterate through colonies
      x <- dat[i,] #Get ith colony data
      coords <- coordMap[ #Want to find the coordinate map info that matches the coordinates of colony i. By the end, should have a single row from coordMap where the x and y coords are within +/- 2 units of the coordinates of colony i
         ((coordMap$Row_Coord-2) <= as.numeric(x[7]) & #Find indices where y coord (row) is greater than coordMap coords -2
             (coordMap$Row_Coord+2) >= as.numeric(x[7])) & #Find indices where y coord is less than coordMap coords +2
            ((coordMap$Col_Coord-2) <= as.numeric(x[6]) & #Find indices where x coord (col) is greater than coordMap coords -2
                (coordMap$Col_Coord+2) >= as.numeric(x[6])),] #Find indices where x coord (col) is less than coordMap coords +2
      if (nrow(coords) == 1) { #Check that a match was found (if no match, coords has no rows)
         dat$Row[i] <- coords$Row #Set row letter of matched coordMap info to data
         dat$Col[i] <- coords$Col #Set column number of match to data
      }
      else { #If no match found, set to NA
         dat$Row[i] <- NA
         dat$Col[i] <- NA
      }
   }
   
   #This section extracts all colony information (isolate, quartet, gene name) from the plate map and inserts it into the colony data
   dat$Well <- str_c(dat$Row,dat$Col) #Set well column (A2, etc.)
   splitName <- strsplit(curPlate,"_")[[1]] #Split file name to extract different plate info
   plateName <- str_c(splitName[1],splitName[2]) #Get plate name (P1A,P2A, etc)
   plateMap <- read_excel(mapName,sheet = str_c(plateName,"_map")) #Read plate map, sheet depends on current plate name
   dat <- left_join(dat, plateMap, by = join_by(Well==Well), unmatched="drop") #Join all colony info by Well. A colony's info is only added if its Well matches a Well in plate map. This always happens unless its well is NA.
   
   #This section adds additional columns to colony data in order for the output files to match previous files made before this code existed. THIS SECTION COULD BE MADE PRETTIER WITH ALTERATIONS TO FILES THAT READ THE DATA OUTPUT BY THIS FILE
   dat$Plate_ID <- plateName #Plate ID is equivalent to plate name P1A
   dat$plate <- splitName[3] #plate is equivalent to plate duplicate number (1 or 2) - plates P1A 1 and P1A 2 are duplicates of eachother
   dat$Array <- splitName[2] #Array is equivalent to array letter (A or B) - P1A and P1B contain the same strains but are randomized differently
   dat <- dat %>% rename(Column = "Col") #Rename col to column
   dat <- dat %>% rename(gene = "Strain") #Rename strain to Gene
   dat <- dat %>% rename(Replicate = "Quartet") #Rename quartet to replicate (this is a colonies position within its quartet -top left:1, top right:2, bot left:3, bot right:4)
   dat <- dat %>% rename(`Colony Area (mm^2)` = "Colony_Area_.mm2.") #Rename colony area column
   
   #This section adds each sorted file to the compilation files
   cond <- strsplit(curPlate,"_")[[1]][4] #Get condition name from file name
   condInd <- which(cond==conditions) #Find condition index
   
   dat$SourcePlate_Name <- as.character(dat$SourcePlate_Name) #ADDED BECAUSE ONE INSTANCE SOURCEPLATENAME WAS DOUBLE NOT CHAR
   
   sortedData[[condInd]] <- bind_rows(sortedData[[condInd]],dat) #Bind the current colony data to the compilation depending on condition name
}

#This section writes the compiled data to a spreadsheet
wb <- createWorkbook(targetDataFilePath) #Initialize workbook
for (condNum in 1:length(conditions)) { #Iterate through conditions
   sortedData[[condNum]] <- sortedData[[condNum]][2:nrow(sortedData[[condNum]]),] #Eliminate the first row of data which is NAs from file initialization
   
   addWorksheet(wb, sheetName = conditions[condNum]) #Add worksheet by condition name
   writeData(wb, sheet = conditions[condNum], sortedData[[condNum]]) #Write to correct worksheet
}

saveWorkbook(wb, file = targetDataFilePath) #Save workbook