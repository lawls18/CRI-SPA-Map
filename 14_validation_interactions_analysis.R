library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)

#The purpose of this code is to perform GxE, GxG, and GxGxE analyses on the MKT1-SAL1 validation data.
#It takes in the processed validation data and performs plate and edge corrections
#It then performs the three GxX analyses separately and saves the results for each.

# FILENAME WITH ALL DATA
liqData_inFile <- 'tables/validation_phenotype_data.tsv'
fileName <- "tables/processed_pixl_data.xlsx"
condNames <- excel_sheets(fileName) #Get conditions from sheet names

# FILENAME TO SAVE GXE MODEL DATA
dataFile <- 'tables/validation_phenotype_allConditions.tsv'
modelFile <- "tables/interaction_models.tsv"

for (i in 1:length(condNames)) { #iterate through conditions
   dat <- read_excel(fileName, sheet = condNames[i]) #Get current condition's data
   dat$Enviro <- condNames[i] #Create environment variable
   
   #Ensure consistent classes for bind_rows
   dat$Isolate <- as.character(dat$Isolate) 
   dat$Replicate <- as.numeric(dat$Replicate)
   dat$plate <- as.numeric(dat$plate)
   dat <- dat %>% rename(colSize = `Colony Area (mm^2)`)
   
   dat$Plate_Dup <- str_c(dat$Plate_ID,"_",dat$plate) #Plate dup is unique ID for each plate
   if (condNames[i] == "YPD37C") { #YPD37C has many "running" colonies. Filter whole plates because it is hard to tell when the running stops/starts, some plates also had all colonies with weird shapes and were removed
      dat <- dat %>% filter(!(Plate_Dup %in% c("P1A_1", "P1A_2", "P2A_2", "P2B_1", "P3A_1", "P3A_2", "P3B_1", "P4A_1", "P4A_2", "P4B_1", "P5B_2", "P8B_2")))
   } else if (condNames[i] == "YPD") { #Two YPD40h plates had the same issue as YPD37C
      dat <- dat %>% filter(!(Plate_Dup %in% c("P3A_1", "P3A_2")))
   }
   
   if (i==1) { #Make new datLong if first condition
      datLong <- dat
   } else { #If not, bind to existing datLong
      datLong <- bind_rows(datLong,dat)
   }
}

datLong <- datLong %>% mutate(IsolateNum = Isolate,  #Rename
                              Isolate = str_c(gene, "_", Isolate), #Add column with combined gene and isolate
                              IsolateType = str_c(gene, "_", Type,IsolateNum),
                              plateName=str_c(Plate_ID,"_",plate)) 
datLong$colSize <- log2(datLong$colSize)
datLong$Column <- str_c("Col",as.character(datLong$Column)) #Make sure col and row are characters for model
datLong$Row <- str_c("Row",as.character(datLong$Row))

datLong$colSizeCorrected <- NA

for (i in 1:length(condNames)) { #Iterate through conditions
   cond <- condNames[i]
   dat <- datLong[datLong$Enviro == cond,] #Get only data for current condition
   
      plateModel <- lmer(colSize ~ (1|gene) + (1|Isolate) + (1|plateName) + (1|Row) + (1|Column), data=dat) #Define plate model
      plateEffects <- data.frame(ranef(plateModel)$plateName) %>% mutate(plateName = rownames(.)) #Get effect for plate duplicate 1 and 2
      rowEffects <- data.frame(ranef(plateModel)$Row) %>% mutate(Row = rownames(.))
      colEffects <- data.frame(ranef(plateModel)$Column) %>% mutate(Column = rownames(.))
      dat <- left_join(dat, plateEffects, by="plateName") #Add plate effect column
      dat$plateEffect <- dat$X.Intercept. #Make plate effect variable
      dat$X.Intercept. <- NULL #Delete X.int variable
      dat <- left_join(dat, rowEffects, by="Row")
      dat$rowEffect <- dat$X.Intercept.
      dat$X.Intercept. <- NULL
      dat <- left_join(dat, colEffects, by="Column")
      dat$colEffect <- dat$X.Intercept.
      dat$X.Intercept. <- NULL
      
      dat$colSizeCorrected <- dat$colSize - dat$plateEffect - dat$rowEffect - dat$colEffect #edge effects for plotting
   datLong$colSizeCorrected[datLong$Enviro == cond] <- dat$colSizeCorrected
}

datLongVal<-datLong[(datLong$Type %in% c("V")),] #Get validation data
datLongVal[datLongVal$gene=="Eng-WT",]$gene <- "wt" #Rename wild type
datLongVal %<>% mutate(IsolateType = paste0(gene, '_', Type, IsolateNum))
datLongVal %<>% select(gene, IsolateType, colSizeCorrected, Enviro) ## restrict to only columns needed in model

liqData <- read_tsv(liqData_inFile, col_names = TRUE, show_col_types = FALSE) %>% 
   mutate(Enviro = 'YPDliq',
          gene = case_when(gene == 'FRG8' ~ 'wt',
                           gene == 'FRG1' ~ 'MKT1-30D', 
                           gene %in% c('FRG2', 'FRG7') ~ 'SAL1-fs',
                           gene == 'FRG3' ~ 'MKT1-30D+SAL1-fs', .default = gene)) %>%
   rename(IsolateType = 'isoType', colSizeCorrected = 'log2growthRatePlateCorrected') %>%
   select(gene, IsolateType, colSizeCorrected, Enviro)

datLongVal <- bind_rows(datLongVal, liqData) %>% rename(Isolate = 'IsolateType')

datLongVal <- datLongVal %>% filter(gene != "RM" & gene != "mkt1KanMXsal1HphMX" & gene != "mkt1::KanMX" & gene != "sal1::KanMX" & gene != "BY4742" & gene != "FRG4" & gene != "FRG5" & gene != "FRG6") #Filter out unused strains
datLongVal$gene <- factor(datLongVal$gene, levels = c("wt","W303","MKT1-30D+SAL1-fs","MKT1-30D","SAL1-fs")) #Order strains

condNames <- c(condNames, 'YPDliq')

# ORF model function
orf_lm <- function(data, ref) {
   lmResults <- bind_rows(lapply(sort(as.character(unique(data$gene[data$gene != ref]))), function(thisGene) { # loop through all genes except WT and create a dataframe of model information
      testDat <- data %>% filter(gene %in% c(thisGene, ref)) # extract all data for gene of interest and WT isolates
      fullModel <- lmer(colSizeCorrected ~ gene + (1|Isolate), data = testDat, REML = FALSE) # create model from corrected colony size and isolate
      res <- c(thisGene, fixef(fullModel)[1], fixef(fullModel)[2], anova(fullModel)[1,'Pr(>F)']) # extract effects from the model
      names(res) <- c('gene', 'geneBaseline', 'geneEffect', 'pValueGene') # name vector containing different model effects
      return(res) # output the effects for each gene
   }))
   lmResults %<>% add_row(gene = ref, pValueGene = '1') # add WT row to results
   lmResults %<>% mutate(geneBaseline = as.numeric(geneBaseline), # convert all model outputs to numerics
                         geneEffect = as.numeric(geneEffect),
                         pValueGene = as.numeric(pValueGene),
                         significance = case_when( # report significance of gene effects based on model p-values
                            (pValueGene < 0.05) ~ 'Nominal',
                            .default = 'NS'))
   return(lmResults)
}

lmGene <- lapply(condNames, function(cond) {
   dat = datLongVal %>% filter(Enviro == cond)
   lm = orf_lm(dat, ref = 'wt')
   return(lm %>% mutate(condition = cond))
}) %>% bind_rows() %>%
   mutate(gene_cond = paste(gene, condition, sep = '_'))

controlCond <- "YPD" #GxE and GxGxE will compare to this condition

# Calculating GxE for all genes except wild type, comparing each condition to control condition
lmGxE <- lapply(condNames[condNames != controlCond], function(cond) { #Iterate through conditions
      dat <- datLongVal[datLongVal$Enviro == cond | datLongVal$Enviro == controlCond,] #Get data for current condition and control condition
      
      dat$Enviro <- factor(dat$Enviro, levels = c(controlCond,cond)) #Order the condition names
      
      lmGxEResults <- data.frame(t(sapply(sort(as.character(unique(dat$gene[dat$gene != "wt"]))), function(thisGene){ #GxE model
         testDat <- dat[dat$gene %in% c(thisGene, "wt"),]
         fullModel <- lmer(colSizeCorrected ~ gene*Enviro + (1|Isolate), data=testDat, REML=FALSE) #Only change is gene*Enviro
         var <- anova(fullModel)
         res <- c(fixef(fullModel)[1], fixef(fullModel)[2], var[1,"Pr(>F)"], var[2,"Pr(>F)"], var[3,"Pr(>F)"]) #New values named below
         names(res) <- c("baseline", "geneEffect", "GpValue","EpValue","GxEpValue") #Name result variables
         res
      })))
      
      lmGxEResults <- rbind(lmGxEResults, wt = c(NA, NA, 1,1,1)) #Add wt row to results
      lmGxEResults$gene <- rownames(lmGxEResults) #Make gene name variable
      rownames(lmGxEResults) <- NULL
      lmGxEResults$condition <- cond

      return(lmGxEResults) #Save GxE data
   }) %>% bind_rows() %>% 
   select(condition, gene, baseline, geneEffect, GpValue, EpValue, GxEpValue) %>%
   mutate(gene_cond = paste(gene, condition, sep = '_'))

#Calculate GxG on each condition
lmGxG <- lapply(condNames, function(cond) {
   dat <- datLongVal[datLongVal$Enviro == cond,] #Only look at current condition
   dat <- dat[dat$gene != "W303",] #Exclude W303 because it includes more variants
   
   dat$mType <- "BY" #Default is BY
   dat$mType[dat$gene == "MKT1-30D" | dat$gene == "MKT1-30D+SAL1-fs"] <- "W303" #These strains have W303 version of MKT1
   
   dat$sType <- "BY" #Default is BY
   dat$sType[dat$gene == "SAL1-fs" | dat$gene == "MKT1-30D+SAL1-fs"] <- "W303" #These strains have W303 version of SAL1
   
   fullModel <- lmer(colSizeCorrected ~ sType*mType + (1|Isolate), data=dat, REML=FALSE) #GxG model
   
   var <- anova(fullModel)
   lmResults <- c(cond, fixef(fullModel)[1], fixef(fullModel)[2], fixef(fullModel)[3], fixef(fullModel)[4],
                  var[1,"Pr(>F)"], var[2,"Pr(>F)"], var[3,"Pr(>F)"]) #New variables for this comparison
   names(lmResults) <- c("condition", "baseline", "sEffect", "mEffect", "sxmEffect",
                         "SpValue", "MpValue", "SxMpValue") #New variable names
   
   return(lmResults) #Save GxG data
}) %>% bind_rows()

# Calculating GxGxE for all conditions except control. Each condition is compared to control
lmGxGxE <- lapply(condNames[condNames != controlCond], function(cond) { #List to save GxGxE data
      dat <- datLongVal[datLongVal$Enviro == cond | datLongVal$Enviro == controlCond,] #Look at both current condition and control condition
      dat <- dat[dat$gene != "W303",] #Once again, exclude W303 because of extra variants
      dat$Enviro <- factor(dat$Enviro, levels = c(controlCond,cond)) #Factorize environments
      
      dat$mType <- "BY" #BY default
      dat$mType[dat$gene == "MKT1-30D" | dat$gene == "MKT1-30D+SAL1-fs"] <- "W303" #Update correct W303 for MKT1
      
      dat$sType <- "BY" #BY default
      dat$sType[dat$gene == "SAL1-fs" | dat$gene == "MKT1-30D+SAL1-fs"] <- "W303" #Update corrent W303 for SAL1
      
      fullModel <- lmer(colSizeCorrected ~ sType*mType*Enviro + (1|Isolate), data=dat, REML=FALSE) #GxGxE model
      
      var <- anova(fullModel)
      lmResults <- c(cond, fixef(fullModel)[1], fixef(fullModel)[2], fixef(fullModel)[3], fixef(fullModel)[4],
                     var[1,"Pr(>F)"], var[2,"Pr(>F)"], var[3,"Pr(>F)"], var[4,"Pr(>F)"], var[5,"Pr(>F)"], var[6,"Pr(>F)"], var[7,"Pr(>F)"]) #Even more new variables
      names(lmResults) <- c("condition", "baseline", "sEffect", "mEffect", "eEffect", #More new variable names
                            "SpValue", "MpValue", "EpValue", "SxMpValue", "SxEpValue", "MxEpValue", "SxMxEpValue")
      
      return(lmResults) #Save GxGxE data
   }) %>% bind_rows()

allModels <- full_join(lmGxGxE %>% select(condition, SxMxEpValue), left_join(lmGxG %>% select(condition, SxMpValue), 
                       left_join(lmGene %>% select(gene_cond, gene, condition, pValueGene), 
                       lmGxE %>% select(gene_cond, GxEpValue), by = 'gene_cond'), by = 'condition'), by = 'condition') %>%
   rename(pValueGxE = 'GxEpValue', pValueGxG = 'SxMpValue', pValueGxGxE = 'SxMxEpValue') %>%
   select(condition, gene, pValueGene, pValueGxE, pValueGxG, pValueGxGxE)

write_tsv(datLongVal, file = dataFile, col_names = TRUE) #Save all model data
write_tsv(allModels, file = modelFile, col_names = TRUE) #Save all model data
