# R version 4.4.3

# load libraries
library(tidyverse) ## v 2.0.0
library(magrittr) ## v 2.0.3
library(readxl) ## v 1.4.3
library(growthcurver) ## v 0.3.1
library(lme4) ## v 1.1-37
library(lmerTest) ## v 3.1-3
library(patchwork) ## v 1.3.0
library(ggforce) ## v 0.4.2
library(extrafont) ## v 0.19

# user-input file paths and info
phenData_inDir <- 'data/phen_data/'
pheno_data_desc_inFile <- 'data/pheno_data_desc.xlsx'
chrXIVLorf_inFile <- 'data/S288C_chrXIVL_ORFs.tsv' ## file containing info on all ORFs found on chrXIV-L in S288C reference

outDir <- '250820_tables/' ## directory for all analysis output files
figDir <- '250820_figures/' ## directory for all figure output files
pheno_model_data_YKOremoved_outFile <- paste0(outDir, 'pheno_model_data_YKOremoved.tsv')
pval_effects_outFile <- paste0(outDir, 'orf_pvalue_effectSizes.tsv')
snp_effects_outFile <- paste0(outDir, 'snp_pVals_effectSizes.tsv')
pheno_plots_YKOremoved_outFile <- paste0(figDir, 'pheno_boxplots_YKOremoved.pdf')
pheno_plots_byRepairType_outFile <- paste0(figDir, 'pheno_boxplots_byRepairType.pdf')
val_plots_outFile <- paste0(figDir, 'val_boxplots.pdf')
MKT1_SAL1_valPlot_outFile <- paste0(figDir, 'MKT1_SAL1_YPDliq_valPlot.pdf')
hygTest_plot_outFile <- paste0(figDir, 'hphmx_test.pdf')
SNP_plot_outFile <- paste0(figDir, 'SNP_plots.pdf')
SNP_plot_zoom_outFile <- paste0(figDir, 'SNP_plots_zoom.pdf')

options(scipen = 999) ## no scientific notation
theme_set(theme_bw(base_size = 10, base_family = 'Arial') + ## set base ggplot fonts
             theme(text = element_text(color = 'black'), axis.text = element_text(color = 'black', size = 8)))

# make output directories
dir.create(outDir, showWarnings = FALSE)
dir.create(figDir, showWarnings = FALSE)

# load in previously-created files
simplified_mappingVar_calls_XIVL <- read_tsv(paste0(outDir, 'CSM_chrXIVLcalls.tsv'), col_names = TRUE, show_col_types = FALSE)
CSMorf_info <- read_tsv(paste0(outDir, 'CSM_orfInfo.tsv'), col_names = TRUE, show_col_types = FALSE)
mappingVars_XIVL <- read_tsv(paste0(outDir, 'mappingVars_XIVL.tsv'), col_names = TRUE, show_col_types = FALSE)
strainInfo <- read_tsv(paste0(outDir, 'strainInfo_wFilters.tsv'), col_names = TRUE, show_col_types = FALSE)
repairTracts_summary <- read_tsv(paste0(outDir, 'repairTracts_summary.tsv'), col_names = TRUE, show_col_types = FALSE)
repairTracts <- read_tsv(paste0(outDir, 'repairTracts.tsv'), col_names = TRUE, show_col_types = FALSE)

# load in all ORFs on chrXIV-L
chrXIVLorf_info <- read_tsv(chrXIVLorf_inFile, show_col_types = FALSE) %>% ## all S288C ORFs on chrXIVL from SGD/alliancemine
   rename(orf = SequenceFeature.secondaryIdentifier, ## rename columns
          gene_symbol = SequenceFeature.symbol, 
          gene_name = SequenceFeature.name, 
          sgd_id = SequenceFeature.primaryIdentifier, 
          start = SequenceFeature.chromosomeLocation.start, 
          stop = SequenceFeature.chromosomeLocation.end, 
          chr = SequenceFeature.chromosomeLocation.locatedOn.primaryIdentifier) %>%
   select(orf, gene_symbol, gene_name, chr, start, stop, sgd_id) ## reorder columns

##### PHENOTYPE DATA PROCESSING #####
# define function to process plate reader data
plateReaderProcessor <- function(file, file_desc){
   
   # load in plate description file
   desc.df <- read_excel(file_desc, col_names = TRUE)
   
   # find OD table in data - must manually input 'theend' in column A of final row of OD table
   OD_loc <- which(read_excel(file, col_names = FALSE, range = 'A1:A700') == '600')
   OD_end <- which(read_excel(file, col_names = FALSE, range = 'A1:A700') == 'theend')
   
   # load in OD table - assumes 96-well plate
   OD <- read_excel(file, col_names = TRUE, range = paste('B', OD_loc + 2,':CU', OD_end, sep=''))
   
   # remove extra data
   rm(OD_loc, OD_end)
   
   # store OD and GFP Time for Grofit
   OD_time <- OD$Time
   
   # define blanks - must contain either the word 'blank' or a 'Y' in Blank column
   
   if('Blank' %in% colnames(desc.df)) { 
      # only run if Blank column exists
      
      contains_blank <- grep('blank', desc.df$Blank, ignore.case = TRUE)
      contains_Y <- grep('Y', desc.df$Blank, ignore.case = TRUE)
      
      blank_num <- unique(c(contains_blank,contains_Y))
      blanks <- desc.df$Wells[blank_num]
      
      if(identical(blanks, character(0))) {warning('Blanks not defined - Use manual entry or change desc.df')}
      
      rm(contains_blank, contains_Y, blank_num)
      
   }  else {
      warning('Blank column not defined - please edit excel sheet')
   }
   
   # define wells 
   wells <- setdiff(desc.df$Wells, blanks) # does not include blanks in wells
   
   removewells <- desc.df$Wells[is.na(desc.df$Desc_1)]
   
   wells <- wells[!is.element(wells, removewells)]
   
   # make a data.frame for growthcurver
   OD.df <- data.frame(OD[,wells])
   
   # time in hours
   OD.df$time <- as.numeric((OD$Time - OD$Time[1])/60/60)
   
   # aggregate all blank columns into one and add them to the df
   OD.df$blank <- rowMeans(OD[,blanks])
   
   # sometimes the plate did not finish the whole 24h; need to remove these time points
   OD.df <- OD.df[sign(OD.df$time) >= 0,]
   
   gc_out <- SummarizeGrowthByPlate(OD.df, bg_correct = 'blank')
   rownames(gc_out) <- gc_out$sample
   gc_out <- gc_out[gc_out$sample != '',]
   gc_out[gc_out$note != '', 2:9] <- NA
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   # Baseline correct with blank data
   
   # Running Well = a mean of all blank wells for each time point
   OD_avg_blankwell <- apply(OD, 1, function(x){
      mean(as.numeric(x[blanks]))
   })
   
   # Grand = avg all time points for the running well
   grdCorr_OD <- mean(OD_avg_blankwell, na.rm = T)
   
   # Blank correct the data
   
   # grand corr
   OD <- OD - grdCorr_OD
   
   # filter tables to contain only relevant wells
   OD <- OD[,colnames(OD) %in% wells]
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   # Output data
   
   # Combined desc row
   catDesc <- apply(desc.df[3:7], 1, function(x){
      x <- x[!is.na(x)] # remove NAs
      x <- str_trim(x)
      cat <- paste(x, collapse = '_') # paste together desc columns
      return(cat)
   })
   
   # Store output data
   ratios.df <- data.frame(
      'Wells' = desc.df$Wells,
      catDesc,
      'growthRate' = gc_out[desc.df$Wells, 'r'],
      desc.df,
      row.names=NULL)
   
   # final return
   return(ratios.df)
}

# make table linking phenotyping data to appropriate description file
hygS_pheno_data_desc <- read_excel(path = pheno_data_desc_inFile, sheet = 'hygS', col_names = TRUE)
hygR_pheno_data_desc <- read_excel(path = pheno_data_desc_inFile, sheet = 'hygR', col_names = TRUE)
val_pheno_data_desc <- read_excel(path = pheno_data_desc_inFile, sheet = 'val', col_names = TRUE)

# process phenotyping data and merge into single table
hygS_pheno_data_processed <- apply(hygS_pheno_data_desc, 1, function(x) {
   data = plateReaderProcessor(file = paste0(phenData_inDir, x[1]), file_desc = paste0(phenData_inDir, x[2]))
   data %<>% mutate(plate = str_sub(str_split_i(x[1], '_', 3), -1, -1), 
                    plate_rep = str_split_i(x[1], '_', 4),
                    type = 'S')
   return(data)
}) %>% bind_rows() %>%
   rename(clone = 'catDesc', background = 'Desc_1', gene = 'Desc_3', isoNum = 'Desc_4') %>% ## rename columns
   filter(is.na(Blank) & clone != '') %>%
   mutate(isoType = str_c(gene, '_', type, isoNum), ## designate isolate with ORF, repair type, and isolate number
          log2growthRate = log2(growthRate), ## take log2 of growth rate
          plateName = str_c(plate, '_', plate_rep)) ## add column for plate name

hygR_pheno_data_processed <- apply(hygR_pheno_data_desc, 1, function(x) {
   data = plateReaderProcessor(file = paste0(phenData_inDir,x[1]), file_desc = paste0(phenData_inDir, x[2]))
   data %<>% mutate(plate = str_sub(str_split_i(x[1], '_', 3), -1, -1), 
                    plate_rep = str_split_i(x[1], '_', 4),
                    type = 'R')
   return(data)
}) %>% bind_rows() %>%
   rename(clone = 'catDesc', background = 'Desc_1', gene = 'Desc_3', isoNum = 'Desc_4') %>% ## rename columns
   filter(is.na(Blank) & clone != '') %>%
   mutate(gene = case_when(gene == 'YNL130W' ~ 'YNL130C', .default = gene), ## fix gene name typo
          isoType = str_c(gene, '_', type, isoNum), ## designate isolate with ORF, repair type, and isolate number
          log2growthRate = log2(growthRate), ## take log2 of growth rate
          plateName = str_c(plate, '_', plate_rep)) ## add column for plate name

val_pheno_data_processed <- apply(val_pheno_data_desc, 1, function(x) {
   data = plateReaderProcessor(file = paste0(phenData_inDir, x[1]), file_desc = paste0(phenData_inDir, x[2]))
   data %<>% mutate(plateName = str_split_i(x[1], '_', 4), 
                    type = 'V')
   return(data)
}) %>% bind_rows() %>%
   rename(clone = 'catDesc', background = 'Desc_1', gene = 'Desc_3', isoNum = 'Desc_4') %>% ## rename columns
   filter(is.na(Blank) & clone != '') %>%
   mutate(isoType = str_c(gene, '_', type, isoNum), ## designate isolate with ORF, repair type, and isolate number
          log2growthRate = log2(growthRate), ## take log2 of growth rate
          gene = case_when(gene == 'FRG7' ~ 'FRG2', gene == 'FRG8' ~ 'wt', .default = gene)) ## combine categories containing the same edits

# convert gene and clone columns to factors
hygS_pheno_data_processed %<>% mutate(gene = factor(gene, levels = c('wt', sort(unique(hygS_pheno_data_processed$gene[hygS_pheno_data_processed$gene != 'wt']), decreasing = TRUE))))
hygR_pheno_data_processed %<>% mutate(gene = factor(gene, levels = c('wt', 'BY', sort(unique(hygR_pheno_data_processed$gene[hygR_pheno_data_processed$gene != 'wt' & hygR_pheno_data_processed$gene != 'BY']), decreasing = TRUE))))
val_pheno_data_processed %<>% mutate(gene = factor(gene, levels = c('wt', sort(unique(val_pheno_data_processed$gene[val_pheno_data_processed$gene != 'wt']), decreasing = TRUE))))

##### LIQUID PHENOTYPING #####
# ORF model function
orf_lm <- function(data, ref) {
   lmResults <- bind_rows(lapply(sort(as.character(unique(data$gene[data$gene != ref]))), function(thisGene) { # loop through all genes except WT and create a dataframe of model information
      testDat <- data %>% filter(gene %in% c(thisGene, ref)) # extract all data for gene of interest and WT isolates
      fullModel <- lmer(log2growthRatePlateCorrected ~ gene + (1|isoType), data = testDat, REML = FALSE) # create model from corrected colony size and isolate
      res <- c(thisGene, fixef(fullModel)[1], fixef(fullModel)[2], anova(fullModel)[1,'Pr(>F)']) # extract effects from the model
      names(res) <- c('gene', 'geneBaseline', 'geneEffect', 'pValueGene') # name vector containing different model effects
      return(res) # output the effects for each gene
   }))
   lmResults %<>% add_row(gene = ref, pValueGene = '1') # add WT row to results
   lmResults %<>% mutate(geneBaseline = as.numeric(geneBaseline), # convert all model outputs to numerics
                         geneEffect = as.numeric(geneEffect),
                         pValueGene = as.numeric(pValueGene),
                         significance = case_when( # report significance of gene effects based on model p-values
                            (pValueGene < (0.05/(nrow(lmResults) - 1))) ~ 'Bonferroni',
                            (pValueGene < 0.05) ~ 'Nominal',
                            .default = 'NS'))
   return(lmResults)
}

# SNP model function
SNP_lm <- function(data, snps) {
   lmResults_SNP <- bind_rows(lapply(snps, function(thisSNP) {
      if(nrow(data %>% filter(get(thisSNP) == 'W303')) > 0) {
         testDat <- data %>% filter((get(thisSNP) == 'W303' | get(thisSNP) == 'BY') & !is.na(get(thisSNP))) %>% select(all_of(thisSNP), 'log2growthRatePlateCorrected', 'isoType')
         colnames(testDat) <- c('SNP', 'log2growthRatePlateCorrected', 'isoType')
         fullModel <- lmer(log2growthRatePlateCorrected ~ SNP + (1|isoType), data = testDat, REML = FALSE)
         res <- c(thisSNP, fixef(fullModel)[1], fixef(fullModel)[2], anova(fullModel)[1,'Pr(>F)'], nrow(data %>% filter(get(thisSNP) == 'W303')), nrow(data %>% filter(get(thisSNP) == 'BY'))) # lmer
         names(res) <- c('SNP', 'SNPbaseline', 'SNPEffect', 'pValueSNP', 'W303_calls', 'BY_calls')
         return(res)
      }
   }))
   lmResults_SNP %<>% mutate(SNPbaseline = as.numeric(SNPbaseline), # convert all model outputs to numerics
                             SNPEffect = as.numeric(SNPEffect),
                             pValueSNP = as.numeric(pValueSNP),
                             significance = case_when( # report significance of gene effects based on model p-values
                                (pValueSNP < (0.05/nrow(lmResults_SNP))) ~ 'Bonferroni',
                                (pValueSNP < 0.05) ~ 'Nominal',
                                .default = 'NS'),
                             W303_calls = as.numeric(W303_calls),
                             BY_calls = as.numeric(BY_calls))
   return(lmResults_SNP)
}

# function to run all ORF and SNP models
phenLinearModels <- function(data, strainList, calls, ref, hyg, orf_info, varList) {
   phenData = left_join(data %>% select(Wells, gene, type, isoType, plateName, growthRate, log2growthRate),
                        strainList, by = 'isoType')
   plateModel = lmer(log2growthRate ~ (1|gene) + (1|isoType) + (1|plateName), data = phenData) # create linear model for colony size
   plateEffects = data.frame(ranef(plateModel)$plateName) %>% mutate(plateName = rownames(.)) # create data frame with plate names and model coefficients
   
   phenData = left_join(phenData, plateEffects, by = 'plateName') # add plate effects column
   phenData %<>% rename(plateEffect = 'X.Intercept.') # rename plate effects column
   
   phenData %<>% mutate(log2growthRatePlateCorrected = log2growthRate - plateEffect) # add column with colony size corrected for plate effects
   
   if (hyg == 'R') {phenData_YKOremoved = phenData %>% filter((aux_removed == FALSE & preSeq_removed == FALSE & YKO_removed == FALSE) | gene %in% c('BY', 'wt'))} else {
      phenData_YKOremoved = phenData %>% filter((aux_removed == FALSE & preSeq_removed == FALSE & YKO_removed == FALSE) | gene == 'wt')}
   
   lmORF_YKOremoved = phenData_YKOremoved %>% orf_lm(data = ., ref = ref) %>% left_join(., orf_info, join_by(gene == orf))
   phenData_YKOremoved = inner_join(phenData_YKOremoved, lmORF_YKOremoved, by = 'gene')
   
   phenData_genos = inner_join(phenData, calls %>% select(!orf), join_by(corrected_isoName == strain))
   phenData_genos = phenData_genos %>%
      filter((aux_removed == FALSE & snp_removed == FALSE & sequenced_illumina == 'yes') | gene == 'wt')
   test_snps = colnames(calls)[which(str_detect(colnames(calls), 'var'))] 
   
   lmSNP = phenData_genos %>% SNP_lm(data = ., snps = test_snps)
   lmSNP = left_join(lmSNP, varList, join_by(SNP == varID))
   
   return(list(phenData_all = phenData, 
               phenData_YKOremoved = phenData_YKOremoved, 
               phenData_genos = phenData_genos, 
               lmORF_YKOremoved = lmORF_YKOremoved,
               lmSNP = lmSNP
   ))}

# transpose genotype call table function
transpose_calls <- function(simp_calls) {
   transpose_calls = as.data.frame(cbind(colnames(simp_calls)[2:ncol(simp_calls)], t(simp_calls[,-1]))) ## transpose calls and add column for strain names
   rownames(transpose_calls) <- NULL ## remove row names
   colnames(transpose_calls) <- c('strain', simp_calls$varID) ## name columns with variant IDs
   transpose_calls %<>% mutate(orf = str_split_i(strain, '_', 1)) ## indicate targeted ORF of each strain
   transpose_calls %<>% arrange(strain) ## sort table by strain name
   return(transpose_calls)
}

# transpose chrXIVL genotype calls
simplified_mappingVar_calls_XIVL[which(simplified_mappingVar_calls_XIVL$varID == 'var7523'), which(simplified_mappingVar_calls_XIVL[which(simplified_mappingVar_calls_XIVL$varID == 'var7523'),] == 'low_qual')] <- 'BY'
simplified_mappingVar_calls_XIVL[which(simplified_mappingVar_calls_XIVL$varID == 'var7523'), c('YNL083W_B8S1_4A7', 'YNL083W_B8S3_4B7', 'YNL083W_B8S4_4C7', 'YNL084C_G6S1_3G10', 'YNL084C_G6S3_3A11', 'YNL085W_C8S4_4B8')] <- 'W303'

transpose_mappingVar_calls_XIVL <- transpose_calls(simplified_mappingVar_calls_XIVL %>% select(varID, starts_with('YNL')))

hygS_models <- phenLinearModels(data = hygS_pheno_data_processed, 
                                strainList = strainInfo, 
                                calls = transpose_mappingVar_calls_XIVL,
                                ref = 'wt',
                                hyg = 'S',
                                orf_info = CSMorf_info,
                                varList = mappingVars_XIVL %>% select(varID, chr, pos, ref, alt))

hygR_models <- phenLinearModels(data = hygR_pheno_data_processed, 
                                strainList = strainInfo, 
                                calls = transpose_mappingVar_calls_XIVL,
                                ref = 'BY',
                                hyg = 'R',
                                orf_info = CSMorf_info,
                                varList = mappingVars_XIVL %>% select(varID, chr, pos, ref, alt))

valLinearModels <- function(data, ref) {
   phenData = data %>% select(Wells, gene, type, isoType, plateName, growthRate, log2growthRate)
   plateModel = lmer(log2growthRate ~ (1|gene) + (1|isoType) + (1|plateName), data = phenData) # create linear model for colony size
   plateEffects = data.frame(ranef(plateModel)$plateName) %>% mutate(plateName = rownames(.)) # create data frame with plate names and model coefficients
   
   phenData = left_join(phenData, plateEffects, by = 'plateName') # add plate effects column
   phenData %<>% rename(plateEffect = 'X.Intercept.') # rename plate effects column
   
   phenData %<>% mutate(log2growthRatePlateCorrected = log2growthRate - plateEffect) # add column with colony size corrected for plate effects
   
   lmORF = phenData %>% orf_lm(data = ., ref = ref)
   
   phenData = inner_join(phenData, lmORF, by = 'gene')
   
   return(list(phenData = phenData, 
               lmORF = lmORF))
}

val_models <- valLinearModels(data = val_pheno_data_processed, ref = 'wt')
val_models_W303ref <- valLinearModels(data = val_pheno_data_processed, ref = 'W303')
val_models_FRG3ref <- valLinearModels(data = val_pheno_data_processed, ref = 'FRG3')

# combine all similar models into single table
allPhenData_YKOremoved <- bind_rows(hygS_models[['phenData_YKOremoved']], hygR_models[['phenData_YKOremoved']] %>% filter(gene != 'wt') %>% mutate(gene = case_when(gene == 'BY' ~ 'wt', .default = gene)))
allPhenData_YKOremoved %<>% mutate(gene = factor(gene, levels = c('wt', sort(unique(allPhenData_YKOremoved$gene[allPhenData_YKOremoved$gene != 'wt']), decreasing = TRUE))),
                                   type = factor(type, levels = c('S', 'R')))

# output modeled data to TSV
write_tsv(allPhenData_YKOremoved %>% select(gene, type, isoType, growthRate, log2growthRate, plateName, plateEffect, log2growthRatePlateCorrected), 
          file = pheno_model_data_YKOremoved_outFile, col_names = TRUE)

##### ORF MODEL ANALYSIS #####
all_pval_effects <- full_join(hygS_models[['phenData_YKOremoved']] %>% group_by(gene, significance) %>% 
                                 summarize(local_YKOremoved_median = median(log2growthRatePlateCorrected, na.rm = T), local_YKOremoved_pVal = mean(pValueGene, na.rm = T)) %>% ungroup() %>% rename(local_YKOremoved_sig = significance),
                              hygR_models[['phenData_YKOremoved']] %>% group_by(gene, significance) %>% 
                                 summarize(distal_YKOremoved_median = median(log2growthRatePlateCorrected, na.rm = T), distal_YKOremoved_pVal = mean(pValueGene, na.rm = T)) %>% ungroup() %>% rename(distal_YKOremoved_sig = significance),
                              by = 'gene')

all_pval_effects %<>%
   mutate(local_YKOremoved_foldChange = 2^(local_YKOremoved_median - all_pval_effects$local_YKOremoved_median[which(all_pval_effects$gene == 'wt')]),
          distal_YKOremoved_foldChange = 2^(distal_YKOremoved_median - all_pval_effects$distal_YKOremoved_median[which(all_pval_effects$gene == 'BY')])) %>%
   select(gene, local_YKOremoved_median, local_YKOremoved_foldChange, local_YKOremoved_pVal, local_YKOremoved_sig,
          distal_YKOremoved_median, distal_YKOremoved_foldChange, distal_YKOremoved_pVal, distal_YKOremoved_sig)

write_tsv(all_pval_effects, pval_effects_outFile, col_names = T)

all_pval_effects %<>% mutate(gene = factor(gene, levels = c('wt', 'BY', sort(unique(all_pval_effects$gene[all_pval_effects$gene != 'wt' & all_pval_effects$gene != 'BY']), decreasing = TRUE))))

##### SNP MODEL ANALYSIS #####
snp_pval_effects <- full_join(hygS_models[['lmSNP']] %>% select(SNP, SNPEffect, pValueSNP, significance, W303_calls, BY_calls) %>% rename(local_SNPEffect = 'SNPEffect', local_pValueSNP = 'pValueSNP', local_significance = 'significance', local_W303_calls = 'W303_calls', local_BY_calls = 'BY_calls'), 
                              hygR_models[['lmSNP']] %>% select(SNP, SNPEffect, pValueSNP, significance, W303_calls, BY_calls) %>% rename(distal_SNPEffect = 'SNPEffect', distal_pValueSNP = 'pValueSNP', distal_significance = 'significance', distal_W303_calls = 'W303_calls', distal_BY_calls = 'BY_calls'), 
                              by = 'SNP') %>%
   inner_join(., mappingVars_XIVL, join_by(SNP == varID)) %>%
   select(SNP, chr, pos, ref, alt, local_BY_calls, local_W303_calls, local_SNPEffect, local_pValueSNP, local_significance, distal_BY_calls, distal_W303_calls, distal_SNPEffect, distal_pValueSNP, distal_significance) %>%
   arrange(as.numeric(str_split_i(SNP, 'r', 2)))

write_tsv(snp_pval_effects, snp_effects_outFile, col_names = T)

#### phenotyping plots ####
model_boxplot <- function(data, ref) {
   strips <- data %>% group_by(gene) %>% summarise(n = n()) %>% ungroup() %>% select(gene) %>%
      mutate(ymin = -Inf, ymax = Inf, xmin = as.numeric(gene) - 0.5, xmax = as.numeric(gene) + 0.5, fill = rep(c('a', 'b'), length.out = nrow(.))) %>%
      pivot_longer(cols = c(ymin, ymax), values_to = 'y', names_to = 'ymin_ymax') %>% select(!ymin_ymax)
   ggplot(data, aes(gene, log2growthRatePlateCorrected)) + #Initialize plot
      geom_ribbon(data = strips, aes(xmin = xmin, xmax = xmax, y = y, group = gene, fill = fill), alpha = 0.4, inherit.aes = F) +
      geom_hline(yintercept = median(unlist(data %>% filter(gene == ref) %>% select(log2growthRatePlateCorrected))), lty = 1, lwd = 0.2) + #Add median line of wt
      geom_boxplot(aes(fill = significance), outlier.shape = NA, lwd = 0.1) +  #Add box plot
      geom_point(aes(group = isoType), alpha = 0.5, size = 0.2, position = position_dodge(width = 0.6)) + #Add all datapoints grouping by isolate
      geom_line(aes(group = isoType), alpha = 0.4, position = position_dodge(width = 0.6), lwd = 0.1) + #Connect isolates by lines
      scale_fill_manual(breaks = c('Bonferroni', 'Nominal', 'NS'), labels = c('Bonferroni', 'Nominal', 'Not Significant'), values = c('Bonferroni'='#FF9999', 'Nominal'='#FFCC99', 'NS'='#E0E0E0', 'a' = 'white', 'b' = 'gray90')) + #Color code for significance
      facet_wrap(facets = vars(type), ncol = 1) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank(), legend.position = 'bottom', legend.title = element_blank(), strip.background = element_blank(), strip.text = element_blank()) + #Adjust gene labels
      ylab(expression(paste('Growth rate [', log[2]*'(doublings per hour)]'))) + xlab('Targeted ORF')
}

model_boxplot_byRepairType <- function(data, ref) {
   strips <- data %>% group_by(gene) %>% summarise(n = n()) %>% ungroup() %>% select(gene) %>%
      mutate(ymin = -Inf, ymax = Inf, xmin = as.numeric(gene) - 0.5, xmax = as.numeric(gene) + 0.5, fill = rep(c('a', 'b'), length.out = nrow(.))) %>%
      pivot_longer(cols = c(ymin, ymax), values_to = 'y', names_to = 'ymin_ymax') %>% select(!ymin_ymax)
   ggplot(data, aes(gene, log2growthRatePlateCorrected)) + #Initialize plot
      geom_ribbon(data = strips, aes(xmin = xmin, xmax = xmax, y = y, group = gene, fill = fill), alpha = 0.4, inherit.aes = F) +
      geom_hline(yintercept = median(unlist(data %>% filter(gene == ref) %>% select(log2growthRatePlateCorrected))), lty = 1, lwd = 0.2) + #Add median line of wt
      geom_boxplot(aes(fill = significance), outlier.shape = NA, lwd = 0.1) +  #Add box plot
      geom_point(aes(group = isoType), alpha = 0.5, size = 0.2, position = position_dodge(width = 0.6)) + #Add all datapoints grouping by isolate
      geom_line(aes(group = isoType), alpha = 0.4, position = position_dodge(width = 0.6), lwd = 0.1) + #Connect isolates by lines
      scale_fill_manual(breaks = c('Bonferroni', 'Nominal', 'NS'), labels = c('Bonferroni', 'Nominal', 'Not Significant'), values = c('Bonferroni'='#FF9999', 'Nominal'='#FFCC99', 'NS'='#E0E0E0', 'a' = 'white', 'b' = 'gray90')) + #Color code for significance
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank(), legend.position = 'bottom', legend.title = element_blank()) + #Adjust gene labels
      ylab(expression(paste('Growth rate [', log[2]*'(doublings per hour)]'))) + xlab('Targeted ORF')
}

model_boxplot_noFacet <- function(data, ref) {
   ggplot(data, aes(gene, log2growthRatePlateCorrected)) + #Initialize plot
      geom_hline(yintercept = median(unlist(data %>% filter(gene == ref) %>% select(log2growthRatePlateCorrected))), lty = 1, lwd = 0.2) + #Add median line of wt
      geom_boxplot(aes(fill = significance), outlier.shape = NA, lwd = 0.1) +  #Add box plot
      geom_point(aes(group = isoType), alpha = 0.5, size = 0.2, position = position_dodge(width = 0.6)) + #Add all datapoints grouping by isolate
      geom_line(aes(group = isoType), alpha = 0.4, position = position_dodge(width = 0.6), lwd = 0.1) + #Connect isolates by lines
      scale_fill_manual(breaks = c('Bonferroni', 'Nominal', 'NS'), labels = c('Bonferroni', 'Nominal', 'Not Significant'), values = c('Bonferroni'='#FF9999', 'Nominal'='#FFCC99', 'NS'='#E0E0E0')) + #Color code for significance
      theme(panel.grid.major.x = element_blank(), legend.position = 'bottom', legend.title = element_blank()) + #Adjust gene labels
      ylab(expression(paste('Growth rate [', log[2]*'(doublings per hour)]'))) + xlab('Strain')
}

# output plot to PDF
ggsave(file = pheno_plots_YKOremoved_outFile, plot = model_boxplot(data = allPhenData_YKOremoved, ref = 'wt'), width = 17.75, units = 'cm')

##### VALIDATION PLOTS #####
MKT1_SAL1_valPlot <- model_boxplot_noFacet(data = val_models[['phenData']] %>% 
                                              filter(gene %in% c('wt', 'FRG1', 'FRG2', 'FRG3', 'W303')) %>%
                                              mutate(gene = factor(gene, levels = c('wt', 'FRG1', 'FRG2', 'FRG3', 'W303')),
                                                     significance = case_when(significance == 'Bonferroni' ~ 'Nominal', .default = significance)), 
                                           ref = 'wt') + 
   geom_text(data = val_models[['phenData']] %>%
                filter(gene %in% c('FRG1', 'FRG2', 'FRG3', 'W303')) %>%
                mutate(gene = factor(gene, levels = c('FRG1', 'FRG2', 'FRG3', 'W303'))), 
             aes(x = gene, y = -0.15, label = ifelse(pValueGene < 0.000001, formatC(pValueGene, format = 'e', digits = 2), round(pValueGene, 4))), size = 3, check_overlap = TRUE) +
   scale_x_discrete(breaks = c('wt', 'FRG1', 'FRG2', 'FRG3', 'W303'), labels = c('wt', 'MKT1-30G', 'SAL1-nfs', 'MKT1-30G\n+ SAL1-nfs', 'W303')) +
   theme(axis.title.x = element_blank(), legend.position = 'none', axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90))

allFRGs_valPlots <- model_boxplot_noFacet(data = val_models[['phenData']] %>% 
                                             filter(gene %in% c('wt', 'RM', 'W303', 'FRG1', 'FRG2', 'FRG3', 'FRG4', 'FRG5', 'FRG6')) %>%
                                             mutate(gene = factor(gene, levels = c('wt', 'RM', 'W303', 'FRG1', 'FRG2', 'FRG3', 'FRG4', 'FRG5', 'FRG6'))), 
                                          ref = 'wt') + 
   geom_text(data = val_models[['phenData']] %>%
                filter(gene %in% c('RM', 'W303', 'FRG1', 'FRG2', 'FRG3', 'FRG4', 'FRG5', 'FRG6')) %>%
                mutate(gene = factor(gene, levels = c('RM', 'W303', 'FRG1', 'FRG2', 'FRG3', 'FRG4', 'FRG5', 'FRG6'))), 
             aes(x = gene, y = -0.15, label = ifelse(pValueGene < 0.000001, formatC(pValueGene, format = 'e', digits = 2), round(pValueGene, 4))), size = 3, check_overlap = TRUE) +
   theme(axis.title.x = element_blank())

# output plot to PDF
ggsave(file = MKT1_SAL1_valPlot_outFile, plot = MKT1_SAL1_valPlot, width = 8.5, height = 8.5, units = 'cm')

##### HphMX TEST PLOT #####
hphmx_plot <- model_boxplot_noFacet(data = hygR_models[['phenData_YKOremoved']] %>% filter(gene %in% c('wt', 'BY')) %>% mutate(gene = factor(gene, levels = c('BY', 'wt'))),
                                    ref = 'BY') + scale_x_discrete(breaks = c('BY', 'wt'), labels = c('CRI-SPA-Map Wildtype', 'BY4742 + HphMX')) +
   theme(axis.title.x = element_blank(), legend.position = 'none')

ggsave(filename = hygTest_plot_outFile, plot = hphmx_plot, width = 8.5, height = 8.5, units = 'cm')

#### SNP plot ####
# function to extract genotype runs for each isolate and combine into single dataframe
geno_runs <- function(simp_calls, min_length = 0, ignore_call = NULL) {
   runs_byIso = apply(simp_calls[, -1], 2, function(iso){ # process each column (isolate) separately, not including varID column
      runs = rle(iso) # extract runs of consecutive calls
      runsDF = data.frame(runs$values, runs$lengths) # convert runs to dataframe format
      colnames(runsDF) <- c('call', 'length') # rename columns
      runsDF$end <- cumsum(runs$lengths) # add column signifying which variant the run ended
      runsDF$start <- runsDF$end - runsDF$length + 1 # add column signifying which variant the run started
      runsDF %<>% mutate(start_var = simp_calls$varID[start], # use varID column of initial data to report start and end variant information
                         end_var = simp_calls$varID[end])
      
      runsDF %<>% filter((length >= min_length & !(call %in% ignore_call)) | end_var == 'CEN14') # filter out short runs and those with undesired calls
      
      runs2 = rle(runsDF$call) # extract runs of desired calls
      runsDF2 = data.frame(runs2$values, runs2$lengths) # convert new runs to dataframe format
      colnames(runsDF2) <- c('call', 'length') # rename columns
      runsDF2$end <- cumsum(runs2$lengths) # add column signifying where new run ends
      runsDF2$start <- runsDF2$end - runsDF2$length + 1 # add column signifying where new run starts
      runsDF2 %<>% mutate(start_var = runsDF$start_var[start], # use varID column of initial data to report start and end variant information
                          end_var = runsDF$end_var[end])
      runsDF2 %<>% select(call, length, start_var, end_var) # extract only desired columns
      
      return(runsDF2)
   })
   
   names(runs_byIso) <- colnames(simp_calls[-1]) # name run dataframes using the column names of initial data (isolates)
   
   allRuns = bind_rows(runs_byIso, .id = 'isolate') # combine all run dataframes into single table, add isolate column
   allRuns <- allRuns %>% inner_join(., strainInfo %>% select(corrected_isoName, orf, repair_type, isoNum, isoType), join_by(isolate == corrected_isoName)) %>%
      inner_join(., mappingVars_XIVL %>% select(varID, pos), join_by(start_var == varID)) %>%
      rename(start_pos = pos) %>%
      inner_join(., mappingVars_XIVL %>% select(varID, pos), join_by(end_var == varID)) %>%  
      rename(end_pos = pos)
   
   return(allRuns)
}

repair_tracts <- geno_runs(simplified_mappingVar_calls_XIVL %>% filter(varID != 'var7325') %>% select(varID, starts_with('YNL')), min_length = 0, ignore_call = c('no_call', 'low_qual')) %>%
   mutate(plot_type = case_when(start_var == end_var ~ 'point', .default = 'line')) %>%
   filter(call == 'W303') %>% 
   inner_join(., bind_rows(hygS_models[['phenData_genos']], hygR_models[['phenData_genos']]) %>% select(isoType, log2growthRatePlateCorrected), by = 'isoType')

snp_plot <- ggplot() + 
   geom_hline(yintercept = median(hygS_models[['phenData_genos']] %>% filter(sequenced_illumina == 'yes') %>% pull(log2growthRatePlateCorrected), na.rm = TRUE), color = 'gray30', lwd = 0.25) +
   geom_hline(yintercept = (-log10(0.000128) - 7.2)/12, linetype = 'dashed', lwd = 0.25, color = 'purple') +
   geom_hline(yintercept = (-log10(0.05/nrow(hygS_models[['lmSNP']])) - 7.2)/12, linetype = 'dashed', lwd = 0.25, color = 'red') +
   geom_rug(data = mappingVars_XIVL %>% filter(varID != 'var7325'), aes(x = pos), color = 'black', linewidth = 0.2, inherit.aes = FALSE) + ## variant tract
   geom_point(data = repair_tracts %>% filter(repair_type == 'S' & plot_type == 'point'), aes(x = start_pos, y = log2growthRatePlateCorrected), color =  '#50C878', shape = 'square', size = 0.1) +
   geom_segment(data = repair_tracts %>% filter(repair_type == 'S'), 
                aes(x = start_pos, xend = end_pos, y = log2growthRatePlateCorrected), color = '#50C878', linewidth = 0.5) + 
   geom_line(data = hygS_models[['lmSNP']], aes(x = pos, y = (-log10(pValueSNP) - 7.2) / 12), color = 'black', linewidth = 0.5) +
   scale_y_continuous((expression(paste('Growth rate [', log[2]*'(doublings per hour)]'))), sec.axis = sec_axis(transform = ~ ((12 * .) + 7.2), name = expression(-log[10]*'(p-value)'))) +
   xlab('Position on chromosome XIV (bp)') +
   theme(axis.title.y.right = element_text(size = 10), legend.position = 'bottom')

snp_plot_zoom <- ggplot() + 
   geom_hline(yintercept = median(hygS_models[['phenData_genos']] %>% filter(sequenced_illumina == 'yes') %>% pull(log2growthRatePlateCorrected), na.rm = TRUE), color = 'gray30', lwd = 0.25) +
   geom_rect(data = left_join(chrXIVLorf_info, hygS_models[['lmORF_YKOremoved']] %>% select(gene, significance), join_by(orf == gene)) %>% 
                mutate(significance = case_when(is.na(significance) ~ 'NT', .default = significance)), 
             aes(xmin = start, xmax = stop, ymin = -0.6, ymax = -0.58, fill = significance), linewidth = 0.3, color = '#E0E0E0') +
   geom_hline(yintercept = (-log10(0.000128) - 7.2)/12, linetype = 'dashed', lwd = 0.25, color = 'purple') +
   geom_hline(yintercept = (-log10(0.05/nrow(hygS_models[['lmSNP']])) - 7.2)/12, linetype = 'dashed', lwd = 0.25, color = 'red') +
   geom_rug(data = mappingVars_XIVL %>% filter(varID != 'var7325'), aes(x = pos), color = 'black', linewidth = 0.2, inherit.aes = FALSE) + ## variant tract
   geom_point(data = repair_tracts %>% filter(repair_type == 'S' & plot_type == 'point'), aes(x = start_pos, y = log2growthRatePlateCorrected), color =  '#50C878', shape = 'square', size = 0.1) +
   geom_segment(data = repair_tracts %>% filter(repair_type == 'S'), 
                aes(x = start_pos, xend = end_pos, y = log2growthRatePlateCorrected), color = '#50C878', linewidth = 0.5) + 
   geom_line(data = hygS_models[['lmSNP']], aes(x = pos, y = (-log10(pValueSNP) - 7.2) / 12), color = 'black', linewidth = 0.5) +
   scale_fill_manual(breaks = c('Bonferroni', 'Nominal', 'NS', 'NT'), labels = c('Bonferroni', 'Nominal', 'Not Significant', 'Not Targeted'), values = c('Bonferroni'='#FF9999', 'Nominal'='#FFCC99', 'NS'='#E0E0E0', 'NT' = 'white')) + #Color code for significance
   scale_y_continuous((expression(paste('Growth rate [', log[2]*'(doublings per hour)]'))), sec.axis = sec_axis(transform = ~ ((12 * .) + 7.2), name = expression(-log[10]*'(p-value)'))) +
   xlab('Position on chromosome XIV (bp)') +
   theme(axis.title.y.right = element_text(size = 10), legend.position = 'none') + 
   coord_cartesian(xlim = c(459000, 477000)) +
   scale_x_continuous(breaks = seq(459000, 477000, by = 9000))

ggsave(filename = SNP_plot_outFile, plot = snp_plot, width = 11.5, height = 9.5, units = 'cm')
ggsave(filename = SNP_plot_zoom_outFile, plot = snp_plot_zoom, width = 7.5, height = 9.5, units = 'cm')
