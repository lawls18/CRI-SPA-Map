# R version 4.5.1

# load libraries
library(tidyverse) ## v 2.0.0
library(magrittr) ## v 2.0.3
library(furrr) ## v 0.3.1
library(lme4) ## v 1.1-37
library(lmerTest) ## v 3.1-3

# user-input file paths and info
outDir <- 'tables/' ## directory for all analysis output files

options(scipen = 999) ## no scientific notation
plan(multisession, workers = 16) ## set parallel parameters
set.seed(18)

# make output directory
dir.create(outDir, showWarnings = FALSE)

# load in previously-created files
local_phenData_genos <- read_tsv(paste0(outDir,'local_phenData_genos.tsv'), col_names = TRUE, show_col_types = FALSE)
distal_phenData_genos <- read_tsv(paste0(outDir,'distal_phenData_genos.tsv'), col_names = TRUE, show_col_types = FALSE)

# define function to permute growth rate across isolates
SNP_perms <- function(data, snps) {
   lmResults_SNP <- bind_cols(lapply(snps, function(thisSNP) {
      if(nrow(data %>% filter(get(thisSNP) == 'W303')) > 0) {
         testDat <- data %>% filter((get(thisSNP) == 'W303' | get(thisSNP) == 'BY') & !is.na(get(thisSNP))) %>% select(all_of(thisSNP), 'perm_growthRate', 'isoType')
         colnames(testDat) <- c('SNP', 'perm_growthRate', 'isoType')
         fullModel <- lmer(perm_growthRate ~ SNP + (1|isoType), data = testDat, REML = FALSE)
         res <- anova(fullModel)[1,'Pr(>F)'] # lmer
         return(res)
      } else { return(NA) }
   }))
   colnames(lmResults_SNP) <- snps
   return(lmResults_SNP)
}

# perform 1000 permutations with local isolates
local_perms <- future_map(1:1000, ## set number of permutations
                          ~ local_phenData_genos %>% ## load data set
                             select(isoType, log2growthRatePlateCorrected, starts_with('var7')) %>% ## select columns with isolate identifier and variant calls
                             mutate(perm_growthRate = sample(log2growthRatePlateCorrected, size = n(), replace = FALSE)) %>% ## sample corrected growth rate without replacement
                             SNP_perms(data = ., snps = colnames(.)[3:ncol(.)]), .progress = TRUE) %>% ## run SNP model on permuted data
   bind_rows() ## merge all permutations into single table

# perform 1000 permutations with distal isolates
distal_perms <- future_map(1:1000, ## set number of permutations
                           ~ distal_phenData_genos %>% ## load data set
                              select(isoType, log2growthRatePlateCorrected, starts_with('var7')) %>% ## select columns with isolate identifier and variant calls
                              mutate(perm_growthRate = sample(log2growthRatePlateCorrected, size = n(), replace = FALSE)) %>% ## sample corrected growth rate without replacement
                              SNP_perms(data = ., snps = colnames(.)[3:ncol(.)]), .progress = TRUE) %>% ## run SNP model on permuted data
   bind_rows() ## merge all permutations into single table

# remove last column - all NAs, does not represent SNP
local_perms %<>% select(!perm_growthRate)
distal_perms %<>% select(!perm_growthRate)

# extract the minimum p-value from each permutation and find the 5% threshold from resulting distribution
quantile(apply(local_perms, 1, function(x) { min(x, na.rm = TRUE)}), probs = c(0.05))
quantile(apply(distal_perms, 1, function(x) { min(x, na.rm = TRUE)}), probs = c(0.05))