# R version 4.5.1

# load libraries
library(tidyverse) ## v 2.0.0
library(magrittr) ## v 2.0.3
library(extrafont) ## v 0.19

# user-input file paths and info
chrXIVLorf_inFile <- 'data/S288C_chrXIVL_ORFs.tsv' ## file containing info on all ORFs found on chrXIV-L in S288C reference

outDir <- 'tables/' ## directory for all analysis output files
figDir <- 'figures/' ## directory for all figure output files
SNP_plot_outFile <- paste0(figDir, 'SNP_plots.pdf')
SNP_plot_zoom_outFile <- paste0(figDir, 'SNP_plots_zoom.pdf')
distal_SNP_plot_outFile <- paste0(figDir, 'distal_SNP_plots.pdf')
distal_SNP_plot_zoom_outFile <- paste0(figDir, 'distal_SNP_plots_zoom.pdf')

options(scipen = 999) ## no scientific notation
theme_set(theme_bw(base_size = 10, base_family = 'Arial') + ## set base ggplot fonts
             theme(text = element_text(color = 'black'), axis.text = element_text(color = 'black', size = 8)))

# make output directories
dir.create(outDir, showWarnings = FALSE)
dir.create(figDir, showWarnings = FALSE)

# load in previously-created files
simplified_mappingVar_calls_XIVL <- read_tsv(paste0(outDir, 'CSM_chrXIVLcalls.tsv'), col_names = TRUE, show_col_types = FALSE)
mappingVars_XIVL <- read_tsv(paste0(outDir, 'mappingVars_XIVL.tsv'), col_names = TRUE, show_col_types = FALSE)
strainInfo <- read_tsv(paste0(outDir, 'strainInfo_wFilters.tsv'), col_names = TRUE, show_col_types = FALSE)
local_phenData_genos <- read_tsv(paste0(outDir, 'local_phenData_genos.tsv'), col_names = TRUE, show_col_types = FALSE)
distal_phenData_genos <- read_tsv(paste0(outDir, 'distal_phenData_genos.tsv'), col_names = TRUE, show_col_types = FALSE)
local_lmSNP <- read_tsv(paste0(outDir, 'local_lmSNP.tsv'), col_names = TRUE, show_col_types = FALSE)
distal_lmSNP <- read_tsv(paste0(outDir, 'distal_lmSNP.tsv'), col_names = TRUE, show_col_types = FALSE)
local_lmORF <- read_tsv(paste0(outDir, 'local_lmORF_YKOremoved.tsv'), col_names = TRUE, show_col_types = FALSE)
distal_lmORF <- read_tsv(paste0(outDir, 'distal_lmORF_YKOremoved.tsv'), col_names = TRUE, show_col_types = FALSE)

# load in all chrXIV-L ORFs
chrXIVLorf_info <- read_tsv(chrXIVLorf_inFile, show_col_types = FALSE) %>% ## all S288C ORFs on chrXIVL from SGD/alliancemine
   rename(orf = SequenceFeature.secondaryIdentifier, ## rename columns
          gene_symbol = SequenceFeature.symbol, 
          gene_name = SequenceFeature.name, 
          sgd_id = SequenceFeature.primaryIdentifier, 
          start = SequenceFeature.chromosomeLocation.start, 
          stop = SequenceFeature.chromosomeLocation.end, 
          chr = SequenceFeature.chromosomeLocation.locatedOn.primaryIdentifier) %>%
   select(orf, gene_symbol, gene_name, chr, start, stop, sgd_id) ## reorder columns

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
transpose_mappingVar_calls_XIVL <- transpose_calls(simplified_mappingVar_calls_XIVL %>% select(varID, starts_with('YNL')))

# manually adjust calls at sal1-1 SNP due to poor calling
transpose_mappingVar_calls_XIVL %<>% mutate(var7523 = case_when(var7523 == 'low_qual'~ 'BY', .default = var7523)) ## change all low quality calls to BY
transpose_mappingVar_calls_XIVL %<>% mutate(var7523 = case_when(strain %in% c('YNL083W_B8S1_4A7', 'YNL083W_B8S3_4B7', 'YNL083W_B8S4_4C7', 'YNL084C_G6S1_3G10', 'YNL084C_G6S3_3A11', 'YNL085W_C8S4_4B8', 'YNL084C_G6R1_3C11', 'YNL084C_G6R2_3D11', 'YNL085W_C8R1_4C8') ~ 'W303', .default = var7523)) ## modify isolates that targeted ORFs near SNP that were determined to be W303 by visual inspection
transpose_mappingVar_calls_XIVL %<>% mutate(var7523 = case_when(orf %in% c('YNL005C', 'YNL008C', 'YNL035C', 'YNL054W', 'YNL055C', 'YNL056W', 'YNL067W', 'YNL070W', 'YNL071W', 'YNL073W', 'YNL074C', 'YNL076W', 'YNL077W', 'YNL078W', 'YNL079C', 'YNL081C', 'YNL082W', 'YNL083W') & str_detect(str_split_i(strain, '_', 2), 'R') ~ 'W303', .default = var7523)) ## change distal isolates with targeted ORFs centromeric to SNP to be W303
transpose_mappingVar_calls_XIVL %<>% mutate(var7523 = case_when(!(orf %in% c('YNL005C', 'YNL008C', 'YNL035C', 'YNL054W', 'YNL055C', 'YNL056W', 'YNL067W', 'YNL070W', 'YNL071W', 'YNL073W', 'YNL074C', 'YNL076W', 'YNL077W', 'YNL078W', 'YNL079C', 'YNL081C', 'YNL082W', 'YNL083W')) & str_detect(str_split_i(strain, '_', 2), 'R') & var7523 == 'no_call' ~ 'BY', .default = var7523)) ## change distal isolates with targeted ORFs telomeric and "no_call" designations to SNP to be BY
simplified_mappingVar_calls_XIVL[which(simplified_mappingVar_calls_XIVL$varID == 'var7523'),] <- c(simplified_mappingVar_calls_XIVL[which(simplified_mappingVar_calls_XIVL$varID == 'var7523'),1:23], transpose_mappingVar_calls_XIVL$var7523) ## apply above adjustments to simplified call table for accurate plotting

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
   inner_join(., bind_rows(local_phenData_genos, distal_phenData_genos) %>% select(isoType, log2growthRatePlateCorrected), by = 'isoType')

snp_plot <- ggplot() + 
   geom_hline(yintercept = median(local_phenData_genos %>% filter(sequenced_illumina == 'yes') %>% pull(log2growthRatePlateCorrected), na.rm = TRUE), color = 'gray30', lwd = 0.25) +
   geom_hline(yintercept = (-log10(0.000207) - 7.2)/12, linetype = 'dashed', lwd = 0.25, color = 'purple') +
   geom_hline(yintercept = (-log10(0.05/nrow(local_lmSNP)) - 7.2)/12, linetype = 'dashed', lwd = 0.25, color = 'red') +
   geom_rug(data = mappingVars_XIVL %>% filter(varID != 'var7325'), aes(x = pos), color = 'black', linewidth = 0.2, inherit.aes = FALSE) + ## variant tract
   geom_point(data = repair_tracts %>% filter(repair_type == 'S' & plot_type == 'point'), aes(x = start_pos, y = log2growthRatePlateCorrected), color =  '#50C878', shape = 'square', size = 0.1) +
   geom_segment(data = repair_tracts %>% filter(repair_type == 'S'), 
                aes(x = start_pos, xend = end_pos, y = log2growthRatePlateCorrected), color = '#50C878', linewidth = 0.5) + 
   geom_line(data = local_lmSNP, aes(x = pos, y = (-log10(pValueSNP) - 7.2) / 12), color = 'black', linewidth = 0.5) +
   scale_y_continuous((expression(paste('Growth rate [', log[2]*'(doublings per hour)]'))), sec.axis = sec_axis(transform = ~ ((12 * .) + 7.2), name = expression(-log[10]*'(p-value)'))) +
   xlab('Position on chromosome XIV (bp)') +
   theme(axis.title.y.right = element_text(size = 10), legend.position = 'bottom')

snp_plot_zoom <- ggplot() + 
   geom_hline(yintercept = median(local_phenData_genos %>% filter(sequenced_illumina == 'yes') %>% pull(log2growthRatePlateCorrected), na.rm = TRUE), color = 'gray30', lwd = 0.25) +
   geom_rect(data = left_join(chrXIVLorf_info, local_lmORF %>% select(gene, significance), join_by(orf == gene)) %>% 
                mutate(significance = case_when(is.na(significance) ~ 'NT', .default = significance)), 
             aes(xmin = start, xmax = stop, ymin = -0.6, ymax = -0.58, fill = significance), linewidth = 0.3, color = '#E0E0E0') +
   geom_hline(yintercept = (-log10(0.000207) - 7.2)/12, linetype = 'dashed', lwd = 0.25, color = 'purple') +
   geom_hline(yintercept = (-log10(0.05/nrow(local_lmSNP)) - 7.2)/12, linetype = 'dashed', lwd = 0.25, color = 'red') +
   geom_rug(data = mappingVars_XIVL %>% filter(varID != 'var7325'), aes(x = pos), color = 'black', linewidth = 0.2, inherit.aes = FALSE) + ## variant tract
   geom_point(data = repair_tracts %>% filter(repair_type == 'S' & plot_type == 'point'), aes(x = start_pos, y = log2growthRatePlateCorrected), color =  '#50C878', shape = 'square', size = 0.1) +
   geom_segment(data = repair_tracts %>% filter(repair_type == 'S'), 
                aes(x = start_pos, xend = end_pos, y = log2growthRatePlateCorrected), color = '#50C878', linewidth = 0.5) + 
   geom_line(data = local_lmSNP, aes(x = pos, y = (-log10(pValueSNP) - 7.2) / 12), color = 'black', linewidth = 0.5) +
   scale_fill_manual(breaks = c('Bonferroni', 'Nominal', 'NS', 'NT'), labels = c('Bonferroni', 'Nominal', 'Not Significant', 'Not Targeted'), values = c('Bonferroni'='#FF9999', 'Nominal'='#FFCC99', 'NS'='#E0E0E0', 'NT' = 'white')) + #Color code for significance
   scale_y_continuous((expression(paste('Growth rate [', log[2]*'(doublings per hour)]'))), sec.axis = sec_axis(transform = ~ ((12 * .) + 7.2), name = expression(-log[10]*'(p-value)'))) +
   xlab('Position on chromosome XIV (bp)') +
   theme(axis.title.y.right = element_text(size = 10), legend.position = 'none') + 
   coord_cartesian(xlim = c(459000, 477000)) +
   scale_x_continuous(breaks = seq(459000, 477000, by = 9000))

ggsave(filename = SNP_plot_outFile, plot = snp_plot, width = 11.5, height = 9.5, units = 'cm')
ggsave(filename = SNP_plot_zoom_outFile, plot = snp_plot_zoom, width = 8, height = 9.5, units = 'cm')

distal_snp_plot <- ggplot() + 
   geom_hline(yintercept = median(distal_phenData_genos %>% filter(sequenced_illumina == 'yes') %>% pull(log2growthRatePlateCorrected), na.rm = TRUE), color = 'gray30', lwd = 0.25) +
   geom_hline(yintercept = (-log10(0.001192) - 12)/(12/0.65), linetype = 'dashed', lwd = 0.25, color = 'purple') +
   geom_hline(yintercept = (-log10(0.05/nrow(distal_lmSNP)) - 12)/(12/0.65), linetype = 'dashed', lwd = 0.25, color = 'red') +
   geom_rug(data = mappingVars_XIVL %>% filter(varID != 'var7325'), aes(x = pos), color = 'black', linewidth = 0.2, inherit.aes = FALSE) + ## variant tract
   geom_point(data = repair_tracts %>% filter(repair_type == 'R'), aes(x = end_pos, y = log2growthRatePlateCorrected), color =  '#50C878', shape = 'square', size = 0.5) +
   geom_line(data = distal_lmSNP, aes(x = pos, y = (-log10(pValueSNP) - 12)/(12/0.65)), color = 'black', linewidth = 0.5) +
   scale_y_continuous((expression(paste('Growth rate [', log[2]*'(doublings per hour)]'))), sec.axis = sec_axis(transform = ~ (((12/0.65) * .) + 12), name = expression(-log[10]*'(p-value)'))) +
   xlab('Position on chromosome XIV (bp)') +
   theme(axis.title.y.right = element_text(size = 10), legend.position = 'bottom')

distal_snp_plot_zoom <- ggplot() + 
   geom_hline(yintercept = median(distal_phenData_genos %>% filter(sequenced_illumina == 'yes') %>% pull(log2growthRatePlateCorrected), na.rm = TRUE), color = 'gray30', lwd = 0.25) +
   geom_rect(data = left_join(chrXIVLorf_info, distal_lmORF %>% select(gene, significance), join_by(orf == gene)) %>% 
                mutate(significance = case_when(is.na(significance) ~ 'NT', .default = significance)), 
             aes(xmin = start, xmax = stop, ymin = -0.6, ymax = -0.58, fill = significance), linewidth = 0.3, color = '#E0E0E0') +
   geom_hline(yintercept = (-log10(0.001192) - 12)/(12/0.65), linetype = 'dashed', lwd = 0.25, color = 'purple') +
   geom_hline(yintercept = (-log10(0.05/nrow(distal_lmSNP)) - 12)/(12/0.65), linetype = 'dashed', lwd = 0.25, color = 'red') +
   geom_rug(data = mappingVars_XIVL %>% filter(varID != 'var7325'), aes(x = pos), color = 'black', linewidth = 0.2, inherit.aes = FALSE) + ## variant tract
   geom_point(data = repair_tracts %>% filter(repair_type == 'R'), aes(x = end_pos, y = log2growthRatePlateCorrected), color =  '#50C878', shape = 'square', size = 0.5) +
   geom_line(data = distal_lmSNP, aes(x = pos, y = (-log10(pValueSNP) - 12)/(12/0.65)), color = 'black', linewidth = 0.5) +
   scale_fill_manual(breaks = c('Bonferroni', 'Nominal', 'NS', 'NT'), labels = c('Bonferroni', 'Nominal', 'Not Significant', 'Not Targeted'), values = c('Bonferroni'='#FF9999', 'Nominal'='#FFCC99', 'NS'='#E0E0E0', 'NT' = 'white')) + #Color code for significance
   scale_y_continuous((expression(paste('Growth rate [', log[2]*'(doublings per hour)]'))), sec.axis = sec_axis(transform = ~ (((12/0.65) * .) + 12), name = expression(-log[10]*'(p-value)'))) +
   xlab('Position on chromosome XIV (bp)') +
   theme(axis.title.y.right = element_text(size = 10), legend.position = 'none') + 
   coord_cartesian(xlim = c(459000, 477000)) +
   scale_x_continuous(breaks = seq(459000, 477000, by = 9000))

ggsave(filename = distal_SNP_plot_outFile, plot = distal_snp_plot, width = 11.5, height = 9.5, units = 'cm')
ggsave(filename = distal_SNP_plot_zoom_outFile, plot = distal_snp_plot_zoom, width = 8, height = 9.5, units = 'cm')