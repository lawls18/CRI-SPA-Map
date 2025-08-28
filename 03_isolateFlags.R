# R version 4.5.1

# load libraries
library(tidyverse) ## v 2.0.0
library(magrittr) ## v 2.0.3
library(readxl) ## v 1.4.3
library(janitor) ## v 2.2.0
library(extrafont) ## v 0.19
library(patchwork) ## v 1.3.0

# user-input file paths and info
strainInfo_inFile <- 'data/CRI_SPA_MAP_INFO.xlsx' ## table with details on each isolate
chrLen_inFile <- 'data/sacCer_chrLen.csv' ## table listing reference chromosome lengths
cov_inDir <- 'data/CNVnator/' ## directory with read coverage files
cov_inFiles <- 'cnv.txt' ## identifier for read coverage files

outDir <- '250820_tables/' ## directory for all analysis output files
figDir <- '250820_figures/' ## directory for all figure output files

strainInfo_wFilters_outFile <- paste0(outDir, 'strainInfo_wFilters.tsv') ## output file for isolate table with filter columns
recomElsewhere_plot_outFile <- paste0(figDir, 'recomElsewhere.pdf') ## output file for SFigure 3
repairTracts_outFile <- paste0(outDir, 'repairTracts.tsv') ## output file for Figure 2A
armplots_outFile <- paste0(figDir, 'armplots.pdf') ## output file for SFigure 2

options(scipen = 999) ## no scientific notation
theme_set(theme_classic(base_size = 10, base_family = 'Arial') + ## set base ggplot fonts
             theme(text = element_text(color = 'black'), axis.text = element_text(color = 'black', size = 8)))

# make output directories
dir.create(outDir, showWarnings = FALSE)
dir.create(figDir, showWarnings = FALSE)

# load in previously-created files
simplified_mappingVar_calls <- read_tsv(paste0(outDir, 'CSM_allCalls.tsv'), col_names = TRUE, show_col_types = FALSE)
simplified_mappingVar_calls_XIVL <- read_tsv(paste0(outDir, 'CSM_chrXIVLcalls.tsv'), col_names = TRUE, show_col_types = FALSE)
CSMorf_info <- read_tsv(paste0(outDir, 'CSM_orfInfo.tsv'), col_names = TRUE, show_col_types = FALSE)
mappingVars <- read_tsv(paste0(outDir, 'mappingVars.tsv'), col_names = TRUE, show_col_types = FALSE)
mappingVars_XIVL <- read_tsv(paste0(outDir, 'mappingVars_XIVL.tsv'), col_names = TRUE, show_col_types = FALSE)

##### BY VARIANTS IN ORF #####
# define function to transpose and reformat simplified call table
transpose_calls <- function(simp_calls) {
   transpose_calls = as.data.frame(cbind(colnames(simp_calls)[2:ncol(simp_calls)], t(simp_calls[,-1]))) ## transpose calls and add column for strain names
   rownames(transpose_calls) <- NULL ## remove row names
   colnames(transpose_calls) <- c('strain', simp_calls$varID) ## name columns with variant IDs
   transpose_calls %<>% mutate(orf = str_split_i(strain, '_', 1)) ## indicate targeted ORF of each strain
   transpose_calls %<>% arrange(strain) ## sort table by strain name
   return(transpose_calls)
}

# transpose genotype call tables
transpose_mappingVar_calls <- transpose_calls(simplified_mappingVar_calls %>% select(varID, starts_with('YNL')))
transpose_mappingVar_calls_XIVL <- transpose_calls(simplified_mappingVar_calls_XIVL %>% select(varID, starts_with('YNL')))

# merge genotype calls and ORF table
orf_vars_validation <- left_join(CSMorf_info, transpose_mappingVar_calls, by = 'orf')

# determine how many of each call was made within the ORF for every isolate
orf_vars_validation$BY_vars <- apply(orf_vars_validation, 1, function(x) { sum(x[unlist(str_split(x['mappingVar_IDs'], ','))] == 'BY')} )
orf_vars_validation$W303_vars <- apply(orf_vars_validation, 1, function(x) { sum(x[unlist(str_split(x['mappingVar_IDs'], ','))] == 'W303')} )
orf_vars_validation$noCall_vars <- apply(orf_vars_validation, 1, function(x) { sum(x[unlist(str_split(x['mappingVar_IDs'], ','))] == 'no_call')} )

# add column to indicate genotype calls within targeted ORF
orf_vars_validation %<>% mutate(orf_vars = case_when(
   BY_vars > 0 & W303_vars == 0 ~ 'BY',
   BY_vars > 0 & W303_vars > 0 ~ 'split',
   W303_vars > 0 & BY_vars == 0 ~ 'W303',
   .default = NA))

# identify isolates with BY calls within the target ORF
BY_orf_strains <- orf_vars_validation %>% 
   filter(orf_vars != 'W303') %>% 
   pull(strain)

##### RECOMBINATION ELSEWHERE ######
# load in file with chromosome lengths from S288C reference
chrLen <- read_csv(chrLen_inFile, show_col_types = FALSE) %>%
   mutate(prev_len = case_when(chr == 'sacCer_chrI' ~ 0, chr == 'sacCer_chrM' ~ NA, .default = lag(cumsum(.$len))), ## add column with cumulative sum of chromsomes before given
          midpoint = len/2, ## add column for center of chromosome
          cum_mid = midpoint + prev_len, ## add column for center of chromosome as cumulative position
          numerals = str_split_i(chr, 'chr', 2)) ## add column to extract roman numerals from chromosome name

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
   allRuns %<>% mutate(
      orf = str_split_i(isolate, '_', 1), # add column extracting target ORF
      iso_id = str_split_i(isolate, '_', 2), # add column extracting isolate ID
      repair_type = case_when( # add column signifying isolate's repair category
         grepl('^[A-Z][0-9]+R[0-9]+', iso_id) ~ 'distal',
         grepl('^[A-Z][0-9]+S[0-9]+', iso_id) ~ 'local',
         grepl('^[A-Z][0-9]+D[0-9]+', iso_id) ~ 'deletion'
      )
   )
   return(allRuns)
}

# identify W303 runs on other chromosome arms
recom_elsewhere <- geno_runs(simplified_mappingVar_calls %>% select(varID, starts_with('YNL')), min_length = 3, ignore_call = c('low_qual', 'no_call'))

recom_elsewhere_isos <- recom_elsewhere %>% 
   filter(call == 'W303' & repair_type != 'deletion' & !(call == 'W303' & (start_var %in% mappingVars_XIVL$varID | end_var %in% mappingVars_XIVL$varID))) %>% ## extract only runs of W303 on other arms and do not include sequenced deletion strains
   pull(isolate) %>% unique()
recom_elsewhere %<>% filter(isolate %in% recom_elsewhere_isos) ## limit table to list all runs of the isolates extracted above

# add coordinate information for runs
recom_elsewhere <- left_join(recom_elsewhere, mappingVars %>% select(varID, chr, pos), join_by(start_var == 'varID')) %>% rename(start_chr = chr, start_pos = pos)
recom_elsewhere <- left_join(recom_elsewhere, mappingVars %>% select(varID, chr, pos), join_by(end_var == 'varID')) %>% rename(end_chr = chr, end_pos = pos)

# add columns for cumulative positions of runs
recom_elsewhere <- recom_elsewhere %>% 
   left_join(., chrLen %>% select(chr, prev_len), join_by(start_chr == chr)) %>%
   mutate(cum_start = start_pos + prev_len) %>% ## add cumulative position for start of run
   select(!prev_len) %>%
   left_join(., chrLen %>% select(chr, prev_len), join_by(end_chr == chr)) %>%
   mutate(cum_end = end_pos + prev_len) %>% ## add cumulative position for end of run
   select(!prev_len) %>%
   mutate(short_iso = str_sub(iso_id, start = -2, end = -1), orf_iso = paste(orf, short_iso, sep = '_')) ## add columns to shorten isolate identifiers

##### COPY NUMBER VARIATIONS #####
# load in files
iso_covFiles <- list.files(path = cov_inDir, pattern = cov_inFiles)
iso_cov <- lapply(iso_covFiles, function(x) {
   read_tsv(paste0(cov_inDir, x), col_names = c('cnv_type', 'cnv_loc', 'cnv_len', 'cov', 'eval1', 'eval2', 'eval3', 'eval4', 'frac_q0'), show_col_types = FALSE) %>%
      select(cnv_type, cnv_loc, cnv_len, cov, eval1) %>% ## extract only necessary columns
      mutate(isolate = str_split_i(x, '\\.', 1), orf = str_split_i(isolate, '_', 1), isolate = str_split_i(isolate, '_', 2)) %>% ## extract isolate information from file name
      separate(col = cnv_loc, into = c('cnv_chr', 'cnv_pos'), sep = ':') %>% ## separate CNV location into chromosome and position
      separate(col = cnv_pos, into = c('cnv_start', 'cnv_end'), sep = '-') %>% ## separate CNV position into start and end
      left_join(., chrLen, join_by(cnv_chr == chr)) %>%
      mutate(cnv_start = as.integer(cnv_start), cnv_end = as.integer(cnv_end),
             cum_start = cnv_start + prev_len,
             cum_end = cnv_end + prev_len) ## add cumulative positions for CNVs
}) %>%
   bind_rows()

aneuploid_isos <- iso_cov %>% filter(cnv_len >= len/2 & eval1 < 0.05/nrow(.)) %>% pull(isolate) ## extract isolates that did not meet assigned thresholds
aneuploid_isos <- c(aneuploid_isos, 'H7R1') ## add isolate identified by manual analysis

##### REPAIR TRACTS #####
# add row for column to each isolate for BY call at chrXIV centromere
simplified_mappingVar_calls_XIVL %<>% add_row(varID = 'CEN14', chr = 'sacCer_chrXIV', pos = 628758)
simplified_mappingVar_calls_XIVL[nrow(simplified_mappingVar_calls_XIVL), which(str_detect(colnames(simplified_mappingVar_calls_XIVL), 'YNL'))] <- 'BY'

# extract chrXIVL repair tracts using runs of genotype calls
repairTracts <- geno_runs(simplified_mappingVar_calls_XIVL %>% select(varID, starts_with('YNL')), min_length = 3, ignore_call = c('no_call', 'low_qual'))

# add coordinate information for run start and end
repairTracts <- left_join(repairTracts, simplified_mappingVar_calls_XIVL %>% select(varID, pos), join_by(start_var == 'varID')) %>% rename(run_start = pos)
repairTracts <- left_join(repairTracts, simplified_mappingVar_calls_XIVL %>% select(varID, pos), join_by(end_var == 'varID')) %>% rename(run_end = pos)

# extend runs to ends of chromosome arm if they start at either of the 2 most telomeric variants
repairTracts %<>% mutate(run_start = case_when(
   start_var == 'var7190' | start_var == 'var7191' ~ 0,
   .default = run_start))

# add ORF coordinate information
repairTracts <- left_join(repairTracts, CSMorf_info %>% select(orf, start, stop), by = 'orf') %>% 
   rename(orf_start = start, orf_end = stop) %>%
   mutate(orf_mid = (orf_start + orf_end)/2, ## add column for midpoint of ORF
          short_iso = str_sub(iso_id, start = -2, end = -1)) ## extract shortened isolate ID

# write repair tract table to tsv file
write_tsv(repairTracts, file = repairTracts_outFile, col_names = TRUE)

##### MOSAIC TRACTS #####
# extract isolates with more than 2 runs of W303 on chrXIVL
chimeric_isos <- repairTracts %>% 
   filter(repair_type != 'deletion') %>% 
   group_by(isolate, call) %>% 
   summarise(n_runs = n()) %>% 
   filter(call == 'W303' & n_runs > 1) %>% 
   pull(isolate)

##### STRAIN REMOVAL #####
snp_removedORFs <- 'YNL079C' ## entire ORF removed from SNP analyses due to YKO strain aneuploidy
tractLength_removedORFS <- c('YNL055C', 'YNL073W', 'YNL133C', 'YNL170W', 'YNL173C', 'YNL220W', 'YNL243W', 'YNL280C', 'YNL284C')

lowHygR_removedORFs <- c('YNL055C', 'YNL170W', 'YNL173C', 'YNL220W', 'YNL243W', 'YNL280C', 'YNL284C') ## YKO strains with low reported HygR
YKO_removedORFs <- c(lowHygR_removedORFs, 'YNL268W', 'YNL073W', 'YNL079C', 'YNL081C', 'YNL133C', 'YNL252C', 'YNL291C') ## YKO strains with suspected issues

# add column to signify whether auxotrophy issues only arose from being LEU+
strainInfo %<>% mutate(onlyLeu = case_when(repair_type == 'S' ~ (rowSums(across(c(lost_in_streaking, hyg_r, g418, sd_cminus_lys, sd_cminus_ura), ~ . == 'yes')) == 0 & sd_cminus_leu == 'yes'),
                                           repair_type == 'R' ~ (rowSums(across(c(lost_in_streaking, g418, sd_cminus_lys, sd_cminus_ura), ~ . == 'yes')) == 0 & sd_cminus_leu == 'yes' & hyg_r == 'yes')))

# add filtering criteria columns
strainInfo <- strainInfo %>%
   rowwise() %>% 
   mutate(aux_removed = (orf == 'YNL268W' | (lost_in_streaking == 'yes' & !(is.na(lost_in_streaking))) | !(unexpected_auxotrophy == 'no' | is.na(unexpected_auxotrophy) | onlyLeu == TRUE)), ## column marking isolates removed due to unexpected auxotrophies
          preSeq_removed = (aux_removed | (hap1_pcr == 'yes' & !is.na(hap1_pcr)) | orf %in% lowHygR_removedORFs), ## column marking isolates with issues detected before sequencing
          YKO_removed = ((inconsistent_fingerprint == 'yes' & !is.na(inconsistent_fingerprint)) | (by_orf_variants == 'yes' & !is.na(by_orf_variants)) | (variants_not_spanning_orf == 'yes' & !is.na(variants_not_spanning_orf)) | orf %in% YKO_removedORFs), ## column marking isolates arising from problematic YKO strains
          CSM_removed = (((het_reads == 'yes' & !is.na(het_reads)) | (aneuploidy == 'yes' & !is.na(aneuploidy)) | ((recombination_elsewhere == 'yes' & !is.na(recombination_elsewhere)) & (hap1_pcr == 'no' | is.na(hap1_pcr)))) & YKO_removed == FALSE), ## column marking isolates with issues not arising from problematic YKO strains
          snp_removed = (aux_removed | 'yes' %in% c(het_reads, aneuploidy, recombination_elsewhere) | orf %in% snp_removedORFs), ## column marking isolates to be removed from SNP analyses
          tract_lengths_removed = (aux_removed | 'yes' %in% c(low_hyg_r, inconsistent_fingerprint, by_orf_variants, variants_not_spanning_orf, mosaic_tract, het_reads) | orf %in% tractLength_removedORFS)) ## column marking isolates to be removed from tract length analyses

# write out strain table with filtering columns to tsv file
write_tsv(strainInfo %>% mutate(mating_type = case_when(is.na(mating_type) ~ '-', .default = mating_type)), file = strainInfo_wFilters_outFile, col_names = TRUE)

#### recombination elsewhere plot ####
# create SFigure 3
recom_elsewhere_plot <- ggplot(data = recom_elsewhere) +
   geom_segment(aes(x = -500, xend = 12071856, y = orf_iso), color = 'black', linewidth = 3) + ## black outline around segment for each isolate
   geom_segment(aes(x = 0, xend = 12071356, y = orf_iso), color = 'gray85', linewidth = 2) + ## gray segment for stretches where allele could be either BY or W303
   geom_segment(aes(x = cum_start, xend = cum_end, y = orf_iso, color = call), linewidth = 2) + ## segments signifying BY or W303 runs for each isolate
   geom_vline(xintercept = cumsum(chrLen$len)[1:15], linewidth = 0.25, linetype = 'dotted') + ## vertical dashed lines showing chromosome boundaries
   annotate('segment', x = 459797 + sum(chrLen$len[1:11]), xend = 468931 + sum(chrLen$len[1:11]), y = 0, yend = length(unique(recom_elsewhere$orf_iso)), color = 'red') + ## vertical line for rDNA locus
   scale_y_discrete(limits = c(sort(unique(recom_elsewhere$orf_iso)))) + ## add space underneath run segments for rDNA label
   scale_x_continuous(breaks = chrLen$cum_mid[1:16], labels = chrLen$numerals[1:16], expand = c(0, 0)) + ## label x-axis with roman numerals and remove padding
   scale_color_manual(values = c('BY' = '#486590', 'W303' = '#50C878'), labels = c('BY', 'W303')) + ## set allele colors
   labs(x = 'Chromosome', y = 'Isolate', color = 'Allele') + ## assign axis and legend labels
   theme(axis.line.y = element_blank(), axis.line.x = element_line(), axis.ticks.y = element_blank(), legend.position = 'bottom') ## remove axis lines and ticks

# output plot to PDF
ggsave(file = recomElsewhere_plot_outFile, plot = recom_elsewhere_plot, width = 11.6, units = 'cm')

#### entire arm variant call plots ####
# create function for plots
armplot <- function(data, vars, strains, type) {
   ggplot(data = pivot_longer(data, cols = !(c(strain, orf)), names_to = 'variant', values_to = 'allele') %>%
             inner_join(., vars %>% select(varID, pos), join_by(variant == 'varID')) %>%
             inner_join(., strains %>% select(corrected_isoName, orf, repair_type, isoType, aux_removed, YKO_removed), join_by(strain == 'corrected_isoName')) %>%
             filter(repair_type == type & aux_removed == FALSE & YKO_removed == FALSE & variant != 'var7325')) +
      geom_point(aes(x = pos, y = isoType, color = allele, size = allele, alpha = allele)) +
      scale_color_manual(values = c('BY' = '#486590', 'W303' = '#50C878'), labels = c('BY', 'W303'), na.value = 'gray') + 
      scale_size_manual(values = c('BY' = 0.5, 'W303' = 2), labels = c('BY', 'W303'), na.value = 0.5) +
      scale_alpha_manual(values = c('BY' = 0.2, 'W303' = 1), labels = c('BY', 'W303'), na.value = 0.2) +
      xlab('Position on Chromosome XIV (bp)') +
      theme(axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(), 
            legend.position = 'none')
}

# create local plot
local_armplot <- armplot(data = transpose_mappingVar_calls_XIVL, vars = mappingVars_XIVL, strains = strainInfo, type = 'S') + 
   ylab('Local Isolates') + 
   theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

# create distal plot
distal_armplot <- armplot(data = transpose_mappingVar_calls_XIVL, vars = mappingVars_XIVL, strains = strainInfo, type = 'R') + 
   geom_rug(data = mappingVars_XIVL, aes(x = pos), color = 'black', linewidth = 0.1, inherit.aes = FALSE) + ## variant tract
   ylab('Distal Isolates')

# combine plots
armplots_combined <- local_armplot / distal_armplot + plot_layout(heights = c(1, 137/288))

# output plots to PDF
ggsave(file = armplots_outFile, plot = armplots_combined, width = 17.75,  units = 'cm')
