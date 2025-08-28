# R version 4.5.1

# load libraries
library(tidyverse) ## v 2.0.0
library(magrittr) ## v 2.0.3
library(readxl) ## v 1.4.3
library(janitor) ## v 2.2.0
library(extrafont) ## v 0.19, must run font_import() followed by loadfonts() if first time using

# user-input file paths and info
strainInfo_inFile <- 'data/CRI_SPA_MAP_INFO.xlsx' ## table with details on each isolate
calledVars_inFile <- 'data/varCalls.txt' ## .txt version of VCF
CSM_orfs_inFile <- 'data/CSM_ORFs_info.tsv' ## table of ORFS targeted in this study
chrXIVLorf_inFile <- 'data/S288C_chrXIVL_ORFs.tsv' ## file containing info on all ORFs found on chrXIV-L in S288C reference

sample_pref <- '/home/albertf/lawle167/005_CSM_250529/alignments/' ## text to remove before sample names in VCF headers
sample_suff <- '_L001_rmdup.bam' ## text to remove after sample names in VCF headers

outDir <- '250820_tables/' ## directory for all analysis output files
figDir <- '250820_figures/' ## directory for all figure output files
mappingVars_outFile <- paste0(outDir, 'mappingVars.tsv') ## output file for high-quality BY/W303 variants
mappingVars_XIVL_outFile <- paste0(outDir, 'mappingVars_XIVL.tsv') ## output file for high-quality BY/W303 variants on chrXIV-L
mitoCalls_outFile <- paste0(outDir, 'CSM_mitoCalls.tsv') ## output file for mitochondrial variant calls
mitoCalls_plot_outFile <- paste0(figDir, 'mitoCalls.pdf') ## output file for SFigure 6
vars_elsewhere_outFile <- paste0(outDir, 'varsElsewhere.tsv') ## output file for background variants
orfInfo_outFile <- paste0(outDir, 'CSM_orfInfo.tsv') ## output file for linking variants to ORFs
chrXIVLcalls_outFile <- paste0(outDir, 'CSM_chrXIVLcalls.tsv') ## output file for chrXIV-L variant calls
allCalls_outFile <- paste0(outDir, 'CSM_allCalls.tsv') ## output file for genome-wide variant calls

options(scipen = 999) ## no scientific notation
theme_set(theme_bw(base_size = 10, base_family = 'Arial') + ## set base ggplot fonts
             theme(text = element_text(color = 'black'), 
                   axis.text = element_text(color = 'black', size = 7), 
                   strip.text = element_text(color = 'black', size = 7)))

# make output directories
dir.create(outDir, showWarnings = FALSE)
dir.create(figDir, showWarnings = FALSE)

##### MAPPING VARIANTS #####
# load input files
calledVars <- read_table(calledVars_inFile, comment = '##', show_col_types = FALSE) 
strainInfo <- bind_rows(read_excel(path = strainInfo_inFile, sheet = 'HygS_ALLStrains', na = c('NA', '')), read_excel(path = strainInfo_inFile, sheet = 'HygR_ALLStrains', na = c('NA', ''))) %>% ## merge HygR and HygS sheets
   clean_names() %>% ## make column names code-friendly
   mutate(corrected_isoName = case_when(!(is.na(renamed_sequencing_id)) ~ renamed_sequencing_id, .default = sequencing_illumina_id_orf_cri_spa_well_id_isolate_seq_plate_well), ## add column with modified isolate name (if present) or original isolate name (if no modification)
          orf = case_when(!is.na(corrected_isoName) ~ str_split_i(corrected_isoName, '_', 1), .default = str_split_i(orf_isolate, '_', 1)), ## extract ORF from isolate name
          isoNum = case_when(!is.na(corrected_isoName) ~ str_sub(str_split_i(corrected_isoName, '_', 2), -1, -1), .default = str_split_i(orf_isolate, '_', 2)), ## extract isolate number from isolate name
          isoType = str_c(orf, '_', repair_type, isoNum)) ## designate isolates with ORF, repair type, and isolate number

# extract table for renaming isolates
name_correction <- strainInfo %>% 
   filter(!is.na(renamed_sequencing_id)) %>% ## extract only isolates that were renamed
   select(sequencing_illumina_id_orf_cri_spa_well_id_isolate_seq_plate_well, renamed_sequencing_id) ## extract only columns with old name and new name

# reformat variant call table
colnames(calledVars) <- c('chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', colnames(calledVars)[10:ncol(calledVars)]) ## rename info columns
calledVars %<>% rename_with(.cols = contains(sample_pref), ~ str_split_i(str_sub(.x, (str_length(sample_pref) + 1), (-(str_length(sample_suff) + 1))), '_S', 1)) ## shorten sample headers
calledVars %<>% rename_with(.cols = name_correction$sequencing_illumina_id_orf_cri_spa_well_id_isolate_seq_plate_well, .fn = ~ name_correction$renamed_sequencing_id) ## fix sample names
calledVars %<>% select(chr, pos, ref, alt, qual, info, format, sort(colnames(.)[10:ncol(.)])) ## extract desired columns and sort samples by name

# add a locus column containing positional and allelic information
calledVars %<>% mutate(locus = paste(chr, pos, ref, alt, sep = '-'))

# extract BY and W303 control isolate calls into separate tables
BYvars <- calledVars %>% select(chr, pos, ref, alt, locus, qual, info, format, starts_with('BYXW'), YFA1752_NA_6F11, YFA1753_NA_6G11, YFA1754_NA_6H11)
W303vars <- calledVars %>% select(chr, pos, ref, alt, locus, qual, info, format, YFA1747_NA_6A12, YFA1748_NA_6B12, YFA1749_NA_6E11, YFA1756_NA_6D11)
tagBYvars <- BYvars %>% select(!(starts_with('BYXW')))

# define function for simplifying variant calls
simplify_calls <- function(data, id, min_qual = 0) {
   simplified_calls <- data %>%
      mutate(across(all_of(grep(paste(id, collapse = '|'), names(data), value = T)), ~ ifelse( ## apply commands to any columns containing any of the strings in the ID variable list
         (as.integer(str_split_i(.x, ':', 5)) >= min_qual & startsWith(.x, '0:')), 'ref', ifelse( ## change all cells that meet GQ threshold and signify reference call to 'ref'
            (as.integer(str_split_i(.x, ':', 5)) >= min_qual & startsWith(.x, '1:')), 'alt', ifelse( ## change all cells that meet GQ threshold and signify alternate call to 'alt'
               startsWith(.x, '.:'), 'no_call', 'low_qual'))))) ## change all cells that meet GQ threshold and signify no call to 'no_call' and remaining to 'low_qual'
}

# simplify variant calls for control isolates
BYvars_simp <- simplify_calls(data = BYvars, id = c('BYXW', 'YFA'))
W303vars_simp <- simplify_calls(data = W303vars, id = 'YFA')
tagBYvars_simp <- simplify_calls(data = BYvars, id = 'YFA')

# count alleles for each variant
BYvars_simp$ref_calls <- apply(BYvars_simp %>% select(starts_with('BYXW'), starts_with('YFA')), 1, function(x) { length(x[x == 'ref']) })
BYvars_simp$alt_calls <- apply(BYvars_simp %>% select(starts_with('BYXW'), starts_with('YFA')), 1, function(x) { length(x[x == 'alt']) })

W303vars_simp$ref_calls <- apply(W303vars_simp %>% select(starts_with('YFA')), 1, function(x) { length(x[x == 'ref']) })
W303vars_simp$alt_calls <- apply(W303vars_simp %>% select(starts_with('YFA')), 1, function(x) { length(x[x == 'alt']) })

tagBYvars_simp$ref_calls <- apply(tagBYvars_simp %>% select(starts_with('YFA')), 1, function(x) { length(x[x == 'ref']) })
tagBYvars_simp$alt_calls <- apply(tagBYvars_simp %>% select(starts_with('YFA')), 1, function(x) { length(x[x == 'alt']) })

# filter isolates based on proportion of alt and high-quality calls
BYvars_simp <- BYvars_simp %>% 
   filter(chr != 'sacCer_chrM') %>% ## remove mitochondrial variants - some controls may be shuffled due to mating
   mutate(prop_called = (ref_calls + alt_calls)/length(colnames(BYvars_simp)[startsWith(colnames(BYvars_simp), 'BYXW') | startsWith(colnames(BYvars_simp), 'YFA')]), ## calculate proportion of isolates with calls that meet GQ threshold
          prop_alt = alt_calls/length(colnames(BYvars_simp)[startsWith(colnames(BYvars_simp), 'BYXW') | startsWith(colnames(BYvars_simp), 'YFA')])) %>% ## calculate proportion of isolates with alt calls that meet GQ threshold
   filter(prop_alt >= 0.75) ## filter for variants that have alt calls in at least 75 percent of control isolates

W303vars_simp <- W303vars_simp %>% 
   mutate(prop_called = (ref_calls + alt_calls)/length(colnames(W303vars_simp)[startsWith(colnames(W303vars_simp), 'YFA')]), ## calculate proportion of isolates with calls that meet GQ threshold
          prop_alt = alt_calls/length(colnames(W303vars_simp)[startsWith(colnames(W303vars_simp), 'YFA')])) %>% ## calculate proportion of isolates with alt calls that meet GQ threshold
   filter(prop_alt >= 0.75) ## filter for variants that have alt calls in at least 75 percent of control isolates

tagBYvars_simp <- tagBYvars_simp %>% 
   filter(chr == 'sacCer_chrM') %>%
   mutate(prop_called = (ref_calls + alt_calls)/length(colnames(tagBYvars_simp)[startsWith(colnames(tagBYvars_simp), 'YFA')]), ## calculate proportion of isolates with calls that meet GQ threshold
          prop_alt = alt_calls/length(colnames(tagBYvars_simp)[startsWith(colnames(tagBYvars_simp), 'YFA')])) %>% ## calculate proportion of isolates with alt calls that meet GQ threshold
   filter(prop_alt >= 0.75) ## filter for variants that have alt calls in at least 75 percent of control isolates

# extract mitochondrial variants from W303 control isolates
mitoVars <- W303vars_simp %>% filter(chr == 'sacCer_chrM') 

# remove mitochondrial variants from nuclear mapping variant list
W303vars_simp %<>% filter(chr != 'sacCer_chrM') 

# load S288C ORF information
chrXIVLorf_info <- read_tsv(chrXIVLorf_inFile, show_col_types = FALSE) %>% ## all S288C ORFs on chrXIVL from SGD/alliancemine
   rename(orf = SequenceFeature.secondaryIdentifier, ## rename columns
          gene_symbol = SequenceFeature.symbol, 
          gene_name = SequenceFeature.name, 
          sgd_id = SequenceFeature.primaryIdentifier, 
          start = SequenceFeature.chromosomeLocation.start, 
          stop = SequenceFeature.chromosomeLocation.end, 
          chr = SequenceFeature.chromosomeLocation.locatedOn.primaryIdentifier) %>%
   select(orf, gene_symbol, gene_name, chr, start, stop, sgd_id) ## reorder columns

# load list of ORFs targeted by CRI-SPA-Map
CSMorf_info <- read_tsv(CSM_orfs_inFile, show_col_types = FALSE)

# identify variants that are not shared between BY and W303
BYonly <- BYvars_simp %>% filter(locus %in% setdiff(BYvars_simp$locus, W303vars_simp$locus)) %>%
   select(chr, pos, ref, alt, locus) %>%
   mutate(alt_strain = 'BY', ref_strain = 'W303') ## add columns specifying which strain is the reference
W303only <- W303vars_simp %>% filter(locus %in% setdiff(W303vars_simp$locus, BYvars_simp$locus)) %>%
   select(chr, pos, ref, alt, locus) %>%
   mutate(alt_strain = 'W303', ref_strain = 'BY') ## add columns specifying which strain is the reference

# combine all distinguishing variants into single table
mappingVars <- bind_rows(BYonly, W303only) %>% 
   mutate(chr_num = case_when(chr == 'sacCer_chrI' ~ 1,
                              chr == 'sacCer_chrII' ~ 2,
                              chr == 'sacCer_chrIII' ~ 3,
                              chr == 'sacCer_chrIV' ~ 4,
                              chr == 'sacCer_chrV' ~ 5,
                              chr == 'sacCer_chrVI' ~ 6,
                              chr == 'sacCer_chrVII' ~ 7,
                              chr == 'sacCer_chrVIII' ~ 8,
                              chr == 'sacCer_chrIX' ~ 9,
                              chr == 'sacCer_chrX' ~ 10,
                              chr == 'sacCer_chrXI' ~ 11,
                              chr == 'sacCer_chrXII' ~ 12,
                              chr == 'sacCer_chrXIII' ~ 13,
                              chr == 'sacCer_chrXIV' ~ 14,
                              chr == 'sacCer_chrXV' ~ 15,
                              chr == 'sacCer_chrXVI' ~ 16)) %>% ## add numeric chromosome number for sorting
   arrange(chr_num, pos) %>% ## sort all variants by chromosome number and position
   mutate(varID = paste0('var', seq(1, (nrow(BYonly) + nrow(W303only))))) %>% ## add unique variant identifier
   select(varID, chr, pos, ref, alt, locus, ref_strain, alt_strain) ## extract only desired columns

# extract variants on chrXIVL
mappingVars_XIVL <- mappingVars %>% filter(chr == 'sacCer_chrXIV' & pos <= 628758)

mappingVars_XIVL <- mappingVars_XIVL %>% 
   rowwise() %>% 
   mutate(orf = case_when(nrow(chrXIVLorf_info[which(chrXIVLorf_info$start <= pos & chrXIVLorf_info$stop >= pos), 'orf']) > 1 ~ as.character(chrXIVLorf_info[which(chrXIVLorf_info$start <= pos & chrXIVLorf_info$stop >= pos & !(is.na(chrXIVLorf_info$gene_name))), 'orf']), ## if variant in more than one ORF, list non-dubious ORF
                          nrow(chrXIVLorf_info[which(chrXIVLorf_info$start <= pos & chrXIVLorf_info$stop >= pos), 'orf']) == 1 ~ as.character(chrXIVLorf_info[which(chrXIVLorf_info$start <= pos & chrXIVLorf_info$stop >= pos), 'orf']), ## if variant within an ORF, report which ORF
                          .default = NA)) ## if variant not within ORF bound, return NA

mappingVars_XIVL %<>% mutate(targeted_orf = case_when(orf %in% CSMorf_info$orf ~ TRUE, ## add column signifying if ORF variant was targeted by CRI-SPA-Map
                                                      is.na(orf) ~ NA, ## NA if no ORF is targeted
                                                      .default = FALSE)) ## false if ORF is not targeted in this study

# output tables to tsv files
write_tsv(mappingVars, mappingVars_outFile, col_names = T) 
write_tsv(mappingVars_XIVL, mappingVars_XIVL_outFile, col_names = T) 

##### MITOCHONDRIAL VARIANTS #####
# filter called variants for known mitochondrial variants - all variants in tagged BY isolates are shared with W303 isolates
calledVars_mito <- calledVars %>% filter(locus %in% setdiff(mitoVars$locus, tagBYvars_simp$locus)) ## identify mitochondrial variants between W303 and unshuffled BY strains
calledVars_mito %<>% mutate(varID = paste0('varM', seq(1, (nrow(calledVars_mito)))), alt_strain = 'W303', ref_strain = 'BY') ## add mitochondrial variant ID and designate reference and alternate alleles

# simplify mitochondrial calls
simplified_mito_calls <- simplify_calls(data = calledVars_mito, id = c('BYXW', 'YFA', 'YNL'), min_qual = 30)
simplified_mito_calls %<>% mutate(across(all_of(c(starts_with('BYXW'), starts_with('YFA'), starts_with('YNL'))), ~ ifelse((.x == 'ref'), 'BY', ifelse((.x == 'alt'), 'W303', .x))))

# define function to transpose and reformat simplified call table
transpose_calls <- function(simp_calls) {
   transpose_calls = as.data.frame(cbind(colnames(simp_calls)[2:ncol(simp_calls)], t(simp_calls[,-1]))) ## transpose calls and add column for strain names
   rownames(transpose_calls) <- NULL ## remove row names
   colnames(transpose_calls) <- c('strain', simp_calls$varID) ## name columns with variant IDs
   transpose_calls %<>% mutate(orf = str_split_i(strain, '_', 1)) ## indicate targeted ORF of each strain
   transpose_calls %<>% arrange(strain) ## sort table by strain name
   return(transpose_calls)
}

# transpose mitochondrial calls
transpose_mito_calls <- simplified_mito_calls %>% select(varID, starts_with('BYXW'), starts_with('YFA'), starts_with('YNL')) %>% transpose_calls()

# write table to tsv files
write_tsv(simplified_mito_calls %>% select(varID, chr, pos, ref, alt, ref_strain, alt_strain, qual, info, format, starts_with('BYXW'), starts_with('YFA'), starts_with('YNL')), mitoCalls_outFile, col_names = T) 

##### VARIANTS ELSEWHERE #####
# simplify all VCF calls using the locus column
simp_all_calls <- simplify_calls(data = calledVars, id = c('BYXW', 'YFA', 'YNL'), min_qual = 30) %>% select(locus, starts_with('BYXW'), starts_with('YFA'), starts_with('YNL'))

# find the proportion of calls across CRI-SPA-Map isolates that were 'alt', low quality, or had no genotype call
simp_all_calls$alt_calls <- apply(simp_all_calls %>% select(!(starts_with('BYXW') | starts_with('YFA') | matches('*_[A-Z][0-9]+D[0-9]+_*'))), 1, function(x) { length(x[x == 'alt']) })
simp_all_calls$prop_no_or_low <- apply(simp_all_calls %>% select(!(starts_with('BYXW') | starts_with('YFA') | matches('*_[A-Z][0-9]+D[0-9]+_*'))), 1, function(x) { (length(x[x == 'no_call']) + length(x[x == 'low_qual']))/(length(x) - 1) })
simp_all_calls %<>% filter(alt_calls >= 3 & (prop_no_or_low < 0.1 | is.na(prop_no_or_low)))

# filter table to only include the CRI-SPA-Map isolates
csm_all_calls <- simp_all_calls %>% select(!(starts_with('BYXW') | starts_with('YFA') | matches('*_[A-Z][0-9]+D[0-9]+_*')))
var_orf_isos <- apply(csm_all_calls, 1, function(x) {
   var = x[['locus']] ## extract variant information
   alt_calls = as.numeric(x['alt_calls']) ## extract how many 'alt' calls were reported in CRI-SPA-Map isolates at this variant
   isos = colnames(csm_all_calls)[which(x == 'alt')] ## extract which isolates had 'alt' calls at this variant
   orfs = str_split_i(isos, '_', 1) ## extract targeted ORF from isolate names
   if(length(which(as.numeric(table(orfs)) >= 4)) > 0) { ## if at least four isolates from the same ORF have 'alt' calls at this variant...
      orfs_4 = paste(names(table(orfs)[which(as.numeric(table(orfs)) >= 4)]), collapse = ',') ## calculate for how many ORFs this is the case
      n_orfs_4 = length(names(table(orfs)[which(as.numeric(table(orfs)) >= 4)])) } else { ## list which ORFs meet this criteria
         orfs_4 = NA
         n_orfs_4 = 0
      }
   if(length(orfs[which(as.numeric(table(orfs)) == 3)]) > 0) { ## if three isolates from the same ORF have 'alt' calls at this variant...
      orfs_3 = paste(names(table(orfs)[which(as.numeric(table(orfs)) == 3)]), collapse = ',') ## calculate for how many ORFs this is the case
      n_orfs_3 = length(names(table(orfs)[which(as.numeric(table(orfs)) == 3)])) } else { ## list which ORFs meet this criteria
         orfs_3 = NA
         n_orfs_3 = 0
      }
   return(data.frame(var, alt_calls, n_orfs_4, orfs_4, n_orfs_3, orfs_3)) ## output row with variant and the ORFs that have at least three isolates with 'alt' calls
}) %>% bind_rows() %>% ## merge rows into table
   mutate(chr = str_split_i(var, '-', 1), ## extract chromosome from variant locus
          pos = as.numeric(str_split_i(var, '-', 2)), ## extract chromosome position from variant locus
          ref = str_split_i(var, '-', 3), ## extract reference allele from variant locus
          alt = str_split_i(var, '-', 4),) ## extract alternate allele from variant locus
rownames(var_orf_isos) <- NULL ## remove row names for table

var_orf_isos <- left_join(var_orf_isos %>% filter(n_orfs_4 > 0 | n_orfs_3 > 0), csm_all_calls %>% select(!alt_calls), join_by(var == locus)) ## extract only variants that have at least one ORF with at least three isolates with 'alt' calls and merge with CRI-SPA-Map isolate variant call table by the locus column
var_orf_isos <- left_join(var_orf_isos, mappingVars %>% select(varID, locus, ref_strain, alt_strain), join_by(var == locus)) ## if the variant is a known BY/W303 variant, add the variant information from mapping variants list
var_orf_isos %<>% select(varID, chr, pos, ref, alt, ref_strain, alt_strain, alt_calls, n_orfs_4, orfs_4, n_orfs_3, orfs_3, starts_with('YNL')) %>% ## extract desired columns
   arrange(desc(alt_calls)) ## sort rows by decreasing amount of 'alt' calls

# output table without mitochondrial variants and known chrXIV-L variants to csv file
write_tsv(var_orf_isos %>% filter(chr != 'sacCer_chrM' & !(!is.na(varID) & chr == 'sacCer_chrXIV')), file = vars_elsewhere_outFile, col_names = TRUE)

##### ORF VARIANTS #####
# filter called variants to only include mapping variants
calledVars_map <- inner_join(mappingVars, calledVars %>% select(locus, qual, info, format, starts_with('BYXW'), starts_with('YFA'), starts_with('YNL')), by = 'locus')

# add chrXIV centromere BY call to all isolates
calledVars_map %<>% add_row(varID = 'CEN14', chr = 'sacCer_chrXIV', pos = 628758, ref_strain = 'BY', .before = min(which((calledVars_map$chr == 'sacCer_chrXIV' & calledVars_map$pos > 628758) | calledVars_map$chr == 'sacCer_chrXV')))

# add chrXIV centromere to mapping variant list
mappingVars %<>% add_row(varID = 'CEN14', chr = 'sacCer_chrXIV', pos = 628758, .before = min(which((mappingVars$chr == 'sacCer_chrXIV' & mappingVars$pos > 628758) | mappingVars$chr == 'sacCer_chrXV')))

# extract chrXIV left arm variant calls at the mapping variants
calledVars_XIVL <- calledVars_map %>% filter(chr == 'sacCer_chrXIV' & pos <= 628758)

# determine how many and which mapping variants are within each ORF
CSMorf_info <- CSMorf_info %>% 
   rowwise() %>% 
   mutate(n_mappingVars = length(which(mappingVars_XIVL$pos >= start & mappingVars_XIVL$pos <= stop))) %>% ## calculate how many mapping variants are within the start and stop coordinates of each targeted ORF
   mutate(mappingVar_IDs = case_when(n_mappingVars > 0 ~ paste(mappingVars_XIVL$varID[which(mappingVars_XIVL$pos >= start & mappingVars_XIVL$pos <= stop)], collapse = ','))) ## extract the variant IDs of the variants in each ORF

# output ORF/variant table to tsv
write_tsv(CSMorf_info, orfInfo_outFile, col_names = T) 

# designate genotype as BY or W303 allele for each strain at each called mapping variant site
simplified_mappingVar_calls <- simplify_calls(data = calledVars_map, id = c('BYXW', 'YFA', 'YNL'), min_qual = 30)
simplified_mappingVar_calls %<>% mutate(across(all_of(starts_with('YNL')), ~ ifelse((.x == 'ref'), ref_strain, 
                                                                                    ifelse((.x == 'alt'), alt_strain, .x))))

# designate chrXIV-L as BY in all CRI-SPA-Map isolates
simplified_mappingVar_calls %<>% mutate(across(all_of(starts_with('YNL')), ~ ifelse((varID == 'CEN14'), 'BY', .x)))

# determine how many isolates have a W303 call at each variant
simplified_mappingVar_calls$W303_calls <- apply(simplified_mappingVar_calls %>% select(starts_with('YNL')), 1, function(x) { length(x[x == 'W303']) })

# extract simplified chrXIVL calls
simplified_mappingVar_calls_XIVL <- simplified_mappingVar_calls %>% 
   filter(chr == 'sacCer_chrXIV' & pos < 628758)

# output simplified calls to tsv files
write_tsv(simplified_mappingVar_calls %>% filter(!(chr == 'sacCer_chrXIV' & pos == 628758)) %>% select(varID, chr, pos, ref, alt, ref_strain, alt_strain, qual, info, format, starts_with('BYXW'), starts_with('YFA'), starts_with('YNL')), allCalls_outFile, col_names = T) 
write_tsv(simplified_mappingVar_calls_XIVL %>% filter(!(chr == 'sacCer_chrXIV' & pos == 628758)) %>% select(varID, chr, pos, ref, alt, ref_strain, alt_strain, qual, info, format, starts_with('BYXW'), starts_with('YFA'), starts_with('YNL')), chrXIVLcalls_outFile, col_names = T)

#### SFIGURE 6 ####
mitoCalls_plot <- ggplot(data = pivot_longer(transpose_mito_calls, cols = !(c(strain, orf)), names_to = 'variant', values_to = 'allele') %>% ## pivot mitochondrial calls
                            mutate(variant = factor(variant, levels = paste0('varM', seq_len(length.out = nrow(calledVars_mito)))), ## add variant ID
                                   short_iso = str_sub(str_split_i(strain, '_', 2), -2, -1), ## extract short isolate identifier from isolate name
                                   orf = case_when(orf == 'BYXW' ~ 'Wildtype', .default = orf)) %>% ## extract ORF from isolate name
                            filter(!str_detect(short_iso, 'D') & short_iso != 'NA') %>% ## remove any sequenced deletion strain calls
                            mutate(short_iso = factor(short_iso, levels = c('R8', 'R7', 'R6', 'R5', 'R4', 'R3', 'R2', 'R1', 'S7', 'S6', 'S5', 'S4', 'S3', 'S2', 'S1')))) + ## set order of isolates
   geom_tile(aes(x = variant, y = short_iso, fill = allele)) + ## make a tile plot
   scale_fill_manual(values = c('BY' = '#486590', 'W303' = '#50C878'), labels = c('BY', 'W303'), na.value = 'white') + ## manually assign allele colors
   facet_wrap(vars(orf), scales = 'free_y', strip.position = 'right') + ## facet plots by targeted ORF
   labs(x = 'Variant', y = 'Isolate', fill = 'Allele') + ## set labels
   guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5)) + ## adjust size of legend
   theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = 'bottom', axis.text.y = element_text(margin = margin(r = 0))) # make theme modifications

# output plot to PDF
ggsave(filename = mitoCalls_plot_outFile, plot = mitoCalls_plot, width = 17.75, units = 'cm')