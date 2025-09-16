# R version 4.5.1

# load libraries
library(tidyverse) ## v 2.0.0
library(magrittr) ## v 2.0.3
library(extrafont)

# user-input file paths and info
aneuploid_inDir <- 'data/aneuploid_coverage/'
coverage_inDir <- 'data/coverage/'
chrLen_inFile <- 'data/sacCer_chrLen.csv'

outDir <- 'tables/' ## directory for all analysis output files
figDir <- 'figures/' ## directory for all figure output files

aneuploid_outFile <- paste0(figDir, 'aneuploid_plot.pdf')
chromDepth_outFile <- paste0(outDir, 'chromDepth.tsv')

options(scipen = 999)
theme_set(theme_bw(base_size = 10, base_family = 'Arial') +
             theme(text = element_text(color = 'black'), axis.text = element_text(color = 'black', size = 8)))

# make output directory
dir.create(figDir, showWarnings = FALSE)

##### ANEUPLOIDY COVERAGE #####
# load in chromosome lengths
chr_len <- read_csv(chrLen_inFile) %>% 
   filter(chr != 'sacCer_chrM') %>% ## remove mitochondrial chromosome
   mutate(cum_len = c(0, cumsum(len)[1:15]), ## add cumulative length column
          midpoint = (len/2) + cum_len, ## add chromosome midpoint column
          numeral = str_split_i(chr, 'chr', 2)) ## extract roman numeral from chromosome name

# define function for determining mean coverages for windows of 5000 bp
determine_mean_covs <- function(file) {
   cov_table <- data.frame(read_table(paste0(aneuploid_inDir, file), col_names = FALSE, show_col_types = FALSE)) 
   colnames(cov_table) <- c('chr', 'pos', 'cov')
   cov_table <- cov_table %>% filter(chr != 'sacCer_chrM') %>% mutate(region = ceiling(pos/5000))
   total_mean_cov <- cov_table %>% pull(cov) %>% mean()
   summ_table <- cov_table %>%
      group_by(chr, region) %>%
      summarise(pos = mean(pos), mean_cov = mean(cov)) %>% 
      mutate(total_mean_cov = total_mean_cov,
             isolate = str_split_i(file, '_covAll.', 1))
   return(summ_table)
}

# load in all coverage data
covFiles <- list.files(aneuploid_inDir)
covDFs <- lapply(covFiles, determine_mean_covs) %>% bind_rows()

# join coverage data with chromosome length data
covDFs <- left_join(covDFs, chr_len, by = 'chr')

covDFs %<>% mutate(cum_pos = pos + cum_len,
                   norm_cov = mean_cov/total_mean_cov,
                   orf = str_split_i(isolate, '_', 1),
                   iso_id = str_split_i(isolate, '_', 2),
                   isoType = str_c(orf, '_', str_sub(iso_id, -2, -1)),
                   aneuploid_chr = case_when(orf %in% c('YNL079C', 'YNL246W') & chr == 'sacCer_chrIX' ~ 'aneuploid',
                                             orf == 'YNL076W' & chr == 'sacCer_chrVIII' ~ 'aneuploid',
                                             orf == 'YNL233W' & chr == 'sacCer_chrI' ~ 'aneuploid', .default = 'euploid'))

covPlot <- ggplot() +
   geom_point(data = covDFs %>% filter(aneuploid_chr == 'aneuploid'), aes(x = cum_pos, y = norm_cov), size = 0.001, color = '#CE6583') + 
   geom_point(data = covDFs %>% filter(aneuploid_chr == 'euploid'), aes(x = cum_pos, y = norm_cov), size = 0.001, color = 'gray') + 
   geom_hline(yintercept = 1, color = 'black', linewidth = 0.25, linetype = 'dashed') +
   geom_hline(yintercept = 2, color = 'black', linewidth = 0.25, linetype = 'dashed') +
   geom_vline(xintercept = chr_len$cum_len[2:16], linewidth = 0.5) +
   ylim(0,3) +
   scale_x_continuous(breaks = chr_len$midpoint[1:16], labels = chr_len$numeral[1:16], expand = c(0.01,0.01)) +
   labs(x = 'Chromosome', y = 'Relative Read Depth') + 
   facet_wrap(~isoType, ncol = 1, nrow = 24, strip.position = 'right') +
   theme(strip.text.y = element_text(angle = 0),
         strip.background = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         legend.position = 'none')

ggsave(aneuploid_outFile, plot = covPlot, height = 10, width = 17.75, units = 'cm')

##### CHROMOSOME COVERAGE #####
# load in files
chrom_covFiles <- list.files(path = coverage_inDir)
chrom_dp <- lapply(chrom_covFiles, function(x) {
   read_table(paste0(coverage_inDir, x), col_names = c('chr', 'start', 'chr_len', 'reads', 'cov_bases', 'coverage', 'mean_dp', 'mean_baseq', 'mean_mapq'), show_col_types = FALSE) %>%
      select(chr, mean_dp) %>% ## extract only necessary columns
      mutate(isolate = str_split_i(x, '_L001_', 1)) %>% ## extract isolate name from file name
      pivot_wider(id_cols = isolate, names_from = chr, values_from = mean_dp) %>% ## create single row for each isolate with mean chromosome depths for the columns
      select(!('#rname'))
}) %>%
   bind_rows() %>% ## combine all isolates into single table
   mutate_at(c(2:ncol(.)), as.numeric) ## convert mean depths to numeric values

chrom_dp <- chrom_dp %>%
   rowwise() %>%
   mutate(mean_nuc_dp = mean(c_across(starts_with('sacCer') & !all_of('sacCer_chrM'))), ## calculate mean depth of nuclear chromosomes for each isolate
          median_nuc_dp = median(c_across(starts_with('sacCer') & !all_of('sacCer_chrM')))) ## calculate median depth of nuclear chromosomes for each isolate

# output table of read depths
write_tsv(chrom_dp, file = chromDepth_outFile, col_names = TRUE)
