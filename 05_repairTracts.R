# R version 4.5.1

# load libraries
library(tidyverse) ## v 2.0.0
library(magrittr) ## v 2.0.3
library(patchwork) ## v 1.3.0
library(ggforce) ## v 0.4.2
library(extrafont) ## v 0.19

# user-input file paths and info
outDir <- '250820_tables/' ## directory for all analysis output files
figDir <- '250820_figures/' ## directory for all figure output files

repairTracts_summary_outFile <- paste0(outDir, 'repairTracts_summary.tsv')
mainRuns_plot_outFile  <- paste0(figDir, 'mainRuns.pdf')
allRuns_plot_outFile  <- paste0(figDir, 'allRuns.pdf')
tractLengths_plot_outFile  <- paste0(figDir, 'tractLengths_waterfall.pdf')

options(scipen = 999) ## no scientific notation
theme_set(theme_bw(base_size = 10, base_family = 'Arial') + ## set base ggplot fonts
             theme(text = element_text(color = 'black'), axis.text = element_text(color = 'black', size = 8)))

# make output directories
dir.create(outDir, showWarnings = FALSE)
dir.create(figDir, showWarnings = FALSE)

# load in previously-created files
repairTracts <- read_tsv(paste0(outDir, 'repairTracts.tsv'), col_names = TRUE, show_col_types = FALSE)
strainInfo <- read_tsv(paste0(outDir, 'strainInfo_wFilters.tsv'), col_names = TRUE, show_col_types = FALSE)
mappingVars_XIVL <- read_tsv(paste0(outDir, 'mappingVars_XIVL.tsv'), col_names = TRUE, show_col_types = FALSE)

##### TRACT LENGTHS #####
# add column with shorter isolate identifier
repairTracts %<>% mutate(isoType = str_c(str_split_i(isolate, '_', 1), '_', short_iso)) 

# reformat and subset tract table for plotting
repairTracts_plot <- inner_join(repairTracts %>% filter(orf_start >= min(mappingVars_XIVL$pos) & orf_end <= max(mappingVars_XIVL$pos)), ## extract runs from isolates with targeted ORFs between the mapping variants on chrXIV-L
                                strainInfo %>% filter(tract_lengths_removed == FALSE) %>% select(isoType, tract_lengths_removed), by = 'isoType') %>% ## join runs with strain info table to filter appropriate isolates
   filter(!(start_var == 'CEN14' & end_var == 'CEN14')) ## remove any runs that begin or end at the centromere

# summarize closest calls to the ORF for each isolate
repairTracts_plot_summary <- repairTracts_plot %>%
   group_by(isolate, orf, repair_type, short_iso, orf_start, orf_end) %>% ## group data by isolate
   summarise(runs = n(), ## determine how many runs each isolate has
             BY_tel = ifelse(length(run_end[call == 'BY']) > 0, min(run_end[call == 'BY']), NA), ## if isolate has more than 1 BY run, extract the minimum value for run end position to find the telomeric BY call closest to ORF
             BY_cen = ifelse(length(run_start[call == 'BY']) > 0, max(run_start[call == 'BY']), NA), ## if isolate has more than 1 BY run, extract the maximum value for run start position to find the centromeric BY call closest to ORF
             W303_tel = ifelse(length(run_start[call == 'W303']) == 1, run_start[call == 'W303'], NA), ## for isolates with only 1 W303 run, extract run start position to find the telomeric W303 call farthest from ORF
             W303_cen = ifelse(length(run_end[call == 'W303']) == 1, run_end[call == 'W303'], NA)) %>% ## for isolates with only 1 W303 run, extract run end position to find the centromeric W303 call farthest from ORF
   ungroup() %>% ## ungroup data
   mutate(dist_BY_5 = case_when( ## determine BY call closest to 5' end of ORF
      (BY_tel <= orf_start & str_detect(orf, 'W')) ~ orf_start - BY_tel,
      (BY_cen >= orf_end & str_detect(orf, 'C')) ~ BY_cen - orf_end,
      .default = NA), ## if closest 5' BY call is within ORF, report NA
      dist_BY_3 = case_when( ## determine BY call closest to 3' end of ORF
         (BY_cen >= orf_end & str_detect(orf, 'W')) ~ BY_cen - orf_end,
         (BY_tel <= orf_start & str_detect(orf, 'C')) ~ orf_start - BY_tel,
         .default = NA), ## if closest 3' BY call is within ORF, report NA
      dist_W303_5 = case_when( ## determine W303 call farthest from 5' end of ORF
         (W303_tel <= orf_start & str_detect(orf, 'W') & runs == 3) ~ orf_start - W303_tel,
         (W303_cen >= orf_end & str_detect(orf, 'C')) ~ W303_cen - orf_end,
         .default = NA), ## if farthest 5' W303 call is within ORF, report NA
      dist_W303_3 = case_when( ## determine W303 call farthest from 3' end of ORF
         (W303_cen >= orf_end & str_detect(orf, 'W')) ~ W303_cen - orf_end,
         (W303_tel <= orf_start & str_detect(orf, 'C') & runs == 3) ~ orf_start - W303_tel,
         .default = NA)) ## if farthest 3' W303 call is within ORF, report NA

# output table of closest calls
write_tsv(repairTracts_plot_summary, file = repairTracts_summary_outFile, col_names = TRUE)

# pivot table for plotting
repairTracts_plot_summary %<>% pivot_longer(cols = starts_with('dist_'), names_to = c('allele', 'bound'), names_prefix = 'dist_', names_sep = '_', values_to = 'dist')

# determine median distances for each allele distance
median(repairTracts_plot_summary %>% filter(allele == 'W303') %>% pull(dist), na.rm = T)
median(repairTracts_plot_summary %>% filter(allele == 'BY') %>% pull(dist), na.rm = T)

#### repair tracts plot - main ####
# extract data for main figure plot
main_run_isos <- repairTracts %>% 
   filter(orf %in% c('YNL303W', 'YNL203C', 'YNL097C', 'YNL077W')) %>% ## use data from 4 ORFs along the chromosome arm for main figure
   mutate(orf = factor(orf, levels = c('YNL303W', 'YNL203C', 'YNL097C', 'YNL077W')), ## change ORF variable to factor for ordering
          short_iso = factor(short_iso, levels = c('R8', 'R7', 'R6', 'R5', 'R4', 'R3', 'R2', 'R1', 'S7', 'S6', 'S5', 'S4', 'S3', 'S2', 'S1'))) ## order isolates on plot by short names

# create plot
main_runs <- ggplot(data = main_run_isos) +
   geom_segment(aes(x = -500, xend = 629258, y = short_iso), color = 'black', linewidth = 3) + ## black chromosome arm outline
   geom_segment(aes(x = 0, xend = 628758, y = short_iso), color = 'gray85', linewidth = 2.5) + ## gray chromosome background for regions that could be BY or W303
   geom_segment(aes(x = run_start, xend = run_end, y = short_iso, color = call), linewidth = 2.5) + ## segments signifying BY or W303 runs
   geom_vline(aes(xintercept = orf_mid), color = 'black', linetype = 'dashed', linewidth = 0.2) + ## vertical dashed line at midpoint of targeted ORF
   geom_rug(data = mappingVars_XIVL %<>% mutate(orf = 'YNL303W'), aes(x = pos), color = 'black', linewidth = 0.1, inherit.aes = FALSE) + ## variant tract
   geom_rug(data = data.frame(pos = 45308, orf = 'YNL303W'), aes(x = pos), color = '#ffe200', linewidth = 0.1, inherit.aes = FALSE) + ## Hyg resistance cassette insertion location
   scale_x_continuous(expand = c(0, 0)) + ## remove padding
   scale_color_manual(values = c('BY' = '#486590', 'W303' = '#50C878')) + ## assign allele colors
   labs(x = 'Chromosome XIV position (bp)', y = 'Isolate', color = 'Allele') + ## assign axis and legend labels
   facet_grid(rows = vars(orf), scales = 'free_y', space = 'free_y') + ## facet by ORF
   theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), legend.position = 'none') ## remove axis lines, ticks, grid, and legend

# output plot to PDF
ggsave(file = mainRuns_plot_outFile, plot = main_runs, width = 8.5, height = 9, units = 'cm')

# #### repair tracts plot - all ####
# all_runs <- ggplot(repairTracts %>% 
#                       mutate(short_iso = factor(short_iso, levels = c('R8', 'R7', 'R6', 'R5', 'R4', 'R3', 'R2', 'R1', 'S7', 'S6', 'S5', 'S4', 'S3', 'S2', 'S1')))) +
#    geom_segment(aes(x = -500, xend = 629258, y = short_iso), color = 'black', linewidth = 5) + ## black chromosome arm outline
#    geom_segment(aes(x = 0, xend = 628758, y = short_iso), color = 'gray85', linewidth = 4) + ## gray chromosome background for regions that could be BY or W303
#    geom_segment(aes(x = run_start, xend = run_end, y = short_iso, color = call), linewidth = 4) + ## segments signifying BY or W303 runs
#    geom_vline(aes(xintercept = orf_mid), color = 'black', linetype = 'dashed', linewidth = 0.5) + ## vertical dashed line at midpoint of targeted ORF
#    annotate('rug', x = mappingVars_XIVL$pos, color = 'black', linewidth = 0.1) + ## variant tract
#    annotate('rug', x = 45308, color = '#ffe200', linewidth = 0.1) + ## Hyg resistance cassette insertion location
#    scale_x_continuous(expand = c(0, 0)) + ## remove padding
#    scale_color_manual(values = c('BY' = '#486590', 'W303' = '#50C878')) + ## assign allele colors
#    labs(x = 'Chromosome XIV position (bp)', y = 'Isolate', color = 'Allele') + ## assign axis and legend labels
#    facet_grid_paginate(facets = vars(orf), ncol = 1, nrow = 4, page = 4, scales = 'free_y', space = 'free') + ## facet by ORF on multiple pages
#    theme(axis.line.y = element_blank(), axis.line.x = element_line(), axis.ticks.y = element_blank(), panel.grid = element_blank()) ## remove axis lines, ticks, and grid
# plot_pages <- n_pages(all_runs) ## output how many pages will be plotted
# 
# # output to multi-page PDF
# pdf(file = allRuns_plot_outFile)
# for (i in 1:plot_pages){
#    print(all_runs + facet_grid_paginate(facets = vars(orf), ncol = 1, nrow = 4, page = i, scales = 'free_y', space = 'free')) + plot_layout(ncol = 1, heights = c(1,15))
# }
# dev.off()

#### tract length distribution plot ####
# create 5' repair tract length waterfall plot
dist5_plot <- ggplot(data = repairTracts_plot_summary %>% filter(bound == 5), aes(x = -dist, color = allele)) + ## filter for 5' data
   geom_vline(data = repairTracts_plot_summary %>% filter(bound == 5 & allele == 'BY'), aes(xintercept = -median(dist, na.rm = TRUE), color = allele), linetype = 'dashed') + ## graph median BY distance as dashed vertical line
   geom_vline(data = repairTracts_plot_summary %>% filter(bound == 5 & allele == 'W303'), aes(xintercept = -median(dist, na.rm = TRUE), color = allele), linetype = 'dashed') + ## graph median W303 distance as dashed vertical line
   stat_ecdf(linewidth = 1) + ## plot cumulative distribution
   labs(x = "5' Distance to ORF (bp)", y = 'Proportion of Isolates') + ## assign axis labels
   coord_cartesian(xlim = c(-15000,0)) + ## set axis limits
   scale_color_manual(values = c('BY' = '#486590', 'W303' = '#50C878')) + ## assign allele colors
   theme(legend.position = 'none') ## remove legend

# create 3' repair tract plot
dist3_plot <- ggplot(data = repairTracts_plot_summary %>% filter(bound == 3), aes(x = dist, color = allele)) + ## filter for 3' data
   geom_vline(data = repairTracts_plot_summary %>% filter(bound == 3 & allele == 'BY'), aes(xintercept = median(dist, na.rm = TRUE), color = allele), linetype = 'dashed') + ## graph median BY distance as dashed vertical line
   geom_vline(data = repairTracts_plot_summary %>% filter(bound == 3 & allele == 'W303'), aes(xintercept = median(dist, na.rm = TRUE), color = allele), linetype = 'dashed') + ## graph median W303 distance as dashed verical line
   stat_ecdf(aes(y = 1 - after_stat(y)), linewidth = 1) + ## plot mirror of cumulative distribution plot
   labs(x = "3' Distance to ORF (bp)") + ## assign axis label
   coord_cartesian(xlim = c(0,15000)) + ## set axis limits
   scale_color_manual(values = c('BY' = '#486590', 'W303' = '#50C878')) + ## assign allele colors
   theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = 'none') ## remove axis title, text, and legend

# position plots
repair_lengths_plot <- dist5_plot + dist3_plot

# output plot to PDF
ggsave(file = tractLengths_plot_outFile, plot = repair_lengths_plot, width = 7.5, height = 4.5, units = 'cm')
