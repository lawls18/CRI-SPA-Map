library(tidyverse)
library(magrittr)
library(readxl)
library(openxlsx)

# LOAD VALIDATION PHENOTYPING DATA
valData <- read_tsv('tables/validation_phenotype_allConditions.tsv', col_names = TRUE, show_col_types = FALSE) %>%
   rename(condition = 'Enviro')

# LOAD MODEL DATA
modelData <- read_tsv('tables/interaction_models.tsv', col_names = TRUE, show_col_types = FALSE) 

# FILE PATH CONTAINING ALL DATA FED TO MODEL (only used to get condition names)
fileName <- "tables/processed_pixl_data.xlsx"
condNames <- excel_sheets(fileName)

theme_set(theme_bw(base_size = 10, base_family = 'Arial') +
             theme(text = element_text(color = 'black'), axis.text = element_text(color = 'black', size = 8)))

# FILENAME TO SAVE VALIDATION BOX PLOTS
valBoxPlotFile <- "figures/validation_plots.pdf"
mainValBoxPlotFile <- "figures/main_validation_plots.pdf"

#Box Plots (with GxGxE)
valData %<>% mutate(gene_cond = paste(gene, condition, sep = '_'))
modelData %<>% mutate(gene_cond = paste(gene, condition, sep = '_'))

valData <- left_join(valData, modelData %>% select(!c(gene, condition)), by = 'gene_cond')

valData %<>% mutate(gene = case_when(
   gene == "MKT1-30D" ~ "MKT1-30G",
   gene == "SAL1-fs" ~ "SAL1-nfs",
   gene == "MKT1-30D+SAL1-fs" ~ "MKT1-30G+SAL1-nfs", .default = gene),
   gene_significance = case_when(pValueGene < 0.05 ~ 'significant', .default = 'NS'))

valData %<>% mutate(gene = factor(gene, levels = c("wt","MKT1-30G","SAL1-nfs","MKT1-30G+SAL1-nfs","W303")),
                    pValueGene_label = case_when(pValueGene < 0.0001 ~ as.character(formatC(pValueGene, format = 'e', digits = 1)), .default = as.character(signif(pValueGene, 2))),
                    pValueGxE_label = case_when(pValueGxE < 0.0001 ~ as.character(formatC(pValueGxE, format = 'e', digits = 1)), .default = as.character(signif(pValueGxE, 2))),
                    pValueGxG_label = case_when(pValueGxG < 0.0001 ~ as.character(formatC(pValueGxG, format = 'e', digits = 1)), .default = as.character(signif(pValueGxG, 2))),
                    pValueGxGxE_label = case_when(pValueGxGxE < 0.0001 ~ as.character(formatC(pValueGxGxE, format = 'e', digits = 1)), .default = as.character(signif(pValueGxGxE, 2))))

yLimits <- valData %>% group_by(condition) %>% 
   summarise(yMax = max(colSizeCorrected),
             yMin = min(colSizeCorrected)) %>%
   ungroup()

valData <- left_join(valData, yLimits, by = 'condition')

valData_main <- valData %>% filter(condition %in% c('YPD', 'LiAc', 'YPDliq'))
valData_main %<>% mutate(condition = factor(condition, levels = c('YPD', 'LiAc', 'YPDliq'), 
                                       labels = c('YPD', 'Lithium Acetate (50 mM)', 'YPD Liquid')))


valData %<>% mutate(condition = factor(condition, levels = c('YPD', 'Caf', 'EtOH', 'LiAc', 'LowG', 'NaCl', 'SDS', 'Tun', 'YNB', 'YPD37C', 'YPG', 'YPDliq'), 
                                             labels = c('YPD', 'Caffeine (2 mg/mL)', 'Ethanol (7.5%)', 'Lithium Acetate (50 mM)', 'Low Glucose (0.25%)', 'Sodium Chloride (750 mM)', 'SDS (0.05%)', 'Tunicamycin (1.6 µM)', 'YNB', 'YPD 37°C', 'YPG', 'YPD Liquid')))

val_plots_main <- ggplot(valData_main, aes(gene, colSizeCorrected)) +
   geom_hline(data = valData_main %>% filter(gene == 'wt') %>% group_by(gene, condition) %>% summarise(median_growth = median(colSizeCorrected, na.rm = TRUE)), aes(yintercept = median_growth), lty = 1, lwd = 0.2) +
   geom_boxplot(aes(fill = gene_significance), outlier.shape = NA, show.legend = TRUE, lwd = 0.1) +  
   geom_point(aes(group=Isolate), alpha=0.5, size = 0.2, position=position_dodge(width=0.6), show.legend = FALSE) +
   geom_line(aes(group=Isolate), alpha=0.4, position=position_dodge(width=0.6), lwd = 0.1, show.legend = FALSE) +
   geom_text(aes(x = 'wt', y = yMax + 0.1, label = ifelse(condition %in% c('YPD', 'YPD Liquid'), 'G:', 'G:\nGxE:')), check_overlap = TRUE, hjust = 0, vjust = 1, size = 6/.pt, family = 'Arial') +
   geom_text(data = valData_main %>% filter(gene != 'wt'), aes(x = gene, y = yMax + 0.1, label = ifelse(condition %in% c('YPD', 'YPD Liquid'), pValueGene_label, str_c(pValueGene_label, '\n', pValueGxE_label))), check_overlap = TRUE, hjust = 0, vjust = 1, size = 6/.pt, family = 'Arial') +
   geom_text(aes(x = 2.5, y = yMin - 0.02, label = ifelse(condition %in% c('YPD', 'YPD Liquid'), str_c('GxG: ', pValueGxG_label), str_c('GxG: ', pValueGxG_label, '\n', 'GxGxE: ', pValueGxGxE_label))), check_overlap = TRUE, hjust = 0, vjust = 0, size = 6/.pt, family = 'Arial') +
   scale_fill_manual(breaks = c('significant', "NS"), labels = c('Nominal', "Not Significant"), values = c("significant"="#FFCC99", "NS"="#E0E0E0"), name = "G pValue Significance", drop = FALSE) + #Color code for significance
   scale_x_discrete(breaks = c('wt', 'MKT1-30G', 'SAL1-nfs', 'MKT1-30G+SAL1-nfs', 'W303'), labels = c('BY', 'MKT1-30G', 'SAL1-nfs', 'MKT1-30G\n+ SAL1-nfs', 'W303')) +
   facet_wrap(facets = vars(condition), ncol = 1, scales = 'free_y') +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = 'none') +
   ylab(expression(paste('Colony size [', log[2]*'(colony area (', mm^2, ')]'))) + xlab("Strain")

ggsave(mainValBoxPlotFile, val_plots_main, height = 16.5, width = 9, units = 'cm')

val_plots <- ggplot(valData, aes(gene, colSizeCorrected)) +
   geom_hline(data = valData %>% filter(gene == 'wt') %>% group_by(gene, condition) %>% summarise(median_growth = median(colSizeCorrected, na.rm = TRUE)), aes(yintercept = median_growth), lty = 1, lwd = 0.2) +
   geom_boxplot(aes(fill = gene_significance), outlier.shape = NA, show.legend = TRUE, lwd = 0.1) +  
   geom_point(aes(group=Isolate), alpha=0.5, size = 0.2, position=position_dodge(width=0.6), show.legend = FALSE) +
   geom_line(aes(group=Isolate), alpha=0.4, position=position_dodge(width=0.6), lwd = 0.1, show.legend = FALSE) +
   geom_text(aes(x = 'wt', y = yMax + 0.1, label = ifelse(condition %in% c('YPD', 'YPD Liquid'), 'G:', 'G:\nGxE:')), check_overlap = TRUE, hjust = 0, vjust = 1, position = position_nudge(x = -0.4), size = 6/.pt, family = 'Arial') +
   geom_text(data = valData %>% filter(gene != 'wt'), aes(x = gene, y = yMax + 0.1, label = ifelse(condition %in% c('YPD', 'YPD Liquid'), pValueGene_label, str_c(pValueGene_label, '\n', pValueGxE_label))), check_overlap = TRUE, hjust = 0, vjust = 1, position = position_nudge(x = -0.4), size = 6/.pt, family = 'Arial') +
   geom_text(aes(x = 3, y = yMin - 0.02, label = ifelse(condition %in% c('YPD', 'YPD Liquid'), str_c('GxG: ', pValueGxG_label), str_c('GxG: ', pValueGxG_label, '\n', 'GxGxE: ', pValueGxGxE_label))), check_overlap = TRUE, hjust = 0, vjust = 0, position = position_nudge(x = -0.5), size = 6/.pt, family = 'Arial') +
   scale_fill_manual(breaks = c('significant', "NS"), labels = c('Nominal', "Not Significant"), values = c("significant"="#FFCC99", "NS"="#E0E0E0"), name = "G pValue Significance", drop = FALSE) + #Color code for significance
   scale_x_discrete(breaks = c('wt', 'MKT1-30G', 'SAL1-nfs', 'MKT1-30G+SAL1-nfs', 'W303'), labels = c('BY', 'MKT1-30G', 'SAL1-nfs', 'MKT1-30G\n+ SAL1-nfs', 'W303')) +
   facet_wrap(facets = vars(condition), ncol = 3, scales = 'free_y') +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = 'none') +
   ylab(expression(paste('Colony size [', log[2]*'(colony area (', mm^2, ')]'))) + xlab("Strain")

ggsave(valBoxPlotFile, val_plots, height = 20, width = 17.75, units = 'cm')
