library(tidyverse)
library(magrittr)
library(readxl)
library(lme4)
library(lmerTest)
library(gridExtra)
library(ggpubr)
library(scales)
library(pheatmap)
library(openxlsx)

# LOAD RDA FILE WITH ALL MODEL DATA
# load("~/Desktop/Aim 1/Phenotyping/Data/modelData_v2_240627.Rda")
load("~/Desktop/CRISPAMap_downloads/sam_files/modelData_ALL_v5_240731.Rda")

# FILE NAME WITH GXE MODEL DATA
load("~/Desktop/CRISPAMap_downloads/sam_files/GxEmodelData_ALL_240730.Rda")

# FILE PATH CONTAINING ALL DATA FED TO MODEL (only used to get condition names)
fileName <- "~/Desktop/CRISPAMap_downloads/sam_files/YNL_Phen_Combined_ALL.xlsx"
condNames <- excel_sheets(fileName)

names(valDatAll) <- condNames
names(lmGxEDatAll) <- condNames
names(lmGxGDatAll) <- condNames
names(lmGxGxEDatAll) <- condNames

condNames <- condNames[!(condNames %in% c('YPD16h', 'YPD24h', 'Cu'))]

valDatAll <- valDatAll[condNames]
lmGxEDatAll <- lmGxEDatAll[condNames]
lmGxGDatAll <- lmGxGDatAll[condNames]
lmGxGxEDatAll <- lmGxGxEDatAll[condNames]

GxE_merged <- bind_rows(lmGxEDatAll[!is.na(lmGxEDatAll)], .id = 'condition') %>% filter(gene != 'wt') %>% select(condition, gene, baseline, geneEffect, GpValue, EpValue, GxEpValue)
rownames(GxE_merged) <- NULL

GxG_merged <- bind_rows(lmGxGDatAll[!is.na(lmGxGDatAll)], .id = 'condition')
GxGxE_merged <- bind_rows(lmGxGxEDatAll[!is.na(lmGxGxEDatAll)], .id = 'condition')

write.xlsx(list('GxE_model' = GxE_merged, 'GxG_model' = GxG_merged, 'GxGxE_model' = GxGxE_merged), file = '~/Desktop/CRI-SPA-Map/250820_tables/interaction_models.xlsx')

theme_set(theme_bw(base_size = 10, base_family = 'Arial') +
             theme(text = element_text(color = 'black'), axis.text = element_text(color = 'black', size = 8)))

# FILENAME TO SAVE VALIDATION BOX PLOTS
valBoxPlotFile <- "~/Desktop/CRI-SPA-Map/250820_figures/Val_Boxplot_240730.pdf"

#Box Plots (with GxGxE)
#Reworked pvals:
for (i in 1:length(condNames)) { #Iterate through conditions
   dat <- valDatAll[[i]] #Get current validation data
   dat <- dat %>% mutate(Significance2 = case_when( #Create significance column with different significance cases
      pValue < 0.05 ~ "Nominal",
      .default = "NS"
   )) 
   valDatAll[[i]] <- dat #Save data with new significance variable
}

plots <- list() #Initialize list to save plots
val_data <- lapply(1:length(condNames), function(i) { #Iterate through all conditions
   datVal <- valDatAll[[i]] #Get validation data
   datVal$gene <- as.character(datVal$gene) #Make sure gene names are strings
   datVal$gene[datVal$gene == "MKT1-30D"] <- "MKT1-30G" #Rename 30D -> 30G
   datVal$gene[datVal$gene == "SAL1-fs"] <- "SAL1-nfs" #Rename fs -> nfs
   datVal$gene[datVal$gene == "MKT1-30D+SAL1-fs"] <- "MKT1-30G+SAL1-nfs" #Rename 30G and nfs
   datVal$gene[datVal$gene == "wt"] <- "Eng-wt" #Rename wild type
   dat <- datVal #Reassign (obsolete)
   
   if (condNames[i] != "YPD40h") { #Check that current condition is not control condition (no GxE or GxGxE for control condition)
      GxEres <- lmGxEDatAll[[i]] #Get current condition's GxE model data
      GxEres$gene <- rownames(GxEres) #Make strain name variable
      GxEres$gene[GxEres$gene == "MKT1-30D"] <- "MKT1-30G" #Do the same strain renames for the model data
      GxEres$gene[GxEres$gene == "SAL1-fs"] <- "SAL1-nfs"
      GxEres$gene[GxEres$gene == "MKT1-30D+SAL1-fs"] <- "MKT1-30G+SAL1-nfs"
      dat <- left_join(dat,GxEres[,c("GxEpValue","gene")], by = "gene") #Add GxE p vals to the validation data
      
      GxGxEres <- lmGxGxEDatAll[[i]] #Get GxGxE model data
      dat$GxGxEpValue <- as.numeric(GxGxEres[length(GxGxEres)])
   } else { 
      GxEpValue <- 1
      GxGxEpValue <- 1}
   
   GxGres <- lmGxGDatAll[[i]] #Control condition (YPD 40h) does have GxG data
   dat$GxGpValue <- as.numeric(GxGres[length(GxGres)]) #Get the GxG p value and convert to scientific notation
   
   dat$Enviro <- condNames[i] #Make environment variable
   
   dat$gene <- factor(dat$gene, levels = c("Eng-wt","MKT1-30G","SAL1-nfs","MKT1-30G+SAL1-nfs","W303")) #Order strain names for plotting
   
   dat %<>% mutate(pValue_label = case_when(pValue < 0.0001 ~ as.character(formatC(pValue, format = 'e', digits = 2)), .default = as.character(signif(pValue, 2))),
                   GxEpValue_label = case_when(GxEpValue < 0.0001 ~ as.character(formatC(GxEpValue, format = 'e', digits = 2)), .default = as.character(signif(GxEpValue, 2))),
                   GxGpValue_label = case_when(GxGpValue < 0.0001 ~ as.character(formatC(GxGpValue, format = 'e', digits = 2)), .default = as.character(signif(GxGpValue, 2))),
                   GxGxEpValue_label = case_when(GxGxEpValue < 0.0001 ~ as.character(formatC(GxGxEpValue, format = 'e', digits = 2)), .default = as.character(signif(GxGxEpValue, 2))),
                   yMax = max(dat$colSizeCorrected), #Get max y value to display GxE p values
                   yMin = min(dat$colSizeCorrected),
                   Replicate = as.numeric(Replicate),
                   plate = as.numeric(plate),
                   IsolateNum = as.numeric(IsolateNum)) #Get min y value to display GxG and GxGxE p values
   return(dat)
}) %>% bind_rows()

val_data %<>% mutate(Enviro = factor(Enviro, levels = c('YPD40h', 'Caf', 'EtOH', 'LiAc', 'LowG', 'NaCl', 'SDS', 'Tun', 'YNB', 'YPD37C', 'YPDliq', 'YPG'), 
                                     labels = c('YPD', 'Caffeine (2 mg/mL)', 'Ethanol (7.5%)', 'Lithium Acetate (50 mM)', 'Low Glucose (0.25%)', 'Sodium Chloride (750 mM)', 'SDS (0.05%)', 'Tunicamycin (1.6 µM)', 'YNB', 'YPD 37°C', 'YPD Liquid', 'YPG')))

val_plots <- ggplot(val_data, aes(gene, colSizeCorrected)) +
   geom_hline(data = val_data %>% filter(gene == 'Eng-wt') %>% group_by(gene, Enviro) %>% summarise(median_growth = median(colSizeCorrected, na.rm = TRUE)), aes(yintercept = median_growth), lty = 1, lwd = 0.2) +
   geom_boxplot(aes(fill = Significance2), outlier.shape = NA, show.legend = TRUE, lwd = 0.1) +  
   geom_point(aes(group=Isolate), alpha=0.5, size = 0.2, position=position_dodge(width=0.6), show.legend = FALSE) +
   geom_line(aes(group=Isolate), alpha=0.4, position=position_dodge(width=0.6), lwd = 0.1, show.legend = FALSE) +
   geom_text(aes(x = 'Eng-wt', y = yMax + 0.1, label = ifelse(Enviro == 'YPD', 'G:', 'G:\nGxE:')), check_overlap = TRUE, vjust = 1, size = 6/.pt, family = 'Arial') +
   geom_text(data = val_data %>% filter(gene != 'Eng-wt'), aes(x = gene, y = yMax + 0.1, label = ifelse(Enviro == 'YPD', pValue_label, str_c(pValue_label, '\n', GxEpValue_label))), check_overlap = TRUE, vjust = 1, size = 6/.pt, family = 'Arial') +
   geom_text(aes(x = 3, y = yMin - 0.02, label = ifelse(Enviro == 'YPD', str_c('GxG: ', GxGpValue_label), str_c('GxG: ', GxGpValue_label, '\n', 'GxGxE: ', GxGxEpValue_label))), check_overlap = TRUE, vjust = 0, size = 6/.pt, family = 'Arial') +
   scale_fill_manual(breaks = c('Nominal', "NS"), labels = c('Nominal', "Not Significant"), values = c("Nominal"="#FFCC99", "NS"="#E0E0E0"), name = "G pValue Significance", drop = FALSE) + #Color code for significance
   scale_x_discrete(breaks = c('Eng-wt', 'MKT1-30G', 'SAL1-nfs', 'MKT1-30G+SAL1-nfs', 'W303'), labels = c('BY', 'MKT1-30G', 'SAL1-nfs', 'MKT1-30G\n+ SAL1-nfs', 'W303')) +
   facet_wrap(facets = vars(Enviro), ncol = 3, scales = 'free_y') +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = 'none') +
   ylab(expression(paste('Colony size [', log[2]*'(colony area (', mm^2, ')]'))) + xlab("Strain")

ggsave(valBoxPlotFile, val_plots, height = 27, width = 17.75, units = 'cm')

LiAc_YPD_data <- val_data %>% filter(as.character(Enviro) %in% c('Lithium Acetate (50 mM)', 'YPD', 'YPD Liquid'))
LiAc_plot <- ggplot(LiAc_YPD_data, aes(gene, colSizeCorrected)) +
   geom_hline(data = LiAc_YPD_data %>% filter(gene == 'Eng-wt') %>% group_by(gene, Enviro) %>% summarise(median_growth = median(colSizeCorrected, na.rm = TRUE)), aes(yintercept = median_growth), lty = 1, lwd = 0.2) +
   geom_boxplot(aes(fill = Significance2), outlier.shape = NA, show.legend = TRUE, lwd = 0.1) +  
   geom_point(aes(group=Isolate), alpha=0.5, size = 0.2, position=position_dodge(width=0.6), show.legend = FALSE) +
   geom_line(aes(group=Isolate), alpha=0.4, position=position_dodge(width=0.6), lwd = 0.1, show.legend = FALSE) +
   geom_text(aes(x = 'Eng-wt', y = yMax + 0.1, label = ifelse(Enviro == 'YPD', 'G:', 'G:\nGxE:')), check_overlap = TRUE, vjust = 1, size = 6/.pt, family = 'Arial') +
   geom_text(data = LiAc_YPD_data %>% filter(gene != 'Eng-wt'), aes(x = gene, y = yMax + 0.1, label = ifelse(Enviro == 'YPD', pValue_label, str_c(pValue_label, '\n', GxEpValue_label))), check_overlap = TRUE, vjust = 1, size = 6/.pt, family = 'Arial') +
   geom_text(aes(x = 3, y = yMin - 0.02, label = ifelse(Enviro == 'YPD', str_c('GxG: ', GxGpValue_label), str_c('GxG: ', GxGpValue_label, '\n', 'GxGxE: ', GxGxEpValue_label))), check_overlap = TRUE, vjust = 0, size = 6/.pt, family = 'Arial') +
   scale_fill_manual(breaks = c('Nominal', "NS"), labels = c('Nominal', "Not Significant"), values = c("Nominal"="#FFCC99", "NS"="#E0E0E0"), name = "G pValue Significance", drop = FALSE) + #Color code for significance
   scale_x_discrete(breaks = c('Eng-wt', 'MKT1-30G', 'SAL1-nfs', 'MKT1-30G+SAL1-nfs', 'W303'), labels = c('BY', 'MKT1-30G', 'SAL1-nfs', 'MKT1-30G\n+ SAL1-nfs', 'W303')) +
   facet_wrap(facets = vars(Enviro), ncol = 2, scales = 'free_y') +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = 'none') +
   ylab(expression(paste('Colony size [', log[2]*'(colony area (', mm^2, ')]'))) + xlab("Strain")
   
ggsave("~/Desktop/CRI-SPA-Map/250820_figures/LiAc_YPD40h_YPDliq.pdf", LiAc_plot, height = 12, width = 17.75, units = 'cm')

