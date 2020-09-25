# ---------------UCSF_VA_NCIRE PROJECT: CyTOF single cell data analysis-----------
#   Jianhong Chen
#   07/06/2020
#
#   Pipeline II-2:
#  High Dim. Analyzing preprocessed with Scaffold Map
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear workspace
rm(list = ls()) # delete all variables in the environment
# Clear console
cat("\014")  # equals to 'control + L'

# data wranggling 
library(tidyverse)
library(reshape2)
library(readxl)

# High Dim Analysis
library(flowCore) # flow cytometry data wranggling & analysis
library(ggcyto) # viz for cyto data
library(limma) # linear model for microarray analysis
library(FlowSOM) # self-organizing map clustering for cytometry data
library(ConsensusClusterPlus) # cluster parameters optimization
library(CATALYST) # wrapper packages in cytometry data analysis
library(Rphenograph) # Phenograph Clustering
## Scafold Map
library(grappolo)
library(vite)
library(panorama)

# visualization
library(ggrepel) # separate text labels on ggplot
library(pheatmap)
library(ComplexHeatmap) # better version of pheatmap

# Color
library(RColorBrewer)
library(randomcoloR)
library(circlize) # color for heatmap

# Dimensional Reduction
library(Rtsne) # t-SNE plot
library(uwot) # UMAP dimensional reduction

# Utilities
library(shiny)

# ------------------------------ Load FCS Data files ----------------------------------

dir_meta = '/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/metadata/LT_old_03.26.2018/'
setwd(dir_meta)

md = read_excel("LT_old_metadata.xlsx") %>%
  transform("condition" = factor(condition),
            "smoker" = factor(smoker),
            "sex" = factor(sex),
            "patient_id" = factor(patient_id),
            "gold" = factor(gold)) %>%
  mutate("sample_id" = str_c("LT", patient_id),# create sample_ids
         "RV_TLC_cate" = as.factor(case_when(
           RV_TLC < mean(RV_TLC, na.rm = TRUE) ~ "Normal_Air",
           RV_TLC >= mean(RV_TLC, na.rm = TRUE) ~ "Abnormal_Air",
           is.na(RV_TLC) ~ "Missing")) )

rownames(md) = md$sample_id

panel = read_excel("panel.xlsx") 

dir_raw = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/gated/LT_old_03.26.2018/"
setwd(dir_raw)
fcs_files = list.files()

fcs_raw = flowCore::read.flowSet(fcs_files, transformation = FALSE,
                                 truncate_max_range = FALSE)

panel_fcs = pData(parameters(fcs_raw[[1]]))
( lineage_markers = panel$marker[panel$lineage == 1] ) 
lineage_markers = str_replace_all(lineage_markers, "_", "-")

# -- Clustering Method III: Scaffold Map
#(1) cluster each file individually

# cluster_fcs_files(fcs_files, num.cores = 2, col.names = lineage_markers ,
#                   num.clusters = 200, asinh.cofactor = 5)


#(2) creating an unsupervised graph
dir_vite = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/scaffold_res/full_text/"
setwd(dir_vite)

input.files = list.files()
G = vite::get_unsupervised_graph_from_files(input.files,
                                             col.names = lineage_markers,
                                             filtering.threshold = 15)










