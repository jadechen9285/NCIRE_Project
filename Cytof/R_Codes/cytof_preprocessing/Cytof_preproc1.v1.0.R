# ---------------UCSF_VA_NCIRE PROJECT: CyTOF single cell data analysis-----------
#   Jianhong Chen
#   05/18/2020
# single cell cytometry data analysis
#     
#  Practice on old Cytof data with Lashimi 
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

# matrix operation
library(matrixStats)

# stat modeling 
library(lme4)
library(multcomp)

# bioconductor 
library(flowCore) # flow cytometry data wranggling & analysis
library(ggcyto) # viz for cyto data
library(limma) # linear model for microarray analysis
library(FlowSOM) # self-organizing map clustering for cytometry data
library(ConsensusClusterPlus) # cluster parameters optimization
library(Rphenograph)  
library(CATALYST) # wrapper packages in cytometry data analysis

# Cytof data analysis GUI
library(premessa)

# visualization
library(ggrepel) # separate text labels on ggplot
library(pheatmap)
library(RColorBrewer)
library(randomcoloR)
library(Rtsne) # t-SNE plot
library(uwot) # UMAP dimensional reduction

# Utilities
library(shiny)


# ------------------- Preprocessing Raw CyTOF fcs data with Premessa ------------------

# -- I. Make sure all fcs files have same panels and same number of columns
paneleditor_GUI() # check alignment of number of columns for raw fcs files

# -- II. Beads Normalization on concatenated files per experiment(barcode)
dir_raw = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/New_Data/LT_Cytof_Blue_07.29.2020/Raw_FCS/"
setwd(dir_raw)
normalizer_GUI()


# -- read raw fcs data
fcs_files = list.files()
fcs_raw20x = flowCore::read.flowSet(fcs_files, transformation = FALSE,
                                   truncate_max_range = FALSE)
pData(parameters(fcs_raw20x[[1]]))


# # -- III. Concatenated normalized files into one set of experiment(per barcode)
dir_norm = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/normalized/LT_set3_07.13.2020/"
setwd(dir_norm)

set.list = list.files()
# set2.list = list.files()[c(11:12, 20:22)]
# set3.list = list.files()[13:19]

concatenate_fcs_files(set.list, output.file = "LT_set3_normed_concatenated_07.13.2020.fcs")
# concatenate_fcs_files(set2.list, output.file = "CO3PD_set2_BALnormed_concatenated.fcs")
# concatenate_fcs_files(set3.list, output.file = "CO3PD_set3_BALnormed_concatenated.fcs")

# -- IV. Debarcoded normalized fcs files into each sample

# debarcoder_GUI()






