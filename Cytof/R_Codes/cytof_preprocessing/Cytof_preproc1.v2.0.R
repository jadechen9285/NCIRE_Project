# ---------------UCSF_VA_NCIRE PROJECT: CyTOF single cell data analysis-----------
#   Jianhong Chen
#   07/30/2020
#  Pre-processing: normalization & debarcoding with CATALYST
#     
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
library(matrixStats) # matrix operation

# stat modeling 
library(lme4)
library(multcomp)

# bioconductor 
library(flowCore) # flow cytometry data wranggling & analysis
library(openCyto)
library(ggcyto) # viz for cyto data
library(limma) # linear model for microarray analysis
library(FlowSOM) # self-organizing map clustering for cytometry data
library(ConsensusClusterPlus) # cluster parameters optimization
library(Rphenograph)  
library(SingleCellExperiment) # for single cell data object
library(CATALYST) # wrapper packages in cytometry data analysis

# Cytof data analysis GUI
library(premessa)

# Viz.
library(RColorBrewer)
library(randomcoloR)
library(GGally) # ggplot pairs plot
library(ggrepel) # separate text labels on ggplot
library(pheatmap) # heatmap
library(ComplexHeatmap) # advanced version of heatmap
library(cowplot) 

# Dimensional Reduction
library(Rtsne) # t-SNE plot
library(uwot) # UMAP dimensional reduction

# Utilities
library(shiny)

##-- set ggplot global setting
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5, size = 25),
             axis.text = element_text(size = 15),
             axis.title = element_text(size = 20),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 15))
jc_color_cate = brewer.pal(n = 9, "Set1")[-3]
jc_color_cont = colorRampPalette(rev(brewer.pal(n=11, name = "RdYlBu")))(100)


# ----------------------------- Load Raw FCS files ---------------------------------
dir_raw = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/LT_Old10x_2018_08.04.2020/Raw_FCS/"
setwd(dir_raw)

fcs_files = list.files()[7:12]

# load raw fcs into single cell data object
sce = prepData(fcs_files)
table(sce$sample_id)
names(int_colData(sce))


# ---------------------------------- Normalization ------------------------------
norm = normCytof(sce, bead = 'dvs', k = 50, 
                 remove_beads = TRUE,
                 overwrite = FALSE)
norm_data = norm$data
# view bead channels
rownames(norm_data)[rowData(norm_data)$bead_ch]


# ---------------------------------- debarcoding --------------------------------
dir_meta = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/LT_Old10x_2018_08.04.2020/metadata/"
setwd(dir_meta)
debarcode_key = read_csv("10x_debarcoder.csv") 

debar_key = as.matrix(debarcode_key[, -1])
rownames(debar_key) = debarcode_key$sample

## -- assignment of preliminary IDs
debar = assignPrelim(norm_data, bc_key = debar_key,
                     assay = 'normexprs', verbose = TRUE)

# view barcode channels
rownames(debar)[rowData(debar)$is_bc]
# view number of events assigned to each barcode population
table(debar$bc_id)

## -- Estimation of separation cutoffs
# instead of one single cutoff, this will estimate sample-specific cutoff
debar2 = estCutoffs(debar)
# view separation cutoff estimates
metadata(debar2)$sep_cutoffs

## -- visualize the separation cutoffs yield
plotYields(debar2, which = 0)

## -- Applying de-convolution parameters
# use the estimated sep_cuttoff
debar3 = applyCutoffs(debar2, assay = 'normexprs', mhl_cutoff = 30,
                     sep_cutoffs = NULL) 
table(debar3$bc_id)

## exclude unassigned events
debar_final = debar3[, debar3$bc_id != 0]
table(debar_final$bc_id)


# --------------------- Write results into csv filess-------------------------

##!!! convert back to "flowset" with one frame per debarcoded ID
# fs = sce2fcs(debar_final, split_by = "bc_id", assay = 'normcounts')

# split check: number of cells per barcode ID
# equals number of cells in each 'flowFrame'
# c(fsApply(fs, nrow)) == table(debar_final$bc_id)

# extra expr data after debarcoding
norm_debarcode_expr = as.data.frame(t(assays(debar_final)$normexprs)) %>%
  mutate('sample_id' = debar_final$bc_id)

## record the number of events at each step
yield_size = data.frame(
  'step' = c('raw', 'normalization', 'debarcode'),
  'num_event' = c(ncol(sce), ncol(norm_data), ncol(debar_final)),
  'note' = c('original', 'remove beads', 'remove unassigned events'),
  'yield' = c(round(ncol(sce)/ncol(sce)*100, 2), 
                 round(ncol(norm_data)/ncol(sce)*100, 2),
                 round(ncol(debar_final)/ncol(sce)*100, 2))
  )


dir_output = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/COPD_BAL_10plex_set1_02.22.2018/processed/"
setwd(dir_output)

## write out fcs files: 
# get sample identifiers
# ids <- fsApply(fs, identifier)
# for (id in ids) {
#   ff <- fs[[id]]                     # subset 'flowFrame'
#   fn <- sprintf("fullPanel_%s.fcs", id) # specify output name that includes ID
#   fn <- file.path(dir_output, fn)         # construct output path
#   write.FCS(ff, fn)                  # write frame to FCS
# }

# write out csv files:
write_excel_csv(norm_debarcode_expr, 
                "fullPanel_norm_debarcode_asinh_expr_08.26.2020.csv")
write_excel_csv(yield_size, "fullPanel_yield_size_08.26.2020.csv")















