# ---------------UCSF_VA_NCIRE PROJECT: CyTOF single cell data analysis-----------
#   Jianhong Chen
#   07/24/2020
# CyTOF fcs files pre-processing II: gating
#     
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear workspace
rm(list = ls()) # delete all variables in the environment
# Clear console
cat("\014")  # equals to 'control + L'


# data wrangling 
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
## Scaffold Map
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


##-- set ggplot global setting
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5, size = 20),
             axis.text = element_text(size = 15),
             axis.title = element_text(size = 15))
jc_color_cate = brewer.pal(n = 9, "Set1")[-3]
jc_color_cont = colorRampPalette(rev(brewer.pal(n=11, name = "RdYlBu")))(100)


# --------------------------------- Load Raw Data---------------------------------

# -- read fcs raw (CD45 gated) data
dir_raw = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/normalized/BAL_02.06.2020/"
setwd(dir_raw)

fcs_files = list.files()
fcs_raw = flowCore::read.flowSet(fcs_files, transformation = FALSE,
                                 truncate_max_range = FALSE)

# extract fcs files columnames & labels
panel_fcs = pData(parameters(fcs_raw[[1]]))
panel_fcs$desc = str_replace_all(panel_fcs$desc, "-", "_")

# -- load metadata for the study
dir_meta = '/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/metadata/BAL_02.06.2020/'
setwd(dir_meta)

# md info: sex: 0 = M, 1 = F; ht: unit in cm; wt: unit in kg
md = read_excel("metadata.xlsx") 

panel = read_excel("panel.xlsx") 
panel$target= str_replace_all(panel$target, "-", "_")
panel$marker = str_replace_all(panel$marker, "-", "_")


# define target markers to study 
# NOTE: no functional markers in this study
( lineage_markers = panel$marker[panel$lineage == 1] )
# spot check
if (all(lineage_markers %in% panel_fcs$desc)) {
  TRUE
} else { # if not, check each indiviual marker 
  sapply(lineage_markers, FUN = function(x){x %in% panel_fcs$desc})
}

## if markers are in different readable format
# once the spot check passed, renamed the markers into more readable version
lineage_markers = map_chr(lineage_markers, .f = function(x){
  str_c(unlist(str_split(x, "_"))[-1], collapse = "_")
})


# --------------------------------- Pre-processing ---------------------------------
# generate sample IDs corresponding to each cell in the 'exprs' matrix
sample_ids = rep(md$sample_id, fsApply(fcs_raw, nrow))

# select required gating channels
gating_markers = c('191Ir_DNA_Intercalator_1', '193Ir_DNA_Intercalator_2', 
                   '198Pt_Live_dead' )

# spot check
all(gating_markers %in% panel_fcs$desc)

## if markers are in different readable format
# once the spot check passed, renamed the markers into more readable version
gating_markers = map_chr(gating_markers, .f = function(x){
  str_c(unlist(str_split(x, "_"))[-1], collapse = "_")
})

# once the spot check passed, renamed the panel markers into more readable version
panel_fcs$desc = map_chr(panel_fcs$desc, .f = function(x){
  str_c(unlist(str_split(x, "_"))[-1], collapse = "_")
})

# Asinh transformed and extracted interested channels
fcs = fsApply(fcs_raw, function(x, cofactor = 5){
  colnames(x) = panel_fcs$desc
  expr = exprs(x)
  expr = asinh(expr[, c(lineage_markers, gating_markers)]/cofactor)
  exprs(x) = expr
  x
})
expr_raw = fsApply(fcs, exprs)
rownames(expr_raw) = sample_ids

## -- **** in order to speed up the analysis of the large raw dataset   ****
# viz only 10% of the original data
#  however the whole analysis still utilize the whole set(mainly in a matrix)

# calculate the 10% sample size per fcs files
p = round(nrow(expr_raw)*0.1 / length(unique(sample_ids)))

## downsample original data
# step 2: create full indices by sample
inds = split(1:length(sample_ids), sample_ids)

# step 3: size of the down-sampling
dsamp_ncells = pmin(table(sample_ids), p) # pmin = parallel min

# step 4: randomly pick indices base on the chosen size
set.seed(78)
dsamp_inds = lapply(names(inds), function(i){
  s = sample(inds[[i]], dsamp_ncells[i], replace = FALSE)
}) %>%
  unlist()

raw_viz = as.data.frame(expr_raw[dsamp_inds, ]) %>%
  mutate("sample_id" = sample_ids[dsamp_inds] )

# ---------------------------------  Gating ---------------------------------------
## -- step 1: filter out live & singlet cells
gating_quantile = colQuantiles(expr_raw[, gating_markers],
                               probs = seq(0, 1, by = 0.05))
LL ='85%'
UL = '100%'
X = "DNA_Intercalator_1"
Y = "DNA_Intercalator_2"

# viz the two DNA-inter channels
## Note: annotate() is much faster than geom_rect()
p =
ggplot(raw_viz,
       aes(x = DNA_Intercalator_1, y = DNA_Intercalator_2)) +
  geom_hex(bins = 128) +
  scale_fill_gradientn("Count", colors = jc_color_cont)

p +
annotate(geom="rect",
         xmin = gating_quantile[X, LL],
         xmax = gating_quantile[X, UL],
         ymin = gating_quantile[Y, LL],
         ymax = gating_quantile[Y, UL],
         colour = "red", fill = NA)

expr_cell = expr_raw %>%
  base::subset(expr_raw[, X] > gating_quantile[X, LL] &
                 expr_raw[, X] < gating_quantile[X, UL] &
                 expr_raw[, Y] > gating_quantile[Y, LL] &
                 expr_raw[, Y] < gating_quantile[Y, UL])

# update the viz dataframe as well
cell_viz = raw_viz %>%
  dplyr::filter(
      DNA_Intercalator_1 > gating_quantile[X, LL] &
        DNA_Intercalator_1 < gating_quantile[X, UL] &
        DNA_Intercalator_2 > gating_quantile[Y, LL] &
        DNA_Intercalator_2 < gating_quantile[Y, UL])


## step 2 gate out dead cells
gating_quantile = colQuantiles(expr_cell[, gating_markers],
                               probs = seq(0, 1, by = 0.05))
LL ='0%'
UL = '80%'
X = "Live_dead"
Y = "DNA_Intercalator_2"

# viz live-dead signal
p =
ggplot(cell_viz,
       aes(x = Live_dead, y = DNA_Intercalator_2)) +
  geom_hex(bins = 128) +
  scale_fill_gradientn("Count", colors = jc_color_cont)

p + annotate(geom="rect",
             xmin = gating_quantile[X, LL],
             xmax =gating_quantile[Y, UL],
             ymin = gating_quantile[Y, '0%'],
             ymax = gating_quantile[Y, '100%'],
             colour = "red", fill = NA)

expr_live = expr_cell %>%
  base::subset(expr_cell[, X] < gating_quantile[X, UL])
# update viz dataframe
live_viz = cell_viz %>%
  dplyr::filter(
    Live_dead < gating_quantile[X, UL] )
  

## -- step 3: Gate out CD45+ cells
gating_quantile = colQuantiles(expr_live[, c(gating_markers, 'CD45')],
                               probs = seq(0, 1, by = 0.05))

LL ='20%'
UL = '100%'
X = "CD45"
Y = "DNA_Intercalator_2"

# viz the CD45 vs DNA intercaltor
p =
  ggplot(live_viz,
         aes(x = CD45, y = DNA_Intercalator_2)) +
  geom_hex(bins = 128) +
  scale_fill_gradientn("Count", colors = jc_color_cont)

p + annotate(geom="rect",
             xmin = gating_quantile[X, LL],
             xmax =gating_quantile[X, UL],
             ymin = gating_quantile[Y, '0%'],
             ymax = gating_quantile[Y, '100%'],
             colour = "red", fill = NA)

# gate out CD45+
expr_CD45 = expr_live %>%
  base::subset(expr_live[, X] > gating_quantile[X, LL])

CD45_viz = live_viz %>%
  dplyr::filter(
    CD45 > gating_quantile[X, LL]
  )

# ----------------------- Check the quality of the Gating----------------------------
expr_final = as.data.frame(expr_CD45) %>%
  dplyr::mutate(
    'sample_id' = rownames(expr_CD45)
  )

gating_size = data.frame('size' = 
                           c(nrow(expr_raw), nrow(expr_cell), 
                           nrow(expr_live), nrow(expr_CD45)),
                         'filter_criteria' = c('raw', 'cell', 'live', 'CD45'))
                           
expr_final_quant = colQuantiles(expr_CD45[, lineage_markers], 
                                probs = seq(0, 1, by = 0.05))

# compute % of non-zero marker expression
p0 = apply(expr_CD45[, lineage_markers], MARGIN = 2, FUN = function(x){
  round(sum(x>0) / length(x) * 100, 2)
})

p1 = apply(expr_CD45[, lineage_markers], MARGIN = 2, FUN = function(x){
  round(sum(x>1) / length(x) * 100, 2)
})

p2 = apply(expr_CD45[, lineage_markers], MARGIN = 2, FUN = function(x){
  round(sum(x>2) / length(x) * 100, 2)
})

ggdf = data.frame(
  'marker' = names(p0),
  'larger0' = p0, 'larger1' = p1, 'larger2' = p2) %>%
  pivot_longer(cols = larger0:larger2, names_to = "criteria", 
               values_to = "percent")

ggplot(ggdf, aes(x = marker, y = percent, fill = criteria)) +
  geom_bar(stat = 'identity', position = position_dodge2()) + 
  geom_text(aes(label=percent),vjust = -1, hjust = 0.5,
            fontface = "bold", size=2.5) +
  scale_fill_manual(values = jc_color_cate) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.05),
        axis.text = element_text(size = 10))

# viz final marker expression qunatiles in heatmap
rownames(expr_final_quant) = str_c(colnames(expr_CD45[, lineage_markers]), 
                                   "(", p0,"%)")
Heatmap(expr_final_quant, name = "markers expression", 
        col = jc_color_cont,
        row_title = "Marker", column_title = "Quantiles %",
        cluster_rows = TRUE, cluster_columns = FALSE,
        border = TRUE, rect_gp = gpar(col = "white", lwd = 2),
        column_names_rot = 45,
        layer_fun = function(j, i, x, y, width, height, fill) {
          # print out intensity text on each cell
          grid.text(sprintf("%.2f", 
                            pindex(expr_final_quant, i, j)), 
                    x, y, gp = gpar(fontsize = 10))}
)

# density of all markers
CD45_viz %>%
  melt(id.var = "sample_id", value.name = "expression",
       variable.name = "antigen") %>%
  ggplot(aes(x = expression, group = sample_id)) +
    geom_density()+
    facet_wrap(~antigen, nrow = 6, scales = "free") +
    ggtitle("Markers Density") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 7), 
          axis.text = element_text(size = 5)) +
    guides(color = guide_legend(ncol = 1))

## zoom in to check on week expressed markers 
rownames(expr_final_quant) = colnames(expr_CD45[, lineage_markers])
test = expr_CD45[, lineage_markers]


expr_final %>%
  dplyr::filter(
    CD45 > expr_final_quant['CD45(100%)', "95%"]
  ) %>%
  melt(id.var = "sample_id", value.name = "expression",
       variable.name = "antigen") %>%
  ggplot(aes(x = expression, group = sample_id)) +
    geom_density()+
    facet_wrap(~antigen, nrow = 6, scales = "free") +
    ggtitle("Markers Density Zoom in CD45{95%, 100%} (5% of the final verision)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 7), 
          axis.text = element_text(size = 5)) +
    guides(color = guide_legend(ncol = 1))

# output gated file into csv file
dir_gate = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/gated/BAL_02.06.2020/"
setwd(dir_gate)

expr_final = expr_final %>%
  dplyr::mutate(
    # add gating procedure note
    'gate_note' = rep('cell: DNA1&2{85%, 100%}; live: Live{0%, 80%}; CD45: CD45{20%, 100%}',
                 nrow(expr_final))
  )

write_csv(expr_final, 'BAL_new_gated_07.29.2020.csv')
write_csv(gating_size, "BAL_new_gated_size_07.29.2020.csv")






