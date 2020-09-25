# ---------------UCSF_VA_NCIRE PROJECT: CyTOF single cell data analysis-----------
#   Jianhong Chen
#   07/24/2020
# CyTOF fcs files pre-processing II: gating version 2.0
#     
# Handle csv files that is preprocessed via CATALYST package
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear workspace
rm(list = ls()) # delete all variables in the environment
# Clear console
cat("\014")  # equals to 'control + L'
gc()


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
theme_update(plot.title = element_text(hjust = 0.5, size = 25),
             axis.text = element_text(size = 15),
             axis.title = element_text(size = 20),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 15))
jc_color_cate = brewer.pal(n = 9, "Set1")[-3]
jc_color_cont = colorRampPalette(rev(brewer.pal(n=11, name = "RdYlBu")))(100)


# --------------------------------- Load Raw Data---------------------------------
dir_raw = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/COPD_BAL_10plex_set1_02.22.2018/processed/"
setwd(dir_raw)

expr_raw = read_csv("COPD_BAL_Feb2018_norm_debarcode_asinh_expr_09.14.2020.csv")
# expr2_raw = read_csv("LT_old_part2_norm_debarcode_asinh_expr_08.05.2020.csv")

# concatenatedd the two parts into one complete expr data
# expr_raw = rbind(expr1_raw, expr2_raw)
# remove '-' from colanmes 
colnames(expr_raw) = str_replace_all(colnames(expr_raw), "-", "_")

# read in yield size file
yield_size = read_csv("COPD_BAL_Feb2018_yield_size_09.14.2020.csv")

# -- load metadata for the study
dir_meta = '/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/COPD_BAL_10plex_set1_02.22.2018/metadata/'
setwd(dir_meta)

panel = read_excel("gating_panel.xlsx")
# remove '-' from labels
panel$target= str_replace_all(panel$target, "-", "_")
panel$marker = str_replace_all(panel$marker, "-", "_")

## fix some of the marker names for the raw expression data
colnames(expr_raw) = str_replace_all(colnames(expr_raw), "-", "_")
colnames(expr_raw)[colnames(expr_raw) == "CD19-56-66b"] = "CD56_19_66b"
colnames(expr_raw)[colnames(expr_raw) == "PD-1"] = "PD_1"
colnames(expr_raw)[colnames(expr_raw) == "PDL-2"] = "PDL2"
colnames(expr_raw)[colnames(expr_raw) == "CCR2"] = "CD192"
colnames(expr_raw)[colnames(expr_raw) == "Galectin9"] = "Galectin_9"
colnames(expr_raw)[colnames(expr_raw) %in% 
                     c("DNAIntercalator1", 'DNAIntercalator2', 'Live_Dead')] = 
  c('DNA_Intercalator_1', 'DNA_Intercalator_2', 'Live_dead')

# define target markers to study 
# NOTE: no functional markers in this study
( lineage_markers = panel$marker[panel$lineage == 1] )
# spot check
if (all(lineage_markers %in% colnames(expr_raw))) {
  TRUE
} else { # if not, check each indiviual marker 
  sapply(lineage_markers, FUN = function(x){x %in% colnames(expr_raw)})
}


# --------------------------------- Pre-processing ---------------------------------
# generate sample IDs corresponding to each cell in the 'exprs' matrix
sample_ids = expr_raw$sample_id

# select required gating channels
gating_markers = c('DNA_Intercalator_1', 'DNA_Intercalator_2', 
                   'Live_dead')

# spot check
all(gating_markers %in% colnames(expr_raw))

## if markers are in different readable format
# once the spot check passed, renamed the markers into more readable version
# gating_markers = map_chr(gating_markers, .f = function(x){
#   str_c(unlist(str_split(x, "_"))[-1], collapse = "_")
# })

# # once the spot check passed, renamed the panel markers into more readable version
# panel_fcs$desc = map_chr(panel_fcs$desc, .f = function(x){
#   str_c(unlist(str_split(x, "_"))[-1], collapse = "_")
# })

## remove useless channels from the raw expr data
expr_mat = expr_raw[, c(lineage_markers, gating_markers)] %>%
  as.matrix()
rownames(expr_mat) = sample_ids

## -- **** in order to speed up the analysis of the large raw dataset   ****
# viz only 10% of the original data
#  however the whole analysis still utilize the whole set(mainly in a matrix)

# calculate the 10% sample size per fcs files
p = round(nrow(expr_mat)*0.1 / length(unique(sample_ids)))

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

raw_viz = as.data.frame(expr_mat[dsamp_inds, ]) %>%
  mutate("sample_id" = sample_ids[dsamp_inds] )

# ---------------------------------  Gating ---------------------------------------

# viz the num of event per samples
expr_mat %>%
  as.data.frame() %>%
  mutate('sample_id' = sample_ids) %>%
  group_by(sample_id) %>%
  tally() %>%
  ungroup() %>%
  ggplot(aes(x = reorder(sample_id, -n), y =n)) + 
  geom_bar(stat = "identity") + 
  ggtitle("Debarcoded Number of Events per Sample") +
  geom_text(aes(label=n),vjust = -1, fontface = "bold", size=2.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title.x = element_blank())



## -- step 1: filter out live & singlet cells
gating_quantile = colQuantiles(expr_mat[, gating_markers],
                               probs = seq(0, 1, by = 0.05))
LL ='90%'
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

expr_cell = expr_mat %>%
  base::subset(expr_mat[, X] > gating_quantile[X, LL] &
                 expr_mat[, X] < gating_quantile[X, UL] &
                 expr_mat[, Y] > gating_quantile[Y, LL] &
                 expr_mat[, Y] < gating_quantile[Y, UL])

# update the viz dataframe as well
cell_viz = raw_viz %>%
  dplyr::filter(
      DNA_Intercalator_1 > gating_quantile[X, LL] &
        DNA_Intercalator_1 < gating_quantile[X, UL] &
        DNA_Intercalator_2 > gating_quantile[Y, LL] &
        DNA_Intercalator_2 < gating_quantile[Y, UL])

# VIZ the num of event after DNA-gating
expr_cell %>%
  as.data.frame() %>%
  mutate('sample_id' = rownames(expr_cell)) %>%
  group_by(sample_id) %>%
  tally() %>%
  ungroup() %>%
  ggplot(aes(x = reorder(sample_id, -n), y =n)) + 
  geom_bar(stat = "identity") + 
  ggtitle("DNAInterC Gating Number of Events per Sample") +
  geom_text(aes(label=n),vjust = -1, fontface = "bold", size=2.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title.x = element_blank())


## step 2 gate out dead cells
gating_quantile = colQuantiles(expr_cell[, gating_markers],
                               probs = seq(0, 1, by = 0.05))
LL ='0%'
UL = '90%'
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
             xmax =gating_quantile[X, UL],
             ymin = gating_quantile[Y, '0%'],
             ymax = gating_quantile[Y, '100%'],
             colour = "red", fill = NA)

expr_live = expr_cell %>%
  base::subset(expr_cell[, X] < gating_quantile[X, UL])
# update viz dataframe
live_viz = cell_viz %>%
  dplyr::filter(
    Live_dead < gating_quantile[X, UL] )
  
# VIZ the num of event after Live-gating
expr_live %>%
  as.data.frame() %>%
  mutate('sample_id' = rownames(expr_live)) %>%
  group_by(sample_id) %>%
  tally() %>%
  ungroup() %>%
  ggplot(aes(x = reorder(sample_id, -n), y =n)) + 
  geom_bar(stat = "identity") + 
  ggtitle("Live/Dead Gating Number of Events per Sample") +
  geom_text(aes(label=n),vjust = -1, fontface = "bold", size=2.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title.x = element_blank())

## -- step 3: Gate out CD45+ cells
gating_quantile = colQuantiles(expr_live[, c(gating_markers, 'CD45')],
                               probs = seq(0, 1, by = 0.05))

LL ='15%'
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
             xmax = gating_quantile[X, UL],
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


# viz the num of event per samples
expr_final %>%
  group_by(sample_id) %>%
  tally() %>%
  ungroup() %>%
  ggplot(aes(x = reorder(sample_id, -n), y =n)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label=n),vjust = -1, fontface = "bold", size=2.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title.x = element_blank())

gating_size = data.frame(
  'step' = c('cell_gate', 'live_gate', 'CD45_gate'),
  'num_event' = c(nrow(expr_cell), nrow(expr_live), nrow(expr_CD45)),
  'note' = c('get actual cell: DNA1&2{90%, 100%}' , 
             'remove dead cell: Live{0%, 90%}',
             'CD45+: CD45{15%, 100%}'),
  'yield' = c(round(nrow(expr_cell)/yield_size$num_event[1]*100, 2),
              round(nrow(expr_live)/yield_size$num_event[1]*100, 2),
              round(nrow(expr_CD45)/yield_size$num_event[1]*100, 2))
  )
                           
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

# viz the marker crietria computed as above into barplot
ggdf = data.frame(
  'marker' = names(p0),
  'larger0' = p0, 'larger1' = p1, 'larger2' = p2) %>%
  pivot_longer(cols = larger0:larger2, names_to = "criteria", 
               values_to = "percent")

ggplot(ggdf, aes(y = reorder(marker, percent),x = percent, fill = criteria)) +
  geom_bar(stat = 'identity', position = position_dodge2()) + 
  geom_text(aes(label=percent),vjust = 0.5, hjust = 0.5,
            fontface = "bold", size=2.5) +
  facet_wrap(~criteria, scales = 'free_y') +
  scale_fill_manual(values = jc_color_cate) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.05),
        axis.text = element_text(size = 10)) +
  labs(x = 'Percentage', y = 'Marker')

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

## boxplot for all markers
CD45_viz %>%
  melt(id.var = "sample_id", value.name = "expression",
       variable.name = "antigen") %>%
  ggplot(aes(x = reorder(antigen, -expression, FUN = median),
             y = expression, group = antigen)) +
  geom_boxplot(outlier.size = 0.5)+
  ggtitle("Markers Boxplot") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 7), 
        axis.text = element_text(size = 15)) +
  labs(x = 'Marker', y = 'Expression')

# output gated file into csv file
dir_output = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/COPD_BAL_10plex_set1_02.22.2018/processed/"
setwd(dir_output)

yield_final = rbind(yield_size, gating_size)

write_excel_csv(expr_final, 'COPD_BAL_Feb2018set1_gated_09.15.2020.csv')
write_excel_csv(yield_final, "COPD_BAL_Feb2018set1_final_yield_09.15.2020.csv")



# ------------- Individual Samples Comparison between Runs -----------------------
## loaded processed & gated data set
# LT_old in 2018
dir_LT_old = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/New_Data/LT_Old10x_2018_08.04.2020/preprocessed_data/"
setwd(dir_LT_old)
LT_old = read_csv('LT_old_gated_08.05.2020.csv')
LT_old$sample_id = str_remove(LT_old$sample_id, ' ')

# LT_new1 in Jul 2020
dir_LT_new1 = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/New_Data/LT_Cytof_Blue_07.29.2020/preprocessed_data/"
setwd(dir_LT_new1)
LT_new1 = read_csv("LT_R2_gated_08.03.2020.csv")
LT_new1$sample_id = str_remove(LT_new1$sample_id, "0")

# LT_new2 in Aug 2020
dir_LT_new2 = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/New_Data/Full_panel_10plex_08.25.2020/preprocessed_data/"
setwd(dir_LT_new2)
LT_new2 = read_csv("fullPanel_gated_08.26.2020.csv")
LT_new2$sample_id = str_remove(LT_new2$sample_id, ' ')

# pre-process the data again so that they match
LT_old_sub = LT_old %>%
  mutate('type' = rep('old', nrow(LT_old))) %>%
  dplyr::filter(sample_id %in% c('LT22','LT46')) %>%
  dplyr::rename(
    'CD192' = CCR2,
    'CD56_19_66b' = CD19_56_66b,
    'PDL2' = PDL_2,
    'Galectin_9' = Galectin9,
    'Live_dead' = Live_Dead,
    'DNA_Intercalator_1' = DNAIntercalator1,
    'DNA_Intercalator_2' = DNAIntercalator2
  ) %>%
  transform(
    sample_id = str_c(sample_id, "_old")
  )

LT_new1_sub = LT_new1 %>%
  mutate(
    'type' = rep('new1', nrow(LT_new1)) ) %>%
  dplyr::filter(sample_id %in% c('LT22', 'LT46')) %>%
  transform(
    sample_id = str_c(sample_id, "_new1")
  )

LT_new2_sub = LT_new2 %>%
  mutate(
    'type' = rep('new2', nrow(LT_new2)) ) %>%
  dplyr::filter(sample_id %in% c('LT22', 'LT46')) %>%
  transform(
    sample_id = str_c(sample_id, "_new2")
  )

## -- combine old and new experiments
# reorder the colnames for both data
LT_old_sub = LT_old_sub[, order(colnames(LT_old_sub))]
LT_new1_sub = LT_new1_sub[, order(colnames(LT_new1_sub))]
LT_new2_sub = LT_new2_sub[, order(colnames(LT_new2_sub))]

# check if both colnames match
all(colnames(LT_old_sub) == colnames(LT_new1_sub))
all(colnames(LT_new1_sub) == colnames(LT_new2_sub))

LT_combine = rbind(LT_old_sub, LT_new1_sub, LT_new2_sub)

# density plot
## make sure each sample get the same number of events during density plot
inds = split(1:length(LT_combine$sample_id), LT_combine$sample_id)

# step 3: size of the down-sampling
dsamp_ncells = pmin(table(LT_combine$sample_id), 5000) # pmin = parallel min

# step 4: randomly pick indices base on the chosen size
set.seed(78)
dsamp_inds = lapply(names(inds), function(i){
  s = sample(inds[[i]], dsamp_ncells[i], replace = FALSE)
}) %>%
  unlist()

# viz down-sample all markers histogram
LT_combine[dsamp_inds, ] %>%
  pivot_longer(cols = c(Axl:PDL2, TLR4:TLR9), names_to = 'markers',
               values_to = 'expression') %>%
  ggplot(aes(x = expression, group = type, fill = type)) +
  stat_bin(bins = 128)+
  facet_wrap(~markers, nrow = 6, scales = "free") +
  ggtitle("Markers Density") +
  scale_fill_manual(values = jc_color_cate) +
  guides(color = guide_legend(ncol = 1))


# viz down-sample all markers boxplot
LT_combine[dsamp_inds, ] %>%
  pivot_longer(cols = c(Axl:PDL2, TLR4:TLR9), names_to = 'markers',
               values_to = 'expression') %>%
  ggplot() +
  geom_boxplot(aes(x = reorder(markers, -expression, FUN = median), 
                   y = expression, color = type, fill = type),
               outlier.size = 0.5,
               position = position_dodge(), alpha = 0.5) +
  ggtitle("Markers Density") +
  scale_color_manual(values = jc_color_cate) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  guides(color = guide_legend(ncol = 1)) +
  labs(x = 'Marker', y = 'Expression')




# output gated file into csv file
dir_output = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/New_Data/Full_panel_10plex_08.25.2020/preprocessed_data/"
setwd(dir_output)

write_csv(LT_combine, "LT_combine_09.08.2020.csv")



 



