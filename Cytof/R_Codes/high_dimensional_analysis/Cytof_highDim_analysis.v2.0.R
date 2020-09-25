# ---------------UCSF_VA_NCIRE PROJECT: CyTOF single cell data analysis---------##
#   Jianhong Chen
#   07/06/2020
#
#   Pipeline II:
#  High Dim. Analyzing preprocessed & gated data
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##

# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear workspace
rm(list = ls()) # delete all variables in the environment
# Clear console
cat("\014")  # equals to 'control + L'


# ------------------------ Load Libraries & Global Setup ----------------------

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

# GUI
library(cytofkit)

# Utilitise
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
#GOLD-0, GOLD-1, GOLD-2, GOLD-3, GOLD-4: 
color_gold = c("#FFDD00","#FFCC00","#FFAA00", "#FF8800", "#FF6600", "#FF4400")
# -- generate heatmap plots: illustrate median marker expression in each cluster
color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", 
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

# ------------------------ Load Raw/Pre-processed Data----------------------
# -- read (CD45 gated) & asinh(cofactor = 5) transfomred data
dir_raw = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/New_Data/Full_panel_10plex_08.25.2020/preprocessed_data/"
setwd(dir_raw)
expr_raw = read_csv("fullPanel_gated_08.26.2020.csv")

# -- load clinical metadata for the study
dir_meta = '/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/New_Data/LT_Old10x_2018_08.04.2020/metadata/'
setwd(dir_meta)

# pft_raw = read_excel("LT_clinical_db_7.10.2020.xlsx", sheet = 1)
# symptom_raw = read_excel("LT_clinical_db_7.10.2020.xlsx", sheet = 4)

panel = read_excel("panel35.xlsx") 
panel$target= str_replace_all(panel$target, "-", "_")
panel$marker = str_replace_all(panel$marker, "-", "_")

# ------------------------------ Processing  Data--------------------------------

# define target markers to study 
# NOTE: no functional markers in this study
( lineage_markers = panel$marker[panel$lineage == 1] )
# make change of the marker names to match with the data
lineage_markers[c(3,13, 22,31)] = c('CD56_19_66b', 'CD192', 
                                    'Galectin_9', 'PDL2')

# spot check
if (all(lineage_markers %in% colnames(expr_raw))) {
  TRUE
} else { # if not, check each indiviual marker 
  sapply(lineage_markers, FUN = function(x){x %in% colnames(expr_raw)})
}

# generate sample IDs corresponding to each cell in the 'exprs' matrix
# sample_ids = expr_raw$sample_id %>%
#   str_remove("HLT") %>%
#   str_remove("LT") %>% 
#   as.numeric() %>%
#   as.factor()
sample_ids = expr_raw$sample_id

# convert raw expression into a matrix data type
expr_mat = expr_raw[, lineage_markers] %>%
  as.matrix()
rownames(expr_mat) =sample_ids
# create an alternative expr version for plotting only (range from 0 to 1)
rng = colQuantiles(expr_mat, probs = c(0.01, 0.99))
expr01 = t((t(expr_mat) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] = 0
expr01[expr01 > 1] = 1
# z-normalization of the raw expression by scale()
expr_z = scale(expr_mat)
expr_z[expr_z < -2] = -2
expr_z[expr_z > 2] = 2

##-- extract useful clinical metadata
# match_id = levels(sample_ids) %>%
#   as.numeric()
#   
# str_subset(colnames(pft_raw), 'pack_year')
# str_subset(colnames(symptom_raw), 'gold')
# 
# pft_sub = data.frame(
#   'sample_id' = as.numeric(pft_raw$`Study ID`)) %>%
#   dplyr::mutate(
#     'age' = as.numeric(pft_raw$Age),
#     'sex' = as.factor(pft_raw$Gender),
#     'gold' = as.factor(pft_raw$`Gold Stage`),
#     'ht_cm' = as.numeric(pft_raw$`Height(cm)`),
#     'wt_kg' = as.numeric(pft_raw$`Weight(kg)`),
#     'pre_fev1' = as.numeric(pft_raw$`Pre_FEV1 (L)`),
#     'pre_fvc' = as.numeric(pft_raw$`Pre_FVC (L)`),
#     'rv' = as.numeric(pft_raw$RV),
#     'tlc' = as.numeric(pft_raw$`TLC(L)`),
#     'smoker' = as.factor(pft_raw$`Current Smoker?`),
#     'copd_smoker' = as.factor(pft_raw$Group)
#   ) %>%
#   dplyr::filter(
#     sample_id %in% match_id
#   ) 

# merged clinical metadata into the marker expression data
# expr_sub = expr_mat %>%
#   as.data.frame() %>%
#   dplyr::mutate(
#     'sample_id' = sample_ids ) %>%
#   merge(pft_sub, by = 'sample_id', all = FALSE)
 
# -- Compute Median intensity per sample per marker
expr_median_sample_tbl = data.frame(sample_id = sample_ids, expr_mat) %>%
  group_by(sample_id) %>%
  summarize_all(median)
# rang 0-1 data is used for visualization only
expr01_median_sample_tbl = data.frame(sample_id = sample_ids, expr01) %>%
  group_by(sample_id) %>%  summarize_all(median)

expr_median_sample = t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) = expr_median_sample_tbl$sample_id

expr01_median_sample = t(expr01_median_sample_tbl[, -1])
colnames(expr01_median_sample) = expr_median_sample_tbl$sample_id

# ------------------------------ EDA -----------------------------------
# num. of events per sample in barplot
expr_sub = expr_raw
expr_sub%>%
  group_by(sample_id) %>%
  tally() %>%
  ungroup() %>%
  ggplot(aes(x = reorder(sample_id, -n),
             y = n) ) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label=n),vjust = -1, fontface = "bold", 
              size=2.5, color = 'red') +
    labs(x = 'sample')

# marker intensity density plot
expr_sub %>%
  pivot_longer(cols = CD68:CD1c, names_to = 'marker', 
               values_to = 'expression') %>%
  ggplot(aes(x = expression, group = sample_id)) +
    geom_density()+
    facet_wrap(~marker, nrow = 6, scales = "free") +
    ggtitle("Markers Density") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 7),
          axis.text = element_text(size = 5)) +
    guides(color = guide_legend(ncol = 1)) +
    scale_color_manual(values = jc_color_cate)

# -- PCA
# projected the original matrixed into PCA space
mds = limma::plotMDS(expr_median_sample, plot = FALSE)
pca.res = prcomp(t(expr_median_sample), scale. = F)
eigen_percent = round(100*pca.res$sdev^2 / sum(pca.res$sdev^2), 1)

ggdf = data.frame(MDS1 = mds$x, MDS2 = mds$y, 
                  PC1 = pca.res$x[, 1], PC2 = pca.res$x[, 2], 
                  PC3 = pca.res$x[, 3], PC4 = pca.res$x[, 4],
                  sample_id = colnames(expr_median_sample))
mm = match(ggdf$sample_id, expr_sub$sample_id)
# ggdf = mutate(ggdf,  "copd_smoker" = expr_sub$copd_smoker[mm])
             
ggplot(ggdf, aes(x = MDS1, y = MDS2))+
  geom_point(size = 2, alpha =0.8) +
  geom_label_repel(aes(label = sample_id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab(str_c("PC1", " (", eigen_percent[1], "%)")) +
  ylab(str_c("PC2", " (", eigen_percent[2], "%)"))


# -- Heatmap viz for median intensity
# rownames(pft_sub) = pft_sub$sample_id
# col_anno = HeatmapAnnotation('copd_smoker' = pft_sub$copd_smoker,
#                              'sex' = pft_sub$sex)
# row_anno = rowAnnotation(df = pft_sub$copd_smoker)

Heatmap(expr_median_sample, name = "sample", 
        col = jc_color_cont,
        row_title = "Marker", 
        column_title = "Marker Median intensity per sample",
        cluster_rows = TRUE, cluster_columns = TRUE,
        border = TRUE, rect_gp = gpar(col = "white", lwd = 2),
        column_names_rot = 90, 
        # top_annotation = col_anno,
        layer_fun = function(j, i, x, y, width, height, fill) {
          # print out intensity text on each cell
          grid.text(sprintf("%.2f", 
                            pindex(expr_median_sample, i, j)), 
                    x, y, gp = gpar(fontsize = 10))}
)


# ------------------------- Dimensional Reduction --------------------------
## downsample original data 
# step 1: check duplicattion
dups = which(!duplicated(expr_mat[, lineage_markers]))

# step 2: create full indices by sample
inds = split(1:length(sample_ids), sample_ids)

# step 3: size of the downsampling
dsamp_ncells = pmin(table(sample_ids), 5000) # pmin = parallel min

# step 4: randomly pick indices base on the chosen size
set.seed(78)
dsamp_inds = lapply(names(inds), function(i){
  s = sample(inds[[i]], dsamp_ncells[i], replace = FALSE)
  intersect(s, dups)
}) %>%
  unlist()

# step 5: subset original expression based on the indices
dsamp_expr = expr_mat[dsamp_inds, lineage_markers]

# -- Run t-SNE 
set.seed(78)
# perplexity resonable range = 5-50 by the authors of t-SNE
tsne_out = Rtsne(dsamp_expr, perplexity = 100, check_duplicates = FALSE,
                 pca = FALSE, verbose = TRUE)

# -- Run UMAP
set.seed(78)
# n_neighbors: small >> local details; large >> global details;
# min_dist: small >> densely packed; large >> more spread out
umap_out = uwot::umap(dsamp_expr, n_neighbors = 20, n_components = 2,
                      learning_rate = 1, min_dist = 0.001,
                      init = "random", n_epochs = 500, scale = FALSE,
                      verbose = TRUE)

# create DR dataframe that contains both tSNE and UMAP
drdf = data.frame(
  # mannual pivot_longer in this case in order to match X and Y exactly
  X = c(tsne_out$Y[, 1],umap_out[, 1]),
  Y = c(tsne_out$Y[, 2],umap_out[, 2]),
  DR = rep(c("tSNE", "UMAP"), c(nrow(tsne_out$Y), nrow(umap_out))),
  rbind(dsamp_expr, dsamp_expr),
  'sample_id' = rep(sample_ids[dsamp_inds], 2))

# check batch effect
ggplot(drdf, aes(x = X, y = Y, color = sample_id)) +
  geom_point(size = 0.8, alpha = 0.8) +
  facet_wrap(~DR, scales = "free") +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
  
# viz one marker
ggplot(drdf, aes(x = X, y = Y, color = CD64)) +
  geom_point(size = 0.8, alpha = 0.8) +
  facet_wrap(~DR, scales = "free") +
  scale_color_gradientn("CD64", colours = jc_color_cont)

# viz all markers in one DR space with z-norm marker expression
data.frame(expr_z[dsamp_inds, ]) %>%
  dplyr::mutate(
    X = dplyr::filter(drdf, DR == 'UMAP')$X,
    Y = dplyr::filter(drdf, DR == 'UMAP')$Y ) %>%
  pivot_longer(cols = CD68:CD1c, names_to = 'marker', 
               values_to = 'value') %>%
  ggplot(aes(x= X, y = Y, color = value)) +
    geom_point(alpha = 0.8, size = 0.8) +
    facet_wrap(~marker) +
    scale_colour_gradientn("Scaled Expr", colors = jc_color_cont)


# # # ------ test with CATALYST UMAp wrapper function
# daf = daFrame(fcs, panel, md, cols_to_use = lineage_markers)
# 
# daf = runDR(daf, "TSNE", cols_to_use = lineage_markers, rows_to_use = 1000)
# daf = runDR(daf, "UMAP", cols_to_use = lineage_markers, rows_to_use = 2000)
# 
# plotDR(daf, "TSNE", color_by = "CD3")
# plotDR(daf, "UMAP", color_by = "CD3")

# ------------------------------ Clustering -------------------------------------

# -- Clustering Method: RphenoGraph using down-sample data
set.seed(78)
phenoGraph.res = Rphenograph(dsamp_expr, k = 40)
modularity(phenoGraph.res)
phenoGraph.cluster = membership(phenoGraph.res)

cldf = drdf %>%
  dplyr::mutate( 
    "pheno_cluster" = rep(factor(phenoGraph.cluster), 2)) 
# add the clinical data back into the clustered expression data
# mm = match(cldf$sample_id, pft_sub$sample_id)
# cldf = data.frame(cldf, pft_sub[mm, -1])
             

# heatmap viz the clustering result:
plot_clustering_heatmap_wrapper = function(expr, expr01,
                                           cell_clustering, color_clusters,
                                           cluster_merging = NULL){
  # --------------------------------------------------------------
  # wrapper function for plotting Hcluster heatmap
  # inputs: @expr: expression sets of the fcs file
  #         @expr01: alternative c(0,1) range expression set of the original one
  #         @cell_clustering: clustered label for each cell
  #         @color_clusters: pre-set color code for plotting
  # --------------------------------------------------------------
  
  # calculate row-wise median expression 
  expr_median = data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_all(median)
  
  expr01_median = data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarise_all(median)
  
  # calculate cluster frequencies
  clustering_table = as.numeric(table(cell_clustering))
  
  # This clustering is based on the markers that were used for the main clustering
  # d_row = dist(expr_median[, colnames(expr)], method = "euclidean")
  # cluster_rows = hclust(d_row, method = "complete")
  
  expr_heat = as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) = expr01_median$cell_clustering
  
  labels_row =  paste0(rownames(expr_heat), " (",
                       round(clustering_table / sum(clustering_table)*100, 2),
                       "%)")
  labels_col = colnames(expr_heat)
  
  rownames(expr_heat) = labels_row
  # Row annotation for the heatmap
  annotation_row = data.frame(cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) = rownames(expr_heat)
  
  color_clusters = color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) = levels(annotation_row$cluster)
  annotation_colors = list(cluster = color_clusters)
  annotation_legend = FALSE
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster 
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }
  
  # colors for the heatmap
  # color = colorRampPalette(rev(brewer.pal(n=9, name = "RdYlBu")))(100)
  color = colorRampPalette(c("gray9","navyblue","lightgoldenrodyellow",
                             "red3","red4"))(n =30)
  # heatmap
  # pheatmap::pheatmap(expr_heat, color = color,
  #                    cluster_cols = TRUE, cluster_rows = TRUE,
  #                    label_col = labels_col, labels_row = labels_row,
  #                    display_numbers = TRUE, number_color = "black",
  #                    fontsize = 8, fontsize_number = 4,
  #                    annotation_row = annotation_row, annotation_colors = annotation_colors,
  #                    annotation_legend = annotation_legend)
  
  # -- Heatmap
  # rownames(condition_anno) = md$sample_id
  # col_anno = HeatmapAnnotation(df = condition_anno)
  # row_anno = rowAnnotation(df = annotation_row)
  
  Heatmap(expr_heat, name = "Scaled Expr",
          col = color,
          row_title = "Cluster", column_title = "Marker Median Expression",
          cluster_rows = TRUE, cluster_columns = TRUE,
          clustering_method_rows = "complete",
          clustering_method_columns = "complete",
          border = TRUE, rect_gp = gpar(col = "white", lwd = 2),
          column_names_rot = 90,
          # right_annotation = row_anno,
          layer_fun = function(j, i, x, y, width, height, fill) {
            # print out intensity text on each cell
            grid.text(sprintf("%.2f",
                              pindex(expr_heat, i, j)),
                      x, y, gp = gpar(fontsize = 6))}
  )
  
}

# phenoGraph cluster heatmap
plot_clustering_heatmap_wrapper(
  expr = cldf[1:length(dsamp_inds), lineage_markers],
  expr01 = expr_z[dsamp_inds, lineage_markers],
  cell_clustering = cldf$pheno_cluster[1:length(dsamp_inds)],
  color_clusters = color_clusters)

# -- Visual Clustering result
ggplot(cldf, aes(x = X, y = Y, color = pheno_cluster)) +
  geom_point(size = 0.8, alpha = 0.8) +
  facet_wrap(~DR, scales = "free") +
  ggtitle("Clustering by PhenoGraph") + 
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

# VIZ cluster with condition in UMAP
cldf %>%
  dplyr::filter(DR == 'UMAP') %>%
  ggplot(aes(x = X, y = Y, color = pheno_cluster)) +
  geom_point(alpha = 0.8, size = 0.8) +
  facet_wrap(~copd_smoker) + 
  ggtitle("Condition: COPD Smoker Group") +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))


# plot all clusters separately
data.frame(expr_z[dsamp_inds, ]) %>%
  dplyr::mutate(
    X = dplyr::filter(cldf, DR == 'UMAP')$X,
    Y = dplyr::filter(cldf, DR == 'UMAP')$Y ,
    'pheno_cluster' = dplyr::filter(cldf, DR == 'UMAP')$pheno_cluster,
    'copd_smoker' = dplyr::filter(cldf, DR == 'UMAP')$pheno_cluster) %>%
  pivot_longer(cols = CD68:CD1c, names_to = 'marker', 
             values_to = 'value') %>%
  ggplot(aes(x = X, y = Y, color = value)) +
  geom_point(size =0.8, alpha = 0.8) +
  facet_grid(pheno_cluster~marker) + 
  scale_color_gradientn("Maker", colors = jc_color_cont)


# ------------------- Annotating Clusters Results ----------------------------------------
# -- load the annotated label table
cluster_annotation = read_excel("marker35_pheno_annotation_08.12.2020.xlsx")

# add the merged & annotated cell cluster into the dataframe
mm = match(cldf$pheno_cluster, cluster_annotation$pheno_cluster)
anno_cluster= cluster_annotation$pheno_annotated[mm]

cldf = mutate(cldf, "cell" = factor(anno_cluster))

# -- viz on the annotated cell clusters
ggplot(cldf, aes(x = X, y = Y, color = cell)) +
  geom_point(size = 0.8) +
  facet_wrap(~DR, scale = 'free') +
  ggtitle("Annotated Cell Clusters") + 
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

# VIZ Annotated cluster on condition in UMAP
cldf %>%
  dplyr::filter(DR == 'UMAP') %>%
  ggplot(aes(x = X, y = Y, color = cell)) +
  geom_point(alpha = 0.8, size = 0.8) +
  facet_wrap(~copd_smoker) + 
  ggtitle("Condition: COPD Smoker Group") +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

# VIZ Annotated cluster on condition in UMAP
cldf %>%
  dplyr::filter(DR == 'UMAP') %>%
  ggplot(aes(x = X, y = Y, color = cell)) +
  geom_point(alpha = 0.8, size = 0.8) +
  facet_wrap(~sample_id) + 
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))


# performed PCA to verify the annotation
anno_expr_median = cldf %>%
  dplyr::select(CD68:CD1c, cell)%>%
  group_by(cell) %>%
  summarise_all(median) %>%
  ungroup()

anno_expr_median_mat = anno_expr_median[, -1] %>%
  as.matrix()
rownames(anno_expr_median_mat) = anno_expr_median$cell

pca_res = prcomp(anno_expr_median_mat, scale. = F)
eigen_percent = round(100*pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
           PC3 = pca_res$x[, 3], PC4 = pca_res$x[, 4], 
           "cell" = rownames(anno_expr_median_mat)) %>%
  ggplot(aes(x = PC1, y = PC3, color = cell))+
    geom_point(size = 2, alpha =0.8) +
    geom_label_repel(aes(label = cell)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    xlab(str_c("PC1", " (", eigen_percent[1], "%)")) +
    ylab(str_c("PC3", " (", eigen_percent[3], "%)")) +
    ggtitle("PhenoGraph Annotation PCA") +
    scale_color_manual(values = color_clusters) +
    theme(legend.position = 'none')


# output the higm dimensional results into csv for further analysis
dir_output = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/New_Data/LT_Old10x_2018_08.04.2020/preprocessed_data/"
setwd(dir_output)

write_csv(cldf, 'LT_old_cldf_08.12.2020.csv')








