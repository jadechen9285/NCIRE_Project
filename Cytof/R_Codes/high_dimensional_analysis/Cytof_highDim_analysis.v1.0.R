# ---------------UCSF_VA_NCIRE PROJECT: CyTOF single cell data analysis-----------
#   Jianhong Chen
#   07/06/2020
#
#   Pipeline II:
#  High Dim. Analyzing preprocessed & gated data
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

# GUI
library(cytofkit)

# Utilitise
library(shiny)

# ------------------------------ Load Pre-processed Data-----------------------------

# -- read fcs raw (CD45 gated) data
dir_raw = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/New_Data/LT_Cytof_Blue_07.29.2020/preprocessed_data/"
setwd(dir_raw)

fcs_files = list.files()
fcs_raw = flowCore::read.flowSet(fcs_files, transformation = FALSE,
                       truncate_max_range = FALSE)

# -- load metadata for the study
dir_meta = '/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/metadata/LT_old_03.26.2018/'
setwd(dir_meta)

# md info: sex: 0 = M, 1 = F; ht: unit in cm; wt: unit in kg
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
panel$target= str_replace_all(panel$target, "-", "_")
panel$marker = str_replace_all(panel$marker, "-", "_")

marker_table = read_excel("Cytof_Markers_table_07.02.2020.xlsx")
marker_table$marker = str_replace_all(marker_table$marker, "-", "_")
marker_table = marker_table[order(marker_table$marker), ]


# ------------------------------ Processing  Data--------------------------------
# define color of conditions
color_conditions = c("#6A3D9A", "#FF7F00") 
names(color_conditions) = levels(md$condition)
color_smokers = c("#66C2A5", "#FC8D62")
names(color_smokers) = levels(md$smoker)
color_LV = brewer.pal(n = 3, "Dark2")
names(color_LV) = levels(md$RV_TLC_cate)

# extract fcs files columnames & labels
panel_fcs = pData(parameters(fcs_raw[[1]]))
panel_fcs$desc = str_replace_all(panel_fcs$desc, "-", "_")

# define target markers to study 
# NOTE: no functional markers in this study
( lineage_markers = panel$marker[panel$lineage == 1] )

# spot check
all(lineage_markers %in% panel_fcs$desc)

# -- Data transformation (arcsinh transformation)
# Why arcsinh transformation?
# Our current mathematical laws don't allow taking the log of a negative number, 
# so different transforms are used that can accommodate negative numbers while 
# also displaying data in a log-like fashion.
# 
# X' = arcsinh(X/cofactor) 
fcs = fsApply(fcs_raw, function(x, cofactor = 5){ 
  colnames(x) = panel_fcs$desc
  expr = exprs(x)
  expr = asinh(expr[, lineage_markers] / cofactor )
  exprs(x) = expr
  x
})

expr = fsApply(fcs, exprs)
# generate sample IDs corresponding to each cell in the 'exprs' matrix
sample_ids = rep(md$sample_id, fsApply(fcs_raw, nrow))
# create an alternative expr version for plotting only (range from 0 to 1)
rng = colQuantiles(expr, probs = c(0.01, 0.99))
expr01 = t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] = 0
expr01[expr01 > 1] = 1

# ------------------------------ EDA & Diagnostic Plots-----------------------------
# num. of events per sample in barplot
data.frame("sample_id" = md$sample_id, "condition" = md$condition,
           "num_event" = fsApply(fcs_raw, nrow) ) %>%
  ggplot(aes(x = reorder(sample_id, -num_event),
             y = num_event, group = condition, fill = condition)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label=num_event),vjust = -1, fontface = "bold", size=2.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.title.x = element_blank())

# marker intensity density plot
ggdf = data.frame(sample_id = sample_ids, expr ) %>%
  melt(id.var = "sample_id", value.name = "expression",
       variable.name = "antigen")

mm = match(ggdf$sample_id, md$sample_id)
ggdf = mutate(ggdf, 
              "condition" = md$condition[mm],
              "smoker" = md$smoker[mm],
              "RV_TLC_cate" = md$RV_TLC_cate[mm])

ggplot(ggdf, aes(x = expression, color = RV_TLC_cate, group = sample_id)) +
  geom_density()+
  facet_wrap(~antigen, nrow = 6, scales = "free") +
  theme_bw() +
  ggtitle("Markers Density") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 7), 
        axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = color_LV)

# -- Complute Median intensity per sample per marker
expr_median_sample_tbl = data.frame(sample_id = sample_ids, expr) %>%
  group_by(sample_id) %>%
  summarize_all(median)
# rang 0-1 data is used for visualization only
expr01_median_sample_tbl = data.frame(sample_id = sample_ids, expr01) %>%
  group_by(sample_id) %>%
  summarize_all(median)

expr_median_sample = t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) = expr_median_sample_tbl$sample_id

expr01_median_sample = t(expr01_median_sample_tbl[, -1])
colnames(expr01_median_sample) = expr_median_sample_tbl$sample_id

# -- PCA
# projected the original matrixed into PCA space
mds = limma::plotMDS(expr_median_sample, plot = FALSE)
pca.res = prcomp(t(expr_median_sample), scale. = F)
eigen_percent = round(100*pca.res$sdev^2 / sum(pca.res$sdev^2), 1)

ggdf = data.frame(MDS1 = mds$x, MDS2 = mds$y, 
                  PC1 = pca.res$x[, 1], PC2 = pca.res$x[, 2], 
                  PC3 = pca.res$x[, 3], PC4 = pca.res$x[, 4],
                  sample_id = colnames(expr_median_sample))
mm = match(ggdf$sample_id, md$sample_id)
ggdf = mutate(ggdf, 
              "condition" = md$condition[mm],
              "smoker" = md$smoker[mm],
              "RV_TLC_cate" = md$RV_TLC_cate[mm] )

ggplot(ggdf, aes(x = MDS1, y = MDS2, color = condition))+
  geom_point(size = 2, alpha =0.8) +
  geom_label_repel(aes(label = sample_id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  xlab(str_c("PC1", " (", eigen_percent[1], "%)")) +
  ylab(str_c("PC2", " (", eigen_percent[2], "%)")) +
  ggtitle("Markers Median Value of each Sample in PCA") +
  scale_color_manual(values = color_conditions) + 
  theme(plot.title = element_text(size = 20, hjust = 0.5))

ggplot(ggdf, aes(x = MDS1, y = MDS2, color = RV_TLC_cate))+
  geom_point(size = 2, alpha =0.8) +
  geom_label_repel(aes(label = sample_id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  xlab(str_c("PC1", " (", eigen_percent[1], "%)")) +
  ylab(str_c("PC2", " (", eigen_percent[2], "%)")) +
  ggtitle("Markers Median Value of each Sample in PCA") +
  scale_color_manual(values = brewer.pal(n = 3, "Dark2")) + 
  theme(plot.title = element_text(size = 20, hjust = 0.5))


                      
# -- Heatmap viz for median intensity

condition_anno = data.frame("RV/TLC" =md$RV_TLC_cate)
# rownames(condition_anno) = md$sample_id
col_anno = HeatmapAnnotation(df = condition_anno)
# row_anno = rowAnnotation(df = condition_anno)

Heatmap(expr01_median_sample, name = "sample01", 
        col = colorRampPalette(rev(brewer.pal(n=9, name = "RdYlBu")))(100),
        row_title = "Marker", column_title = "Marker Median intensity per sample",
        cluster_rows = TRUE, cluster_columns = TRUE,
        border = TRUE, rect_gp = gpar(col = "white", lwd = 2),
        column_names_rot = 45, top_annotation = col_anno,
        layer_fun = function(j, i, x, y, width, height, fill) {
          # print out intensity text on each cell
          grid.text(sprintf("%.2f", 
                            pindex(expr01_median_sample, i, j)), 
                    x, y, gp = gpar(fontsize = 10))}
)


# -- Dimensional Reduction: tSNE & UMAP viz

## downsample original data 
# step 1: check duplicattion
dups = which(!duplicated(expr[, lineage_markers]))

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
dsamp_expr = expr[dsamp_inds, lineage_markers]

# -- Run t-SNE 
set.seed(78)
# perplexity resonable range = 5-50 by the authors of t-SNE
tsne_out = Rtsne(dsamp_expr, perplexity = 100, check_duplicates = FALSE,
                 pca = FALSE, verbose = TRUE)

# -- Run UMAP
set.seed(78)
# n_neighbors: small >> local details; large >> global details;
# min_dist: small >> densely packed; large >> more spread out
umap_out = uwot::umap(dsamp_expr, n_neighbors = 15, n_components = 2,
                      learning_rate = 1, min_dist = 0.01,
                      init = "random", n_epochs = 500, scale = FALSE,
                      verbose = TRUE)

# create DR dataframe that contains both tSNE and UMAP
drdf = data.frame(X = c(tsne_out$Y[, 1],umap_out[, 1]),
                  Y = c(tsne_out$Y[, 2],umap_out[, 2]),
                  DR = rep(c("tSNE", "UMAP"), c(nrow(tsne_out$Y), nrow(umap_out))),
                  dsamp_expr) %>%
  mutate('sample_id' = rep(sample_ids[dsamp_inds], 2))
mm = match(drdf$sample_id, md$sample_id)
drdf = mutate(drdf,
             "condition" = md$condition[mm],
             "smoker" = md$smoker[mm],
             "patient_id" = md$patient_id[mm],
             "gold" = md$gold[mm],
             "sex" = md$sex[mm],
             "RV_TLC_cate" = md$RV_TLC_cate[mm])

# check batch effect
ggplot(drdf, aes(x = X, y = Y, color = patient_id)) +
  geom_point(size = 0.8, alpha = 0.8) +
  facet_wrap(~DR, scales = "free") +
  theme_bw() +
  scale_color_manual(values = brewer.pal(n=11, name = "Spectral"))

ggplot(drdf, aes(x = X, y = Y, color = CD3)) +
  geom_point(size = 0.8, alpha = 0.8) +
  facet_wrap(~DR, scales = "free") +
  theme_bw() +  
  scale_color_gradientn("CD3", colours = colorRampPalette(rev(
    brewer.pal(n=11, name = "Spectral")))(50))

# plot tSNE color by smoker status
ggplot(drdf, aes(x = X, y = Y, color = condition)) +
  geom_point(size = 0.8, alpha = 0.5) +
  facet_wrap(~DR, scales = "free") +
  theme_bw() +
  scale_color_manual(values = color_conditions)

# # # ------ test with CATALYST UMAp wrapper function
# daf = daFrame(fcs, panel, md, cols_to_use = lineage_markers)
# 
# daf = runDR(daf, "TSNE", cols_to_use = lineage_markers, rows_to_use = 1000)
# daf = runDR(daf, "UMAP", cols_to_use = lineage_markers, rows_to_use = 2000)
# 
# plotDR(daf, "TSNE", color_by = "CD3")
# plotDR(daf, "UMAP", color_by = "CD3")

# ------------------------------ Clustering -------------------------------------

# -- Cluster method I: FloSOM & metaclustering 
set.seed(78)
# (1) Read input fcs files
fsom = FlowSOM::ReadInput(fcs, transform = FALSE, scale = FALSE)
# (2) build self-organizing map(SOM)
som = FlowSOM::BuildSOM(fsom,colsToUse = lineage_markers)
# (3) build minimum spanning tree (MST)
mst = FlowSOM::BuildMST(som, silent = FALSE, tSNE = TRUE)
PlotStars(mst, main = "Spanning Tree Visualization of the Nodes")
# (4) Metaclustering into 20 clusters with ConsensusClusterPlus
codes = som$map$codes
plot_outdir = "consensus_plot"
nmc = 20

# ConsensusClustering SOM grids (default: 10x10)
mc = ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100,
                          pItem = 0.9, pFeature = 1, title = plot_outdir,
                          plot= "png", clusterAlg = "hc", 
                          innerLinkage = "average",finalLinkage = "average", 
                          distance = "euclidean", seed = 78)

# get cluster is for each cell
code_clustering1 = mc[[nmc]]$consensusClass
cell_clustering1 = code_clustering1[som$map$mapping[, 1]]

# -- Clustering Method II: RphenoGraph using resample data along tSNE viz.
set.seed(78)
phenoGraph.res = Rphenograph(dsamp_expr, k = 40)
modularity(phenoGraph.res[[2]])
phenoGraph.cluster = membership(phenoGraph.res[[2]])
color_pheno = distinctColorPalette(k = length(unique(phenoGraph.cluster)))

cldf = drdf %>%
  mutate("phenoGraph_cluster" = rep(factor(phenoGraph.cluster), 2),
         "som_cluster" = rep(factor(cell_clustering1[dsamp_inds], 
                                    levels = 1:nmc), 2))

# -- generate heatmap plots: illustrate median marker expression in each cluster
color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", 
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

# Complexheat viz the clustering result:
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
  
  # calculate the median expression
  expr_median = data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_all(median)
  expr01_median = data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarise_all(median)
  
  # calculate cluster frequencies
  clustering_table = as.numeric(table(cell_clustering))
  
  # This clustering is based on the markers that were used for the main clustering
  d = dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows = hclust(d, "average")
  
  expr_heat = as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) = expr01_median$cell_clustering
  
  labels_row =  paste0(rownames(expr_heat), " (",
                       round(clustering_table / sum(clustering_table)*100, 2),
                       "%)")
  labels_col = colnames(expr_heat)
  
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
  color = colorRampPalette(rev(brewer.pal(n=9, name = "RdYlBu")))(100)
  
  pheatmap::pheatmap(expr_heat, color = color,
           cluster_cols = TRUE, cluster_rows = cluster_rows,
           label_col = labels_col, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 4,
           annotation_row = annotation_row, annotation_colors = annotation_colors,
           annotation_legend = annotation_legend)
}

# FlowSOM cluster heatmap
plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers], 
                                expr01 = expr01[, lineage_markers],
                                cell_clustering = cell_clustering1,
                                color_clusters = color_clusters)

# phenoGraph cluster heatmap
plot_clustering_heatmap_wrapper(expr = expr[dsamp_inds, lineage_markers],
                                expr01 = expr01[dsamp_inds, lineage_markers],
                                cell_clustering = 
                                  cldf$phenoGraph_cluster[1:length(dsamp_inds)],
                                color_clusters = color_pheno)

# # -- Visual Clustering result with UMAP
##!!! from this point: Use UMAP and phenoGraph only
cldf = dplyr::filter(cldf, DR == "UMAP") %>%
  dplyr::rename('UMAP1' = X, 'UMAP2' = Y  )
pheno_p =
  ggplot(cldf, aes(x = UMAP1, y = UMAP2, color = phenoGraph_cluster)) +
  geom_point(size = 0.8) +
  theme_bw() +
  ggtitle("Clustering by PhenoGraph") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

pheno_p + facet_wrap(~smoker) + 
  ggtitle("Clustering by PhenoGraph by Smoker")
pheno_p + facet_wrap(~condition) + 
  ggtitle("Clustering by PhenoGraph by Condition")
pheno_p + facet_wrap(~RV_TLC_cate) +
  ggtitle("Clustering by PhenoGraph by Lung Volume")
pheno_p + facet_wrap(~sample_id) +
  ggtitle("Clustering by PhenoGraph by Sample ID")

# plot t-SNE color by flowsom metaclustering
som_p = ggplot(cldf, aes(x = UMAP1, y = UMAP2, color = som_cluster)) +
  geom_point(size = 0.8) +
  theme_bw() +
  ggtitle("FlowSOM MetaClustering") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

som_p + facet_wrap(~smoker)

# plot all markers intensity(facet)
cldf %>%
  pivot_longer(cols = CD68:CD1c, names_to = "marker",
               values_to = "intensity") %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = intensity)) +
  geom_point(size = 0.8) +
  facet_wrap(~marker) +
  theme_bw() +
  scale_color_gradientn("Maker", colours = colorRampPalette(rev(
    brewer.pal(n=11, name = "Spectral")))(50))

# plot all clusters separately
cldf %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = CD206)) +
  geom_point(size = 0.8) +
  facet_wrap(~phenoGraph_cluster, nrow = 3) +
  theme_bw() +
  scale_color_gradientn("Maker", colours = colorRampPalette(rev(
    brewer.pal(n=11, name = "Spectral")))(50))

  
# ------------------- Annotating Clusters Results ----------------------------------------
# -- load the annotated label table
cluster_merging_md = read_excel("LT_manual_annotation_07.02.2020.xlsx")

# add the merged & annotated cell cluster into the dataframe
mm = match(phenoGraph.cluster, cluster_merging_md$original_cluster)
merge_cluster = cluster_merging_md$new_cluster[mm]

mm = match(cldf$phenoGraph_cluster, cluster_merging_md$original_cluster)
merge_cluster_sub = cluster_merging_md$new_cluster[mm]

cldf = mutate(cldf, "cell" = factor(merge_cluster_sub))

# -- viz on the annotated cell clusters
plot_clustering_heatmap_wrapper(expr = expr[dsamp_inds, lineage_markers],
                                expr01 = expr01[dsamp_inds, lineage_markers],
                                cell_clustering = cldf$phenoGraph_cluster,
                                color_clusters = color_clusters,
                                cluster_merging = cluster_merging_md)

plot_clustering_heatmap_wrapper(expr = expr[dsamp_inds, lineage_markers],
                                expr01 = expr01[dsamp_inds, lineage_markers],
                                cell_clustering = merge_cluster,
                                color_clusters = color_clusters)

cell_p = 
ggplot(cldf, aes(x = UMAP1, y = UMAP2, color = cell)) +
  geom_point(size = 0.8) +
  theme_bw() +
  ggtitle("Annotated Cell Clusters") + 
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))

cell_p + facet_wrap(~condition)
cell_p + facet_wrap(~smoker)
cell_p + facet_wrap(~RV_TLC_cate)
cell_p + facet_wrap(~sample_id)

# output the higm dimensional results into csv for further analysis
dir_proc = "/Users/jianhongchen/Desktop/Career/UCSF_NCIRE_VA_Research_Statistician/Projects/06_CyTOF/Data/processed/LT_old_set_analyzed_07.09.2020/"
setwd(dir_proc)

write_csv(cldf, 'LT_old_processed_07.09.2020.csv')








