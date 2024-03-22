rm(list = ls())
RNGversion("3.5.3")

library(readxl)
library(HDCytoData)
library(ggplot2)
library(CATALYST)
library(FlowSOM)
library(cowplot)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(scater)
library(diffcyt)

### 5) Data Import ####

url <- "https://zenodo.org/records/10039274/files"
md <- "PBMC8_metadata.xlsx"
download.file(file.path(url, md), destfile = md, mode = "wb")
md <- read_excel(md)
head(data.frame(md))

fs <- Bodenmiller_BCR_XL_flowSet()

panel <- "PBMC8_panel_v3.xlsx"
download.file(file.path(url, panel), destfile = panel, mode = "wb")
panel <- read_excel(panel)
head(data.frame(panel))
all(panel$fcs_colname %in% colnames(fs))

#### 7) Data organization #####

md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))
md$sample_id <- factor(md$sample_id, 
                       levels = md$sample_id[order(md$condition)])

sce <- prepData(fs, panel, md, features = panel$fcs_colname)

##### 8) Diagnostic plots #######

p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 6
p

n_cells(sce) 
plotCounts(sce, group_by = "sample_id", color_by = "condition")

pbMDS(sce, color_by = "condition", label_by = "sample_id")

plotExprHeatmap(sce, scale = "last",
                hm_pal = rev(hcl.colors(10, "YlGnBu")))

### 9) Marker ranking based on the non-redundancy score ####

plotNRS(sce, features = "type", color_by = "condition")

##### 10) Cell population identification #######
set.seed(1234)
sce <- cluster(sce, features = "type",
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)

plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta20", 
                bars = TRUE, perc = TRUE)

plotClusterExprs(sce, k = "meta20", features = "type")

plotMultiHeatmap(sce, 
                 hm1 = "type", hm2 = "pS6", k = "meta20", 
                 row_anno = FALSE, bars = TRUE, perc = TRUE)

sce <- runDR(sce, "TSNE", cells = 500, features = "type")
sce <- runDR(sce, "UMAP", cells = 1e3, features = "type")
plotDR(sce, "UMAP", color_by = "CD4")

p1 <- plotDR(sce, "TSNE", color_by = "meta20") + 
  theme(legend.position = "none")
p2 <- plotDR(sce, "UMAP", color_by = "meta20")
lgd <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")
plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))

plotDR(sce, "UMAP", color_by = "meta20", facet_by = "sample_id")

plotDR(sce, "UMAP", color_by = "meta20", facet_by = "condition")

plotCodes(sce, k = "meta20") #plotCodes not showing class labels

plotMultiHeatmap(sce, 
                 hm1 = "type", hm2 = "pS6", k = "som100", m = "meta20", 
                 row_anno = FALSE, col_anno = FALSE, bars = TRUE, perc = TRUE)

merging_table1 <- "PBMC8_cluster_merging1.xlsx"
download.file(file.path(url, merging_table1), 
              destfile = merging_table1, mode = "wb")
merging_table1 <- read_excel(merging_table1)
head(data.frame(merging_table1))

merging_table1$new_cluster <- factor(merging_table1$new_cluster, 
                                     levels = c("B-cells IgM+", "B-cells IgM-", "CD4 T-cells",
                                                "CD8 T-cells", "DC", "NK cells", "monocytes", "surface-"))

sce <- mergeClusters(sce, k = "meta20", 
                     table = merging_table1, id = "merging1")
plotDR(sce, "UMAP", color_by = "merging1")
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta20", m = "merging1")

plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "merging1")

merging_table2 <- "PBMC8_cluster_merging2_v3.xlsx"
download.file(file.path(url, merging_table2), 
              destfile = merging_table2, mode = "wb")
merging_table2 <- read_excel(merging_table2)
data.frame(merging_table2)

# convert to factor with merged clusters in desired order
merging_table2$new_cluster <- factor(
  merging_table2$new_cluster, 
  levels = levels(merging_table1$new_cluster))

# apply manual merging
sce <- mergeClusters(sce, k = "meta12", 
                     table = merging_table2, id = "merging2")
table(manual = cluster_codes(sce)[cluster_ids(sce), "merging2"],
      algorithm = cluster_codes(sce)[cluster_ids(sce), "meta8"] )

plot_grid(labels = c("A", "B"),
          plotDR(sce, "UMAP", color_by = "merging2"),
          plotDR(sce, "UMAP", color_by = "meta8"))

###### 11) DIF ANALYSIS #########

FDR_cutoff <- 0.05
plotAbundances(sce, k = "merging1", by = "sample_id")
plotAbundances(sce, k = "merging1", by = "cluster_id", shape_by = "patient_id")

ei <- metadata(sce)$experiment_info
(da_formula1 <- createFormula(ei, 
                              cols_fixed = "condition", 
                              cols_random = "sample_id"))

(da_formula2 <- createFormula(ei, 
                              cols_fixed = "condition", 
                              cols_random = c("sample_id", "patient_id")))
contrast <- createContrast(c(0, 1))

da_res1 <- diffcyt(sce, 
                   formula = da_formula1, contrast = contrast,
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "merging1", verbose = FALSE)

da_res2 <- diffcyt(sce, 
                   formula = da_formula2, contrast = contrast,
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "merging1", verbose = FALSE)
names(da_res1)

rowData(da_res1$res) 

table(rowData(da_res1$res)$p_adj < FDR_cutoff)

table(rowData(da_res2$res)$p_adj < FDR_cutoff)

topTable(da_res2, show_props = TRUE, format_vals = TRUE, digits = 2)

plotDiffHeatmap(sce, rowData(da_res2$res), all = TRUE, fdr = FDR_cutoff)

p <- plotMedExprs(sce, k = "merging1", 
                  facet_by = "cluster_id", shape_by = "patient_id")
p$facet$params$ncol <- 2
p

ds_formula1 <- createFormula(ei, cols_fixed = "condition")
ds_formula2 <- createFormula(ei, 
                             cols_fixed = "condition", cols_random = "patient_id")

ds_res1 <- diffcyt(sce, 
                   formula = ds_formula1, contrast = contrast,
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                   clustering_to_use = "merging1", verbose = FALSE)
table(rowData(ds_res1$res)$p_adj < FDR_cutoff)

ds_res2 <- diffcyt(sce, 
                   formula = ds_formula2, contrast = contrast,
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                   clustering_to_use = "merging1", verbose = FALSE)
table(rowData(ds_res2$res)$p_adj < FDR_cutoff)

topTable(ds_res2, top_n = 5, order_by = "cluster_id", 
         show_meds = TRUE, format_vals = TRUE, digits = 3)

plotDiffHeatmap(sce, rowData(ds_res2$res), top_n = 50, fdr = FDR_cutoff)

sce <- mergeClusters(sce, k = "meta20", id = "merging_all",
                     table = data.frame(old_cluster = seq_len(20), new_cluster = "all"))

p <- plotMedExprs(sce, features = "state", shape_by = "patient_id")
p$facet$params$ncol <- 3
p

# fit linear model
ds_res3 <- diffcyt(sce, 
                   formula = ds_formula1, contrast = contrast,
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                   clustering_to_use = "merging_all", verbose = FALSE)

# fit linear mixed model with patient ID as random effect
ds_res4 <- diffcyt(sce, 
                   formula = ds_formula2, contrast = contrast,
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                   clustering_to_use = "merging_all", verbose = FALSE)

table(rowData(ds_res3$res)$p_adj < FDR_cutoff)

table(rowData(ds_res4$res)$p_adj < FDR_cutoff)

topTable(ds_res4, top_n = 5, order_by = "p_adj",
         show_meds = TRUE, format_vals = TRUE, digits = 3)

plotDiffHeatmap(sce, rowData(ds_res4$res), all = TRUE)
