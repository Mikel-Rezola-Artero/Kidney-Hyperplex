####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(ggplot2)
library(scales)
library(dplyr)
library(pheatmap)
library(ggpubr)
library(tidyr)

####2. LOAD GLOMERULI DATA####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(getActiveDocumentContext()$path),
                             "/Only glomeruli/Processed_Glomeruli.rds"))


####3. FIGURE 2D-F####

#Complement deposits in Glomeruli across disease and cell types
LN <- DotPlot(Seurat.obj[, Seurat.obj$Disease == "LN"],
              features = c("C1q.Pos", "MBL.Pos", "C4d.Pos", 
                           "C3c.Pos", "C3d.Pos", "C9.Pos"),
              scale = FALSE) + 
  labs(title = "Complement Deposits in Lupus Nephritis") + 
  scale_color_gradient(low = "white", high = "blue", 
                       limits = c(0, 1), oob = squish) +
  scale_size(limits = c(15, 100),range = c(0.1,25))  + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))


DN <- DotPlot(Seurat.obj[, Seurat.obj$Disease == "Diabetes"],
              features = c("C1q.Pos", "MBL.Pos", "C4d.Pos", 
                           "C3c.Pos", "C3d.Pos", "C9.Pos"),
              scale = FALSE) + 
  labs(title = "Complement Deposits in Diabetes Nephritis") + 
  scale_color_gradient(low = "white", high = "blue", 
                       limits = c(0, 1), oob = squish) +
  scale_size(limits = c(15, 100),range = c(0.1,25))  + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))


CCE <- DotPlot(Seurat.obj[, Seurat.obj$Disease == "CCE"],
               features = c("C1q.Pos", "MBL.Pos", "C4d.Pos", 
                            "C3c.Pos", "C3d.Pos", "C9.Pos"),
               scale = FALSE) + 
  labs(title = "Complement Deposits in CCE") + 
  scale_color_gradient(low = "white", high = "blue", 
                       limits = c(0, 1), oob = squish) +
  scale_size(limits = c(15, 100),range = c(0,25))  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )

#Plot FIGURE 2D-F
LN / DN / CCE


####4. FIGURE 2G####

#Filter only Lupus slides' data
Lupus <- Seurat.obj[, Seurat.obj$Disease == "LN"]

#Filter Complement & Immune marker positivity, patient and glomeruli metadata 
cluster_counts <-  data.frame("Glomeruli" = Lupus$Glomeruli,
                              "Patient" = Lupus$slide,
                              "C1q" = Lupus$C1q.Pos,
                              "C3c" = Lupus$C3c.Pos,
                              "C3d" = Lupus$C3d.Pos,
                              "MBL" = Lupus$MBL.Pos,
                              "C4d" = Lupus$C4d.Pos,
                              "C9" = Lupus$C9.Pos,
                              "CD68" = Lupus$CD68.Pos,
                              "C5aR1" = Lupus$C5aR1.Pos,
                              "CD11b" = Lupus$CD11b.Pos,
                              "MPO" = Lupus$MPO.Pos,
                              "CD3" = Lupus$CD3.Pos,
                              "CD20" = Lupus$CD20.Pos)

#For each Glomeruli annotation and Patient pair,
#calculate frequency of positive cells
freq_df <- cluster_counts %>%
  group_by(Glomeruli,Patient) %>%
  summarise(
    C1q = mean(C1q),
    MBL = mean(MBL),
    C3c = mean(C3c),
    C3d = mean(C3d),
    C4d = mean(C4d),
    C9 = mean(C9),
    CD68 = mean(CD68),
    C5aR1 = mean(C5aR1),
    CD11b = mean(CD11b),
    MPO = mean(MPO),
    CD3 = mean(CD3),
    CD20 = mean(CD20)
  )

mat <- freq_df[,3:14]
#note that all features are in the same scale from 0 to 1, no scaling needed

#Plot Heatmap (FIGURE 2G)
heatmap.a <-pheatmap(mat,
                     cluster_rows = T,
                     cluster_cols = F,
                     cellwidth = 80,
                     cellheight = 4,  
                     clustering_distance_rows = "euclidean",
                     clustering_method = "complete",
                     scale = "none",
                     cutree_rows = 2,
                     show_rownames = F,
                     border_color = NA,
                     fontsize_col = 20,
                     breaks = seq(0,1, by =0.01),
                     color = colorRampPalette(c("lightgrey","yellow","red","brown"))(length(seq(0,1, by =0.01))),
                     angle_col = 315)


####5. FIGURE 2H####

#Extract clusters, run pairwise wilcox and plot volcano plot

#Extract clusters and merge with markers:
glomeruli.clust <- cbind(mat, 
                         cluster = cutree(heatmap.a$tree_row, k = 2))      


#Run wilcox test for each marker
results.wilcox <- compare_means(
  formula = value ~ cluster,
  data = glomeruli.clust %>%
    pivot_longer(cols = all_of(names(glomeruli.clust)[1:12]),
                 names_to = "marker", values_to = "value"),
  group.by = "marker",
  method = "wilcox.test",
  p.adjust.method = "bonferroni"
)

#Calculate median difference for each comparison
diff_table <- glomeruli.clust %>%
  group_by(cluster) %>%
  summarize(across(all_of(names(glomeruli.clust)[1:12]), median, na.rm = TRUE)) %>%
  pivot_longer(-cluster, names_to = "marker", values_to = "median") %>%
  pivot_wider(names_from = cluster, values_from = median, names_prefix = "cluster_") %>%
  mutate(
    median_diff = cluster_1 - cluster_2,
  )

#Merge adj.p.value and median_diff info
results <- data.frame(marker = results.wilcox$marker,
                      neg_log10_p = -log10(results.wilcox$p.adj),
                      median_diff = diff_table$median_diff,
                      significant = results.wilcox$p.adj < 0.05)

#Get top 5 markers by significance
top5 <- results[order(results$neg_log10_p,decreasing = T),][1:5,]

#Define colours
results <- results %>%
  mutate(
    category = case_when(
      marker %in% top5$marker ~ "top5",
      significant ~ "significant",
      TRUE ~ "not_significant"
    )
  )

#Volcano plot (FIGURE 2H)
ggplot(results, aes(x = median_diff, y = neg_log10_p)) +
  geom_point(aes(color = category, alpha = significant), size = 5) +
  geom_text(data = top5, aes(label = marker), 
            vjust = 1,hjust= -.2, size = 5,
            fontface = "bold") +
  scale_color_manual(
    values = c(
      "top5" = "red",
      "significant" = "blue",
      "not_significant" = "gray"
    )
  ) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  labs(
    x = "Median Diff. (Cluster 1 - Cluster 2)",
    y = "-log10(adj.p.value)"
  ) +
  theme_classic() + NoLegend() +
  xlim(0, 0.5)


####6. Session Information#####

sessionInfo()

# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3;  LAPACK version 3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
# [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Paris
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tidyr_1.3.1        ggpubr_0.6.0       pheatmap_1.0.12    dplyr_1.1.4       
# [5] ggplot2_3.5.2      Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0          
# [9] rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3         
# [4] rlang_1.1.6            magrittr_2.0.3         RcppAnnoy_0.0.22      
# [7] matrixStats_1.5.0      ggridges_0.5.4         compiler_4.3.2        
# [10] spatstat.geom_3.4-1    png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         stringr_1.5.1          crayon_1.5.3          
# [16] pkgconfig_2.0.3        fastmap_1.2.0          backports_1.5.0       
# [19] labeling_0.4.3         utf8_1.2.6             promises_1.3.3        
# [22] purrr_1.0.4            jsonlite_2.0.0         goftest_1.2-3         
# [25] later_1.4.2            spatstat.utils_3.1-4   broom_1.0.8           
# [28] irlba_2.3.5.1          parallel_4.3.2         cluster_2.1.6         
# [31] R6_2.6.1               ica_1.0-3              stringi_1.8.7         
# [34] RColorBrewer_1.1-3     spatstat.data_3.1-6    reticulate_1.42.0     
# [37] car_3.1-3              parallelly_1.45.0      spatstat.univar_3.1-3 
# [40] lmtest_0.9-40          scattermore_1.2        Rcpp_1.0.14           
# [43] tensor_1.5             future.apply_1.11.0    zoo_1.8-12            
# [46] sctransform_0.4.2      httpuv_1.6.16          Matrix_1.6-5          
# [49] splines_4.3.2          igraph_2.1.4           tidyselect_1.2.1      
# [52] dichromat_2.0-0.1      abind_1.4-8            codetools_0.2-19      
# [55] spatstat.random_3.4-1  miniUI_0.1.1.1         spatstat.explore_3.4-3
# [58] listenv_0.9.0          lattice_0.22-7         tibble_3.2.1          
# [61] plyr_1.8.9             withr_3.0.2            shiny_1.10.0          
# [64] ROCR_1.0-11            Rtsne_0.17             future_1.58.0         
# [67] fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7       
# [70] fitdistrplus_1.1-11    pillar_1.10.2          carData_3.0-5         
# [73] KernSmooth_2.23-26     plotly_4.10.4          generics_0.1.3        
# [76] RcppHNSW_0.6.0         scales_1.4.0           globals_0.18.0        
# [79] xtable_1.8-4           glue_1.8.0             lazyeval_0.2.2        
# [82] tools_4.3.2            data.table_1.17.4      RSpectra_0.16-1       
# [85] ggsignif_0.6.4         RANN_2.6.1             dotCall64_1.2         
# [88] cowplot_1.1.3          grid_4.3.2             nlme_3.1-168          
# [91] patchwork_1.3.0        Formula_1.2-5          cli_3.6.5             
# [94] spatstat.sparse_3.1-0  spam_2.11-1            viridisLite_0.4.2     
# [97] uwot_0.2.3             gtable_0.3.6           rstatix_0.7.2         
# [100] digest_0.6.37          progressr_0.15.1       ggrepel_0.9.6         
# [103] htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1     
# [106] lifecycle_1.0.4        httr_1.4.7             mime_0.13             
# [109] MASS_7.3-60 
