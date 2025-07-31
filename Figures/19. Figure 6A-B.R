####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(readxl)
library(dplyr)
library(pheatmap)
library(ggplot2)

####2. IMPORT DATA####

#Set directory
setwd(dirname(getActiveDocumentContext()$path))

#Load Processed data
Seurat.obj <- readRDS("Only skin/Processed_Skin.rds")

#Load data file with artery annotation data
Skin <- read_xlsx("Data/CCE_Skin_Arteries.xlsx")

#Add name-corrected artery annotation data to Seurat object
Skin$`Analysis Region` <- factor(Skin$`Analysis Region`)
levels(Skin$`Analysis Region`) <-  c("Papillary",
                                     "Reticular Deep",
                                     "Reticular Sup.")
Skin$`Analysis Region` <- factor(Skin$`Analysis Region`,levels = c("Reticular Deep",
                                                                   "Reticular Sup.",
                                                                   "Papillary"))
Seurat.obj$Arteries <- Skin$`Analysis Region`


####3. FIGURE 6A####

#Plot Immune sub-cluster immune frequencies by artery type

#Extract immune marker and artery information
cluster_counts <-  data.frame("Artery" = Seurat.obj$Arteries,
                              "CD68" = Seurat.obj$CD68.Pos,
                              "CD3" = Seurat.obj$CD3.Pos,
                              "CD20" = Seurat.obj$CD20.Pos)

#Calculate frequencies by artery type
freq_df <- cluster_counts %>%
  group_by(Artery) %>%
  summarise(
    CD68 = mean(CD68),
    CD3 = mean(CD3),
    CD20 = mean(CD20)
  )

#Re-order data for heatmap:
mat <- freq_df[,2:4]
rownames(mat) <- freq_df$Artery

#Plot Heatmap (FIGURE 6A)
pheatmap(mat,
         cluster_rows = F,
         cluster_cols = F,
         cellwidth = 80,
         cellheight = 60,
         display_numbers = TRUE,
         number_format = "%.2f",
         scale = "none",
         show_rownames = T,
         border_color = NA,
         breaks = seq(0,0.2, by =0.01),
         color = colorRampPalette(c("white","lightyellow","orange","red",
                                    "brown"))(length(seq(0,0.2, by =0.01))),
         angle_col = 315)


####4. FIGURE 6B####

#Plot Associations between immune cells in CCE skin arteries and C5aR1 and CD11b

#Extract immune marker data from Seurat object
Immune_counts <-  data.frame("Clusters" = Idents(Seurat.obj),
                            "CD3" = Seurat.obj$CD4.Pos,
                            "CD68" = Seurat.obj$CD68.Pos,
                            "C5aR1" = Seurat.obj$C5aR1.Pos,
                            "CD11b" = Seurat.obj$CD11b.Pos)

#Filter only immune cells:
Immune_counts <- Immune_counts[ Immune_counts$Clusters == "T cells" |
                                  Immune_counts$Clusters == "Macrophages", ]

#Run fisher test for p.value and O.R. calculation
#Set comparisons
marker_pairs <- list(
  c("C5aR1", "CD3"),c("C5aR1", "CD68"),
  c("CD11b", "CD3"),c("CD11b", "CD68")
)

#Run fisher test for each pair
fisher_results <- lapply(marker_pairs, function(pair) {
  
  #Extract counts per pair
  x <- Immune_counts[[pair[1]]]
  y <- Immune_counts[[pair[2]]]
  
  #Do counts table
  tab <- table(x, y)
  
  #Run Fisher's test on table
  test <- fisher.test(tab)
  
  #Store results
  data.frame(
    Marker1 = pair[1],
    Marker2 = pair[2],
    OR = test$estimate,
    p.value = test$p.value
  )
})

#Combine results into one data frame and adjust p.value
fisher_df <- do.call(rbind, fisher_results)
fisher_df$p.adj <- p.adjust(fisher_df$p.value, method = "bonferroni")

#Filter only relevant comparisons for the plot
plot_df <- fisher_df %>%
  filter(Marker1 %in% c("C5aR1", "CD11b") & 
           Marker2 %in% c("CD3", "CD68"))


#Create bubble plot
ggplot(plot_df, aes(x = Marker1, y = Marker2)) +
  geom_point(aes(size = OR, fill = p.adj, alpha = log2(OR+1)), shape = 21, color = "black") +
  scale_size_continuous(range = c(3, 12)) +
  scale_fill_gradient(low = "red", high = "white", name = "Adjusted\np-value") +
  theme_minimal() +
  labs(
    x = "Complement Markers",
    y = "Immune Markers",
    size = "Odds Ratio",
    title = "Complement Markers vs Immune Markers (CCE)"
  ) +
  theme(
    axis.text.x = element_text(angle = 0,size = 10, 
                               hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold",size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

####5. Session Information####

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
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Paris
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.5.2      pheatmap_1.0.12    dplyr_1.1.4        readxl_1.4.5      
# [5] Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0           rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3         
# [4] rlang_1.1.6            magrittr_2.0.3         RcppAnnoy_0.0.22      
# [7] matrixStats_1.5.0      ggridges_0.5.4         compiler_4.3.2        
# [10] spatstat.geom_3.4-1    png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3       
# [16] fastmap_1.2.0          labeling_0.4.3         promises_1.3.3        
# [19] purrr_1.0.4            jsonlite_2.0.0         goftest_1.2-3         
# [22] later_1.4.2            spatstat.utils_3.1-4   irlba_2.3.5.1         
# [25] parallel_4.3.2         cluster_2.1.6          R6_2.6.1              
# [28] ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3    
# [31] spatstat.data_3.1-6    reticulate_1.42.0      parallelly_1.45.0     
# [34] spatstat.univar_3.1-3  cellranger_1.1.0       lmtest_0.9-40         
# [37] scattermore_1.2        Rcpp_1.0.14            tensor_1.5            
# [40] future.apply_1.11.0    zoo_1.8-12             sctransform_0.4.2     
# [43] httpuv_1.6.16          Matrix_1.6-5           splines_4.3.2         
# [46] igraph_2.1.4           tidyselect_1.2.1       dichromat_2.0-0.1     
# [49] abind_1.4-8            codetools_0.2-19       spatstat.random_3.4-1 
# [52] miniUI_0.1.1.1         spatstat.explore_3.4-3 listenv_0.9.0         
# [55] lattice_0.22-7         tibble_3.2.1           plyr_1.8.9            
# [58] withr_3.0.2            shiny_1.10.0           ROCR_1.0-11           
# [61] Rtsne_0.17             future_1.58.0          fastDummies_1.7.5     
# [64] survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.1-11   
# [67] pillar_1.10.2          KernSmooth_2.23-26     plotly_4.10.4         
# [70] generics_0.1.3         RcppHNSW_0.6.0         scales_1.4.0          
# [73] globals_0.18.0         xtable_1.8-4           glue_1.8.0            
# [76] lazyeval_0.2.2         tools_4.3.2            data.table_1.17.4     
# [79] RSpectra_0.16-1        RANN_2.6.1             dotCall64_1.2         
# [82] cowplot_1.1.3          grid_4.3.2             tidyr_1.3.1           
# [85] nlme_3.1-168           patchwork_1.3.0        cli_3.6.5             
# [88] spatstat.sparse_3.1-0  spam_2.11-1            viridisLite_0.4.2     
# [91] uwot_0.2.3             gtable_0.3.6           digest_0.6.37         
# [94] progressr_0.15.1       ggrepel_0.9.6          htmlwidgets_1.6.4     
# [97] farver_2.1.2           htmltools_0.5.8.1      lifecycle_1.0.4       
# [100] httr_1.4.7             mime_0.13              MASS_7.3-60           