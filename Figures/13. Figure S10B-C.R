####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(ggplot2)
library(scales)

####2. LOAD ARTERY DATA####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(getActiveDocumentContext()$path),
                             "/Only arteries/Processed_Artery.rds"))

####3. FIGURE S10B####

#Plot annotated clusters as DimPlot (FIGURE S10B)
DimPlot(object = Seurat.obj,repel = T,
        label = T,label.box = T,label.size = 5,
        cols = c("gold","grey","red","violet"),
        pt.size = 1.5,reduction = "umap") + NoLegend()


####4. FIGURE S10C####

#Check distribution of key cell markers as DotPlot (FIGURE S10C)
DotPlot(Seurat.obj,
        c("aSMA",
          "CD31","VCAM1",
          "CD45","CD68","C5aR1","CD11b","CD163","CD3","CD20"
        ),
        scale = T,
        col.min = 0.6,
        col.max = 1.5,
        dot.min = F)+ scale_size(range = c(15, 15)) + guides(size = "none") +
  scale_color_gradient(low = "white",
                       high = "blue",
                       limits = c(0.6, 1.5),
                       oob = squish) + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())


####5. Session Information#####

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
#   [1] scales_1.4.0       ggplot2_3.5.2      Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0          
# [6] rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.6           
# [5] magrittr_2.0.3         RcppAnnoy_0.0.22       matrixStats_1.5.0      ggridges_0.5.4        
# [9] compiler_4.3.2         spatstat.geom_3.4-1    png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3        fastmap_1.2.0         
# [17] promises_1.3.3         purrr_1.0.4            jsonlite_2.0.0         goftest_1.2-3         
# [21] later_1.4.2            spatstat.utils_3.1-4   irlba_2.3.5.1          parallel_4.3.2        
# [25] cluster_2.1.6          R6_2.6.1               ica_1.0-3              stringi_1.8.7         
# [29] RColorBrewer_1.1-3     spatstat.data_3.1-6    reticulate_1.42.0      parallelly_1.45.0     
# [33] spatstat.univar_3.1-3  lmtest_0.9-40          scattermore_1.2        Rcpp_1.0.14           
# [37] tensor_1.5             future.apply_1.11.0    zoo_1.8-12             sctransform_0.4.2     
# [41] httpuv_1.6.16          Matrix_1.6-5           splines_4.3.2          igraph_2.1.4          
# [45] tidyselect_1.2.1       dichromat_2.0-0.1      abind_1.4-8            codetools_0.2-19      
# [49] spatstat.random_3.4-1  miniUI_0.1.1.1         spatstat.explore_3.4-3 listenv_0.9.0         
# [53] lattice_0.22-7         tibble_3.2.1           plyr_1.8.9             withr_3.0.2           
# [57] shiny_1.10.0           ROCR_1.0-11            Rtsne_0.17             future_1.58.0         
# [61] fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.1-11   
# [65] pillar_1.10.2          KernSmooth_2.23-26     plotly_4.10.4          generics_0.1.3        
# [69] RcppHNSW_0.6.0         globals_0.18.0         xtable_1.8-4           glue_1.8.0            
# [73] lazyeval_0.2.2         tools_4.3.2            data.table_1.17.4      RSpectra_0.16-1       
# [77] RANN_2.6.1             dotCall64_1.2          cowplot_1.1.3          grid_4.3.2            
# [81] tidyr_1.3.1            nlme_3.1-168           patchwork_1.3.0        cli_3.6.5             
# [85] spatstat.sparse_3.1-0  spam_2.11-1            viridisLite_0.4.2      dplyr_1.1.4           
# [89] uwot_0.2.3             gtable_0.3.6           digest_0.6.37          progressr_0.15.1      
# [93] ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1     
# [97] lifecycle_1.0.4        httr_1.4.7             mime_0.13              MASS_7.3-60   
