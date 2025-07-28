####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(ggplot2)
library(scales)

####2. LOAD GLOMERULI DATA####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(getActiveDocumentContext()$path),
                             "/Only arteries/Processed_Artery.rds"))


####3. FIGURE 3D-F####

#Complement deposits in Arteries across disease and cell types
LN <- DotPlot(Seurat.obj[, Seurat.obj$Disease == "LN"],
              features = c("C1q.Pos", "MBL.Pos", "C4d.Pos", 
                           "C3c.Pos", "C3d.Pos", "C9.Pos"),
              scale = FALSE) + 
  labs(title = "Complement Deposits in Lupus Nephritis") + 
  scale_color_gradient(low = "white", high = "blue", 
                       limits = c(0, 1), oob = squish) +
  scale_size(limits = c(5, 100),range = c(0.1,25))  + 
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
  scale_size(limits = c(5, 100),range = c(0.1,25))  + 
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
  scale_size(limits = c(5, 100),range = c(0,25))  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )

#PLOT FIGURE FIGURE 3D-F
LN / DN / CCE


####4. Session Information####

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
# [17] labeling_0.4.3         promises_1.3.3         purrr_1.0.4            jsonlite_2.0.0        
# [21] goftest_1.2-3          later_1.4.2            spatstat.utils_3.1-4   irlba_2.3.5.1         
# [25] parallel_4.3.2         cluster_2.1.6          R6_2.6.1               ica_1.0-3             
# [29] stringi_1.8.7          RColorBrewer_1.1-3     spatstat.data_3.1-6    reticulate_1.42.0     
# [33] parallelly_1.45.0      spatstat.univar_3.1-3  lmtest_0.9-40          scattermore_1.2       
# [37] Rcpp_1.0.14            tensor_1.5             future.apply_1.11.0    zoo_1.8-12            
# [41] sctransform_0.4.2      httpuv_1.6.16          Matrix_1.6-5           splines_4.3.2         
# [45] igraph_2.1.4           tidyselect_1.2.1       dichromat_2.0-0.1      abind_1.4-8           
# [49] codetools_0.2-19       spatstat.random_3.4-1  miniUI_0.1.1.1         spatstat.explore_3.4-3
# [53] listenv_0.9.0          lattice_0.22-7         tibble_3.2.1           plyr_1.8.9            
# [57] withr_3.0.2            shiny_1.10.0           ROCR_1.0-11            Rtsne_0.17            
# [61] future_1.58.0          fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7       
# [65] fitdistrplus_1.1-11    pillar_1.10.2          KernSmooth_2.23-26     plotly_4.10.4         
# [69] generics_0.1.3         RcppHNSW_0.6.0         globals_0.18.0         xtable_1.8-4          
# [73] glue_1.8.0             lazyeval_0.2.2         tools_4.3.2            data.table_1.17.4     
# [77] RSpectra_0.16-1        RANN_2.6.1             dotCall64_1.2          cowplot_1.1.3         
# [81] grid_4.3.2             tidyr_1.3.1            nlme_3.1-168           patchwork_1.3.0       
# [85] cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-1            viridisLite_0.4.2     
# [89] dplyr_1.1.4            uwot_0.2.3             gtable_0.3.6           digest_0.6.37         
# [93] progressr_0.15.1       ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2          
# [97] htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7             mime_0.13             
# [101] MASS_7.3-60     