####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(ggplot2)
library(scales)

####2. LOAD PROCESSED TONSIL DATA####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(getActiveDocumentContext()$path),
                             "/Data/Processed_Tonsils.rds"))


####3. FIGURE S3E ####

DimPlot(object = Seurat.obj,repel = T,
        label = T,label.box = T,label.size = 5,
        cols = c("red","gold","mediumseagreen","skyblue","purple","pink"),
        pt.size = 2,reduction = "umap") + NoLegend()

####4. FIGURE S3F ####

FeaturePlot(object = Seurat.obj,
            pt.size = 1, repel = T,
            reduction = "umap",label = T,
            features = c("CD31","CD20","CD3",
                         "CD4","CD8","CD68"),
            min.cutoff = "q40",max.cutoff = "q99",
            ncol = 3)

####5. FIGURE S3G ####

DotPlot(Seurat.obj,
        c("CD3",
          "CD20","Ki67",
          "CD31",
          "CD68","C1q","C5aR1"
        ),
        scale = T,
        col.min = 0.1,
        col.max = 1,
        dot.min = F) + 
  scale_size(range = c(15, 15)) + 
  guides(size = "none") +
  scale_color_gradient(low = "white",
                       high = "blue",
                       limits = c(0.1, 1),
                       oob = squish)


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
#   [1] scales_1.4.0       ggplot2_3.5.2      Seurat_5.3.0       SeuratObject_5.1.0
# [5] sp_2.2-0           rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.3        
# [4] spatstat.utils_3.1-4   farver_2.1.2           rmarkdown_2.29        
# [7] vctrs_0.6.5            ROCR_1.0-11            spatstat.explore_3.4-3
# [10] base64enc_0.1-3        htmltools_0.5.8.1      Formula_1.2-5         
# [13] sctransform_0.4.2      parallelly_1.45.0      KernSmooth_2.23-26    
# [16] phenoptr_0.3.2         htmlwidgets_1.6.4      ica_1.0-3             
# [19] plyr_1.8.9             plotly_4.10.4          zoo_1.8-12            
# [22] igraph_2.1.4           mime_0.13              lifecycle_1.0.4       
# [25] iterators_1.0.14       pkgconfig_2.0.3        Matrix_1.6-5          
# [28] R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.1-11   
# [31] future_1.58.0          shiny_1.10.0           digest_0.6.37         
# [34] fdrtool_1.2.18         colorspace_2.1-1       patchwork_1.3.0       
# [37] tensor_1.5             RSpectra_0.16-1        irlba_2.3.5.1         
# [40] Hmisc_5.2-3            labeling_0.4.3         progressr_0.15.1      
# [43] spatstat.sparse_3.1-0  httr_1.4.7             polyclip_1.10-7       
# [46] abind_1.4-8            compiler_4.3.2         withr_3.0.2           
# [49] glasso_1.11            htmlTable_2.4.3        backports_1.5.0       
# [52] psych_2.5.3            fastDummies_1.7.5      MASS_7.3-60           
# [55] corpcor_1.6.10         gtools_3.9.5           tools_4.3.2           
# [58] pbivnorm_0.6.0         foreign_0.8-90         lmtest_0.9-40         
# [61] httpuv_1.6.16          future.apply_1.11.0    goftest_1.2-3         
# [64] nnet_7.3-19            glue_1.8.0             quadprog_1.5-8        
# [67] nlme_3.1-168           promises_1.3.3         grid_4.3.2            
# [70] checkmate_2.3.2        Rtsne_0.17             cluster_2.1.6         
# [73] reshape2_1.4.4         generics_0.1.3         gtable_0.3.6          
# [76] spatstat.data_3.1-6    tidyr_1.3.1            data.table_1.17.4     
# [79] spatstat.geom_3.4-1    RcppAnnoy_0.0.22       ggrepel_0.9.6         
# [82] RANN_2.6.1             foreach_1.5.2          pillar_1.10.2         
# [85] stringr_1.5.1          spam_2.11-1            RcppHNSW_0.6.0        
# [88] later_1.4.2            splines_4.3.2          dplyr_1.1.4           
# [91] lattice_0.22-7         deldir_1.0-9           survival_3.8-3        
# [94] tidyselect_1.2.1       miniUI_0.1.1.1         pbapply_1.7-2         
# [97] knitr_1.50             gridExtra_2.3          scattermore_1.2       
# [100] stats4_4.3.2           xfun_0.52              qgraph_1.9.8          
# [103] matrixStats_1.5.0      stringi_1.8.7          lazyeval_0.2.2        
# [106] yaml_2.3.10            evaluate_1.0.3         codetools_0.2-19      
# [109] tibble_3.2.1           cli_3.6.5              uwot_0.2.3            
# [112] rpart_4.1.24           xtable_1.8-4           reticulate_1.42.0     
# [115] lavaan_0.6-19          dichromat_2.0-0.1      Rcpp_1.0.14           
# [118] spatstat.random_3.4-1  globals_0.18.0         png_0.1-8             
# [121] spatstat.univar_3.1-3  parallel_4.3.2         dotCall64_1.2         
# [124] jpeg_0.1-10            listenv_0.9.0          viridisLite_0.4.2     
# [127] ggridges_0.5.4         purrr_1.0.4            rlang_1.1.6           
# [130] cowplot_1.1.3          mnormt_2.1.1 
