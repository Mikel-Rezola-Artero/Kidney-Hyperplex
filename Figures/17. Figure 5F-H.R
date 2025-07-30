####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(readxl)
library(ggplot2)
library(ggbreak)

####2. IMPORT DATA####

#Set directory
setwd(dirname(getActiveDocumentContext()$path))

#Load Processed data
Seurat.obj <- readRDS("Only skin/Processed_Skin.rds")

#Load original data file with coordinate data
Skin <- read_xlsx("Data/CCE_Skin_Arteries.xlsx")


####3. FIGURE 5F####

#Plot skin arteries by cell cluster and artery type

#Ensure columns are numeric and calculate centroids
Skin$x <- (as.numeric(Skin$XMin) + as.numeric(Skin$XMax))/2
Skin$y <- (as.numeric(Skin$YMin) + as.numeric(Skin$YMax))/2

#Add cell type information
Skin$Clusters <- Idents(Seurat.obj)

#Fix names for aesthetics and transform into factors
Skin$`Analysis Region` <- factor(Skin$`Analysis Region`)
levels(Skin$`Analysis Region`) <-  c("Papillary Dermal Arteries",
                                     "Deep Reticular Dermal Arteries",
                                     "Superficial Reticular Dermal Arteries")
Skin$`Analysis Region` <- factor(Skin$`Analysis Region`,levels = c("Papillary Dermal Arteries",
                                                                   "Superficial Reticular Dermal Arteries",
                                                                   "Deep Reticular Dermal Arteries"))
Skin$Clusters <- factor(Skin$Clusters)

# Plot all arteries (FIGURE 5F)
ggplot(Skin, aes(x = x, y = y)) +
  geom_point(aes(color = Clusters), alpha = 0.5, size = 0.25) +
  stat_ellipse(aes(group = `Analysis Region`,
                   linetype = `Analysis Region`),alpha = 0.5) +
  ggtitle("Skin Arteries in CCE") + 
  scale_y_break(c(2500,7000)) + ylim(0,13800)+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.3,face = "bold"),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank()) +
  
  scale_color_manual(values = c("orange", "brown", "grey", "black", "blue")) + 
  guides(color = guide_legend(override.aes = list(size = 3)))


####4. FIGURE 5G####

#Plot C' deposits spatially

#Loop through all Complement markers to make plots

plot_list <- list()

#Markers
Complement <- c("C4d.Pos","C3d.Pos","C9.Pos","VCAM.1.Pos")

#Marker colors
Colors <- c("blue","red","black","brown")



for (i in 1:4) {
  
  Comp <- Complement[i]
  pos_color <- Colors[i]
  
  df <- Skin
  
  #Ensure columns are numeric
  df$x <- (as.numeric(df$XMin) + as.numeric(df$XMax))/2
  df$y <- (as.numeric(df$YMin) + as.numeric(df$YMax))/2
  
  df$Comp <- Seurat.obj@meta.data[, Comp]
  
  
  #Plot generation
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(color = as.numeric(Comp)), 
               alpha = 0.4, 
               size = 0.3) +
    theme_classic() +
    labs(title = gsub(".Pos","",Comp)) +
    scale_colour_gradientn(colours = c("#F2F2F2", pos_color)) + 
    guides(color = guide_legend(override.aes = list(size = 3))) + 
    ggbreak::scale_y_break(c(2500,7000)) + ylim(0,13800)+
    theme(legend.position="none",
          plot.title = element_text(
            hjust = 0.5,      
            face = "bold",
            size = 16 ),
          axis.text.y = element_text(size=7),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.text.y.right = element_blank())
  
  plot_list[[i]] <- p 
  
}

#Plot Complement deposits spatially (FIGURE 5G)
(plot_list[[1]] | plot_list[[2]]) /
  (plot_list[[3]] | plot_list[[4]])


####5. FIGURE 5H####

#C' deposits and VCAM-1 positivity

#Add artery type data to Seurat object
Seurat.obj$`Analysis Region` <- Skin$`Analysis Region`

#Rename and re-order levels for aesthetics
levels(Seurat.obj$`Analysis Region`) <- c("Papillary",
                                          "Reticular Sup.",
                                          "Reticular Deep")
Seurat.obj$`Analysis Region` <- factor(Seurat.obj$`Analysis Region`,
                                       levels=c("Reticular Deep",
                                                "Reticular Sup.",
                                                "Papillary"))
#Add Idents as metadata column
Seurat.obj$Clusters <- Idents(Seurat.obj)


#Plot complement and VCAM-1 positivity frequency per artery type (FIGURE 5K)
DotPlot(Seurat.obj,
        group.by = c("Analysis Region"),
        c("C4d.Pos","C3d.Pos","C9.Pos","VCAM.1.Pos"),
        scale = F) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

####6. Session Information####

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
#   [1] ggbreak_0.1.4      readxl_1.4.5       ggplot2_3.5.2      Seurat_5.3.0      
# [5] SeuratObject_5.1.0 sp_2.2-0           rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3         
# [4] rlang_1.1.6            magrittr_2.0.3         RcppAnnoy_0.0.22      
# [7] spatstat.geom_3.4-1    matrixStats_1.5.0      ggridges_0.5.4        
# [10] compiler_4.3.2         png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3       
# [16] fastmap_1.2.0          labeling_0.4.3         promises_1.3.3        
# [19] purrr_1.0.4            aplot_0.2.6            jsonlite_2.0.0        
# [22] goftest_1.2-3          later_1.4.2            spatstat.utils_3.1-4  
# [25] irlba_2.3.5.1          parallel_4.3.2         cluster_2.1.6         
# [28] R6_2.6.1               ica_1.0-3              stringi_1.8.7         
# [31] RColorBrewer_1.1-3     spatstat.data_3.1-6    reticulate_1.42.0     
# [34] parallelly_1.45.0      spatstat.univar_3.1-3  cellranger_1.1.0      
# [37] lmtest_0.9-40          scattermore_1.2        Rcpp_1.0.14           
# [40] tensor_1.5             future.apply_1.11.0    zoo_1.8-12            
# [43] sctransform_0.4.2      httpuv_1.6.16          Matrix_1.6-5          
# [46] splines_4.3.2          igraph_2.1.4           tidyselect_1.2.1      
# [49] dichromat_2.0-0.1      abind_1.4-8            spatstat.random_3.4-1 
# [52] codetools_0.2-19       miniUI_0.1.1.1         spatstat.explore_3.4-3
# [55] listenv_0.9.0          lattice_0.22-7         tibble_3.2.1          
# [58] plyr_1.8.9             withr_3.0.2            shiny_1.10.0          
# [61] ROCR_1.0-11            Rtsne_0.17             gridGraphics_0.5-1    
# [64] future_1.58.0          fastDummies_1.7.5      survival_3.8-3        
# [67] polyclip_1.10-7        fitdistrplus_1.1-11    pillar_1.10.2         
# [70] KernSmooth_2.23-26     ggfun_0.1.8            plotly_4.10.4         
# [73] generics_0.1.3         RcppHNSW_0.6.0         scales_1.4.0          
# [76] globals_0.18.0         xtable_1.8-4           glue_1.8.0            
# [79] lazyeval_0.2.2         tools_4.3.2            data.table_1.17.4     
# [82] RSpectra_0.16-1        RANN_2.6.1             fs_1.6.6              
# [85] dotCall64_1.2          cowplot_1.1.3          grid_4.3.2            
# [88] tidyr_1.3.1            nlme_3.1-168           patchwork_1.3.0       
# [91] cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-1           
# [94] viridisLite_0.4.2      dplyr_1.1.4            uwot_0.2.3            
# [97] gtable_0.3.6           yulab.utils_0.2.0      digest_0.6.37         
# [100] progressr_0.15.1       ggrepel_0.9.6          ggplotify_0.1.2       
# [103] htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1     
# [106] lifecycle_1.0.4        httr_1.4.7             mime_0.13             
# [109] MASS_7.3-60           
