####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(ggplot2)

####2. LOAD TONSIL DATA####

#Set directory
setwd(dirname(getActiveDocumentContext()$path))

#Load Processed Data
Seurat.obj <- readRDS("Data/Processed_Tonsils.rds")

#Load original data files with coordinate information
files <- list.files(path = "Data/",#Access folder with data
                    pattern = "Tonsil[1-4]-Square.*\\.csv$",#Select relevant files 
                    full.names = TRUE)#Take the names of the files for analysis
tonsil_list <- lapply(files, read.csv2,sep = ",", check.names = F)#Import data as list
names(tonsil_list) <- paste0("Tonsil",1:4)#Name list elements

####3. FIGURE 1B ####

#Modify graphics to show 4 graphs at a time
plot_list <- list()

#Loop through all Tonsil slides
for (name in names(tonsil_list)) {
  df <- tonsil_list[[name]]
  
  #Ensure columns are numeric
  df$x <- (as.numeric(df$XMin) + as.numeric(df$XMax))/2
  df$y <- (as.numeric(df$YMin) + as.numeric(df$YMax))/2
  
  df$clusters <- Idents(Seurat.obj)[Seurat.obj$slide == name]
  
  
  # Plot
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(color = as.factor(clusters)), alpha = 0.5, size = 0.5) +
    theme_classic() +
    labs(title = name) +
    scale_color_manual(values = c("blue", "yellow", "red", "brown", "pink", "black")) + 
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  plot_list[[name]] <- p
  
}

#QC to check if cluster distribution makes sense according to prior tonsil knowledge
(plot_list[[1]] | plot_list[[2]]) /
  (plot_list[[3]] | plot_list[[4]])# Combine 4 plots into a 2x2 grid


#Plot Tonsil 3 clusters spatially (FIGURE 1B)
plot_list[[3]] + scale_y_reverse() #reverse y axis to match image


####4. FIGURE 1C ####

#Loop through all Tonsil slides for each C deposit
Complement <- c("C1q.Pos","C4d.Pos","C3d.Pos","C9.Pos")
Colors <- c("brown","blue","red","black")

for (i in 1:4) {
  
  Comp <- Complement[i]
  pos_color <- Colors[i]
  
  df <- tonsil_list[["Tonsil3"]]
  
  #Ensure columns are numeric
  df$x <- (as.numeric(df$XMin) + as.numeric(df$XMax))/2
  df$y <- (as.numeric(df$YMin) + as.numeric(df$YMax))/2
  
  df$Comp <- Seurat.obj@meta.data[Seurat.obj$slide == "Tonsil3", Comp]
  
  
  # Plot
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(color = as.numeric(Comp)), alpha = as.numeric(Comp), size = 1) +
    theme_classic() +
    labs(title = gsub(".Pos","",Comp)) +
    scale_colour_gradientn(colours = c("lightgrey", pos_color)) + 
    guides(color = guide_legend(override.aes = list(size = 3))) + 
    theme(legend.position="none",
          plot.title = element_text(
            hjust = 0.5,      
            face = "bold",
            size = 16 )) + scale_y_reverse()
  
  plot_list[[i]] <- p 
  
}

#Plot Tonsil 3 clusters spatially (FIGURE 1C)
(plot_list[[1]] | plot_list[[2]]) /
  (plot_list[[3]] | plot_list[[4]])# Combine 4 plots into a 2x2 grid


####5. FIGURE 1D ####

#Obtain Tonsil3 cluster data with coordinates
Tonsil3 <- tonsil_list[["Tonsil3"]]
df <- data.frame(x = (as.numeric(Tonsil3$XMin) + as.numeric(Tonsil3$XMax))/2,#x axis
                 y = (as.numeric(Tonsil3$YMin) + as.numeric(Tonsil3$YMax))/2,#y axis
                 Phenotype = Idents(Seurat.obj)[Seurat.obj$slide == "Tonsil3"],#Idents
                 IDs = colnames(Seurat.obj)[Seurat.obj$slide == "Tonsil3"]#Cell IDs
)

#Redefine cell labels and remove endothelial and T cells
pheno <- as.character(df$Phenotype)
pheno[Tonsil3$`C9 - Cy5 Positive Classification` == 1 &
        Tonsil3$`C4d - Cy5 Positive Classification` == 1] <- "C4d+/C9+ cell"
pheno[df$Phenotype == "Myeloid" &
        Tonsil3$`C5aR1 Positive Classification` == 1] <- "Myeloid C5aR1+"
pheno[df$Phenotype == "Myeloid" &
        Tonsil3$`C5aR1 - Cy5 Positive Classification` == 0] <- "Myeloid C5aR1-"
df$Phenotype <- pheno
df <- df[df$Phenotype != "Endothelial" & 
           df$Phenotype != "CD4/CD8 T cell",]

#Plot C9+ cells, C5aR1+ myeloid cells and C5aR1- myeloid cells in Tonsil 3
ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color = as.factor(Phenotype)),
             alpha = 0.4, size = 0.75) +
  theme_classic() + theme() + 
  labs(title="Tonsil 3") + 
  scale_color_manual(values = c("orange","yellow","lightblue",
                                "red","brown","black"),
                     guide = guide_legend(override.aes = list(size = 4))) + 
  scale_y_reverse() #FIGURE 1D

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
#   [1] ggplot2_3.5.2      Seurat_5.3.0       rstudioapi_0.17.1  SeuratObject_5.1.0 sp_2.2-0          
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
# [69] RcppHNSW_0.6.0         scales_1.4.0           globals_0.18.0         xtable_1.8-4          
# [73] glue_1.8.0             lazyeval_0.2.2         tools_4.3.2            data.table_1.17.4     
# [77] RSpectra_0.16-1        RANN_2.6.1             dotCall64_1.2          cowplot_1.1.3         
# [81] grid_4.3.2             tidyr_1.3.1            nlme_3.1-168           patchwork_1.3.0       
# [85] cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-1            viridisLite_0.4.2     
# [89] dplyr_1.1.4            uwot_0.2.3             gtable_0.3.6           digest_0.6.37         
# [93] progressr_0.15.1       ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2          
# [97] htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7             mime_0.13             
# [101] MASS_7.3-60   
