####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(phenoptr)
library(ggplot2)
library(dplyr)
library(jpeg)

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


####3. FIGURE 1F####

#Get tonsil data
Tonsil <- tonsil_list[["Tonsil3"]]

#Calculate centroids
df <- data.frame(
  x = (as.numeric(Tonsil$XMin) + as.numeric(Tonsil$XMax)) / 2,
  y = (as.numeric(Tonsil$YMin) + as.numeric(Tonsil$YMax)) / 2
)

#Add Seurat clusters
df$clusters <- Idents(Seurat.obj)[Seurat.obj$slide == "Tonsil3"]

#Prepare for distance calculation
cds <- df %>%
  mutate(
    `Cell X Position` = x,
    `Cell Y Position` = y,
    Phenotype = as.character(clusters),
    `Cell ID` = row_number()
  ) %>%
  select(`Cell X Position`, `Cell Y Position`, Phenotype, `Cell ID`) %>%
  as_tibble()

#Detect correct marker columns for phenotype reassignment
C4d_col   <- grep("^C4d.*Positive Classification$", colnames(Tonsil), value = TRUE)[1]
C9_col    <- grep("^C9.*Positive Classification$", colnames(Tonsil), value = TRUE)[1]
C5aR1_col <- grep("^C5aR1.*Positive Classification$", colnames(Tonsil), value = TRUE)[1]

#Reassign phenotypes
pheno <- cds$Phenotype
pheno[Tonsil[[C9_col]] == 1 & Tonsil[[C4d_col]] == 1] <- "C4d+/C9+ cell"
pheno[cds$Phenotype == "Myeloid" & Tonsil[[C5aR1_col]] == 1] <- "Myeloid C5aR1+"
pheno[cds$Phenotype == "Myeloid" & Tonsil[[C5aR1_col]] == 0] <- "Myeloid C5aR1-"

cds$Phenotype <- pheno

#Filter relevant cells
C9.cells <- cds %>% 
  filter(Phenotype!='Endothelial' &
           Phenotype!='CD4/CD8 T cell')

#Compute per-cell nearest neighbor distances
dist.C9 <- find_nearest_distance(C9.cells)
gc()
C9_with_distance <- bind_cols(C9.cells, dist.C9)

#Plot FIGURE 1F
ggplot(C9_with_distance %>% filter(Phenotype != "C4d+/C9+ cell"), 
       aes(`Distance to C4d+/C9+ cell`, color = Phenotype)) +
  geom_density(size = 1) +
  scale_y_continuous("Cell probability within each cluster") +
  scale_color_manual(values = c("red","yellow","orange","grey","black")) +
  ggtitle(paste("Tonsil 3 Distance from C4d+/C9+ cells"))


####4. FIGURE 1G####

#Re-align Tonsil3 data for image overlapping
csd_with_distance2 <- C9_with_distance
csd_with_distance2$`Cell X Position` <- C9_with_distance$`Cell X Position`- min(C9_with_distance$`Cell X Position`)
csd_with_distance2$`Cell Y Position` <- C9_with_distance$`Cell Y Position`- min(C9_with_distance$`Cell Y Position`)

#Filter to just Myeloid C5aR1+ cells & C4d+/C9+ cells
Myeloid_cells = csd_with_distance2 %>% filter(select_rows(C9_with_distance, "Myeloid C5aR1+"))
Comp_cells = csd_with_distance2 %>% filter(select_rows(C9_with_distance, "C4d+/C9+ cell"))

#For each C4d+/C9+ cell, join with the data for the nearest Myeloid C5aR1+ cell
Comp_to_Myeloid = Comp_cells %>% left_join(Myeloid_cells, by=c('Cell ID Myeloid C5aR1+'='Cell ID'),
                                           suffix=c('', '.B'))

#Read a background image and make a base plot (DAPI image of slide)
background_path = 
  system.file("~/Idris stuff/Integrated complementomics/CCE Lupus DN/Tonsil IFA multiplex/Analysis Paper/Data/Tonsil 3.jpg", package='phenoptr')
background = readJPEG("~/Idris stuff/Integrated complementomics/CCE Lupus DN/Tonsil IFA multiplex/Analysis Paper/Data/Tonsil 3.jpg") %>% as.raster()
xlim = c(0, 12989.5)#max X axis - min X axis for "csd_with_distance2"
ylim = c(0, 10061)#max Y axis - min Y axis "csd_with_distance2"
base_plot = ggplot(mapping=aes(`Cell X Position`, `Cell Y Position`)) %>% 
  phenoptr:::add_scales_and_background(background, xlim, ylim, scale_color= "black") +
  labs(x='Cell X Position', y='Cell Y Position') +
  scale_color_manual('Phenotype', 
                     values=c('Myeloid C5aR1+'='red','C4d+/C9+ cell'='lightgreen' ),
                     guide = guide_legend(override.aes = list(size = 4)))

#Add distances to nearest neighbors as lines and points representing cells
base_plot + geom_segment(data=Comp_to_Myeloid,
                         aes(xend=`Cell X Position.B`, yend=`Cell Y Position.B`),
                         color='white') +
  geom_point(data=Myeloid_cells, aes(color='Myeloid C5aR1+'), size=0.1) +
  geom_point(data=Comp_cells, aes(color='C4d+/C9+ cell'), size=0.1) +
  labs(title='Nearest Myeloid C5aR1+ to each C4d+/C9+ cell') #FIGURE 1G


####5. EXTRA####

#Function to calculate distances between reassigned Seurat clusters 
analyze_tonsil_distances <- function(tonsil_list, Seurat.obj) {
  
  results <- list()
  
  for (tonsil_name in names(tonsil_list)) {
    message(paste("Processing:", tonsil_name))
    
    #Get tonsil data
    Tonsil <- tonsil_list[[tonsil_name]]
    
    #Calculate centroids
    df <- data.frame(
      x = (as.numeric(Tonsil$XMin) + as.numeric(Tonsil$XMax)) / 2,
      y = (as.numeric(Tonsil$YMin) + as.numeric(Tonsil$YMax)) / 2
    )
    
    #Add Seurat clusters
    df$clusters <- Idents(Seurat.obj)[Seurat.obj$slide == tonsil_name]
    
    #Prepare for distance calculation
    cds <- df %>%
      mutate(
        `Cell X Position` = x,
        `Cell Y Position` = y,
        Phenotype = as.character(clusters),
        `Cell ID` = row_number()
      ) %>%
      select(`Cell X Position`, `Cell Y Position`, Phenotype, `Cell ID`) %>%
      as_tibble()
    
    #Detect correct marker columns for phenotype reassignment
    C4d_col   <- grep("^C4d.*Positive Classification$", colnames(Tonsil), value = TRUE)[1]
    C9_col    <- grep("^C9.*Positive Classification$", colnames(Tonsil), value = TRUE)[1]
    C5aR1_col <- grep("^C5aR1.*Positive Classification$", colnames(Tonsil), value = TRUE)[1]
    
    #Reassign phenotypes
    pheno <- cds$Phenotype
    pheno[Tonsil[[C9_col]] == 1 & Tonsil[[C4d_col]] == 1] <- "C4d+/C9+ cell"
    pheno[cds$Phenotype == "Myeloid" & Tonsil[[C5aR1_col]] == 1] <- "Myeloid C5aR1+"
    pheno[cds$Phenotype == "Myeloid" & Tonsil[[C5aR1_col]] == 0] <- "Myeloid C5aR1-"
    
    cds$Phenotype <- pheno
    
    #Filter relevant cells
    C9.cells <- cds %>% 
      filter(Phenotype %in% c("Myeloid C5aR1+", "Myeloid C5aR1-", "C4d+/C9+ cell"))
    
    #Compute per-cell nearest neighbor distances
    dist.C9 <- find_nearest_distance(C9.cells)
    gc()
    C9_with_distance <- bind_cols(C9.cells, dist.C9)
    
    #Plot
    p <- ggplot(C9_with_distance %>% filter(Phenotype != "C4d+/C9+ cell"), 
                aes(`Distance to C4d+/C9+ cell`, color = Phenotype)) +
      geom_density(size = 1) +
      scale_y_continuous("Cell probability within each cluster") +
      scale_color_manual(values = c("grey","black")) +
      ggtitle(paste("Distance distribution:", tonsil_name))
    
    #Compute medians
    med_C5aR1_neg <- median(dist.C9$`Distance to C4d+/C9+ cell`[C9_with_distance$Phenotype == "Myeloid C5aR1-"])
    med_C5aR1_pos <- median(dist.C9$`Distance to C4d+/C9+ cell`[C9_with_distance$Phenotype == "Myeloid C5aR1+"])
    
    results[[tonsil_name]] <- list(
      plot = p,
      medians = list(C5aR1_minus = med_C5aR1_neg, C5aR1_plus = med_C5aR1_pos)
    )
  }
  
  return(results)
}

#Calculate distances and make plots
res <- analyze_tonsil_distances(tonsil_list, Seurat.obj)

#Plot all distances between Myeloid C5aR1+/- cells and C4d+/C9+ cells
library(patchwork)
library(cowplot)

(res$Tonsil1$plot|res$Tonsil2$plot) /
  (res$Tonsil3$plot|res$Tonsil4$plot)


#Compare median distances to nd C4d+/C9+ cells between Myeloid C5aR1+/- cells

#Extract median distances
med_dists <- do.call(rbind, lapply(names(res), function(name) {
  med <- res[[name]]$medians
  data.frame(
    Tonsil = name,
    Median_C5aR1_minus = med$C5aR1_minus,
    Median_C5aR1_plus  = med$C5aR1_plus,
    stringsAsFactors = FALSE
  )
}))

#Transform distances into log scale
med_dists$Median_C5aR1_minus <- log2(med_dists$Median_C5aR1_minus)
med_dists$Median_C5aR1_plus <- log2(med_dists$Median_C5aR1_plus)

#Plot median distances
library(ggpubr)
library(tidyr)
ggpaired(med_dists,
         cond1 = "Median_C5aR1_minus",
         cond2 = "Median_C5aR1_plus",
         line.color = "gray",
         line.size = 0.5,
         palette = "jco") +
  stat_compare_means(paired = TRUE, method = "t.test") +
  labs(title = "Paired Comparison of Median Distances",
       y = "Median Distance")


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
#   [1] jpeg_0.1-10        cowplot_1.1.3      patchwork_1.3.0    tidyr_1.3.1       
# [5] ggpubr_0.6.0       dplyr_1.1.4        phenoptr_0.3.2     ggplot2_3.5.2     
# [9] Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0           rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.3         spatstat.utils_3.1-4  
# [5] farver_2.1.2           rmarkdown_2.29         vctrs_0.6.5            ROCR_1.0-11           
# [9] spatstat.explore_3.4-3 base64enc_0.1-3        rstatix_0.7.2          htmltools_0.5.8.1     
# [13] broom_1.0.8            Formula_1.2-5          sctransform_0.4.2      parallelly_1.45.0     
# [17] KernSmooth_2.23-26     htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9            
# [21] plotly_4.10.4          zoo_1.8-12             igraph_2.1.4           mime_0.13             
# [25] lifecycle_1.0.4        iterators_1.0.14       pkgconfig_2.0.3        Matrix_1.6-5          
# [29] R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.1-11    future_1.58.0         
# [33] shiny_1.10.0           digest_0.6.37          fdrtool_1.2.18         colorspace_2.1-1      
# [37] tensor_1.5             RSpectra_0.16-1        irlba_2.3.5.1          Hmisc_5.2-3           
# [41] labeling_0.4.3         progressr_0.15.1       spatstat.sparse_3.1-0  httr_1.4.7            
# [45] polyclip_1.10-7        abind_1.4-8            compiler_4.3.2         withr_3.0.2           
# [49] glasso_1.11            htmlTable_2.4.3        backports_1.5.0        carData_3.0-5         
# [53] psych_2.5.3            fastDummies_1.7.5      ggsignif_0.6.4         MASS_7.3-60           
# [57] corpcor_1.6.10         ggsci_3.2.0            gtools_3.9.5           tools_4.3.2           
# [61] pbivnorm_0.6.0         foreign_0.8-90         lmtest_0.9-40          httpuv_1.6.16         
# [65] future.apply_1.11.0    goftest_1.2-3          nnet_7.3-19            glue_1.8.0            
# [69] quadprog_1.5-8         nlme_3.1-168           promises_1.3.3         grid_4.3.2            
# [73] checkmate_2.3.2        Rtsne_0.17             cluster_2.1.6          reshape2_1.4.4        
# [77] generics_0.1.3         gtable_0.3.6           spatstat.data_3.1-6    data.table_1.17.4     
# [81] car_3.1-3              utf8_1.2.6             spatstat.geom_3.4-1    RcppAnnoy_0.0.22      
# [85] ggrepel_0.9.6          RANN_2.6.1             foreach_1.5.2          pillar_1.10.2         
# [89] stringr_1.5.1          spam_2.11-1            RcppHNSW_0.6.0         later_1.4.2           
# [93] splines_4.3.2          lattice_0.22-7         deldir_1.0-9           survival_3.8-3        
# [97] tidyselect_1.2.1       miniUI_0.1.1.1         pbapply_1.7-2          knitr_1.50            
# [101] gridExtra_2.3          scattermore_1.2        stats4_4.3.2           xfun_0.52             
# [105] qgraph_1.9.8           matrixStats_1.5.0      stringi_1.8.7          lazyeval_0.2.2        
# [109] evaluate_1.0.3         codetools_0.2-19       tibble_3.2.1           cli_3.6.5             
# [113] uwot_0.2.3             rpart_4.1.24           xtable_1.8-4           reticulate_1.42.0     
# [117] lavaan_0.6-19          dichromat_2.0-0.1      Rcpp_1.0.14            spatstat.random_3.4-1 
# [121] globals_0.18.0         png_0.1-8              spatstat.univar_3.1-3  parallel_4.3.2        
# [125] dotCall64_1.2          listenv_0.9.0          viridisLite_0.4.2      scales_1.4.0          
# [129] ggridges_0.5.4         crayon_1.5.3           purrr_1.0.4            rlang_1.1.6           
# [133] mnormt_2.1.1  
