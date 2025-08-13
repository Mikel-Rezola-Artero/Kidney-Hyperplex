####1. LOAD REQUIRED LIBRARIES####
library(Seurat)
library(rstudioapi)
library(ggplot2)

####2. IMPORT & PROCESS snRNA-seq LIVER DATA####

#Set data location & select file with snRNAseq data
setwd(dirname(getActiveDocumentContext()$path))
data.dir <- "Liver single nuclei/snRNAseq Data"

#Single-cell transcriptomic RNA-seq data obtained from 
#https://www.livercellatlas.org/download.php 
#(Liver Cell Atlas: Human > All liver cells).
#Data contains counts for all the cells before filtering and annotation matrix
#containes cell IDs post-quality filtering (used to filter counts data here)
#for more details check PMID: 35021063 

#Load the liver dataset
Liver <- Read10X(data.dir,
                 gene.column = 1)

#Create Seurat object
Liver <- CreateSeuratObject(counts = Liver,
                            min.cells = 3,
                            min.features = 200) 

#Load annotation data and UMAP embeddings
UMAP.dir <- "Liver single nuclei/annot_humanAll.csv"
annot <- read.csv2(UMAP.dir, sep = ",")

#Filter out low quality cells using author annotation labels
Liver_filtered <- subset(Liver, cells = Cells(Liver)[Cells(Liver) %in% annot$cell])

#Re-order annotations according to cell order in Seurat object
annot.ord <- annot[match(Cells(Liver_filtered), annot$cell),]

#Check samples containing healthy and diseased samples
table(annot.ord$patient,annot.ord$diet)

#Process and add UMAP embeddings data
umap_embeddings <- as.matrix(annot.ord[,1:2]) #Take UMAP embeddings
umap_embeddings <- apply(umap_embeddings,2,as.numeric) #Make sure they are numeric
rownames(umap_embeddings) <- Cells(Liver_filtered)  # Add cell names
colnames(umap_embeddings) <- c("UMAP_1", "UMAP_2")  # Set UMAP dimension names
Liver_filtered[["umap"]] <- CreateDimReducObject(embeddings = umap_embeddings,
                                                 key = "UMAP_",
                                                 assay = DefaultAssay(Liver_filtered)) #Add to Seurat

#Add metadata to Seurat object 
Liver_filtered <- AddMetaData(Liver_filtered, annot.ord$diet, col.name = "Condition")
Liver_filtered <- AddMetaData(Liver_filtered, annot.ord$annot, col.name = "Cells")

#Check cells in diseased and healthy liver
DimPlot(Liver_filtered, reduction = "umap",group.by = "Cells",
        split.by = "Condition",label = T)

#Filter only healthy liver samples
Liver_filtered <- subset(Liver_filtered, Condition == "Lean")

#Normalize data 
Liver_filtered <- NormalizeData(Liver_filtered)

#Save processed data
#saveRDS(Liver_filtered, file = paste0("Liver single nuclei/",
#                                      "Liver Cells Processed.rds"))


####3. FIGURE S14C####

DimPlot(Liver_filtered, reduction = "umap",
        group.by = "Cells",label = T,label.size = 4,repel = T) #FIGURE S14C


####4. FIGURE S14D####

DotPlot(Liver_filtered,features =c("CD68","C1QA","C1QB","C1QC", "C5AR1","ITGAM",
                                   "CYP2E1","C1S","C1R","C3","C5","C6","C7","C9",
                                   "CFH","CFB"),
        dot.scale = 10,col.max = 1,group.by = "Cells") +
  scale_size(range = c(1, 8)) + RotatedAxis() + 
  scale_color_gradientn(colors = c("white", "lightblue","blue"))# FIGURE S14C


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
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8   
# [6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Paris
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.5.2      rstudioapi_0.17.1  Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.6            magrittr_2.0.3         RcppAnnoy_0.0.22      
# [7] matrixStats_1.5.0      ggridges_0.5.4         compiler_4.3.2         spatstat.geom_3.4-1    png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3        fastmap_1.2.0          promises_1.3.3         purrr_1.0.4           
# [19] jsonlite_2.0.0         goftest_1.2-3          later_1.4.2            spatstat.utils_3.1-4   irlba_2.3.5.1          parallel_4.3.2        
# [25] cluster_2.1.6          R6_2.6.1               ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3     spatstat.data_3.1-6   
# [31] reticulate_1.42.0      parallelly_1.45.0      spatstat.univar_3.1-3  lmtest_0.9-40          scattermore_1.2        Rcpp_1.1.0            
# [37] tensor_1.5             future.apply_1.11.0    zoo_1.8-12             sctransform_0.4.2      httpuv_1.6.16          Matrix_1.6-5          
# [43] splines_4.3.2          igraph_2.1.4           tidyselect_1.2.1       dichromat_2.0-0.1      abind_1.4-8            codetools_0.2-19      
# [49] spatstat.random_3.4-1  miniUI_0.1.1.1         spatstat.explore_3.4-3 listenv_0.9.0          lattice_0.22-7         tibble_3.3.0          
# [55] plyr_1.8.9             withr_3.0.2            shiny_1.10.0           ROCR_1.0-11            Rtsne_0.17             future_1.58.0         
# [61] fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.1-11    pillar_1.11.0          KernSmooth_2.23-26    
# [67] plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0         scales_1.4.0           globals_0.18.0         xtable_1.8-4          
# [73] glue_1.8.0             lazyeval_0.2.2         tools_4.3.2            data.table_1.17.4      RSpectra_0.16-1        RANN_2.6.1            
# [79] dotCall64_1.2          cowplot_1.1.3          grid_4.3.2             tidyr_1.3.1            nlme_3.1-168           patchwork_1.3.0       
# [85] cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-1            viridisLite_0.4.2      dplyr_1.1.4            uwot_0.2.3            
# [91] gtable_0.3.6           digest_0.6.37          progressr_0.15.1       ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2          
# [97] htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7             mime_0.13              MASS_7.3-60           
# > sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3;  LAPACK version 3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Paris
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.5.2      rstudioapi_0.17.1  Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.6            magrittr_2.0.3        
# [6] RcppAnnoy_0.0.22       matrixStats_1.5.0      ggridges_0.5.4         compiler_4.3.2         spatstat.geom_3.4-1   
# [11] png_0.1-8              vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3       
# [16] fastmap_1.2.0          promises_1.3.3         purrr_1.0.4            jsonlite_2.0.0         goftest_1.2-3         
# [21] later_1.4.2            spatstat.utils_3.1-4   irlba_2.3.5.1          parallel_4.3.2         cluster_2.1.6         
# [26] R6_2.6.1               ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3     spatstat.data_3.1-6   
# [31] reticulate_1.42.0      parallelly_1.45.0      spatstat.univar_3.1-3  lmtest_0.9-40          scattermore_1.2       
# [36] Rcpp_1.1.0             tensor_1.5             future.apply_1.11.0    zoo_1.8-12             sctransform_0.4.2     
# [41] httpuv_1.6.16          Matrix_1.6-5           splines_4.3.2          igraph_2.1.4           tidyselect_1.2.1      
# [46] dichromat_2.0-0.1      abind_1.4-8            codetools_0.2-19       spatstat.random_3.4-1  miniUI_0.1.1.1        
# [51] spatstat.explore_3.4-3 listenv_0.9.0          lattice_0.22-7         tibble_3.3.0           plyr_1.8.9            
# [56] withr_3.0.2            shiny_1.10.0           ROCR_1.0-11            Rtsne_0.17             future_1.58.0         
# [61] fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.1-11    pillar_1.11.0         
# [66] KernSmooth_2.23-26     plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0         scales_1.4.0          
# [71] globals_0.18.0         xtable_1.8-4           glue_1.8.0             lazyeval_0.2.2         tools_4.3.2           
# [76] data.table_1.17.4      RSpectra_0.16-1        RANN_2.6.1             dotCall64_1.2          cowplot_1.1.3         
# [81] grid_4.3.2             tidyr_1.3.1            nlme_3.1-168           patchwork_1.3.0        cli_3.6.5             
# [86] spatstat.sparse_3.1-0  spam_2.11-1            viridisLite_0.4.2      dplyr_1.1.4            uwot_0.2.3            
# [91] gtable_0.3.6           digest_0.6.37          progressr_0.15.1       ggrepel_0.9.6          htmlwidgets_1.6.4     
# [96] farver_2.1.2           htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7             mime_0.13             
# [101] MASS_7.3-60           