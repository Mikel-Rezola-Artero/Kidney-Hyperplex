####1. LOAD LIBRARIES####
library(rstudioapi)
library(readxl)
library(ggplot2)
library(dplyr)
library(Seurat)

####2. IMPORT DATA####

#Set directory to script location & obtain file name
setwd(dirname(getActiveDocumentContext()$path))
files <- list.files(path = "Data/", pattern = "Skin", full.names = TRUE)

#Create list to store the data
Skin_IF <- list()

#Load CCE Skin data
Skin_IF$CCE <- read_xlsx(files)

#Check if any data are missing
any(is.na(Skin_IF$CCE))


####3. Q.C.: CHECK OF AREAS SELECTED FROM THE SKIN (ONLY ARTERIES)####

#Filter Arteries only
Art_IF <- lapply(Skin_IF, function(df) {
  df[grep("arteries", df$`Analysis Region`, value = FALSE), ]
})

#Plot Skin Arteries only for Q.C.
plot_data <- do.call(rbind, lapply(names(Art_IF), function(name) {
  df <- Art_IF[[name]][,c(2,5,7)]}))
ggplot(plot_data, aes(x = as.numeric(`XMin`), y = as.numeric(`YMin`))) +
  geom_point(aes(color = `Analysis Region`, 
                 shape = `Analysis Region`), 
             size = 1, alpha = 0.3) +
  labs(x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() 
#We clearly see artery structures distributed spatially in the skin with
#two regions: Papillary Dermis = Close to epidermis, Reticular Dermis = Deep Skin
#In the Reticular Dermis we have two types of arteries: "Prof" (Deep Arteries) and "Sup" (Superficial arteries)


####4. Q.C.: COMMON Abs/MARKERS ACROSS DATASETS FOR CLUSTERING####

#Check raw column names
colnames(Art_IF$CCE)#No clear Ab naming

#Take Cell Intensity Columns
Art_intensity <- lapply(Art_IF, function(df) {
  df %>% select(contains("Cell Intensity"))
})

#Rename elements to create variables for the intensity object
names(Art_intensity) <- paste0(names(Art_intensity),".Int")

#Correct Cell intensity column names
colnames(Art_intensity$CCE.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Art_intensity$CCE.Int))
colnames(Art_intensity$CCE.Int)[c(18)] <- c("VCAM-1")

#Check corrected column names
Reduce(intersect, lapply(Art_intensity, colnames))
#Ab names have been corrected


####5. PROCESS FLUORESCENCE INTENSITIES FOR MXSAMPLE OBJECT CREATION####

# Convert all columns to numeric
Art_intensity[["CCE.Int"]][] <- lapply(Art_intensity[["CCE.Int"]], as.numeric)

# Remove DAPI from Intensity Data
Art_intensity[["CCE.Int"]] <- Art_intensity[["CCE.Int"]][, -1]


####6. PROCESS SUPERVISED CLASSIFICATION DATA FOR MXSAMPLE OBJECT CREATION####

#Take Cell Classification Columns with manual classification labels
Art_classification <- lapply(Art_IF, function(df) {
  df %>% select(contains("Positive Classification"))
})

#Rename elements to create variables for each Cell Classification object
names(Art_classification) <- paste0(names(Art_classification),".meta")

#Correct classification column names
colnames(Art_classification$CCE.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Art_classification$CCE.meta))
colnames(Art_classification$CCE.meta)[c(17)] <- c("VCAM-1")

#Check corrected column names
Reduce(intersect, lapply(Art_classification, colnames))

#Common marker names in supervised and non-supervised data
Markers <- intersect(
  Reduce(intersect, lapply(Art_intensity, colnames)),#Intensity
  Reduce(intersect, lapply(Art_classification, colnames))#Classification
)


####7. COMBINE OBJECTS FOR SEURAT CLUSTERING####

#Fuse Ab intensity data

#Extract selected Ab intensities
selected_intensities <- lapply(Art_intensity, function(df) {
  df[, Markers, drop = FALSE]
})

#Transform intesity list into df for selected Abs
Merged.Intensities <- do.call(rbind, selected_intensities)


#Fuse Supervised Classifications

#Extract slected Ab supervised labels and add add Image and Slide columns
selected_meta <- Map(function(df, name, idx) {
  df_selected <- df[, Markers, drop = FALSE]#extract and order column markers
  df_selected$slide <- name#add slide name information
  df_selected$image <- paste0("image", idx)#add image number information
  return(df_selected)
}, 
Art_classification,#df
gsub(".meta","",names(Art_classification)),#name
seq_along(Art_classification)#idx
)


#Transform Supervised Classifications list into df for selected Abs
Merged.Meta <- do.call(rbind, selected_meta)

# Append ".Pos" suffix to meta Ab data
marker_cols <- 1:(ncol(Merged.Meta) - 2)  # exclude slide and image columns
colnames(Merged.Meta)[marker_cols] <- paste0(colnames(Merged.Meta)[marker_cols], ".Pos")


#Create df for mxsample object creation

#Check if rows are correctly aligned
table(gsub("Int.","",rownames(Merged.Intensities)) == 
        gsub("meta.","",rownames(Merged.Meta)))


#Combine intensity and meta data
mx_sample <- cbind(Merged.Intensities,Merged.Meta)

#Add Artery Region Annotations
art.info <- lapply(Art_IF, function (df) {
  arts <- df$`Analysis Region`
})
mx_sample$Arteries <- unlist(art.info)


####8. CREATE SEURAT OBJECT AND IDENTIFY CELL TYPES#####

#Fluorescence Intensity data
seurat_int <- data.frame(mx_sample[,c(1:30)])#extract intensity data
seurat_int[is.na(seurat_int)]#check for missing data
#Add nuclear size data
seurat_int$Nuc.size <- unlist(lapply(Art_IF, function(df) df[["Nucleus Area (µm²)"]]))

#Create Seurat object for clustering
Seurat_mxnorm <- CreateSeuratObject(counts = t(seurat_int),
                                    meta.data = data.frame(cbind(Merged.Meta,mx_sample[c(63)]),
                                                           row.names = rownames(seurat_int)))

#Normalize "data" slot with pseudolog transformed IF values
Seurat_mxnorm@assays$RNA$data <- log1p(Seurat_mxnorm@assays$RNA$counts)

#Set variable features to all continuous data
VariableFeatures(Seurat_mxnorm) <- rownames(Seurat_mxnorm@assays$RNA$counts)

#Scale data prior to dimensional reduction
Seurat_mxnorm <- ScaleData(Seurat_mxnorm,scale.max = 5)

#Run PCA with cell marker Ab intensities
Seurat_mxnorm <- RunPCA(Seurat_mxnorm, 
                        assay = "RNA",
                        features = c("C5aR1",#Myeloid mainly
                                     "CD68",#Macrophage
                                     "CD163",#M2 macrophage
                                     "aSMA",#vSMC
                                     "CD45",#immune cells
                                     "CD3",#T cells
                                     "VCAM.1",#Activated Endothelial cells
                                     "CD34",#Fibro/Endo
                                     "CD11b",#Phagocytes
                                     "CD4",#CD4 T cells
                                     "CD8",#CD8 T cells
                                     "CD31",#Endothelial
                                     "vWF",#Endothelial
                                     "Nuc.size"#All
                        ),
                        npcs = 14,
                        approx = FALSE)

#Run UMAP for aesthetics
Seurat_mxnorm <- RunUMAP(Seurat_mxnorm, dims = 1:14, 
                         n.components = 2)

#Find Clusters using the Dimensional Reduction Markers
Seurat_mxnorm <- FindNeighbors(Seurat_mxnorm, dims = 1:14)
Seurat_mxnorm <- FindClusters(Seurat_mxnorm,resolution = 0.4)

#Check clusters
DimPlot(object = Seurat_mxnorm,repel = T,
        label = T,label.box = T,label.size = 5,
        pt.size = 2,reduction = "umap") + NoLegend()

#Check distribution of key cell markers as DotPlot
DotPlot(Seurat_mxnorm,
        c("aSMA",#vSMC
          "CD31",#Endothelial
          "vWF",#Endothelial
          "VCAM.1",#Activated Endothelial cells
          "CD34",#Fibro/Endo/
          "C5aR1",#Myeloid mainly
          "CD68",#Macrophage
          "CD163",#M2 macrophage
          "CD11b",#Phagocytes
          "CD3",#T cells
          "CD4",#CD4 T cells
          "CD8"#CD8 T cells
        ),
        scale = T,
        col.min = 0.5,
        col.max = 2,cols = c("white","blue"),
        dot.min = F)

#Check Main IF markers
FeaturePlot(object = Seurat_mxnorm,ncol = 2,
            pt.size = 1, repel = T,label.size = 3,
            reduction = "umap",label = T,
            features = c("CD3",#T cells
                         "CD68",#CD4 T cells
                         "CD31",#CD8 T cells
                         "aSMA"
            ),min.cutoff = "q50",max.cutoff = "q99"
)

#Check Main supervised labels
FeaturePlot(object = Seurat_mxnorm,ncol = 2,
            pt.size = 1, repel = T,label.size = 3,
            reduction = "umap",label = T,
            features = c("CD3.Pos",#T cells
                         "CD68.Pos",#CD4 T cells
                         "CD31.Pos",#CD8 T cells
                         "aSMA.Pos"
            )
)

#Cluster 0 = aSMA
#Cluster 1 = Endothelial
#Cluster 2 = Unknown
#Cluster 3 = Macrophages
#Cluster 4 = T cells

#Rename clusters according to Cell Markers
new_cluster_names <- c("vSMC",
                       "Endothelial",
                       "Unknown",
                       "Macrophages",
                       "T cells")

names(new_cluster_names) <- levels(Seurat_mxnorm)
Seurat_mxnorm <- RenameIdents(Seurat_mxnorm, 
                              new_cluster_names)

#Store Processed Seurat Object with labels
saveRDS(Seurat_mxnorm, file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),
                                     "/Only skin/Processed_Skin.rds"))

# If curious plot clusters as 3D
# Seurat_mxnorm$clusters <- Idents(Seurat_mxnorm)
# Seurat_mxnorm <- RunUMAP(Seurat_mxnorm, dims = 1:6, n.components = 2)
# SCP::CellDimPlot3D(srt = Seurat_mxnorm,
#                    group.by = "clusters")


####9. Session Information#####

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
#   [1] future_1.58.0      Seurat_5.3.0       SeuratObject_5.1.0
# [4] sp_2.2-0           dplyr_1.1.4        ggplot2_3.5.2     
# [7] readxl_1.4.5       rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3         
# [4] rlang_1.1.6            magrittr_2.0.3         RcppAnnoy_0.0.22      
# [7] spatstat.geom_3.4-1    matrixStats_1.5.0      ggridges_0.5.4        
# [10] compiler_4.3.2         png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         stringr_1.5.1          crayon_1.5.3          
# [16] pkgconfig_2.0.3        fastmap_1.2.0          labeling_0.4.3        
# [19] promises_1.3.3         purrr_1.0.4            aplot_0.2.6           
# [22] jsonlite_2.0.0         goftest_1.2-3          later_1.4.2           
# [25] spatstat.utils_3.1-4   irlba_2.3.5.1          parallel_4.3.2        
# [28] cluster_2.1.6          R6_2.6.1               ica_1.0-3             
# [31] spatstat.data_3.1-6    stringi_1.8.7          RColorBrewer_1.1-3    
# [34] reticulate_1.42.0      spatstat.univar_3.1-3  parallelly_1.45.0     
# [37] lmtest_0.9-40          scattermore_1.2        cellranger_1.1.0      
# [40] Rcpp_1.0.14            tensor_1.5             future.apply_1.11.0   
# [43] zoo_1.8-12             sctransform_0.4.2      httpuv_1.6.16         
# [46] Matrix_1.6-5           splines_4.3.2          igraph_2.1.4          
# [49] tidyselect_1.2.1       abind_1.4-8            dichromat_2.0-0.1     
# [52] spatstat.random_3.4-1  codetools_0.2-19       miniUI_0.1.1.1        
# [55] spatstat.explore_3.4-3 listenv_0.9.0          lattice_0.22-7        
# [58] tibble_3.2.1           plyr_1.8.9             shiny_1.10.0          
# [61] withr_3.0.2            ROCR_1.0-11            Rtsne_0.17            
# [64] gridGraphics_0.5-1     fastDummies_1.7.5      survival_3.8-3        
# [67] polyclip_1.10-7        fitdistrplus_1.1-11    pillar_1.10.2         
# [70] KernSmooth_2.23-26     ggfun_0.1.8            plotly_4.10.4         
# [73] generics_0.1.3         RcppHNSW_0.6.0         scales_1.4.0          
# [76] globals_0.18.0         xtable_1.8-4           glue_1.8.0            
# [79] lazyeval_0.2.2         tools_4.3.2            data.table_1.17.4     
# [82] RSpectra_0.16-1        RANN_2.6.1             fs_1.6.6              
# [85] dotCall64_1.2          cowplot_1.1.3          grid_4.3.2            
# [88] tidyr_1.3.1            nlme_3.1-168           patchwork_1.3.0       
# [91] cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-1           
# [94] viridisLite_0.4.2      uwot_0.2.3             gtable_0.3.6          
# [97] yulab.utils_0.2.0      ggbreak_0.1.4          digest_0.6.37         
# [100] progressr_0.15.1       ggrepel_0.9.6          ggplotify_0.1.2       
# [103] htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1     
# [106] lifecycle_1.0.4        httr_1.4.7             mime_0.13             
# [109] MASS_7.3-60           
