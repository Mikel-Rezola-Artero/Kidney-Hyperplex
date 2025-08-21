####1. LOAD LIBRARIES####
library(rstudioapi)
library(dplyr)
#Load reticulate and specify conda env (with "scikit-image") BEFORE loading mxnorm
library(reticulate)
use_condaenv("image-analysis-env", required = TRUE)
library(mxnorm)
library(plotly)
library(Seurat)

####2. IMPORT DATA####

#Set directory to script location & obtain files names
setwd(dirname(getActiveDocumentContext()$path))
files <- list.files(path = "Data/",#Access folder with data
                    pattern = "Tonsil[1-4]-Square.*\\.csv$",#Select relevant files 
                    full.names = TRUE)#Take the names of the files for analysis

#Import cell-segmented data from HALO Analysis
tonsil_list <- lapply(files, read.csv2,sep = ",", check.names = F)#Import data as list
names(tonsil_list) <- paste0("Tonsil",1:4)#Name list elements

####3. Q.C.: CHECK OF AREAS SELECTED FROM EACH IMAGE####

#Plot cells for quick annotation quality check

#Modify graphics to show 4 graphs at a time
par(mfrow = c(2, 2))

#Loop through all Tonsil slides
for (name in names(tonsil_list)) {
  df <- tonsil_list[[name]]
  
  #Ensure columns are numeric
  x <- as.numeric(df$XMin)
  y <- as.numeric(df$YMin)
  
  # Plot
  plot(x, y, col = "blue", main = name,
       xlab = "XMin", ylab = "YMin", pch = 20, cex = 0.6)
  
  #Remove the x and y vectors from environment
  rm(list = c("x","y"))
}

#There are some "white gaps", which are the blood vessels

#Set graphics to default
par(mfrow = c(1, 1))  

####4. Q.C.: COMMON Abs/MARKERS ACROSS DATASETS FOR CLUSTERING####

#Check raw column names
Reduce(intersect, lapply(tonsil_list, colnames))#No consensus in Ab naming

#Take Cell Intensity Columns
tonsil_intensity <- lapply(tonsil_list, function(df) {
  df %>% select(contains("Cell Intensity"))
})

#Rename elements to create variables for each intensity object
names(tonsil_intensity) <- paste0(names(tonsil_intensity),".Int")

#Assign each dataset to a variable
list2env(tonsil_intensity, envir = .GlobalEnv)#Assign each dataset to a variable

#Correct Cell intensity column names manually (no common patterns/spacers)
colnames(Tonsil1.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Tonsil1.Int))
colnames(Tonsil1.Int)#new names check Tonsil 1

colnames(Tonsil2.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Tonsil2.Int))
colnames(Tonsil2.Int) <- gsub("Abcam|8D6", "", colnames(Tonsil2.Int))
colnames(Tonsil2.Int)#new names check Tonsil 2

colnames(Tonsil3.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Tonsil3.Int))
colnames(Tonsil3.Int) <- gsub("Abcam", "", colnames(Tonsil3.Int))
colnames(Tonsil3.Int)#new names check Tonsil 3

colnames(Tonsil4.Int) <- gsub("[_ ].*Cell Intensity", "", colnames(Tonsil4.Int))
colnames(Tonsil4.Int) <- gsub("Abcam|-IB|\\(rb\\)", "", colnames(Tonsil4.Int))
colnames(Tonsil4.Int)[3] <- "FH"
colnames(Tonsil4.Int)#new names check Tonsil 4

#Add the objects with corrected column names back to the intensities list
for (i in 1:4) {
  name <- names(tonsil_intensity)[i]
  tonsil_intensity[[name]] <- get(name)
  rm(list = name)
}

#Check corrected column names
Reduce(intersect, lapply(tonsil_intensity, colnames))

#Markers for Cell Intensity: 
#Complement: "C1q" "C3c" "C3d" "C5aR1" "C1r" "C1s" "C4d" "FH" "FB" "C9"
#Immune Myeloid: "CD45" "CD68" "CD163" "CD11b" "MPO"
#Immune Lympho: "CD45" "CD20" "CD3" "CD8" "CD4" "IgG" 
#Stroma: "CD31" "aSMA"
#Cell State: "Ki67"

####5. PROCESS FLUORESCENCE INTENSITIES FOR MXNORM OBJECT CREATION####

for (i in 1:4) {
  int_name <- paste0("Tonsil", i, ".Int")
  orig_name <- paste0("Tonsil", i)
  
  # Convert all columns to numeric
  tonsil_intensity[[int_name]][] <- lapply(tonsil_intensity[[int_name]], as.numeric)
  
  # Remove DAPI from Intensity Data
  tonsil_intensity[[int_name]] <- tonsil_intensity[[int_name]][, -1]
}

####6. PROCESS SUPERVISED CLASSIFICATION DATA FOR MXNORM OBJECT CREATION####

#Take Cell Classification Columns with manual classification labels
tonsil_classification <- lapply(tonsil_list, function(df) {
  df %>% select(contains("Positive Classification"))
})

#Rename elements to create variables for each Cell Classification object
names(tonsil_classification) <- paste0(names(tonsil_classification),".meta")

#Assign each dataset to a variable
list2env(tonsil_classification, envir = .GlobalEnv)#Assign each dataset to a variable

#Correct classification column names manually (no common patterns/spacers)
colnames(Tonsil1.meta) <- gsub("[-_ ].*Classification", "", colnames(Tonsil1.meta))
colnames(Tonsil1.meta)#new names check Tonsil 1

colnames(Tonsil2.meta) <- gsub("[-_ ].*Classification", "", colnames(Tonsil2.meta))
colnames(Tonsil2.meta) <- gsub("Abcam|8D6", "", colnames(Tonsil2.meta))
colnames(Tonsil2.meta)#new names check Tonsil 2

colnames(Tonsil3.meta) <- gsub("[-_ ].*Classification", "", colnames(Tonsil3.meta))
colnames(Tonsil3.meta) <- gsub("Abcam", "", colnames(Tonsil3.meta))
colnames(Tonsil3.meta)#new names check Tonsil 3

colnames(Tonsil4.meta) <- gsub("[_ ].*Classification", "", colnames(Tonsil4.meta))
colnames(Tonsil4.meta) <- gsub("Abcam|-IB|\\(rb\\)", "", colnames(Tonsil4.meta))
colnames(Tonsil4.meta)[3] <- "FH"
colnames(Tonsil4.meta)#new names check Tonsil 4


#Add the objects with corrected column names back to the classification list
for (i in 1:4) {
  name <- names(tonsil_classification)[i]
  tonsil_classification[[name]] <- get(name)
  rm(list = name)
}

#Check corrected column names
Reduce(intersect, lapply(tonsil_classification, colnames))

#Markers for Cell Metadata: 
#Complement: "C1q" "C3c" "C3d" "C5aR1" "C1r" "C1s" "C4d" "FH" "FB" "C9"
#Immune Myeloid: "CD45" "CD68" "CD163" "CD11b" "MPO"
#Immune Lympho: "CD45" "CD20" "CD3" "CD8" "CD4" "IgG" 
#Stroma: "CD31" "aSMA"
#Cell State: "Ki67"

#Store common marker names in supervised and non-supervised data
Markers <- intersect(
  Reduce(intersect, lapply(tonsil_intensity, colnames)),#Intensity
  Reduce(intersect, lapply(tonsil_classification, colnames))#Classification
)

####7. COMBINED OBJECT FOR MXNORM BATCH CORRECTION####

#Fuse Ab intensity data

#Extract common Ab intensities across slides
selected_intensities <- lapply(tonsil_intensity, function(df) {
  df[, Markers, drop = FALSE]
})
#Check if column names and order are correct in each slide
name_list <- lapply(selected_intensities, colnames)
col_consistency <- sapply(name_list, function(x) identical(x, name_list[[1]]))
col_consistency#Are all markers accross slides in the correct order?
#Merge intensity data for common Abs
Merged.Intensities <- do.call(rbind, selected_intensities)

#Fuse Image Supervised Classifications

#Extract common Ab supervised labels across slides and add add Image and Slide columns
selected_meta <- Map(function(df, name, idx) {
  df_selected <- df[, Markers, drop = FALSE]#extract and order column markers
  df_selected$slide <- name#add slide name information
  df_selected$image <- paste0("image", idx)#add image number information
  return(df_selected)
}, 
tonsil_classification,#df
gsub(".meta","",names(tonsil_classification)),#name
seq_along(tonsil_classification)#idx
)
#Check if column names and order are correct in each slide
name_list <- lapply(selected_meta, colnames)
col_consistency <- sapply(name_list, function(x) identical(x, name_list[[1]]))
col_consistency#Are all markers accross slides in the correct order?
#Merge Supervised Classifications data for common Abs
Merged.Meta <- do.call(rbind, selected_meta)
# Append ".Pos" suffix to meta Ab data
marker_cols <- 1:(ncol(Merged.Meta) - 2)  # exclude slide and image columns
colnames(Merged.Meta)[marker_cols] <- paste0(colnames(Merged.Meta)[marker_cols], ".Pos")

#Create df for mxnorm object creation

#Check if rows are correctly aligned
table(gsub("Int.","",rownames(Merged.Intensities)) == 
        gsub("meta.","",rownames(Merged.Meta)))
#Combine intensity and meta data
mx_sample <- cbind(Merged.Intensities,Merged.Meta)

####8. PERFORM BATCH CORRECTION WITH MXNORM####

#Create mxnorm object
mx_norm = mx_dataset(data=mx_sample,
                     slide_id="slide",
                     image_id="image",
                     marker_cols=colnames(mx_sample)[1:23],
                     metadata_cols=grep("Pos",colnames(mx_sample), value = TRUE))
summary(mx_norm)#Check summary data

#Correct scale intensity batch effect using "log10_mean_divide"
mx_norm = mx_normalize(mx_data = mx_norm,
                       transform = "log10_mean_divide",
                       method="None")
summary(mx_norm)#Check summary data post batch correction

#Run an Otsu discordance score analysis to determine how well our normalization method performs
mx_norm = run_otsu_discordance(mx_norm,
                               table="both",
                               threshold_override = NULL,
                               plot_out = FALSE)
summary(mx_norm)#Check summary data including Otsu discordance

#Visualize the results of the Otsu miss-classification analysis: raw vs corrected
ggplotly(plot_mx_discordance(mx_norm))#MPO seems not completely corrected

#Visualize the densities of the marker values normalized vs raw
density_plot<- plot_mx_density(mx_norm)
# ggsave(filename = file.path("Outputs/", "density_plot_wide.png"),
#        plot = density_plot,
#        width = 35,  # Increase width to fit facets
#        height = 6,
#        dpi = 300)#Check png file for plot output

#Check % of positive cells for every marker according to supervised labeling
(apply(mx_norm$norm_data[grep("Pos",colnames(mx_sample), value = TRUE)],2,sum)/
    nrow(mx_norm$norm_data[grep("Pos",colnames(mx_sample), value = TRUE)]))*100 
#1% > cells classified as positive for MPO & FH; almost no positive cells for aSMA

#Do UMAP & check how the batch correction integration did
mx_norm = run_reduce_umap(mx_norm,
                          table="both",
                          marker_list = c("CD68","CD20","CD3","CD31",
                                           "aSMA","CD45"),
                          downsample_pct = 0.2,
                          metadata_col = colnames(mx_norm$norm_data)[c(26:48)])

#Q.C., Plot UMAP with raw and batch corrected data
plot_mx_umap(mx_norm,metadata_col = "slide")

#Check if main cell markers present logical distributions in batch corrected data
plot_mx_umap(mx_norm,metadata_col = "CD31.Pos")
plot_mx_umap(mx_norm,metadata_col = "CD20.Pos")
plot_mx_umap(mx_norm,metadata_col = "CD3.Pos")
plot_mx_umap(mx_norm,metadata_col = "CD4.Pos")
plot_mx_umap(mx_norm,metadata_col = "CD8.Pos")
plot_mx_umap(mx_norm,metadata_col = "CD68.Pos")

#We see good segregation of key markers, we can use batch corrected data for clustering

####9. CREATE SEURAT OBJECT AND IDENTIFY CELL TYPES#####

#Batch corrected data
seurat_int <- data.frame(mx_norm$norm_data[,c(3:25)])#extract batch corrected intensity data
seurat_int[is.na(seurat_int)]#check for missing data
seurat_int$Nuc.size <- unlist(lapply(tonsil_list, function(df) df[["Nucleus Area (µm²)"]]))#add nuclear size data
#Non-Batch corrected data
seurat_raw_int <- data.frame(mx_norm$data[,c(3:25)])#extract batch corrected intensity data
seurat_raw_int[is.na(seurat_raw_int)]#check for missing data
seurat_raw_int$Nuc.size <- unlist(lapply(tonsil_list, function(df) df[["Nucleus Area (µm²)"]]))#add nuclear size data

#Create Seurat objects which we will use for clustering
Seurat_raw <- CreateSeuratObject(counts = t(seurat_raw_int),
                                 meta.data = data.frame(Merged.Meta,
                                                        row.names = rownames(seurat_raw_int)))
Seurat_mxnorm <- CreateSeuratObject(counts = t(seurat_int),
                                    meta.data = data.frame(Merged.Meta,
                                                           row.names = rownames(seurat_int)))
#Create "data" slot with the batch corrected/normalized values
Seurat_raw@assays$RNA$data <- Seurat_raw@assays$RNA$counts
Seurat_mxnorm@assays$RNA$data <- Seurat_mxnorm@assays$RNA$counts

#Set variable features to all continuous data
VariableFeatures(Seurat_raw) <- rownames(Seurat_raw@assays$RNA$counts)
VariableFeatures(Seurat_mxnorm) <- rownames(Seurat_mxnorm@assays$RNA$counts)

#Scale data prior to dimensional reduction
Seurat_raw <- ScaleData(Seurat_raw,scale.max = 5)
Seurat_mxnorm <- ScaleData(Seurat_mxnorm,scale.max = 5)

#Run PCA using the features determined
Seurat_raw <- RunPCA(Seurat_raw, 
                     assay = "RNA",
                     features = c("CD68","CD20","CD3","CD31",
                                  "aSMA","CD45","Nuc.size"),
                     npcs = 6,
                     approx = FALSE)
Seurat_mxnorm <- RunPCA(Seurat_mxnorm, 
                        assay = "RNA",
                        features = c("CD68","CD20","CD3","CD31",
                                     "aSMA","CD45","Nuc.size"),
                        npcs = 6, 
                        approx = FALSE)

#Determine the ‘dimensionality’ of the dataset: choosing the number of PCs to used for further analysis.
ElbowPlot(Seurat_raw) | ElbowPlot(Seurat_mxnorm) #we will use the 6 PCs

#Run UMAP
Seurat_raw <- RunUMAP(Seurat_raw, dims = 1:6, n.components = 2)
Seurat_mxnorm <- RunUMAP(Seurat_mxnorm, dims = 1:6, n.components = 2)

#Check batch effect (Figure S2C-D)
p1 <- DimPlot(object = Seurat_raw,group.by = "slide",
              label = T,label.box = T,label.size = 5,raster = T,
              repel = T, pt.size = 2, reduction = "umap")+ NoLegend()
p2 <- DimPlot(object = Seurat_mxnorm,group.by = "slide",
              label = T,label.box = T,label.size = 5,raster = T,
              repel = T, pt.size = 2, reduction = "umap")+ NoLegend()
p1 | p2 #Figure S3C-D

#Find Clusters in the batch corrected data
Seurat_mxnorm <- FindNeighbors(Seurat_mxnorm,dims = 1:6) #Compute SNN graph for clustering
Seurat_mxnorm <- FindClusters(Seurat_mxnorm,resolution = 0.2) #Compute clusters on the SNN

#Check clusters
DimPlot(object = Seurat_mxnorm,repel = T,
        label = T,label.box = T,label.size = 5,
        cols = c("red","gold","mediumseagreen","skyblue","purple","pink"),
        pt.size = 1.5,reduction = "umap",raster = T) + NoLegend()

#Check distribution of key cell markers as DotPlot
DotPlot(Seurat_mxnorm,
        c("CD3","CD4","CD8",#adaptive T cells
          "CD20","Ki67",#adaptive replicating B cells
          "Nuc.size",#all
          "CD31",#endothelial
          "CD68","C1q","C5aR1"#Myeloid
        ),
        scale = T,
        col.min = 0,
        col.max = 2,cols = c("white","blue"),
        dot.min = F) 

#Check distribution of key cell markers as FeaturePlot
FeaturePlot(object = Seurat_mxnorm,raster = T,
            pt.size = 1, repel = T,ncol = 3,label.size = 3,
            reduction = "umap",label = T,
            features = c("CD3",#adaptive T cells
                         "CD20","Ki67",#adaptive replicating B cells
                         "Nuc.size",#all
                         "CD31",#endothelial
                         "CD68"),
            min.cutoff = "q50",max.cutoff = "q99")

#Check distribution of supervised labels
FeaturePlot(object = Seurat_mxnorm,ncol = 3,
            pt.size = 1,label = T,raster = T,
            reduction = "umap", repel = T,
            features = c("CD31.Pos", "CD20.Pos",
                         "CD3.Pos","CD68.Pos",
                         "Ki67.Pos"))

#Cluster 0 = T cells
#Cluster 1 = B cell Ki67 Pos 1
#Cluster 2 = B cell Ki67 Pos 2
#Cluster 3 = Endothelial Cell
#Cluster 4 = B cell Ki67 Neg
#Cluster 5 = Myeloid Cell

#Rename clusters according to Cell Markers
new_cluster_names <- c("CD4/CD8 T cell",
                       "B cell Ki67 High 1",
                       "B cell Ki67 High 2",
                       "Endothelial",
                       "B cell Ki67 Neg",
                       "Myeloid")
names(new_cluster_names) <- levels(Seurat_mxnorm)
Seurat_mxnorm <- RenameIdents(Seurat_mxnorm, 
                              new_cluster_names)

#Store Processed Seurat Object with labels
saveRDS(Seurat_mxnorm, file = paste0(dirname(getActiveDocumentContext()$path),
                                     "/Data/Processed_Tonsils.rds"))

####10. Session Information#####

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
#   [1] future_1.58.0      Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0          
# [5] plotly_4.10.4      ggplot2_3.5.2      mxnorm_1.0.3       dplyr_1.1.4       
# [9] rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.3        
# [4] spatstat.utils_3.1-4   farver_2.1.2           nloptr_2.2.1          
# [7] rmarkdown_2.29         vctrs_0.6.5            ROCR_1.0-11           
# [10] spatstat.explore_3.4-3 minqa_1.2.8            base64enc_0.1-3       
# [13] htmltools_0.5.8.1      Formula_1.2-5          sctransform_0.4.2     
# [16] parallelly_1.45.0      KernSmooth_2.23-26     phenoptr_0.3.2        
# [19] htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9            
# [22] zoo_1.8-12             igraph_2.1.4           mime_0.13             
# [25] lifecycle_1.0.4        iterators_1.0.14       pkgconfig_2.0.3       
# [28] Matrix_1.6-5           R6_2.6.1               fastmap_1.2.0         
# [31] rbibutils_2.3          fitdistrplus_1.1-11    shiny_1.10.0          
# [34] digest_0.6.37          fdrtool_1.2.18         colorspace_2.1-1      
# [37] patchwork_1.3.0        tensor_1.5             RSpectra_0.16-1       
# [40] irlba_2.3.5.1          Hmisc_5.2-3            labeling_0.4.3        
# [43] progressr_0.15.1       spatstat.sparse_3.1-0  polyclip_1.10-7       
# [46] httr_1.4.7             abind_1.4-8            compiler_4.3.2        
# [49] withr_3.0.2            glasso_1.11            htmlTable_2.4.3       
# [52] backports_1.5.0        psych_2.5.3            fastDummies_1.7.5     
# [55] MASS_7.3-60            corpcor_1.6.10         gtools_3.9.5          
# [58] tools_4.3.2            pbivnorm_0.6.0         foreign_0.8-90        
# [61] lmtest_0.9-40          httpuv_1.6.16          future.apply_1.11.0   
# [64] goftest_1.2-3          nnet_7.3-19            glue_1.8.0            
# [67] quadprog_1.5-8         nlme_3.1-168           promises_1.3.3        
# [70] grid_4.3.2             checkmate_2.3.2        Rtsne_0.17            
# [73] cluster_2.1.6          reshape2_1.4.4         generics_0.1.3        
# [76] spatstat.data_3.1-6    gtable_0.3.6           tidyr_1.3.1           
# [79] data.table_1.17.4      spatstat.geom_3.4-1    RcppAnnoy_0.0.22      
# [82] ggrepel_0.9.6          RANN_2.6.1             foreach_1.5.2         
# [85] pillar_1.10.2          stringr_1.5.1          spam_2.11-1           
# [88] RcppHNSW_0.6.0         later_1.4.2            splines_4.3.2         
# [91] lattice_0.22-7         deldir_1.0-9           survival_3.8-3        
# [94] tidyselect_1.2.1       miniUI_0.1.1.1         pbapply_1.7-2         
# [97] knitr_1.50             reformulas_0.4.1       gridExtra_2.3         
# [100] scattermore_1.2        stats4_4.3.2           xfun_0.52             
# [103] qgraph_1.9.8           matrixStats_1.5.0      stringi_1.8.7         
# [106] lazyeval_0.2.2         yaml_2.3.10            boot_1.3-31           
# [109] evaluate_1.0.3         codetools_0.2-19       tibble_3.2.1          
# [112] cli_3.6.5              uwot_0.2.3             rpart_4.1.24          
# [115] xtable_1.8-4           reticulate_1.42.0      Rdpack_2.6.4          
# [118] lavaan_0.6-19          dichromat_2.0-0.1      Rcpp_1.0.14           
# [121] spatstat.random_3.4-1  globals_0.18.0         png_0.1-8             
# [124] spatstat.univar_3.1-3  parallel_4.3.2         dotCall64_1.2         
# [127] jpeg_0.1-10            lme4_1.1-37            listenv_0.9.0         
# [130] viridisLite_0.4.2      scales_1.4.0           ggridges_0.5.4        
# [133] purrr_1.0.4            rlang_1.1.6            cowplot_1.1.3         
# [136] mnormt_2.1.1          
