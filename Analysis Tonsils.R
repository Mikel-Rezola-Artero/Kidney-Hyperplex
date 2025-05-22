####1. LOAD LIBRARIES####
library(rstudioapi)
library(dplyr)
library(mxnorm)
library(readxl)
library(writexl)
library(ggplot2)
library(Seurat)
library(Matrix)
library(ggpubr)
library(tibble)
library(phenoptr)
library(plotly)
library(factoextra)
library(cluster)
library(corrplot)
library(dbscan)
library(psych)
library(Matrix)
library(corrplot)
library(qgraph)


####2. IMPORT DATA####

#Set directory to script location & obtain files names
setwd(dirname(getActiveDocumentContext()$path))
files <- list.files(path = "Data/",#Access folder with data
                    pattern = "Tonsil[1-4]-Square.*\\.csv$",#Select relevant files 
                    full.names = TRUE)#Take the names of the files for analysis

#Import cell-level data from HALO Analysis
tonsil_list <- lapply(files, read.csv2,sep = ",", check.names = F)#Import data as list
names(tonsil_list) <- paste0("Tonsil",1:4)#Name list elements by environment variable name

####3. Q.C.: CHECK OF AREAS SELECTED FROM EACH IMAGE####

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

#There are some "white gaps", which are the arteries/veins

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
#Immune Lympho: "CD20" "CD3" "CD8" "CD4" "IgG" 
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
#Immune Lympho: "CD20" "CD3" "CD8" "CD4" "IgG" 
#Stroma: "CD31" "aSMA"
#Cell State: "Ki67"

#Common marker names in supervised and non-supervised data
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
#Merge intesity datasets for common Abs
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
#Merge Supervised Classifications datasets for common Abs
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

#Visualize the results of the Otsu miss-classification analysis: raw vs corrected
ggplotly(plot_mx_discordance(mx_norm))#MPO seems not completely corrected

#Visualize the densities of the marker values normalized vs raw
density_plot<- plot_mx_density(mx_norm)
ggsave(filename = file.path("Outputs/", "density_plot_wide.png"),
       plot = density_plot,
       width = 35,  # Increase width to fit facets
       height = 6,
       dpi = 300)#Check png file for plot output

#Check % of positive cells for every marker according to supervised labeling
(apply(mx_norm$norm_data[grep("Pos",colnames(mx_sample), value = TRUE)],2,sum)/
    nrow(mx_norm$norm_data[grep("Pos",colnames(mx_sample), value = TRUE)]))*100 
#1% > cells classified as positive for MPO which could explain the lack of batch correction

#Do UMAP & check how the batch correction integration did
mx_norm = run_reduce_umap(mx_norm,
                          table="both",
                          marker_list = colnames(mx_norm$norm_data)[c(7,11,16,
                                                                      20,24,25)],
                          downsample_pct = 0.2,
                          metadata_col = colnames(mx_norm$norm_data)[c(26:48)])

#Q.C., Plot UMAP with raw and batch corrected data
plot_mx_umap(mx_norm,metadata_col = "slide")

#Check if main cell markers present logical distributions in batch corrected data
# plot_mx_umap(mx_norm,metadata_col = "CD31.Pos")
# plot_mx_umap(mx_norm,metadata_col = "CD20.Pos")
# plot_mx_umap(mx_norm,metadata_col = "CD3.Pos")
# plot_mx_umap(mx_norm,metadata_col = "CD4.Pos")
# plot_mx_umap(mx_norm,metadata_col = "CD8.Pos")
# plot_mx_umap(mx_norm,metadata_col = "CD68.Pos")

#Careful!!! Batch normalization is not suitable to compare Intensities across 
#slides because over corrects markers expected to be different across conditions.
#This is why for comparison between Intensity Marker between slides for a cell
#cluster needs to be done with normalized Intensity (non-batch effect corrected) 
#adjusted by batch effect as a covariate
#For example by Logistic Regression with a coefficient for the marker and
#coefficient for batch predicting cell cluster for example

#Store Non-Batch corrected data for potential comparisons across slides
write_xlsx(x = mx_norm$data, 
           path = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),
                         "/Outputs/RawIFdata.xlsx"))

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
                     npcs = 6,# = 1 - n features 
                     approx = FALSE)
Seurat_mxnorm <- RunPCA(Seurat_mxnorm, 
                        assay = "RNA",
                        features = c("CD68","CD20","CD3","CD31",
                                     "aSMA","CD45","Nuc.size"),
                        npcs = 6,# = 1 - n features 
                        approx = FALSE)

#Determine the ‘dimensionality’ of the dataset: choosing the number of PCs to used for further analysis.
ElbowPlot(Seurat_raw) | ElbowPlot(Seurat_mxnorm) #we will use the 6 PCs

#Run UMAP
Seurat_raw <- RunUMAP(Seurat_raw, dims = 1:6, n.components = 2)
Seurat_mxnorm <- RunUMAP(Seurat_mxnorm, dims = 1:6, n.components = 2)

#Check batch effect
p1 <- DimPlot(object = Seurat_raw,group.by = "slide",
              label = T,label.box = T,label.size = 5,raster = T,
              repel = T, pt.size = 2, reduction = "umap")+ NoLegend()
p2 <- DimPlot(object = Seurat_mxnorm,group.by = "slide",
              label = T,label.box = T,label.size = 5,raster = T,
              repel = T, pt.size = 2, reduction = "umap")+ NoLegend()
p1 | p2 #Plot Batch Correction Effect

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
          "CD68","C1q","C5aR1","CD163"#Myeloid
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
saveRDS(Seurat_mxnorm, file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),
                                      "/Data/Processed_Tonsils.rds"))

####10. MAKE NON-SPATIAL FIGURES USING CLUSTER INFORMATION####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),
                                 "/Data/Processed_Tonsils.rds"))
#FIGURE: Plot clusters
p1 <- DimPlot(object = Seurat.obj,repel = T,
              label = T,label.box = T,label.size = 5,
              cols = c("red","gold","mediumseagreen","skyblue","purple","pink"),
              pt.size = 2,reduction = "umap") + NoLegend()
p1

#FIGURE: Plot Marker intensities as FeaturePlot
FeaturePlot(object = Seurat.obj,
            pt.size = 1, repel = T,
            reduction = "umap",label = T,
            features = c("CD31","CD20","CD3",
                         "CD4","CD8","CD68"),
            min.cutoff = "q40",max.cutoff = "q99",
            ncol = 3)

#FIGURE: Plot Marker intensities as DotPlot
DotPlot(Seurat.obj,
        c("CD3",#adaptive T cells
          "CD20","Ki67",#adaptive replicating B cells
          #"Nuc.size",#all
          "CD31",#endothelial
          "CD68","C1q","C5aR1"#Myeloid
        ),
        scale = T,
        col.min = 0.1,
        col.max = 1,
        dot.min = F) + 
  scale_size(range = c(15, 15)) + guides(size = "none") +
  scale_color_gradient(low = "white",
                       high = "blue",
                       limits = c(0.1, 1),
                       oob = scales::squish)

# #Only after DimRed we can add to data slot our pseudolog10 raw IF data for C deposits (optional)
# rawdat <- read_xlsx(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/Outputs/RawIFdata.xlsx"))
# 
# Seurat.obj@assays$RNA$data <- as(Matrix(as.matrix(t(log10(data.frame(rawdat[,c(3:25)],
#                                                                          row.names = colnames(Seurat.obj@assays$RNA$data))+1))), 
#                                         sparse = TRUE), "dgCMatrix")
# #This data is only suitable to compare within slides
# # Check non-batch corrected C staining as VlnPlot
# VlnPlot(Seurat.obj[,Seurat.obj$slide == "Tonsil1.meta"],
#          c("C1q",
#             "C3c",
#            "C3d",
#            "C4d",
#            "C9"),
#          pt.size = 0,
#          slot = "data") & theme(text = element_text(size = 10))
# # Check batch corrected C staining as VlnPlot
# VlnPlot(Seurat.obj[,Seurat.obj$slide == "Tonsil1.meta"],
#         c("C1q",
#           "C3c",
#           "C3d",
#           "C4d",
#           "C9"),
#         pt.size = 0,
#         slot = "count") & theme(text = element_text(size = 10))

#No clear visualization of Complement deposits when using continuous variables,
#use metadata to determine positive and negative cells (Supervised Classification)


#EXAMPLE OF MARKER COMPARISON: Plot Ki67 by cluster and compare "B cell Ki67 High 1" vs all
my_comparisons <- list( c("B cell Ki67 High 1", "CD4/CD8 T cell"),
                        c("B cell Ki67 High 1", "B cell Ki67 High 2"),
                        c("B cell Ki67 High 1", "B cell Ki67 Neg"),
                        c("B cell Ki67 High 1", "Endothelial"),
                        c("B cell Ki67 High 1", "Myeloid"))
VlnPlot(Seurat.obj,
        "Ki67",slot = "counts",
        pt.size = 0) + 
  geom_boxplot(
    width = 0.2,        # Adjust boxplot width
    outlier.shape = NA, # Hide outlier points for a cleaner look
    fill = "white",     # Boxplot fill color
    color = "black"     # Boxplot border color
  ) + labs(y = "Intensity") +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test", 
                     p.adjust.method = "bonferroni") +
  ylim(0,2) + guides(fill = "none")

#The problem of this type of test with very high number of points per group is 
#that everything is significant even with small Fold Changes:
table(Idents(Seurat.obj))
#The group with less cells has 22749 cells!!!

####11. MAKE SPATIAL FIGURES USING CLUSTER INFORMATION####

#Load original data files
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
files <- list.files(path = "Data/",#Access folder with data
                    pattern = "Tonsil[1-4]-Square.*\\.csv$",#Select relevant files 
                    full.names = TRUE)#Take the names of the files for analysis
tonsil_list <- lapply(files, read.csv2,sep = ",", check.names = F)#Import data as list
names(tonsil_list) <- paste0("Tonsil",1:4)#Name list elements by environment variable name

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

#FIGURE: Plot Tonsil clusters spatially
# Combine 4 plots into a 2x2 grid
(plot_list[[1]] | plot_list[[2]]) /
  (plot_list[[3]] | plot_list[[4]])


#FIGURE: Plot Tonsil 4 clusters spatially
plot_list[[4]]


#Plot Tonsil 4 C deposits spatially
#Loop through all Tonsil slides
Complement <- c("C1q.Pos","C4d.Pos","C3d.Pos","C9.Pos")
Colors <- c("brown","blue","red","black")

for (i in 1:4) {
  
  Comp <- Complement[i]
  pos_color <- Colors[i]
  
  df <- tonsil_list[["Tonsil4"]]
  
  #Ensure columns are numeric
  df$x <- (as.numeric(df$XMin) + as.numeric(df$XMax))/2
  df$y <- (as.numeric(df$YMin) + as.numeric(df$YMax))/2
  
  df$Comp <- Seurat.obj@meta.data[Seurat.obj$slide == "Tonsil4", Comp]
  
  
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
            size = 16 ))
  
  plot_list[[i]] <- p
  
}

#FIGURE: Plot Tonsil 4 clusters spatially
# Combine 4 plots into a 2x2 grid
(plot_list[[1]] | plot_list[[2]]) /
  (plot_list[[3]] | plot_list[[4]])

####12. INFER FOLLICLES WITH MACHINE LEARNING & CORRELATE DEPOSITS IN TONSIL 4####

#Extract Tonsil 4 data
Tonsil4 <- tonsil_list[["Tonsil4"]]
df <- data.frame(x = (as.numeric(Tonsil4$XMin) + as.numeric(Tonsil4$XMax))/2,#x axis
                 y = (as.numeric(Tonsil4$YMin) + as.numeric(Tonsil4$YMax))/2,#y axis
                 Phenotype = Idents(Seurat.obj)[Seurat.obj$slide == "Tonsil4"],#Idents
                 IDs = colnames(Seurat.obj)[Seurat.obj$slide == "Tonsil4"]#Cell IDs
)
#Extract only replicating B cells
unique(df$Phenotype)#we don't include B cell Ki67 Neg because they are in "mantle zone"
df.clust <- df[df$Phenotype == "B cell Ki67 High 1" |
                 df$Phenotype == "B cell Ki67 High 2",
               1:2]
rownames(df.clust) <- df$IDs[df$Phenotype == "B cell Ki67 High 1" |
                               df$Phenotype == "B cell Ki67 High 2"]

#Set the graphics to display two plots
par(mfrow = c(1,2))

#How many follicles?
plot(x = as.numeric(df.clust$x),#x axis
     y = as.numeric(df.clust$y),#y axis
     pch = 16,cex = .5,
     main = "B cells Ki67 Positive")
#We see around 10 follicles

#Run hdbscan to automatically detect follicles
#We inspect the 100-NN distance plot to decide minPts.
#We want minimum 100 B cells to form a cluster/follicle
#We run clustering with hdbscan with minPts = 100 B cells:
db <- hdbscan(df.clust,minPts = 100)
num_clusters <- max(db$cluster)#unclassified cells are defined as "0" 
cluster_colors <- c("white", rainbow(num_clusters))
plot(x = as.numeric(df.clust$x),#x axis
     y = as.numeric(df.clust$y),#y axis
     pch = as.numeric(db$cluster),cex = .5,
     col = cluster_colors[db$cluster + 1], xlab = "x",ylab = "y",
     main = "Automatic Follicle Annotation with Hierarchical DBSCAN")
#The algorithm detects 10 follicles

#Set the graphics back to normal
par(mfrow = c(1,1))


#Correlation between complement deposits across follicles

#Process data
corr.data <- data.frame(Tonsil4[df$Phenotype == "B cell Ki67 High 1" |
                                  df$Phenotype == "B cell Ki67 High 2",
                                c(13,19,23,29,39,49,59,63,67,71,121)]) #extract deposits info on follicular B cells
corr.data$cluster <- db$cluster
corr.data <- corr.data %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))#calculate frequency of positivity by follicle
#Correct column names
colnames(corr.data) <- gsub("[_.].*Positive.Classification", "", colnames(corr.data))
colnames(corr.data) <- gsub("Abcam", "", colnames(corr.data))
colnames(corr.data)[3] <- "FH"
#Remove unclassified cells
corr.data <- corr.data[!corr.data$cluster == 0,-1]

#Plot correlations
M1 <- cor(corr.data,
          use = "pairwise.complete.obs",
          method = "spearman")#calculate correlations
testRes <- corr.test(corr.data,
                     use = "pairwise",
                     method="spearman",
                     adjust="fdr", 
                     alpha=.05,
                     ci=TRUE,
                     minlength=5,
                     normal=TRUE)#calculate p-values
corrplot(M1,p.mat = as.matrix(forceSymmetric(testRes$p,uplo = "U")), insig = "blank" ,
         sig.level = 0.05,
         order = 'hclust',
         tl.cex = 1.25,
         addrect = 1,
         type = 'upper',
         tl.srt = 50,
         mar=c(0,0,0,0))#plot as corrplot

#Plot as correlation network
colnames(M1) <- sub("_Pos","",colnames(M1))
rownames(M1) <- colnames(M1)
qgraph::qgraph(M1,
               layout="spring",esize = 40,vsize=15,
               color = c("purple","green","brown","red",
                         "lightblue","skyblue"),
               borders = T,
               threshold = 0.4,#show only correlations above 0.4
               label.cex = 1.2,
               labels = colnames(M1),
               graph = "cor",
               edge.width = 0.5,
               node.resolution = 100,
               repulsion = 0.8,
               label.scale=F,
               curved = TRUE,
               graph = "association",
               posCol = c("brown"),negCol =c("blue"))






####13. PLOT DEPOSITS AND CLUSTERS IN TONSIL 4####

#FIGURE: Plot C9+ cells, C5aR1+ myeloid cells and C5aR1- myeloid cells in Tonsil4
Tonsil4 <- tonsil_list[["Tonsil4"]]
df <- data.frame(x = (as.numeric(Tonsil4$XMin) + as.numeric(Tonsil4$XMax))/2,#x axis
                 y = (as.numeric(Tonsil4$YMin) + as.numeric(Tonsil4$YMax))/2,#y axis
                 Phenotype = Idents(Seurat.obj)[Seurat.obj$slide == "Tonsil4"],#Idents
                 IDs = colnames(Seurat.obj)[Seurat.obj$slide == "Tonsil4"]#Cell IDs
)
pheno <- as.character(df$Phenotype)
pheno[Tonsil4$`C9 - Cy5 Positive Classification` == 1 &
        Tonsil4$`C4d - Cy5 Positive Classification` == 1] <- "C4d+/C9+ cell"
pheno[df$Phenotype == "Myeloid" &
        Tonsil4$`C5aR1 - Cy5 Positive Classification` == 1] <- "Myeloid C5aR1+"
pheno[df$Phenotype == "Myeloid" &
        Tonsil4$`C5aR1 - Cy5 Positive Classification` == 0] <- "Myeloid C5aR1-"
df$Phenotype <- pheno
df <- df[df$Phenotype != "Endothelial" & 
           df$Phenotype != "CD4/CD8 T cell",]
ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color = as.factor(Phenotype)),
             alpha = 0.4, size = 0.75) +
  theme_classic() + theme() + 
  labs(title="Tonsil 4") + 
  scale_color_manual(values = c("orange","yellow","lightblue","red","brown","black"),# Apply custom colors
                     guide = guide_legend(override.aes = list(size = 4)))

####14. SPATIAL DISTANCE ANALYSIS####

####Calculate distances between Seurat clusters in Tonsil4
#Calculate cell centroids
Tonsil4 <- tonsil_list[["Tonsil4"]]
df <- data.frame(x = (as.numeric(Tonsil4$XMin) + as.numeric(Tonsil4$XMax))/2,#x axis
                 y = (as.numeric(Tonsil4$YMin) + as.numeric(Tonsil4$YMax))/2#y axis
)
#Add Seurat clusters
df$clusters <- Idents(Seurat.obj)[Seurat.obj$slide == "Tonsil4"]
#Create object for distance calculations
cds <- df
# cds <- cds[1:3]
colnames(cds) <- c("Cell X Position","Cell Y Position","Phenotype")
cds$"Cell ID" <- 1:nrow(cds)
cds$Phenotype <- as.character(cds$Phenotype)
cds <- as_tibble(cds)
#Calculate distances between cell cluster based on Nearest Neighbor distance
#This takes time and a lot of RAM!!!!!
distances <- find_nearest_distance(cds
                                   #[cds$`Cell X Position` < 15000 & #You can subset an area of the image for this
                                   #     cds$`Cell Y Position` > 12500,]#You can subset an area of the image for this
)
gc()#to liberate memory

#Plot distances between Myeloid cells and the other clusters in Tonsil 4
csd_with_distance <- bind_cols(cds
                               #[cds$`Cell X Position` < 15000 & #You can subset an area of the image for this
                               #       cds$`Cell Y Position` > 12500,]#You can subset an area of the image for this
                               , distances)
csd_with_distance %>% group_by(Phenotype) %>% 
  select(Phenotype, starts_with('Distance to')) %>% 
  summarize_all(~round(mean(.), 1))
ggplot(csd_with_distance[csd_with_distance$Phenotype != "Myeloid",], aes(`Distance to Myeloid`, color=Phenotype)) +
  geom_density(size=1) + 
  scale_color_manual(values = c("orange","yellow","blue","brown","red"))   # Apply custom colors

#For more dist measurements, follow tutorial in:https://akoyabio.github.io/phenoptr/articles/computing_distances.html

#FIGURE: Plot distances between C9+ cells and C5aR1+ and C5aR1- Myeloid cells in Tonsil4
pheno <- unlist(cds$Phenotype)
pheno[Tonsil4$`C9 - Cy5 Positive Classification` == 1 &
        Tonsil4$`C4d - Cy5 Positive Classification` == 1] <- "C4d+/C9+ cell"
pheno[cds$Phenotype == "Myeloid" &
        Tonsil4$`C5aR1 - Cy5 Positive Classification` == 1] <- "Myeloid C5aR1+"
pheno[cds$Phenotype == "Myeloid" &
        Tonsil4$`C5aR1 - Cy5 Positive Classification` == 0] <- "Myeloid C5aR1-"
C9.cells <- cds
C9.cells$Phenotype <- pheno
C9.cells <- C9.cells %>% filter(Phenotype!='Endothelial' &
                                  Phenotype!='CD4/CD8 T cell')
dist.C9 <- find_nearest_distance(C9.cells)#takes less time because we filtered out cells
C9_with_distance <- bind_cols(C9.cells, dist.C9)
ggplot(C9_with_distance[C9_with_distance$Phenotype != "C4d+/C9+ cell",], aes(`Distance to C4d+/C9+ cell`, color=Phenotype)) +
  geom_density(size=1) + ggplot2::scale_y_continuous("Cell probability within each cluster") +
  scale_color_manual(values = c("red","yellow","orange","grey","black"))   # Apply custom colors

#Mean distance from B cells in C5aR1+ and C5aR1- cells
median(dist.C9$`Distance to C4d+/C9+ cell`[C9_with_distance$Phenotype == "Myeloid C5aR1-"])
median(dist.C9$`Distance to C4d+/C9+ cell`[C9_with_distance$Phenotype == "Myeloid C5aR1+"])
median(dist.C9$`Distance to C4d+/C9+ cell`[C9_with_distance$Phenotype == "B cell Ki67 High 1"])
median(dist.C9$`Distance to C4d+/C9+ cell`[C9_with_distance$Phenotype == "B cell Ki67 High 2"])
median(dist.C9$`Distance to C4d+/C9+ cell`[C9_with_distance$Phenotype == "B cell Ki67 Neg"])

#Plot as histogram
# ggplot(C9_with_distance[C9_with_distance$Phenotype != "C4d+/C9+ cell",],
#       aes(x = `Distance to C4d+/C9+ cell`, fill = Phenotype)) +
#  geom_histogram(aes(y = after_stat(count) / tapply(after_stat(count), after_stat(fill), sum)[after_stat(fill)]),
#                 bins = 30, alpha = 0.5, position = "identity") +
#  scale_y_continuous("Cell percentage within each cluster", labels = scales::percent_format()) +
#  scale_fill_manual(values = c("red","yellow","orange","grey","black")) +
#  theme_minimal()


####Illustrate calculated distances in Tonsil4
#Re-align Tonsil4 data for image overlapping
csd_with_distance2 <- C9_with_distance
csd_with_distance2$`Cell X Position` <- C9_with_distance$`Cell X Position`- min(C9_with_distance$`Cell X Position`)
csd_with_distance2$`Cell Y Position` <- C9_with_distance$`Cell Y Position`- min(C9_with_distance$`Cell Y Position`)
# Filter to just Myeloid C5aR1+ cells
Myeloid_cells = csd_with_distance2 %>% filter(select_rows(C9_with_distance, "Myeloid C5aR1+"))
Comp_cells = csd_with_distance2 %>% filter(select_rows(C9_with_distance, "C4d+/C9+ cell"))
# For each C4d+/C9+ cell, join with the data for the nearest Myeloid C5aR1+ cell
Comp_to_Myeloid = Comp_cells %>% left_join(Myeloid_cells, by=c('Cell ID Myeloid C5aR1+'='Cell ID'),
                                      suffix=c('', '.B'))
# Read a background image and make a base plot (DAPI image of slide)
background_path = 
  system.file("~/Idris stuff/Integrated complementomics/CCE Lupus DN/Tonsil IFA multiplex/Analysis Paper/Data/Tonsil 4.jpg", package='phenoptr')
background = jpeg::readJPEG("~/Idris stuff/Integrated complementomics/CCE Lupus DN/Tonsil IFA multiplex/Analysis Paper/Data/Tonsil 4.jpg") %>% as.raster()
xlim = c(0, 12989.5)#max X axis - min X axis for "csd_with_distance2"
ylim = c(0, 10061.5)#max Y axis - min Y axis "csd_with_distance2"
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
  labs(title='Nearest Myeloid C5aR1+ to each C4d+/C9+ cell')


# ####RECOVER CELL TYPE LABELS FOR HALO VALIDATION####
# Tonsil.pheno <- data.frame("slide" = Seurat_pseudo2$slide)
# Tonsil.pheno$"B cell"[Idents(Seurat_pseudo2) == "B cell"] <- 1
# Tonsil.pheno$"B cell"[Idents(Seurat_pseudo2) != "B cell"] <- 0
# Tonsil.pheno$"T cell"[Idents(Seurat_pseudo2) == "CD4/CD8 T cell"] <- 1
# Tonsil.pheno$"T cell"[Idents(Seurat_pseudo2) != "CD4/CD8 T cell"] <- 0
# Tonsil.pheno$"Myeloid"[Idents(Seurat_pseudo2) == "Macro/Mono/Myeloid"] <- 1
# Tonsil.pheno$"Myeloid"[Idents(Seurat_pseudo2) != "Macro/Mono/Myeloid"] <- 0
# Tonsil.pheno$"Endothelial"[Idents(Seurat_pseudo2) == "Endothelial"] <- 1
# Tonsil.pheno$"Endothelial"[Idents(Seurat_pseudo2) != "Endothelial"] <- 0
# library(writexl)
# write_xlsx(Tonsil.pheno[Tonsil.pheno$slide == "Tonsil1",],"Tonsil1.Labels.xlsx")
# write_xlsx(Tonsil.pheno[Tonsil.pheno$slide == "Tonsil2",],"Tonsil2.Labels.xlsx")
# write_xlsx(Tonsil.pheno[Tonsil.pheno$slide == "Tonsil3",],"Tonsil3.Labels.xlsx")
# write_xlsx(Tonsil.pheno[Tonsil.pheno$slide == "Tonsil4",],"Tonsil4.Labels.xlsx")
