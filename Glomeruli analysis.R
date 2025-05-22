####1. LOAD LIBRARIES####
library(rstudioapi)
library(ggplot2)
library(dplyr)
library(mxnorm)
library(plotly)
library(writexl)
library(Seurat)
library(tidyr)
library(ggpubr)
library(scales)
library(corrplot)
library(psych)
library(Matrix)

# library(readxl)
# library(Matrix)
# library(tibble)
# library(phenoptr)
# library(factoextra)
# library(cluster)

# library(dbscan)

# library(corrplot)
# library(qgraph)

####2. IMPORT DATA####

#Set directory to script location & obtain files names
setwd(dirname(getActiveDocumentContext()$path))
files <- list.files(path = "Data/", pattern = "\\.csv$", full.names = TRUE)

#Create list to store the data
Kidney_IF <- list()

#Load DN data
Kidney_IF$Diab1 <- read.csv2(grep("41_DNv2.csv",files, value = TRUE),
                             sep = ",", check.names = F)

Kidney_IF$Diab2 <- read.csv2(grep("_DN.csv",files, value = TRUE)[2],
                   sep = ",", check.names = F)

Kidney_IF$Diab3 <- read.csv2(grep("_DN.csv",files, value = TRUE)[3],
                   sep = ",", check.names = F)

#Load CCE data
Kidney_IF$Chol1 <- read.csv2(grep("5_CCEv2.csv",files, value = TRUE),
                   sep = ",", check.names = F)
#Fix incomplete columns
Kidney_IF[["Chol1"]]$'CD163-50 - TRITC Cell Intensity'[is.na(Kidney_IF[["Chol1"]]$'CD163-50 - TRITC Cell Intensity')] <- Kidney_IF[["Chol1"]]$`CD163-50 Cell Intensity`[Kidney_IF[["Chol1"]]$'CD163-50 - TRITC Cell Intensity' == ""]
Kidney_IF[["Chol1"]]$`CD163-50 - TRITC Positive Classification`[is.na(Kidney_IF[["Chol1"]]$`CD163-50 - TRITC Positive Classification`)] <- Kidney_IF[["Chol1"]]$`CD163-50 Positive Classification`[is.na(Kidney_IF[["Chol1"]]$`CD163-50 - TRITC Positive Classification`)]
Kidney_IF$Chol1 <- Kidney_IF$Chol1[1:131]#keep first 131 columns

Kidney_IF$Chol2 <- read.csv2(grep("_CCEv3.csv",files, value = TRUE),
                   sep = ",", check.names = F)

Kidney_IF$Chol3 <- read.csv2(grep("_CCE.csv",files, value = TRUE)[3],
                   sep = ",", check.names = F,
                   colClasses = c(rep(NA, 79),
                                  rep("NULL", 4),#remove columns 80-83
                                  rep(NA, 34),
                                  rep("NULL", 10)#remove columns 118-127
                                  ))

#Load Lupus data
Kidney_IF$Lupus1 <- read.csv2(grep("_Lupus.csv",files, value = TRUE)[1],
                    sep = ",", check.names = F,
                    colClasses = c(rep(NA, 162),rep("NULL", 2))#keep first 162 columns
                    )

Kidney_IF$Lupus2 <- read.csv2(grep("_Lupus.csv",files, value = TRUE)[2],
                    sep = ",", check.names = F,
                    colClasses = c(rep(NA, 141),rep("NULL", 8))#keep first 141 columns
                    )

Kidney_IF$Lupus3 <- read.csv2(grep("_Lupus.csv",files, value = TRUE)[3],
                    sep = ",", check.names = F)

Kidney_IF$Lupus4 <- read.csv2(grep("_Lupus.csv",files, value = TRUE)[4],
                    sep = ",", check.names = F)

Kidney_IF$Lupus5 <- read.csv2(grep("_Lupus.csv",files, value = TRUE)[5],
                    sep = ",", check.names = F,
                    colClasses = c(rep(NA, 162),rep("NULL", 5))#keep first 162 columns
                    )

Kidney_IF$Lupus6 <- read.csv2(grep("_Lupus.csv",files, value = TRUE)[6],
                    sep = ",", check.names = F)

Kidney_IF$Lupus7 <- read.csv2(grep("_Lupusv2.csv",files, value = TRUE),
                    sep = ",", check.names = F,
                    colClasses = c(rep(NA, 156),rep("NULL", 2))#keep first 156 columns
                    )

Kidney_IF$Lupus8 <- read.csv2(grep("_Lupus.csv",files, value = TRUE)[8],
                    sep = ",", check.names = F,
                    colClasses = c(rep(NA, 158),rep("NULL", 2))#keep first 158 columns
                    )

Kidney_IF$Lupus9 <- read.csv2(grep("_Lupus.csv",files, value = TRUE)[9],
                    sep = ",", check.names = F)

####3. FILTER GLOMERULI & CHECK AREAS SELECTED FROM EACH IMAGE####

#Filter Glomeruli cells
Glom_IF <- lapply(Kidney_IF, function(df) {
 df[grep("Glomeruli", df$`Analysis Region`, value = FALSE), ]
  })

#Plot Glomeruli for Q.C.
plot_data <- do.call(rbind, lapply(names(Glom_IF), function(name) {
  df <- Glom_IF[[name]][,c(5,7)]
  cbind(df, Dataset = name)}))
ggplot(plot_data, aes(x = as.numeric(`XMin`), y = as.numeric(`YMin`))) +
  geom_point(aes(color = Dataset), size = 0.01, alpha = 0.1) +
  facet_wrap(~Dataset, ncol = 3) +
  # coord_equal() +
  labs(x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")
#We clearly see glomerular structures distributed spatially in our filtered data

####4. Q.C.: COMMON Abs/MARKERS ACROSS DATASETS FOR CLUSTERING####

#Check raw column names
Reduce(intersect, lapply(Glom_IF, colnames))#No consensus in Ab naming

#Take Cell Intensity Columns
Glom_intensity <- lapply(Glom_IF, function(df) {
  df %>% select(contains("Cell Intensity"))
})

#Rename elements to create variables for each intensity object
names(Glom_intensity) <- paste0(names(Glom_intensity),".Int")

#Assign each dataset to a variable
list2env(Glom_intensity, envir = .GlobalEnv)#Assign each dataset to a variable

#Correct Cell intensity column names manually (no common patterns/spacers)


#DN
colnames(Diab1.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Diab1.Int))
colnames(Diab1.Int) <- gsub("Rb|Ms", "", colnames(Diab1.Int))
colnames(Diab1.Int)[c(16,19)] <- c("Megalin","VCAM-1")

colnames(Diab2.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Diab2.Int))
colnames(Diab2.Int)[c(25,17)] <- c("ITGA8","VCAM-1")

colnames(Diab3.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Diab3.Int))
colnames(Diab3.Int)[c(25,17)] <- c("ITGA8","VCAM-1")


#CCE
colnames(Chol1.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Chol1.Int))
colnames(Chol1.Int) <- gsub("Rb", "", colnames(Chol1.Int))
colnames(Chol1.Int)[c(15,17,24,20)] <- c("C9","Megalin","Synaptopodin","VCAM-1")

colnames(Chol2.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Chol2.Int))
colnames(Chol2.Int) <- gsub("Rb", "", colnames(Chol2.Int))
colnames(Chol2.Int)[c(14,16,19)] <- c("C9","Megalin","VCAM-1")

colnames(Chol3.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Chol3.Int))
colnames(Chol3.Int) <- gsub("Rb", "", colnames(Chol3.Int))
colnames(Chol3.Int)[c(15,17,14,18)] <- c("C9","Megalin","Synaptopodin","VCAM-1")


#Lupus
colnames(Lupus1.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Lupus1.Int))
colnames(Lupus1.Int)[c(25,17)] <- c("ITGA8","VCAM-1")

colnames(Lupus2.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Lupus2.Int))
colnames(Lupus2.Int) <- gsub("Rb", "", colnames(Lupus2.Int))
colnames(Lupus2.Int)[c(9,11,15)] <- c("Megalin","C9","VCAM-1")

colnames(Lupus3.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Lupus3.Int))
colnames(Lupus3.Int)[c(25,17)] <- c("ITGA8","VCAM-1")

colnames(Lupus4.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Lupus4.Int))
colnames(Lupus4.Int)[c(25,17)] <- c("ITGA8","VCAM-1")

colnames(Lupus5.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Lupus5.Int))
colnames(Lupus5.Int) <- gsub("Cy", "", colnames(Lupus5.Int))
colnames(Lupus5.Int)[c(17,25,28,29,33)] <- c("VCAM-1","IgA","Megalin","CD62P","IgM")

colnames(Lupus6.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Lupus6.Int))
colnames(Lupus6.Int)[c(25,17)] <- c("ITGA8","VCAM-1")

colnames(Lupus7.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Lupus7.Int))
colnames(Lupus7.Int) <- gsub("Cy", "", colnames(Lupus7.Int))
colnames(Lupus7.Int)[c(17,18,21,25,26,28,33)] <- c("Synaptopodin","ITGA8","VCAM-1","MUM1","CD34","Megalin","IgM")

colnames(Lupus8.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Lupus8.Int))
colnames(Lupus8.Int) <- gsub("Cy", "", colnames(Lupus8.Int))
colnames(Lupus8.Int)[c(17,18,21,25,26,28,33)] <- c("Synaptopodin","ITGA8","VCAM-1","MUM1","CD34","Megalin","IgM")

colnames(Lupus9.Int) <- gsub("[-_ ].*Cell Intensity", "", colnames(Lupus9.Int))
colnames(Lupus9.Int) <- gsub("Cy", "", colnames(Lupus9.Int))
colnames(Lupus9.Int)[c(17,18,21,25,26,28,33)] <- c("Synaptopodin","ITGA8","VCAM-1","MUM1","CD34","Megalin","IgM")



#Add the objects with corrected column names back to the intensities list
for (i in 1:15) {
  name <- names(Glom_intensity)[i]
  Glom_intensity[[name]] <- get(name)
  rm(list = name)
}

#Check corrected column names
Reduce(intersect, lapply(Glom_intensity, colnames))

#Markers for Cell Intensity: 
#Complement: "C1q" "MBL" "C3c" "C3d" "C5aR1" "C1r" "C1s" "C4d" "FH" "FB" "C9"
#Immune Myeloid: "CD45" "CD68" "CD163" "CD11b" "MPO"
#Immune Lympho: "CD20" "CD3" "CD8" "CD4" "IgG" 
#Stroma: "CD31" "aSMA" "Megalin"
#Cell State: "Ki67" "VCAM-1"

####5. PROCESS FLUORESCENCE INTENSITIES FOR MXNORM OBJECT CREATION####

for (i in 1:15) {
  int_name <- names(Glom_intensity)[i]
  orig_name <- names(Glom_IF)[i]
  
  # Convert all columns to numeric
  Glom_intensity[[int_name]][] <- lapply(Glom_intensity[[int_name]], as.numeric)
  
  # Remove DAPI from Intensity Data
  Glom_intensity[[int_name]] <- Glom_intensity[[int_name]][, -1]
}

####6. PROCESS SUPERVISED CLASSIFICATION DATA FOR MXNORM OBJECT CREATION####

#Take Cell Classification Columns with manual classification labels
Glom_classification <- lapply(Glom_IF, function(df) {
  df %>% select(contains("Positive Classification"))
})

#Rename elements to create variables for each Cell Classification object
names(Glom_classification) <- paste0(names(Glom_classification),".meta")

#Assign each dataset to a variable
list2env(Glom_classification, envir = .GlobalEnv)#Assign each dataset to a variable


#Correct classification column names manually (no common patterns/spacers)


#DN
colnames(Diab1.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Diab1.meta))
colnames(Diab1.meta) <- gsub("Rb|Ms", "", colnames(Diab1.meta))
colnames(Diab1.meta)[c(16,19)] <- c("Megalin","VCAM-1")

colnames(Diab2.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Diab2.meta))
colnames(Diab2.meta)[c(24,16)] <- c("ITGA8","VCAM-1")

colnames(Diab3.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Diab3.meta))
colnames(Diab3.meta)[c(24,16)] <- c("ITGA8","VCAM-1")


#CCE
colnames(Chol1.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Chol1.meta))
colnames(Chol1.meta) <- gsub("Rb", "", colnames(Chol1.meta))
colnames(Chol1.meta)[c(15,17,24,20)] <- c("C9","Megalin","Synaptopodin","VCAM-1")

colnames(Chol2.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Chol2.meta))
colnames(Chol2.meta) <- gsub("Rb", "", colnames(Chol2.meta))
colnames(Chol2.meta)[c(14,16,19)] <- c("C9","Megalin","VCAM-1")

colnames(Chol3.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Chol3.meta))
colnames(Chol3.meta) <- gsub("Rb", "", colnames(Chol3.meta))
colnames(Chol3.meta)[c(15,17,14,19)] <- c("C9","Megalin","Synaptopodin","VCAM-1")


#Lupus
colnames(Lupus1.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Lupus1.meta))
colnames(Lupus1.meta)[c(24,16)] <- c("ITGA8","VCAM-1")

colnames(Lupus2.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Lupus2.meta))
colnames(Lupus2.meta) <- gsub("Rb", "", colnames(Lupus2.meta))
colnames(Lupus2.meta)[c(9,11,15)] <- c("Megalin","C9","VCAM-1")

colnames(Lupus3.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Lupus3.meta))
colnames(Lupus3.meta)[c(24,16)] <- c("ITGA8","VCAM-1")

colnames(Lupus4.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Lupus4.meta))
colnames(Lupus4.meta)[c(24,16)] <- c("ITGA8","VCAM-1")

colnames(Lupus5.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Lupus5.meta))
colnames(Lupus5.meta) <- gsub("Cy", "", colnames(Lupus5.meta))
colnames(Lupus5.meta)[c(16,23,26,30)] <- c("VCAM-1","IgA","Megalin","IgM")

colnames(Lupus6.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Lupus6.meta))
colnames(Lupus6.meta)[c(24,16)] <- c("ITGA8","VCAM-1")

colnames(Lupus7.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Lupus7.meta))
colnames(Lupus7.meta) <- gsub("Cy", "", colnames(Lupus7.meta))
colnames(Lupus7.meta)[c(16,17,20,24,25,27,32)] <- c("Synaptopodin","ITGA8","VCAM-1","MUM1","CD34","Megalin","IgM")

colnames(Lupus8.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Lupus8.meta))
colnames(Lupus8.meta) <- gsub("Cy", "", colnames(Lupus8.meta))
colnames(Lupus8.meta)[c(16,17,20,24,25,27,32)] <- c("Synaptopodin","ITGA8","VCAM-1","MUM1","CD34","Megalin","IgM")

colnames(Lupus9.meta) <- gsub("[-_ ].*Positive Classification", "", colnames(Lupus9.meta))
colnames(Lupus9.meta) <- gsub("Cy", "", colnames(Lupus9.meta))
colnames(Lupus9.meta)[c(16,17,20,24,25,27,32)] <- c("Synaptopodin","ITGA8","VCAM-1","MUM1","CD34","Megalin","IgM")


#Add the objects with corrected column names back to the classification list
for (i in 1:15) {
  name <- names(Glom_classification)[i]
  Glom_classification[[name]] <- get(name)
  rm(list = name)
}

#Check corrected column names
Reduce(intersect, lapply(Glom_classification, colnames))

#Markers for Cell Classification (NA for IgG, CD4 & FH classifiers): 
#Complement: "C1q" "MBL" "C3c" "C3d" "C5aR1" "C1r" "C1s" "C4d" "C9"
#Immune Myeloid: "CD45" "CD68" "CD163" "CD11b" "MPO"
#Immune Lympho: "CD20" "CD3" "CD8" 
#Stroma: "CD31" "aSMA" "Megalin"
#Cell State: "Ki67" "VCAM-1"

#Common marker names in supervised and non-supervised data
Markers <- intersect(
  Reduce(intersect, lapply(Glom_intensity, colnames)),#Intensity
  Reduce(intersect, lapply(Glom_classification, colnames))#Classification
)

####7. COMBINED OBJECT FOR MXNORM BATCH CORRECTION####

#Fuse Ab intensity data

#Extract common Ab intensities across slides
selected_intensities <- lapply(Glom_intensity, function(df) {
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
Glom_classification,#df
gsub(".meta","",names(Glom_classification)),#name
seq_along(Glom_classification)#idx
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

#Add Glomeruli Region Annotations
glom.info <- lapply(Glom_IF, function (df) {
  gloms <- df$`Analysis Region`
})
mx_sample$Glomeruli <- unlist(glom.info)

#Add Disease Annotations
disease.info <- unlist(
  lapply(names(Glom_IF), function(df_name) {
    # Extract disease name by removing trailing numbers
    disease <- gsub("\\d+$", "", df_name)
    # Repeat for all rows in current dataframe
    rep(disease, nrow(Glom_IF[[df_name]]))
  })
)
disease.info[disease.info == "Diab"] <- "Diabetes"
disease.info[disease.info == "Chol"] <- "CCE"
disease.info[disease.info == "Lupus"] <- "LN"
mx_sample$Disease <- disease.info


####8. PERFORM BATCH CORRECTION WITH MXNORM####

#Modify VCAM-1 name to make it compatible with mxnorm
colnames(mx_sample) <- gsub("-","",colnames(mx_sample))

#Create mxnorm object
mx_norm = mx_dataset(data=mx_sample,
                     slide_id="slide",
                     image_id="image",
                     marker_cols=colnames(mx_sample)[1:21],
                     metadata_cols=c(grep("Pos",colnames(mx_sample),value = TRUE),
                                          "Glomeruli", "Disease")) 
                                     
summary(mx_norm)#Check summary data

#Correct scale intensity batch effect using "log10_mean_divide"
mx_norm = mx_normalize(mx_data = mx_norm,
                       transform = "log10_mean_divide",
                       method="None")
summary(mx_norm)#Check summary data post batch correction
#The test stat has diminished by ~1/2

#Run an Otsu discordance score analysis to determine how well our normalization method performs
mx_norm = run_otsu_discordance(mx_norm,
                               table="both",
                               threshold_override = NULL,
                               plot_out = FALSE)

#Visualize the results of the Otsu miss-classification analysis: raw vs corrected
ggplotly(plot_mx_discordance(mx_norm))#some markers are not completely corrected

#Check % of positive cells for every marker according to supervised labeling
(apply(mx_norm$norm_data[grep("Pos",colnames(mx_sample), value = TRUE)],2,sum)/
    nrow(mx_norm$norm_data[grep("Pos",colnames(mx_sample), value = TRUE)]))*100 
#The markers with lack of batch correction correspond with very low positivity frequency (< 1%)


#Visualize the densities of the marker values normalized vs raw
density_plot<- plot_mx_density(mx_norm)
ggsave(filename = file.path("Only glomeruli/", "density_plot_wide.png"),
       plot = density_plot,
       width = 35,  # Increase width to fit facets
       height = 6,
       dpi = 300)#Check png file for plot output


#Do UMAP & check how the batch correction integration did
mx_norm = run_reduce_umap(mx_norm,
                          table="both",
                          marker_list = colnames(mx_norm$norm_data)[c(6,15,17,18,19)],
                          downsample_pct = 0.2,
                          metadata_col = colnames(mx_norm$norm_data)[c(24:44)])

#Q.C., Plot UMAP with raw and batch corrected data
plot_mx_umap(mx_norm,metadata_col = "slide")

#Check if main cell markers present logical distributions in batch corrected data
# plot_mx_umap(mx_norm,metadata_col = "CD31.Pos")
# plot_mx_umap(mx_norm,metadata_col = "CD45.Pos")
# plot_mx_umap(mx_norm,metadata_col = "VCAM1.Pos")
# plot_mx_umap(mx_norm,metadata_col = "aSMA.Pos")
# plot_mx_umap(mx_norm,metadata_col = "CD68.Pos")

#Careful!!! Batch normalization is not suitable to compare Intensities across 
#slides because over corrects markers expected to be different across conditions.
#This is why for comparison between Intensity Marker between slides for a cell
#cluster needs to be done with normalized Intensity (non-batch effect corrected) 
#adjusted by batch effect as a covariate
#For example by Logistic Regression with a coefficient for the marker and
#coefficient for batch predicting cell cluster for example

#Store Non-Batch corrected data for potential comparisons across slides
write_xlsx(x = mx_norm$data, path = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),
                                           "/Only glomeruli/RawIFdataGlomeruli.xlsx"))


####9. CREATE SEURAT OBJECT AND IDENTIFY CELL TYPES#####

#Batch corrected data
seurat_int <- data.frame(mx_norm$norm_data[,c(3:23)])#extract batch corrected intensity data
seurat_int[is.na(seurat_int)]#check for missing data
seurat_int$Nuc.size <- unlist(lapply(Glom_IF, function(df) df[["Nucleus Area (µm²)"]]))#add nuclear size data

#Create Seurat object for clustering
Seurat_mxnorm <- CreateSeuratObject(counts = t(seurat_int),
                                    meta.data = data.frame(cbind(Merged.Meta,mx_sample[c(45,46)]),
                                                           row.names = rownames(seurat_int)))

#Create "data" slot with the batch corrected/normalized values
Seurat_mxnorm@assays$RNA$data <- Seurat_mxnorm@assays$RNA$counts

#Set variable features to all continuous data
VariableFeatures(Seurat_mxnorm) <- rownames(Seurat_mxnorm@assays$RNA$counts)

#Scale data prior to dimensional reduction
Seurat_mxnorm <- ScaleData(Seurat_mxnorm,scale.max = 5)

#Run PCA using the features determined
Seurat_mxnorm <- RunPCA(Seurat_mxnorm, 
                        assay = "RNA",
                        features = c("CD31",#Endothelial
                                     "CD45",#Immune
                                     "aSMA",#Mesangium
                                     "VCAM1",#PECs
                                     "Nuc.size"#All
                        ),
                        npcs = 5,
                        approx = FALSE)

#Run UMAP
Seurat_mxnorm <- RunUMAP(Seurat_mxnorm, dims = 1:5, n.components = 2)

##Find Clusters using the Dimensional Reduction Markers
Seurat_mxnorm <- FindNeighbors(Seurat_mxnorm, dims = 1:5)
Seurat_mxnorm <- FindClusters(Seurat_mxnorm,resolution = 0.2)

#Check clusters
DimPlot(object = Seurat_mxnorm,repel = T,
        label = T,label.box = T,label.size = 5,
        pt.size = 2,reduction = "umap") + NoLegend()

#Check distribution of key cell markers as DotPlot
DotPlot(Seurat_mxnorm,
        c("CD31",#Endothelial
          "CD45",#Immune
          "CD68","C5aR1",#Myeloid
          "VCAM1",#PECs
          "Nuc.size",#All
          "aSMA"#Mesangium
        ),
        scale = T,
        col.min = 0.7,
        col.max = 2,cols = c("white","blue"),
        dot.min = F) 

#Check distribution of key cell markers
FeaturePlot(object = Seurat_mxnorm,raster = T,
            pt.size = 2, repel = T,ncol = 3,label.size = 3,
            reduction = "umap",label = T,
            features = c("CD31",#Endothelial
                         "CD45",#Immune
                         "VCAM1",#PECs
                         "CD68",#Myeloid
                         "aSMA",#Mesangium
                         "CD3"#T cells
            ),min.cutoff = "q50",max.cutoff = "q100")

#Check distribution of key cell markers as supervised labels
FeaturePlot(object = Seurat_mxnorm,raster = T,
            pt.size = 2, repel = T,ncol = 3,label.size = 3,
            reduction = "umap",label = T,
            features = c("CD31.Pos",#Endothelial
                         "CD45.Pos",#Immune
                         "VCAM.1.Pos",#PECs
                         "CD68.Pos",#Myeloid
                         "aSMA.Pos",#Mesangium
                         "CD3.Pos"#T cells
            ))

#Cluster 0 = PECs/PT-Conv
#Cluster 1 = Endothelial
#Cluster 2 = Unknown 
#Cluster 3 = Unknown
#Cluster 4 = Immune
#Cluster 5 = Mesangial
#Cluster 6 = Unknown

#Rename clusters according to Cell Markers
new_cluster_names <- c("PECs/PT-Conv",
                       "Endothelial",
                       "Unknown",
                       "Unknown",
                       "Immune",
                       "Mesangial",
                       "Unknown")
names(new_cluster_names) <- levels(Seurat_mxnorm)
Seurat_mxnorm <- RenameIdents(Seurat_mxnorm, 
                              new_cluster_names)

#Store Processed Seurat Object with labels
saveRDS(Seurat_mxnorm, file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),
                                      "/Only glomeruli/Processed_Glomeruli.rds"))

# If curious plot clusters as 3D
# Seurat_mxnorm <- RunUMAP(Seurat_mxnorm, dims = 1:6, n.components = 3)
# SCP::CellDimPlot3D(srt = Seurat_mxnorm,
#                    group.by = "RNA_snn_res.0.2"
# )


####10. MAKE NON-SPATIAL FIGURES USING CLUSTER INFORMATION####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),
                             "/Only glomeruli/Processed_Glomeruli.rds"))
#SUPPLEMENTARY FIGURE
#Plot clusters
DimPlot(object = Seurat.obj,repel = T,
        label = T,label.box = T,label.size = 5,
        cols = c("orange","red","grey","violet","yellow"),
        pt.size = 1.5,reduction = "umap") + NoLegend()

# #Check if integration worked
# DimPlot(object = Seurat.obj,repel = T,
#         label = T,label.box = T,label.size = 5,group.by = "slide",
#         pt.size = 1,reduction = "umap")
# Yes it did!!

#SUPPLEMENTARY FIGURE
#Plot Markers for each cluster
#Check distribution of key cell markers as DotPlot
DotPlot(Seurat.obj,
        c("VCAM1","Megalin",
          "CD31",
          "CD45","CD68","C5aR1","CD11b",
          "aSMA"
        ),
        scale = T,
        col.min = 0.7,
        col.max = 1.5,
        dot.min = F)+ scale_size(range = c(15, 15)) + guides(size = "none") +
  scale_color_gradient(low = "white",
                       high = "blue",
                       limits = c(0.7, 1.5),
                       oob = scales::squish) + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())


#SUPPLEMENTARY FIGURE
#Plot Immune cluster frequencies by glomeruli and disease

#Extract data from Seurat object
cluster_counts <-  data.frame("Cells" = Idents(Seurat.obj),
                              "Glomeruli" = Seurat.obj$Glomeruli,
                              "Disease" = Seurat.obj$Disease
)

#Calculate immune cell frequency per glomeruli
Immune_freq <- cluster_counts %>%
  group_by(Disease, Glomeruli, Cells) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>% #ask
  group_by(Disease, Glomeruli) %>%
  mutate(Frequency = Cell_Count / sum(Cell_Count)) %>%
  ungroup() %>% filter(Cells == "Immune")

#Set disease order
Immune_freq$Disease <- factor(Immune_freq$Disease, levels = c("LN","Diabetes","CCE"))

#Check by Kruskal Wallis if all markers' variability is explained by Disease
compare_means(
  Frequency ~ Disease,
  data = Immune_freq,
  method = "kruskal.test",
  p.adjust.method = "bonferroni")
#Answer: yes

#Plot
ggplot(Immune_freq, aes(x = Disease, y = Frequency, fill = Disease)) +
  geom_jitter(size = 0.4) + geom_violin(width = 0.7) + 
  geom_boxplot(width = 0.5,outlier.size = 0.2,alpha = 0.7) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values =c("LN" = "red","Diabetes" ="green","CCE" = "magenta")) +
  theme_classic() +
  labs(
    title = "Immune Cell Frequency per Glomeruli",
    x = "Disease Group",
    y = "Cell Frequency (%)",
    fill = "Disease"
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("LN", "Diabetes"), c("LN", "CCE")),
    p.adjust.method = "bonferroni",
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    hide.ns = TRUE
  ) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,face = "bold"),
        legend.position = "none")

#SUPPLEMENTARY FIGURE
#Plot Immune marker positivity in Immune cluster per disease
DotPlot(Seurat.obj[, Idents(Seurat.obj) == "Immune"],
        features = c("CD68.Pos", "CD163.Pos", "CD11b.Pos", "C5aR1.Pos",
                     "MPO.Pos", "CD3.Pos", "CD20.Pos"),
        group.by = "Disease", 
        scale = FALSE) + 
  scale_size(range = c(2, 12)) +  
  scale_color_gradient(low = "white", high = "blue", oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  guides(
    size = guide_legend(title = "Percentage Positive"),
    color = guide_colorbar(title = "Average Expression")
  )

#SUPPLEMENTARY FIGURE
#Plot Associations between immune cells in Lupus and C5aR1 and CD11b

#Obtain Lupus Immune cell data
Lupus.Imm <- Seurat.obj[, Seurat.obj$Disease == "LN" & 
                          Idents(Seurat.obj) == "Immune"]

#Extract data from Seurat object
Imune_counts <-  data.frame("CD20" = Lupus.Imm$CD20.Pos,
                              "CD3" = Lupus.Imm$CD3.Pos,
                              "CD68" = Lupus.Imm$CD68.Pos,
                              "C5aR1" = Lupus.Imm$C5aR1.Pos,
                              "CD11b" = Lupus.Imm$CD11b.Pos)

#Run fisher test for p.value and O.R. calculation
#Set comparisons
marker_pairs <- list(
  c("C5aR1", "CD3"),c("C5aR1", "CD20"),c("C5aR1", "CD68"),
  c("CD11b", "CD3"),c("CD11b", "CD20"),c("CD11b", "CD68")
)
#Run fisher test for each pair
fisher_results <- lapply(marker_pairs, function(pair) {
  
  #Extract counts per pair
  x <- Imune_counts[[pair[1]]]
  y <- Imune_counts[[pair[2]]]
  
  #Do counts table
  tab <- table(x, y)
  
  #Run Fisher's test on table
  test <- fisher.test(tab)
  
  #Store results
  data.frame(
    Marker1 = pair[1],
    Marker2 = pair[2],
    OR = test$estimate,
    p.value = test$p.value
  )
})

#Combine results into one data frame and adjust p.value
fisher_df <- do.call(rbind, fisher_results)
fisher_df$p.adj <- p.adjust(fisher_df$p.value, method = "bonferroni")

#Filter only relevant comparisons for the plot
plot_df <- fisher_df %>%
  filter(Marker1 %in% c("C5aR1", "CD11b") & 
           Marker2 %in% c("CD3", "CD20", "CD68"))

#Create bubble plot
ggplot(plot_df, aes(x = Marker1, y = Marker2)) +
  geom_point(aes(size = OR, fill = p.adj, alpha = log2(OR+1)), shape = 21, color = "black") +
  scale_size_continuous(range = c(3, 12)) +
  scale_fill_gradient(low = "red", high = "white", name = "Adjusted\np-value") +
  theme_minimal() +
  labs(
    x = "Complement Markers",
    y = "Immune Markers",
    size = "Odds Ratio",
    title = "Complement Markers vs Immune Markers (Lupus)"
  ) +
  theme(
    axis.text.x = element_text(angle = 0,size = 10, 
                               hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold",size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



#MAIN FIGURE
#Complement deposits in Glomeruli across disease and cell types
LN <- DotPlot(Seurat.obj[, Seurat.obj$Disease == "LN"],
              features = c("C1q.Pos", "MBL.Pos", "C4d.Pos", 
                           "C3c.Pos", "C3d.Pos", "C9.Pos"),
              scale = FALSE) + 
  labs(title = "Complement Deposits in Lupus Nephritis") + 
  scale_color_gradient(low = "white", high = "blue", 
                       limits = c(0, 1), oob = scales::squish) +
  scale_size(limits = c(15, 100),range = c(0.1,25))  + 
  theme(legend.position = "none",
        axis.title.y = element_blank(),
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
                       limits = c(0, 1), oob = scales::squish) +
  scale_size(limits = c(15, 100),range = c(0.1,25))  + 
  theme(legend.position = "none",
        axis.title.y = element_blank(),
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
                       limits = c(0, 1), oob = scales::squish) +
  scale_size(limits = c(15, 100),range = c(0,25))  + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
        )

#PLOT FIGURE 
LN / DN / CCE

# #OPTION 2 ALTERNATIVE PLOT C deposits by Glomeruli
# cluster_counts <-  data.frame("Cells" = Idents(Seurat.obj),
#                               "Glomeruli" = Seurat.obj$Glomeruli,
#                               "Disease" = Seurat.obj$Disease,
#                               "C1q" = Seurat.obj$C1q.Pos,
#                               "MBL" = Seurat.obj$MBL.Pos,
#                               "C4d" = Seurat.obj$C4d.Pos,
#                               "C3c" = Seurat.obj$C3c.Pos,
#                               "C3d" = Seurat.obj$C3d.Pos,
#                               "C9" = Seurat.obj$C9.Pos)
# 
# #Summarize by Glomeruli & Cell Type
# positive_df <- cluster_counts %>%
#   group_by(Cells, Disease, Glomeruli) %>%
#   summarise(
#     C3c = mean(C3c == 1),
#     C3d = mean(C3d == 1),
#     C4d = mean(C4d == 1),
#     C9 = mean(C9 == 1),
#     C1q = mean(C1q == 1),
#     MBL = mean(MBL == 1)
#   ) %>%
#   pivot_longer(cols = c(-Cells, -Glomeruli, -Disease),
#                names_to = "Marker", values_to = "Frequency")
# 
# #Re-order for plotting purposes
# positive_df$Cells <- factor(positive_df$Cells,levels = c("Endothelial","Mesangial","Immune","PECs/PT-Conv","Unknown"))
# positive_df$Marker <- factor(positive_df$Marker,levels = c("C1q","MBL","C4d","C3c","C3d","C9"))
# positive_df$Disease <- factor(positive_df$Disease, levels = c("LN","Diabetes","CCE"))
# 
# # Plot positive frequencies complement by Glomeruli & Cell Type
# ggplot(positive_df, aes(x = Cells, y = Frequency, fill = Marker)) +
#   geom_jitter(size = 0.01) + geom_violin(width = 0.7) + geom_boxplot(width = 0.5,outlier.size = 0.2,alpha = 0.7) + facet_grid(Disease ~ Marker,switch = "y") +
#   scale_y_continuous(labels = scales::percent_format()) +
#   scale_fill_manual(values =c("C1q" = "lightblue","MBL" ="green","C3c" = "orange", "C3d" = "red", "C4d" = "violet", "C9" = "royalblue")) +
#   theme_classic() +
#   labs(
#     title = "Complement Deposits on Glomeruli",
#     x = "Cell Type",
#     y = "Positive Cell Frequency (%)",
#     fill = "Marker"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         title = element_text(face = "bold"),
#         plot.title = element_text(hjust = 0.5,face = "bold"))


#SUPPLEMENTARY FIGURE
#Show C deposits on for Endothelial Cells, Mesangial Cells & "PECs/PT-Conv" at Glomeruli level
Endo.Mes.PECs <- Seurat.obj[ ,
                             (Idents(Seurat.obj) == "Endothelial" |
                                Idents(Seurat.obj) == "Mesangial" |
                                Idents(Seurat.obj) == "PECs/PT-Conv") ]
cluster_counts <-  data.frame("Cells" = Idents(Endo.Mes.PECs),
                              "Glomeruli" = Endo.Mes.PECs$Glomeruli,
                              "C1q" = Endo.Mes.PECs$C1q.Pos,
#                             "MBL" = Endo.Mes.PECs$MBL.Pos, #Very Low positivity
                              "C4d" = Endo.Mes.PECs$C4d.Pos,
                              "C3c" = Endo.Mes.PECs$C3c.Pos,
                              "C3d" = Endo.Mes.PECs$C3d.Pos,
                              "C9" = Endo.Mes.PECs$C9.Pos,
                              "Patient" = Endo.Mes.PECs$slide,
                              "Disease" = Endo.Mes.PECs$Disease
)

#Summarize by Glomeruli & Cell Type
positive_df <- cluster_counts %>%
  group_by(Cells, Disease, Glomeruli) %>%
  summarise(
    C3c = mean(C3c == 1),
    C3d = mean(C3d == 1),
    C4d = mean(C4d == 1),
    C9 = mean(C9 == 1),
    C1q = mean(C1q == 1),
  )  %>%
  pivot_longer(cols = c(-Cells, -Glomeruli, -Disease),
               names_to = "Marker", values_to = "Frequency")

#Order factor for plot aesthethics
positive_df$Cells <- factor(positive_df$Cells,levels = c("Endothelial","Mesangial","PECs/PT-Conv"))
positive_df$Marker <- factor(positive_df$Marker,levels = c("C1q","C4d","C3c","C3d","C9"))
positive_df$Disease <- factor(positive_df$Disease, levels = c("LN","Diabetes","CCE"))


#Check by Kruskal Wallis if all markers' variability is explained by Disease
compare_means(
  Frequency ~ Disease,
  data = positive_df,
  group.by = c("Cells", "Marker"),
  method = "kruskal.test",
  p.adjust.method = "bonferroni"
)
#Answer: yes

#Calculate multiple comparisons using wicox test with Bonferroni Correction
comparisons_df <- compare_means(
  Frequency ~ Disease,
  data = positive_df,
  group.by = c("Cells", "Marker"),
  ref.group = "LN",
  method = "wilcox.test",
  p.adjust.method = "bonferroni"
)

#Set y coordinate for p.values
comparisons_df$y.position <- 1.05
comparisons_df$y.position[seq(1, nrow(comparisons_df), by = 2)] <- 1.15

#Plot Comparisons across disease
ggplot(positive_df, aes(x = Disease, y = Frequency, fill = Marker)) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.5) +
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.size = 0.2, alpha = 0.7) +
  stat_pvalue_manual(
    comparisons_df,
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    hide.ns = TRUE) +
  facet_grid(Cells ~ Marker, switch = "y") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.2)) +
  scale_fill_manual(values = c(
    "C1q" = "lightblue", "MBL" = "green", "C3c" = "orange",
    "C3d" = "red", "C4d" = "violet", "C9" = "royalblue"
  )) +
  theme_classic() +
  labs(
    title = "Complement Deposits on Glomeruli across Diseases",
    x = "Disease Group",
    y = "Positive Cell Frequency (%)",
    fill = "Marker"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) + scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.2),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)  # Ends at 100%
  )


######GLOMERULI SUBTYPES IN LUPUS#######
Lupus <- Seurat.obj[, Seurat.obj$Disease == "LN"]
cluster_counts <-  data.frame("Cells" = Idents(Lupus),
                              "Glomeruli" = Lupus$Glomeruli,
                              "C1q" = Lupus$C1q.Pos,
                              "C3c" = Lupus$C3c.Pos,
                              "C3d" = Lupus$C3d.Pos,
                              "MBL" = Lupus$MBL.Pos,
                              "C4d" = Lupus$C4d.Pos,
                              "C9" = Lupus$C9.Pos,
                              "CD68" = Lupus$CD68.Pos,
                              "C5aR1" = Lupus$C5aR1.Pos,
                              "CD11b" = Lupus$CD11b.Pos,
                              "CD3" = Lupus$CD3.Pos,
                              "CD20" = Lupus$CD20.Pos)


freq_df <- cluster_counts %>%
  group_by(Cells, Glomeruli) %>%
  summarise(
    C1q = mean(C1q == 1),
    MBL = mean(MBL == 1),
    C3c = mean(C3c == 1),
    C3d = mean(C3d == 1),
    C4d = mean(C4d == 1),
    C9 = mean(C9 == 1),
    CD68 = mean(CD68 == 1),
    C5aR1 = mean(C5aR1 == 1),
    CD11b = mean(CD11b == 1),
    CD3 = mean(CD3 == 1),
    CD20 = mean(CD20 == 1)
  )

colnames(freq_df)
mat <- freq_df[,3:13]

library(pheatmap)
heatmap.a <-pheatmap(mat,
                     cluster_rows = T,
                     cluster_cols = F,
                     cellwidth = 80,
                     cellheight = 1,  
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "correlation",
                     clustering_method = "ward.D",
                     scale = "none",
                     cutree_rows = 4,
                     show_rownames = F,
                     fontsize_col = 20,
                     breaks = seq(0,1, by =0.01),
                     color = colorRampPalette(c("#f2f7f7","yellow","orange","red","brown"))(length(seq(0,1, by =0.01))),
                     angle_col = 315)

####CORRELATIONS ACROSS GLOMERULI IN LUPUS####
#Correlations for Lupus glomeruli
library(corrplot)
library(psych)
library(Matrix)
freq_df

#Lupus
cor_matrix <- cor(mat,use = "pairwise.complete.obs", method = "spearman")
testRes <- corr.test(mat,
                     use = "pairwise",
                     method="spearman",
                     adjust="bonferroni", 
                     alpha=.05,
                     ci=TRUE,
                     minlength=5,
                     normal=TRUE)
corrplot(cor_matrix,
         p.mat = as.matrix(forceSymmetric(testRes$p,uplo = "U")),
         insig = "blank",order = "hclust",
         hclust.method = "ward.D2",
         sig.level = 0.05,
         tl.cex = 1,
         addrect = 1,
         tl.col = "black",
         type = 'upper',tl.srt = 45,
         mar=c(0,0,0,0),is.corr = T)


#Calculate make all p-vals adjusted p-vals
testRes <- as.matrix(forceSymmetric(testRes$p,uplo = "U")) > 0.05
for(i in 1:nrow(cor_matrix)){
  for(j in 1:nrow(cor_matrix)){
    if(testRes[i,j] > 0.05){
      cor_matrix[i,j] <- 0.001
    }
  }
}
colnames(cor_matrix)
#Plot only signifcant correlations with rho > 0.3
library(qgraph)
colnames(cor_matrix) <- sub("_Pos","",colnames(cor_matrix))
rownames(cor_matrix) <- colnames(cor_matrix)
qgraph(cor_matrix,
       layout="spring",esize = 40,vsize=10,
       color = c("green","royalblue","orange","violet","red",
                 "magenta","yellow","lightyellow","lightblue","magenta","gray","pink"),
       borders = T,threshold =0.3,
       label.cex = 1.2,
       labels = colnames(cor_matrix),
       graph = "cor",
       edge.width = 0.35,
       node.resolution = 100,
       repulsion = 0.7,
       label.scale=F,
       graph = "association",
       posCol = c("brown"),negCol =c("blue"))


####(OPTIONAL) EXTRA FOR SUMMARY ABSTRACT######
#Plot curved heatmap for summary abstract
#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),
                             "/Only glomeruli/Processed_Glomeruli.rds"))

#Complement deposits and activation markers by disease
cluster_counts <-  data.frame("Cells" = Idents(Seurat.obj),
                              "C3c" = Seurat.obj$C3c.Pos,
                              "C3d" = Seurat.obj$C3d.Pos,
                              "MBL" = Seurat.obj$MBL.Pos,
                              "C1q" = Seurat.obj$C1q.Pos,
                              "C4d" = Seurat.obj$C4d.Pos,
                              "C9" = Seurat.obj$C9.Pos,
                              "Patient" = Seurat.obj$slide,
                              "Disease" = Seurat.obj$Disease
)
#Plot data across Disease
positive_df <- cluster_counts %>%
    group_by(Cells, Disease) %>%
    summarise(
      C3c = mean(C3c == 1),
      C3d = mean(C3d == 1),
      C4d = mean(C4d == 1),
      C9 = mean(C9 == 1),
      C1q = mean(C1q == 1),
      MBL = mean(MBL == 1)
    ) %>%
    pivot_longer(cols = c(-Cells, -Disease),
                 names_to = "Marker", values_to = "Frequency")

#Re-order for plotting purposes
positive_df$Cells <- factor(positive_df$Cells,levels = c("Endothelial","Mesangial","Immune","PECs/PT-Conv","Unknown"))
positive_df$Marker <- factor(positive_df$Marker,levels = c("C9","C3c",
                                                           "C3d","C4d",
                                                           "MBL","C1q"))

ggplot(positive_df, aes(x = interaction(Cells,Disease), y = Marker, fill = Frequency)) +
  geom_tile(color = "black") +  # Heatmap tiles with white grid lines
  scale_fill_gradient2(low = "white", high = "red",midpoint = 0.1,limits = c(0,0.9)) +  # Heatmap color scale
  coord_radial(start = pi/6, end = 2*pi/3, inner.radius = 0.3)+  # Circular transformation
  theme_minimal() + 
  labs(x = "Marker",
       y = "Disease & Cell Type",
       fill = "Positive Cell Frequency (%)",
       title = "Circular Heatmap of Complement Deposits") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 10, face = "bold"),
        panel.grid = element_blank())+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),  # Tilt x-axis labels
    axis.text.y = element_text(angle = 45, hjust = 1, vjust = 1, size = 15, face = "bold"),  # Tilt y-axis labels
    panel.grid = element_blank()  # Remove background grid
  )  # Remove background grid
#####






####MAKE NON-SPATIAL FIGURES USING CLUSTER INFORMATION####
####UNDER DEVELOPMENT: Spatial Dimplot test...#####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
files <- list.files("Data/")#Take the names of the files for analysis
Chol2 <- read.csv2(paste0("Data/",
                          grep("_CCEv2.csv",files, value = TRUE)[1]),sep = ",", check.names = F)
df <- data.frame(x = (as.numeric(Chol2$XMin[grep("Artery",Chol2$`Analysis Region`,value =F)]) + as.numeric(Chol2$XMax[grep("Artery",Chol2$`Analysis Region`,value =F)]))/2,#x axis
                 y = (as.numeric(Chol2$YMin[grep("Artery",Chol2$`Analysis Region`,value =F)]) + as.numeric(Chol2$YMax[grep("Artery",Chol2$`Analysis Region`,value =F)]))/2#y axis
)
unique(Seurat_pseudo2$slide)
df$clusters <- Idents(Seurat_pseudo2)[Seurat_pseudo2$slide == "CCE2"]
ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color = as.factor(clusters)),
             alpha = 0.5, size = 0.75) +
  theme_classic() + theme() + 
  labs(title= "") +
  scale_color_manual(values = c("yellow","blue","red","black"))   # Apply custom colors
