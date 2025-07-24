####1. LOAD LIBRARIES####
library(rstudioapi)
library(ggplot2)
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
files <- list.files(path = "Data/", full.names = TRUE)

#Create list and import cell-segmented data from HALO Analysis
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

#Filter only Glomeruli cells based on manual annotation by nephrologist
Glom_IF <- lapply(Kidney_IF, function(df) {
  df[grep("Glomeruli", df$`Analysis Region`, value = FALSE), ]
})


####3. Q.C.: CHECK OF AREAS SELECTED FROM EACH IMAGE####

#Plot Glomeruli cells for quick annotation quality check
plot_data <- do.call(rbind, lapply(names(Glom_IF), function(name) {
  df <- Glom_IF[[name]][,c(5,7)]
  cbind(df, Dataset = name)}))

ggplot(plot_data, aes(x = as.numeric(`XMin`), y = as.numeric(`YMin`))) +
  geom_point(aes(color = Dataset), size = 0.01, alpha = 0.1) +
  facet_wrap(~Dataset, ncol = 3) +
  labs(x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")
#We clearly see glomerular structures distributed spatially in our filtered data

#From previous EDA (not shown) we detected one strange glomeruli with 
#abnormaly high infiltration of T cells
df <- Glom_IF[["Lupus8"]][,c(2,5,7,84)]

df$`CD3-Eq13 - Cy5 Positive Classification` <- factor(df$`CD3-Eq13 - Cy5 Positive Classification`)
ggplot(df, aes(x = as.numeric(`XMin`), y = as.numeric(`YMin`))) +
  geom_point(aes(color = `CD3-Eq13 - Cy5 Positive Classification`, 
                 shape = `Analysis Region`), size = 0.5, alpha = 0.4) +
  labs(x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal()
#Note "green" glomerulus indicating high levels of CD3 Pos cells in Glom 2

#We will filter out of the analysis one glomerulus 2 for patient Lupus8 
#as high frequencies of CD3+ T cells suggests TLS/Tubulo-Int compartment
#annotation "contamination" or sclerotic glomeruli

#Filter out Glomeruli 2 from Lupus8
Glom_IF[["Lupus8"]] <- Glom_IF[["Lupus8"]] %>%
  filter(!`Analysis Region` == "Glomeruli 2")


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
#Complement: "C1q" "MBL" "C3c" "C3d" "C5aR1" "C1r" "C1s" "C4d" "C9"
#Immune Myeloid: "CD45" "CD68" "CD163" "CD11b" "MPO"
#Immune Lympho:  "CD45" "CD20" "CD3" 
#Stroma: "CD31" "aSMA" "Megalin"
#Cell State: "Ki67" "VCAM-1"


####5. PROCESS FLUORESCENCE INTENSITIES FOR MXNORM OBJECT CREATION####

for (i in 1:15) {
  int_name <- names(Glom_intensity)[i]
  orig_name <- names(Glom_IF)[i]
  
  #Convert all columns to numeric
  Glom_intensity[[int_name]][] <- lapply(Glom_intensity[[int_name]], as.numeric)
  
  #Remove DAPI from Intensity Data
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
#Immune Lympho:  "CD45" "CD20" "CD3" 
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
    
    #Extract disease name by removing trailing numbers
    disease <- gsub("\\d+$", "", df_name)
    
    #Repeat for all rows in current dataframe
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
#The test stat has diminished by >1/2

#Run an Otsu discordance score analysis to determine how well our normalization method performs
mx_norm = run_otsu_discordance(mx_norm,
                               table="both",
                               threshold_override = NULL,
                               plot_out = FALSE)
summary(mx_norm)#Check summary data including Otsu discordance

#Visualize the results of the Otsu miss-classification analysis: raw vs corrected
ggplotly(plot_mx_discordance(mx_norm))#some markers are not completely corrected

#Check % of positive cells for every marker according to supervised labeling
(apply(mx_norm$norm_data[grep("Pos",colnames(mx_sample), value = TRUE)],2,sum)/
    nrow(mx_norm$norm_data[grep("Pos",colnames(mx_sample), value = TRUE)]))*100 
#The markers with low of batch correction correspond with very low positivity frequency (< 1%)


#Visualize the densities of the marker values normalized vs raw
density_plot <- plot_mx_density(mx_norm)
# ggsave(filename = file.path("Only glomeruli/", "density_plot_wide.png"),
#        plot = density_plot,
#        width = 35,  # Increase width to fit facets
#        height = 6,
#        dpi = 300)#Check png file for plot output


#Do UMAP & check how the batch correction integration did
mx_norm = run_reduce_umap(mx_norm,
                          table="both",
                          marker_list = colnames(mx_norm$norm_data)[c(15,17,18,19)],
                          downsample_pct = 0.2,
                          metadata_col = colnames(mx_norm$norm_data)[c(24:44)])

#Q.C., Plot UMAP with raw and batch corrected data
plot_mx_umap(mx_norm,metadata_col = "slide")

#Check if main cell markers present logical distributions in batch corrected data
plot_mx_umap(mx_norm,metadata_col = "CD31.Pos")
plot_mx_umap(mx_norm,metadata_col = "CD45.Pos")
plot_mx_umap(mx_norm,metadata_col = "VCAM1.Pos")
plot_mx_umap(mx_norm,metadata_col = "aSMA.Pos")

#We see good segregation of key markers, we can use batch corrected data for clustering


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
                                     "VCAM1",#PECs & PT_VCAM1
                                     "Nuc.size"#All
                        ),
                        npcs = 5,
                        approx = FALSE)

#Run UMAP
Seurat_mxnorm <- RunUMAP(Seurat_mxnorm, dims = 1:5, 
                         n.components = 2)

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
          "VCAM1",#PECs & PT_VCAM1
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
                         "VCAM1",#PECs & PT_VCAM1
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

#Cluster 0 = Unknown
#Cluster 1 = Endothelial
#Cluster 2 = Mesangial
#Cluster 3 = PECs/PT_VCAM1
#Cluster 4 = Unknown
#Cluster 5 = Immune
#Cluster 6 = Unknown

#Rename clusters according to Cell Markers
new_cluster_names <- c("Unknown",
                       "Endothelial",
                       "Mesangial",
                       "PECs/PT_VCAM1",
                       "Unknown",
                       "Immune",
                       "Unknown")

names(new_cluster_names) <- levels(Seurat_mxnorm)
Seurat_mxnorm <- RenameIdents(Seurat_mxnorm, 
                              new_cluster_names)

#Store Processed Seurat Object with labels
saveRDS(Seurat_mxnorm, file = paste0(dirname(getActiveDocumentContext()$path),
                                     "/Only glomeruli/Processed_Glomeruli.rds"))

# If curious plot clusters as 3D
# Seurat_mxnorm$clusters <- Idents(Seurat_mxnorm)
# Seurat_mxnorm <- RunUMAP(Seurat_mxnorm, dims = 1:6, n.components = 2)
# library(SCP)
# CellDimPlot3D(srt = Seurat_mxnorm,
#               group.by = "clusters")

####10. Session Information#####

sessionInfo()

# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /home/mikel/miniconda3/envs/image-analysis-env/lib/libmkl_rt.so.2;  LAPACK version 3.10.1
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
#   [1] future_1.58.0      Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0           plotly_4.10.4     
# [6] mxnorm_1.0.3       reticulate_1.42.0  dplyr_1.1.4        ggplot2_3.5.2      rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.3         spatstat.utils_3.1-4  
# [5] SuppDists_1.1-9.9      farver_2.1.2           nloptr_2.2.1           vctrs_0.6.5           
# [9] ROCR_1.0-11            spatstat.explore_3.4-3 minqa_1.2.8            htmltools_0.5.8.1     
# [13] pROC_1.18.5            caret_7.0-1            sctransform_0.4.2      parallelly_1.45.0     
# [17] KernSmooth_2.23-26     htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9            
# [21] zoo_1.8-12             lubridate_1.9.4        igraph_2.1.4           mime_0.13             
# [25] lifecycle_1.0.4        iterators_1.0.14       pkgconfig_2.0.3        Matrix_1.6-5          
# [29] R6_2.6.1               fastmap_1.2.0          rbibutils_2.3          fitdistrplus_1.1-11   
# [33] shiny_1.10.0           digest_0.6.37          patchwork_1.3.0        tensor_1.5            
# [37] RSpectra_0.16-1        irlba_2.3.5.1          crosstalk_1.2.1        labeling_0.4.3        
# [41] progressr_0.15.1       spatstat.sparse_3.1-0  timechange_0.3.0       polyclip_1.10-7       
# [45] abind_1.4-8            httr_1.4.7             compiler_4.3.2         proxy_0.4-27          
# [49] withr_3.0.2            fastDummies_1.7.5      MASS_7.3-60            lava_1.8.1            
# [53] rappdirs_0.3.3         ModelMetrics_1.2.2.2   tools_4.3.2            lmtest_0.9-40         
# [57] httpuv_1.6.16          future.apply_1.11.0    goftest_1.2-3          nnet_7.3-19           
# [61] glue_1.8.0             nlme_3.1-168           promises_1.3.3         grid_4.3.2            
# [65] Rtsne_0.17             cluster_2.1.6          reshape2_1.4.4         generics_0.1.3        
# [69] recipes_1.3.1          spatstat.data_3.1-6    gtable_0.3.6           class_7.3-22          
# [73] tidyr_1.3.1            data.table_1.17.4      utf8_1.2.6             spatstat.geom_3.4-1   
# [77] RcppAnnoy_0.0.22       ggrepel_0.9.6          RANN_2.6.1             foreach_1.5.2         
# [81] pillar_1.10.2          stringr_1.5.1          spam_2.11-1            RcppHNSW_0.6.0        
# [85] later_1.4.2            splines_4.3.2          lattice_0.22-7         deldir_1.0-9          
# [89] survival_3.8-3         tidyselect_1.2.1       miniUI_0.1.1.1         pbapply_1.7-2         
# [93] gridExtra_2.3          reformulas_0.4.1       scattermore_1.2        stats4_4.3.2          
# [97] hardhat_1.4.1          timeDate_4041.110      matrixStats_1.5.0      stringi_1.8.7         
# [101] lazyeval_0.2.2         yaml_2.3.10            boot_1.3-31            codetools_0.2-19      
# [105] kSamples_1.2-10        tibble_3.2.1           cli_3.6.5              uwot_0.2.3            
# [109] rpart_4.1.24           xtable_1.8-4           Rdpack_2.6.4           dichromat_2.0-0.1     
# [113] Rcpp_1.0.14            spatstat.random_3.4-1  globals_0.18.0         png_0.1-8             
# [117] spatstat.univar_3.1-3  parallel_4.3.2         gower_1.0.2            dotCall64_1.2         
# [121] lme4_1.1-37            listenv_0.9.0          viridisLite_0.4.2      ipred_0.9-15          
# [125] scales_1.4.0           prodlim_2025.04.28     e1071_1.7-13           ggridges_0.5.4        
# [129] crayon_1.5.3           purrr_1.0.4            rlang_1.1.6            cowplot_1.1.3    