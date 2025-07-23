####1. LOAD LIBRARIES####
library(rstudioapi)
library(dbscan)
library(dplyr)
library(qgraph)

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


####3. FIGURE S4B####

#Extract Tonsil 3 data
Tonsil3 <- tonsil_list[["Tonsil3"]]
df <- data.frame(x = (as.numeric(Tonsil3$XMin) + as.numeric(Tonsil3$XMax))/2,#x axis
                 y = (as.numeric(Tonsil3$YMin) + as.numeric(Tonsil3$YMax))/2,#y axis
                 Phenotype = Idents(Seurat.obj)[Seurat.obj$slide == "Tonsil3"],#Idents
                 IDs = colnames(Seurat.obj)[Seurat.obj$slide == "Tonsil3"]#Cell IDs
)

#Extract only replicating B cells
df.clust <- df[df$Phenotype == "B cell Ki67 High 1" |
                 df$Phenotype == "B cell Ki67 High 2",
               1:2]
rownames(df.clust) <- df$IDs[df$Phenotype == "B cell Ki67 High 1" |
                               df$Phenotype == "B cell Ki67 High 2"]

#Set the graphics to display two plots
par(mfrow = c(1,2))

#How many follicles?
plot(x = as.numeric(df.clust$x),#x axis
     y = 1/as.numeric(df.clust$y),#y axis
     xlab = "x", ylab= "y",
     pch = 16,cex = .5,
     main = "B cells Ki67 Positive")
#We see around 8 follicles

#Run hdbscan to automatically detect follicles
#We inspect the 100-NN distance plot to decide minPts.
#We want minimum 100 B cells to form a cluster/follicle
#We run clustering with hdbscan with minPts = 100 B cells:
db <- hdbscan(df.clust,minPts = 200)
num_clusters <- max(db$cluster)#unclassified cells are defined as "0" 
cluster_colors <- c("white", rainbow(num_clusters))
plot(x = as.numeric(df.clust$x),#x axis
     y = 1/as.numeric(df.clust$y),#y axis
     pch = as.numeric(db$cluster),cex = .5,
     col = cluster_colors[db$cluster + 1], xlab = "x",ylab = "y",
     main = "Follicle Annotation with Hierarchical DBSCAN")
#The algorithm detects 7 follicles

#Set the graphics back to normal
par(mfrow = c(1,1))


####4. FIGURE S4C####

#Correlation between complement deposits across follicles

#Process data
corr.data <- data.frame(Tonsil3[df$Phenotype == "B cell Ki67 High 1" |
                                  df$Phenotype == "B cell Ki67 High 2",
                                c(23,29,41,51,69,73)]) #extract deposits info on follicular B cells
corr.data$cluster <- db$cluster
corr.data <- corr.data %>%
  group_by(cluster) %>%
  summarise(across(everything(),
                   mean, na.rm = TRUE)
            )#calculate frequency of positivity by follicle

#Correct column names
colnames(corr.data) <- gsub("[_.].*Positive.Classification", "", colnames(corr.data))
colnames(corr.data) <- gsub("Abcam", "", colnames(corr.data))

#Remove unclassified cells
corr.data <- corr.data[!corr.data$cluster == 0,-1]

#Calculate correlations
M1 <- cor(corr.data,
          use = "pairwise.complete.obs",
          method = "spearman")

#Plot as correlation network (FIGURE S4C)
colnames(M1) <- sub("_Pos","",colnames(M1))
rownames(M1) <- colnames(M1)
qgraph(M1,
       layout="spring",
       esize = 40,
       vsize=15,
       color = c("purple","green","brown","red",
                 "lightblue","skyblue"),
       borders = T,
       threshold = 0.45,#show only correlations above 0.4
       label.cex = 1.2,
       labels = colnames(M1),
       graph = "cor",
       edge.width = 0.5,
       node.resolution = 100,
       repulsion = 0.8,
       label.scale=F,
       graph = "association",
       posCol = c("brown"),negCol =c("blue"))


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
#   [1] qgraph_1.9.8       dplyr_1.1.4        dbscan_1.2.2       SeuratObject_5.1.0 sp_2.2-0           rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] dotCall64_1.2       spam_2.11-1         gtable_0.3.6        xfun_0.52           ggplot2_3.5.2      
# [6] htmlwidgets_1.6.4   psych_2.5.3         lattice_0.22-7      quadprog_1.5-8      vctrs_0.6.5        
# [11] tools_4.3.2         generics_0.1.3      stats4_4.3.2        parallel_4.3.2      tibble_3.2.1       
# [16] cluster_2.1.6       phenoptr_0.3.2      pkgconfig_2.0.3     Matrix_1.6-5        data.table_1.17.4  
# [21] checkmate_2.3.2     RColorBrewer_1.1-3  lifecycle_1.0.4     compiler_4.3.2      farver_2.1.2       
# [26] stringr_1.5.1       mnormt_2.1.1        codetools_0.2-19    htmltools_0.5.8.1   glasso_1.11        
# [31] fdrtool_1.2.18      htmlTable_2.4.3     Formula_1.2-5       crayon_1.5.3        pillar_1.10.2      
# [36] Hmisc_5.2-3         iterators_1.0.14    rpart_4.1.24        abind_1.4-8         foreach_1.5.2      
# [41] parallelly_1.45.0   nlme_3.1-168        lavaan_0.6-19       gtools_3.9.5        tidyselect_1.2.1   
# [46] digest_0.6.37       future_1.58.0       stringi_1.8.7       listenv_0.9.0       reshape2_1.4.4     
# [51] fastmap_1.2.0       grid_4.3.2          colorspace_2.1-1    cli_3.6.5           magrittr_2.0.3     
# [56] base64enc_0.1-3     dichromat_2.0-0.1   future.apply_1.11.0 pbivnorm_0.6.0      withr_3.0.2        
# [61] foreign_0.8-90      corpcor_1.6.10      scales_1.4.0        backports_1.5.0     rmarkdown_2.29     
# [66] globals_0.18.0      jpeg_0.1-10         igraph_2.1.4        nnet_7.3-19         gridExtra_2.3      
# [71] progressr_0.15.1    png_0.1-8           pbapply_1.7-2       evaluate_1.0.3      knitr_1.50         
# [76] rlang_1.1.6         Rcpp_1.0.14         glue_1.8.0          R6_2.6.1            plyr_1.8.9  
