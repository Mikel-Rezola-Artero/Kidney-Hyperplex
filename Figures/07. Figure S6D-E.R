####1. LOAD REQUIRED LIBRARIES####
library(rstudioapi)
library(tidyverse)
library(Matrix)
library(Seurat)

####2. IMPORT scRNA-seq GLOMERULI DATA & PROCESS RPKM MATRIX####

#Set data location & select file with RPKMs
setwd(dirname(getActiveDocumentContext()$path))
files <- list.files(path = "Data/", full.names = TRUE)

#Load the glomeruli dataset as dataframe
Glomeruli <- read.csv(file=files,sep="\t")

#Smart-seq2 single-cell RNA-seq data obtained from microdissected Glomeruli
#obtained from GSE160048, data consists on RPKMs from quality filtered cells
#fro more details check PMID: 33837218

#Clean gene names
Glomeruli <- separate(Glomeruli, X, into = c("HGNC", "accession", "dup_tag"), sep = "\\|", fill = "right") #separate HGNC from accession and dup tag
rownames(Glomeruli) <- make.names(Glomeruli$HGNC,unique = T)#name rows by HGNC symbol
#meta <- Glomeruli[, c("HGNC", "accession", "dup_tag")] #In case they are needed later on
Glomeruli <- Glomeruli[,-c(1:3)]#remove non-cell columns

#Have a quick look of the data matrix
Glomeruli[1:4,1:5]

#The dataset has been already quality corrected:
#Low-quality cells were excluded when they failed to meet the following criteria:
# (1) ≥50,000 sequence reads; 
# (2) ≥40% of reads uniquely aligned to the genome; 
# (3) ≥40% of these reads mapping to RefSeq annotated exons; 
# (4) <10% of uniquely mapped reads from ERCC spike-ins;
# (5) ≥500 genes with RPKM ≥1;
#Additionally, doublets detected by Scrublet were further removed.


####3. CREATE SEURAT OBJECT & PROCESS DATA####

#Create the Seurat Object
Glomeruli <- CreateSeuratObject(counts = Matrix(data.matrix(Glomeruli), 
                                                sparse = T),
                                project = "Glomeruli")
#NOTE THAT DATA ARE RPKM NOT COUNTS!!!

#Transform data (data is RPKM already depth-normalized, just pseudolog data)
Glomeruli@assays$RNA$data <- log1p(Glomeruli@assays$RNA$counts)#Transform normalized RPKM data
range(Glomeruli@assays$RNA@layers$counts)#Check non-transformed
range(Glomeruli@assays$RNA@layers$data)#Check transformed data


#Identify highly variable features using "FindVariableFeatures" and 
#the "dispersion" method (does not assume raw counts)
Glomeruli <- FindVariableFeatures(Glomeruli,
                                  selection.method = "dispersion",#we use dispersion as we don't have counts
                                  nfeatures = 2000)

#Check top variable features
LabelPoints(plot = VariableFeaturePlot(Glomeruli),
            points = head(VariableFeatures(Glomeruli), 10),
            repel = TRUE, xnudge = 0, ynudge = 0)

#Calculate Cell Cycle Score (to remove any effects due to cell cycle)
Glomeruli <- CellCycleScoring(object = Glomeruli, g2m.features = cc.genes$g2m.genes,
                              s.features = cc.genes$s.genes)

#Scale Variable Features for Dim Red & regress out cell cycle & nFeatures
Glomeruli <- ScaleData(Glomeruli,
                       vars.to.regress = c("nFeature_RNA",
                                           "S.Score",
                                           "G2M.Score"))#scale on 2000 identified feature variables with removal of Cell cycle and mitochondria/RBC gene effects

#Run PCA for 50 dimensions on the identified variable features
Glomeruli <- RunPCA(Glomeruli,npcs = 50)

#Check how many PCs to use downstream
ElbowPlot(Glomeruli, ndims=50, reduction = "pca")#Seems that the first 12 PCs are enough

# #Determine the number of informative PCs by JackStraw
# Glomeruli <- JackStraw(Glomeruli,num.replicate = 100,dims = 50)
# Glomeruli <- ScoreJackStraw(Glomeruli, dims =1:50)
# JackStrawPlot(Glomeruli, dims = 1:50)

#Cluster Cells using the first 12 PCs
Glomeruli <- FindNeighbors(Glomeruli,dims = 1:12)
Glomeruli <- FindClusters(Glomeruli,resolution = 0.05)#we get 5 big clusters

#Non-linear dimensional reduction (UMAP) to plot clusters
Glomeruli <- RunUMAP(Glomeruli, dims =1:12)
DimPlot(Glomeruli, reduction = "umap",
        label = TRUE, repel = TRUE)

#See major glomeruli/kidney markers

#As FeaturePlot
FeaturePlot(Glomeruli,features = c("NPHS1",#Podocytes
                                   "CDH5",#Endothelial Cells
                                   "PTPRC",#Immune cells
                                   "PDGFRA",#Mesangial & Mesangial-like cells
                                   "PDGFRB",#Mesangial & Mesangial-like cells
                                   "ACTA2",#Mesangial & Mesangial-like cells
                                   "VCAM1",#Tubular
                                   "LRP2"#PECs & PT 
),
label = T,ncol=4,pt.size = 2,alpha = 0.5)

#As ViolinPlot
VlnPlot(Glomeruli,features = c("NPHS1",#Podocytes
                               "CDH5",#Endothelial Cells
                               "PTPRC",#Immune cells
                               "PDGFRA",#Mesangial & Mesangial-like cells
                               "PDGFRB",#Mesangial & Mesangial-like cells
                               "ACTA2",#Mesangial & Mesangial-like cells
                               "VCAM1",#Tubular
                               "LRP2"#PECs & PT 
),ncol=4)

#Calculate DEGs between clusters to confirm differential expression of marker genes
All.markers <- FindAllMarkers(Glomeruli, min.pct = 0.5,#at least half cells express gene
                              only.pos = TRUE, logfc.threshold = 0.5,
                              slot = "data",test.use = "wilcox_limma")

#Check marker enrichment
All.markers[All.markers$gene %in% c("NPHS1",#Podocytes
                                    "CDH5",#Endothelial Cells
                                    "PTPRC",#Immune cells
                                    "PDGFRA",#Mesangial & Mesangial-like cells
                                    "PDGFRB",#Mesangial & Mesangial-like cells
                                    "ACTA2",#Mesangial & Mesangial-like cells
                                    "VCAM1",#Tubular
                                    "LRP2"#PECs & PT 
                                    ),
            c(2:7)]


#Conclusion of cluster labelling
#Cluster 0 = Immune cells (PTPRC expression)
#Cluster 1 = Endothelial cells (CDH5 expression)
#Cluster 2 = Tubular (LRP2 expression)
#Cluster 3 = Mesangial & Mesangial-like (PDGFRA/B & ACTA2 expression)
#Cluster 4 = Podocytes (NPHS1 expression)

#Annotate clusters
Cell.types <- Glomeruli@meta.data$seurat_clusters
levels(Cell.types) <- c("Immune",
                        "Endothelial",
                        "Tubular",
                        "Mesangial",
                        "Podocyte")#rename levels

#Add to metadata
Glomeruli$`Cell type` <- Cell.types

#Check for comparison with original author cell clusters
#Check cell numbers in each cluster:
sum(Glomeruli$`Cell type` == "Immune")#319 Immune cells (3 B cells + 263 MNP + 58 T & NK cells = 324) 5 cells diff 
sum(Glomeruli$`Cell type` == "Endothelial")#246 endothelial cells (241 GEC) 5 cells diff 
sum(Glomeruli$`Cell type` == "Tubular")#117 Tubular Compartment cells (12 CD+cTAL + 98 PTC + 6 DLT = 116) 1 cell diff
sum(Glomeruli$`Cell type` == "Mesangial")#57 Mesangial cells (58 MLCs) 1 cell diff
sum(Glomeruli$`Cell type` == "Podocyte")#27 Podocytes (same as in original paper)
#Number of cells per cluster is almost identical to original papers (less than 2% difference per cluster)


####4. FIGURE S6D####

#Set new labels
Idents(Glomeruli) <- "Cell type"

#Plot annotated clusters
DimPlot(Glomeruli, reduction = "umap",
        label = TRUE, repel = TRUE) + NoLegend()


####5. FIGURE S6E####

#Plot Markers per cluster
FeaturePlot(Glomeruli,
            features = c( "PTPRC",#Immune cells
                          "PECAM1",#Endothelial Cells
                          "SYNPO",#Podocyte
                          "PDGFRA",#Mesangial & Mesangial-like cells
                          "PDGFRB",#Mesangial & Mesangial-like cells
                          "ACTA2",#Mesangial & Mesangial-like cells
                          "LRP2",#Tubular
                          "CUBN"#Tubular
            ),
            label = T, 
            cols = c("lightgrey","lightyellow", "yellow",
                     "gold","orange", "red","brown"),
            ncol=2,
            pt.size = 2,
            alpha = 0.5,
            min.cutoff = "q1",
            repel = TRUE)

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
#   [1] future_1.58.0      Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0          
# [5] Matrix_1.6-5       lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1     
# [9] dplyr_1.1.4        purrr_1.0.4        readr_2.1.4        tidyr_1.3.1       
# [13] tibble_3.2.1       ggplot2_3.5.2      tidyverse_2.0.0    rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3         
# [4] rlang_1.1.6            magrittr_2.0.3         RcppAnnoy_0.0.22      
# [7] spatstat.geom_3.4-1    matrixStats_1.5.0      ggridges_0.5.4        
# [10] compiler_4.3.2         png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         crayon_1.5.3           pkgconfig_2.0.3       
# [16] fastmap_1.2.0          labeling_0.4.3         promises_1.3.3        
# [19] tzdb_0.5.0             ggbeeswarm_0.7.2       jsonlite_2.0.0        
# [22] goftest_1.2-3          later_1.4.2            spatstat.utils_3.1-4  
# [25] irlba_2.3.5.1          parallel_4.3.2         cluster_2.1.6         
# [28] R6_2.6.1               ica_1.0-3              stringi_1.8.7         
# [31] RColorBrewer_1.1-3     spatstat.data_3.1-6    limma_3.58.1          
# [34] reticulate_1.42.0      parallelly_1.45.0      spatstat.univar_3.1-3 
# [37] lmtest_0.9-40          scattermore_1.2        Rcpp_1.0.14           
# [40] tensor_1.5             future.apply_1.11.0    zoo_1.8-12            
# [43] sctransform_0.4.2      httpuv_1.6.16          splines_4.3.2         
# [46] igraph_2.1.4           timechange_0.3.0       tidyselect_1.2.1      
# [49] dichromat_2.0-0.1      abind_1.4-8            spatstat.random_3.4-1 
# [52] codetools_0.2-19       miniUI_0.1.1.1         spatstat.explore_3.4-3
# [55] listenv_0.9.0          lattice_0.22-7         plyr_1.8.9            
# [58] shiny_1.10.0           withr_3.0.2            ROCR_1.0-11           
# [61] ggrastr_1.0.2          Rtsne_0.17             fastDummies_1.7.5     
# [64] survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.1-11   
# [67] pillar_1.10.2          KernSmooth_2.23-26     plotly_4.10.4         
# [70] generics_0.1.3         RcppHNSW_0.6.0         hms_1.1.3             
# [73] scales_1.4.0           globals_0.18.0         xtable_1.8-4          
# [76] glue_1.8.0             lazyeval_0.2.2         tools_4.3.2           
# [79] data.table_1.17.4      RSpectra_0.16-1        RANN_2.6.1            
# [82] dotCall64_1.2          cowplot_1.1.3          grid_4.3.2            
# [85] nlme_3.1-168           patchwork_1.3.0        presto_1.0.0          
# [88] beeswarm_0.4.0         vipor_0.4.7            cli_3.6.5             
# [91] spatstat.sparse_3.1-0  spam_2.11-1            viridisLite_0.4.2     
# [94] uwot_0.2.3             gtable_0.3.6           digest_0.6.37         
# [97] progressr_0.15.1       ggrepel_0.9.6          htmlwidgets_1.6.4     
# [100] farver_2.1.2           htmltools_0.5.8.1      lifecycle_1.0.4       
# [103] httr_1.4.7             statmod_1.5.0          mime_0.13             
# [106] MASS_7.3-60 
