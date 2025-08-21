####1. LOAD REQUIRED LIBRARIES####
library(HCATonsilData)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(Seurat)

####2. IMPORT TONSIL ATLAS DATA####
#We take the data from the Tonsil Atlas, discovery cohort:
#https://pmc.ncbi.nlm.nih.gov/articles/PMC10869140/
#The HCATonsilData package provides access to its 5 main types of assays: 
#RNA, ATAC, Multiome, CITE-seq and Spatial, 
#https://github.com/massonix/HCATonsilData/blob/devel/vignettes/HCATonsilData.Rmd


#We obtain the `SingleCellExperiment` object with gene expression (RNA)
#data as follows:
sce <- HCATonsilData(assayType = "RNA", cellType = "All")
table(sce$assay)#Check the methods used for generating the RNA data

#This object consists of 377,988 cells profiled with scRNA-seq (3P)
#and 84,364 cells profiled with multiome, 
#for a total of 462,352 cells ,
#(37,378 genes were quantified across all of these).

#We also have metadata as "colData":
#Here's a brief explanation of all the variables in the colData slot of the
#`SingleCellExperiment` objects:
#     barcode: the cell barcode. Combination of [GEM well](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/glossary) and cellranger 10X barcode.
#     donor_id: the donor-specific identifier, which can be used to retrieve donor-level metadata form the table above.
#     gem_id: we gave a unique hashtag to each [GEM well](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/glossary) (10X Chip Channel) in our dataset. This allows to traceback all the metadata for a given cell.       
#     library_name: each GEM well can give rise to multiple Illumina libraries. For example, one Multiome GEM well will give rise to 2 illumina libraries (ATAC and RNA).
#     assay: 3P (scRNA-seq) or multiome
#     sex
#     age
#     age_group: kid, young adult, old adult
#     hospital: hospital where the tonsils where obtained [Hospital Clinic](https://www.clinicbarcelona.org/), [CIMA](https://cima.cun.es/), or Newcastle.
#     cohort_type: discovery or validation cohort
#     cause_for_tonsillectomy: condition that led to the tonsillectomy (e.g. tonsillitis, sleep apnea, etc.)
#     is_hashed: whether we used [cell hashing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1) or not.
#     preservation: whether we processed the sample fresh or frozen.
#     nCount_RNA, nFeature_RNA: number of counts and features (genes) detected per cell.
#     pct_mt, pct_ribosomal: percentage of counts that map to mitochondrial (^MT) or ribosomal (^RPS) genes.
#     pDNN_hashing, pDNN_scrublet, pDNN_union: proportion of doublet nearest neighbors (pDNN) using different doublet annotations.
#     S.Score, G2M.Score Phase, CC.Difference: outputs of the [CellCycleScoring](https://rdrr.io/cran/Seurat/man/CellCycleScoring.html) from `Seurat`.
#     scrublet_doublet_scores: doublet score obtained with [Scrublet](https://www.cell.com/cell-systems/fulltext/S2405-4712(18)30474-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471218304745%3Fshowall%3Dtrue)
#     scrublet_predicted_doublet: doublet annotation obtained with Scirpy (doublet or singlet)
#     doublet_score_scDblFinder: doublet score obtained with [scDblFinder](https://f1000research.com/articles/10-979), which was run in the validation cohort following [the most recent best practices for single-cell analysis](https://www.nature.com/articles/s41576-023-00586-w).
#     annotation_level_1: annotation used at our level 1 clustering (9 main compartments).
#     annotation_level_1_probability: annotation confidence for the level 1 annotation (relevant for validation cohort, as it implied KNN annotation)
#     annotation_figure_1: annotation used in the figure 1 of our manuscript. This annotation consisted of grouping the final subtypes into main cell types that were distinguishable in the UMAP.
#     annotation_20220215, annotation_20220619, annotation_20230508: time-stamped annotation for cell types.
#     annotation_20230508_probability: annotation confidence for the final annotation (relevant for validation cohort, as it implied KNN annotation)
#     UMAP_1_level_1, UMAP_2_level_1: UMAP1 coordinates of the figure 1B of the article.
#     annotation_20220215: see above.
#     UMAP_1_20220215, UMAP_2_20220215: UMAP coordinates used in figures of the preprint for each cell type.

#Let's do a DimPlot using ggplot2
#First let's extract the UMAP coordinates and the cell labels
umap_coords <- reducedDim(sce, "UMAP")
cell_types <- sce$annotation_level_1

#Now let's combine UMAP coordinates and the cell type metadata into a data frame
plot_data <- data.frame(UMAP1 = umap_coords[, 1], 
                        UMAP2 = umap_coords[, 2], 
                        CellType = as.character(cell_types))

####3. FIGURE S3H ####

#Modify the labels to fit the manuscript
plot_data$CellType[plot_data$CellType == "PDC"] <-  "pDCs"
plot_data$CellType[plot_data$CellType == "GCBC"] <-  "GC B cells"
plot_data$CellType[plot_data$CellType == "NBC_MBC"] <-  "Naïve & Memory B cells"
plot_data$CellType[plot_data$CellType == "preBC"] <-  "Pre-B cells"
plot_data$CellType[plot_data$CellType == "PC"] <-  "Plasma cells"
plot_data$CellType[plot_data$CellType == "myeloid"] <-  "Myeloid cells"
plot_data$CellType[plot_data$CellType == "Cytotoxic"] <-  "Cytotoxic cells"
plot_data$CellType[plot_data$CellType == "epithelial"] <-  "Epithelial cells"
plot_data$CellType[plot_data$CellType == "FDC"] <-  "Follicular DCs"

#Filter out non-labelled cells
plot_data <- plot_data[!is.na(cell_types),]

#Calculate the median centroids to get the coordinates 
#for the cell type labels on the graph
centroids_median <- plot_data %>%
  group_by(CellType) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))

# Generate the DimPlot using ggplot2 without Pre-B and Pre-T cells (FIGURE S2H)
ggplot(plot_data[plot_data$CellType != "Pre-B cells" & 
                   plot_data$CellType != "preTC",], 
       aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.4, alpha = 0.8) +
  labs(x = "UMAP_1", y = "UMAP_2") + 
  theme_classic() +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position="none") + 
  geom_text(data = centroids_median[centroids_median$CellType != "Pre-B cells" & 
                                      centroids_median$CellType != "preTC",],
            aes(label = CellType),
            col = "black",
            size = 4,
            vjust = -0.5,
            hjust =0.75)# FIGURE S3H


####4. FIGURE S3I ####

#Extract the counts from the sce object for the cells with labels 
#and the genes of interest
counts <- logcounts(sce)#extract normalized counts
counts.seurat <- counts[c("C1QA","C1QB","C1QC","C5AR1",
                          "CD3E", "CD4","CD8A","PECAM1","MKI67",
                          "MS4A1","CD68","CD163"),]#filter by gene

#Transform counts into Seurat compatible object
counts.seurat <- as(counts.seurat, "dgCMatrix")#transform to dgcMatrix
counts.seurat <- counts.seurat[,!is.na(cell_types)]#filter out non-labelled cells
counts.seurat <- counts.seurat[,plot_data$CellType != "Pre-B cells" & 
                                 plot_data$CellType != "preTC"]#filter out Pre-B and T cells

# Create a Seurat object with metadata
Tonsil <- CreateSeuratObject(counts = counts.seurat,
                             min.cells = 0,
                             min.features = 0)
Tonsil@assays$RNA$data <- counts.seurat #add the normalized counts to the data slot
Tonsil$Cells <- plot_data$CellType[plot_data$CellType != "Pre-B cells" & plot_data$CellType != "preTC"]

#Do the DotPlot (FIGURE S3I)
DotPlot(Tonsil,features = c("C1QA","C1QB","C1QC","C5AR1",
                            "CD3E", "CD4","CD8A","MKI67",
                            "MS4A1","CD68","CD163"),
        dot.scale = 10,col.max = 1,group.by = "Cells") +
  scale_size(range = c(1, 8)) + RotatedAxis() + 
  scale_color_gradientn(colors = c("white", "lightblue","blue"))# FIGURE S3I

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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Seurat_5.3.0                SeuratObject_5.1.0          sp_2.2-0                   
# [4] ggplot2_3.5.2               dplyr_1.1.4                 SingleCellExperiment_1.24.0
# [7] SummarizedExperiment_1.32.0 Biobase_2.62.0              GenomicRanges_1.54.1       
# [10] GenomeInfoDb_1.38.8         IRanges_2.36.0              MatrixGenerics_1.14.0      
# [13] matrixStats_1.5.0           S4Vectors_0.40.2            BiocGenerics_0.52.0        
# [16] rhdf5_2.46.1                HCATonsilData_1.0.0        
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22              splines_4.3.2                
# [3] later_1.4.2                   bitops_1.0-7                 
# [5] filelock_1.0.2                tibble_3.2.1                 
# [7] polyclip_1.10-7               fastDummies_1.7.5            
# [9] lifecycle_1.0.4               globals_0.18.0               
# [11] lattice_0.22-7                MASS_7.3-60                  
# [13] magrittr_2.0.3                plotly_4.10.4                
# [15] rmarkdown_2.29                yaml_2.3.10                  
# [17] httpuv_1.6.16                 sctransform_0.4.2            
# [19] spam_2.11-1                   spatstat.sparse_3.1-0        
# [21] reticulate_1.42.0             cowplot_1.1.3                
# [23] pbapply_1.7-2                 DBI_1.2.3                    
# [25] RColorBrewer_1.1-3            abind_1.4-8                  
# [27] zlibbioc_1.48.2               Rtsne_0.17                   
# [29] purrr_1.0.4                   RCurl_1.98-1.17              
# [31] rappdirs_0.3.3                GenomeInfoDbData_1.2.11      
# [33] ggrepel_0.9.6                 irlba_2.3.5.1                
# [35] listenv_0.9.0                 spatstat.utils_3.1-4         
# [37] goftest_1.2-3                 RSpectra_0.16-1              
# [39] spatstat.random_3.4-1         fitdistrplus_1.1-11          
# [41] parallelly_1.45.0             codetools_0.2-19             
# [43] DelayedArray_0.32.0           tidyselect_1.2.1             
# [45] farver_2.1.2                  phenoptr_0.3.2               
# [47] BiocFileCache_2.10.2          base64enc_0.1-3              
# [49] spatstat.explore_3.4-3        jsonlite_2.0.0               
# [51] progressr_0.15.1              ggridges_0.5.4               
# [53] survival_3.8-3                iterators_1.0.14             
# [55] foreach_1.5.2                 tools_4.3.2                  
# [57] ica_1.0-3                     Rcpp_1.0.14                  
# [59] glue_1.8.0                    gridExtra_2.3                
# [61] SparseArray_1.6.2             xfun_0.52                    
# [63] HDF5Array_1.30.1              withr_3.0.2                  
# [65] BiocManager_1.30.26           fastmap_1.2.0                
# [67] rhdf5filters_1.14.1           digest_0.6.37                
# [69] R6_2.6.1                      mime_0.13                    
# [71] scattermore_1.2               tensor_1.5                   
# [73] jpeg_0.1-10                   dichromat_2.0-0.1            
# [75] spatstat.data_3.1-6           RSQLite_2.4.1                
# [77] tidyr_1.3.1                   generics_0.1.3               
# [79] data.table_1.17.4             httr_1.4.7                   
# [81] htmlwidgets_1.6.4             S4Arrays_1.6.0               
# [83] uwot_0.2.3                    pkgconfig_2.0.3              
# [85] gtable_0.3.6                  blob_1.2.4                   
# [87] lmtest_0.9-40                 XVector_0.42.0               
# [89] htmltools_0.5.8.1             dotCall64_1.2                
# [91] scales_1.4.0                  png_0.1-8                    
# [93] SpatialExperiment_1.12.0      spatstat.univar_3.1-3        
# [95] knitr_1.50                    rstudioapi_0.17.1            
# [97] reshape2_1.4.4                rjson_0.2.21                 
# [99] nlme_3.1-168                  curl_6.3.0                   
# [101] cachem_1.1.0                  zoo_1.8-12                   
# [103] stringr_1.5.1                 BiocVersion_3.18.1           
# [105] KernSmooth_2.23-26            parallel_4.3.2               
# [107] miniUI_0.1.1.1                AnnotationDbi_1.64.1         
# [109] pillar_1.10.2                 grid_4.3.2                   
# [111] vctrs_0.6.5                   RANN_2.6.1                   
# [113] promises_1.3.3                dbplyr_2.5.0                 
# [115] xtable_1.8-4                  cluster_2.1.6                
# [117] evaluate_1.0.3                magick_2.8.7                 
# [119] cli_3.6.5                     compiler_4.3.2               
# [121] rlang_1.1.6                   crayon_1.5.3                 
# [123] future.apply_1.11.0           labeling_0.4.3               
# [125] plyr_1.8.9                    stringi_1.8.7                
# [127] deldir_1.0-9                  viridisLite_0.4.2            
# [129] Biostrings_2.70.3             lazyeval_0.2.2               
# [131] spatstat.geom_3.4-1           Matrix_1.6-5                 
# [133] ExperimentHub_2.10.0          RcppHNSW_0.6.0               
# [135] patchwork_1.3.0               bit64_4.6.0-1                
# [137] future_1.58.0                 Rhdf5lib_1.24.2              
# [139] KEGGREST_1.42.0               shiny_1.10.0                 
# [141] interactiveDisplayBase_1.40.0 AnnotationHub_3.10.1         
# [143] ROCR_1.0-11                   igraph_2.1.4                 
# [145] memoise_2.0.1                 bit_4.6.0   
