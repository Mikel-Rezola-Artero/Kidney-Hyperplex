####1. LOAD LIBRARIES####
library(rstudioapi)
library(biomaRt)
library(DESeq2)
library(tximport)
library(dplyr)
library(purrr)
library(magrittr)
library(mixOmics)
library(msigdbr)
library(GSVA)
library(pheatmap)

####2. PREPARE COUNTS MATRIX FROM RAW DATA####

#Set working directory to script location
setwd(dirname(getActiveDocumentContext()$path))

#Folder containing HTSeq count files
count_dir <- "PBMC_Avacopan_Patient"

#Get all results files
files <- list.files(count_dir,full.names = T)

##Load "gene.results" RSEM files
names(files) <- c("Before_1","Before_2","Before_3",
                  "During_1","During_2","During_3",
                  "After_1","After_2","After_3")
txi.rsem <- tximport(files, 
                     type = "rsem", txIn = FALSE,
                     txOut = FALSE)
txi.rsem$length <- ifelse(txi.rsem$length==0,1,txi.rsem$length) # Setting the length to 1, to avoid error for the dds object

#Make metadata
Metadata <- data.frame(
  Treatment_TimePoint = factor(c(rep("Before",3),rep("During",3),rep("After",3)),
                               levels = c("Before", "During","After"))
)

#Make dds object
dds <- DESeqDataSetFromTximport(txi.rsem, Metadata, ~ Treatment_TimePoint)

#Extract counts
counts <- counts(dds, normalized=FALSE)
head(counts)

#Set Ensembl database that matches GTF used for read alignments
ensembl.database <- useEnsembl(biomart = "ensembl",
                               dataset = "hsapiens_gene_ensembl", 
                               version = 115, 
                               host = "https://www.ensembl.org")

#Obtain HGNC symbols for each ENSEMBL ID and fix symbols
HGNC.symbols <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id'), 
                      filters = 'ensembl_gene_id', 
                      values = row.names(counts), 
                      mart = ensembl.database,
                      useCache = FALSE)
HGNC.symbols %<>% mutate(across(everything(), na_if, ""))
rownames(HGNC.symbols) <- make.names(HGNC.symbols$ensembl_gene_id, unique=TRUE)

#Merge HGNC symbols with counts matrix
dataMerged_normal <- merge(HGNC.symbols["hgnc_symbol"],as.data.frame(counts),by="row.names",all.x=TRUE)

#Remove genes with missing HGNC
dataMerged_normal <- dataMerged_normal[!is.na(dataMerged_normal$hgnc_symbol),-1]

#Add HGNC symbols as row names and remove missing data
colnames(dataMerged_normal)[1] <- "GeneSymbol"
dataMerged_normal$GeneSymbol[which(!complete.cases(dataMerged_normal))]#check if any genes have NA values
dataMerged_normal <- dataMerged_normal[complete.cases(dataMerged_normal),]#remove genes with missing data
rownames(dataMerged_normal) <- make.names(dataMerged_normal$GeneSymbol,unique = T)#add HGNC as rownames
dataMerged_normal <-dataMerged_normal[,-1]#remove HGNC column

#Create final counts matrix for analysis
counts <- dataMerged_normal[!(apply(dataMerged_normal,1,sum) == 0),]#remove genes with 0 counts

#Save merged processed Count Matrix
write.csv(counts, file = "counts_HGNC.csv", row.names = TRUE)

####3. EXPLORE PATIENT DATA ACROSS AVACOPAN TREATMENT ON PCA PLOT####

#Load prepared counts
counts <- read.csv(file = "counts_HGNC.csv")
rownames(counts) <- counts$X #set rownames as HGNC symbols
counts <- counts[,-1] #remove "X" column with HGNC symbols

#Check counts
head(counts)

#Make metadata
Metadata <- data.frame(
  Treatment_TimePoint = factor(gsub("*_.","",colnames(counts)),
                               levels = c("Before", "During","After"))
)

#Create DESeq2 object with metadata
dds <- DESeqDataSetFromMatrix(countData = counts, colData   = Metadata,
                              design    = ~ Treatment_TimePoint)

#Calculate size factors
dds <- estimateSizeFactors(dds)

#Filter low abundance transcripts (minimum 10 counts in at least 3 samples)
idx <- rowSums(counts(dds, normalized=TRUE) >= 10) >= 3
dds <- dds[idx, ]

#Apply variance stabilizing transformation (VST)
vst_data <- vst(dds, blind = TRUE)

#Select the top 500 most variable genes across samples
ntop <- 500
rv <- rowVars(assay(vst_data))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t(assay(vst_data)[select, ])

#Run PCA 9 dimensions to decide if 2 PCs are enough
tune.pca.multi <- tune.pca(mat, ncomp = 9, scale = TRUE)
plot(tune.pca.multi)
#we see that the first 2 PCs explain ~70% of the variance
sum(tune.pca.multi$prop_expl_var[["X"]][1:2])

#Run PCA on 2 components
pca.final <- pca(mat, ncomp = 2, scale = TRUE)

#Plot PCA (FIGURE S21A)
plotIndiv(pca.final,
          group = Metadata$Treatment_TimePoint,
          ellipse = T,
          ellipse.level = 0.75,
          ind.names = FALSE,
          style = "lattice",
          legend = TRUE,
          title = 'PCA (PC1 vs PC2)')


####4. CALCULATE GSVA SCORES IN HALLMARK PATHWAYS ACROSS AVACOPAN TREATMENT####

#Get all Hallmark genesets (H) for a species human
all <- msigdbr(species = "Homo sapiens", collection = "H")
x <- all[, c("gs_name", "gene_symbol")]
colnames(x) <- c("GO", "GENE")
h_list <- split(
  x$GENE,  # genes
  x$GO     # GO terms
)

#Perform GSVA using Hallmark Genesets
gsva_results.h <- gsva(
  assay(vst_data),
  h_list,
  method = "gsva",
  kcdf = "Gaussian",
  mx.diff = TRUE,
  verbose = FALSE)

#Order samples by annotation
sample_order <- order(Metadata$Treatment_TimePoint)

#Select the 15 top most variable pathways according to GSVA scores
pathway_var  <- apply(gsva_results.h, 1, sd, na.rm = TRUE)
top_pathways <- names(sort(pathway_var, decreasing = TRUE))[1:15]
gsva_top <- gsva_results.h[top_pathways, sample_order, drop = FALSE]
gsva_top #Check

#Check significantly variable pathways across treatment timepoints
p.vals <- vector()
for (i in 1:nrow(gsva_results.h)) {
  genesets <- gsva_results.h[i,]
  res <- aov(genesets ~ Metadata$Treatment_TimePoint, data = data.frame(gsva_results.h))
  res <- summary(res)
  p.vals[i] <- res[[1]][5]$`Pr(>F)`[1]
  names(p.vals)[i] <- rownames(gsva_results.h)[i]
}
p.vals <- p.adjust(p.vals,method = "fdr")

#Filter by significance
signif_genesets <- names(p.vals)[p.vals < 0.05]
gsva_top <- gsva_results.h[signif_genesets, sample_order, drop = FALSE]

#Create heatmap annotations
anno_col <- data.frame(
  Treatment_TimePoint = Metadata$Treatment_TimePoint[sample_order],
  row.names = colnames(gsva_top) # Forces an absolute, perfect match
)

#Make colors for annotations
make_palette <- function(x) {
  x <- as.factor(x)
  lev <- levels(x)
  setNames(grDevices::hcl.colors(length(lev), palette = "Dark 3"), lev)
}
anno_colors <- list(
  Treatment_TimePoint = make_palette(anno_col$Treatment_TimePoint)
)

#Define heatmap palette
heat_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)

#Plot FIGURE S21B
pheatmap(
  gsva_top,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  annotation_names_col = TRUE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  cellwidth = 25,
  cellheight = 15,
  breaks = seq(-2.5,2.5, by =0.01),
  color = colorRampPalette(c(rep("blue",1),
                             rep("white",1),
                             rep("red",1)))(length(seq(-2.5,2.5, by =0.01))),
  border_color = NA,
  angle_col = "315"
)


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
#   [1] pheatmap_1.0.12             GSVA_1.50.5                 msigdbr_24.1.0             
# [4] mixOmics_6.26.0             ggplot2_3.5.2               lattice_0.22-7             
# [7] MASS_7.3-60                 magrittr_2.0.3              purrr_1.0.4                
# [10] dplyr_1.1.4                 tximport_1.30.0             DESeq2_1.42.1              
# [13] SummarizedExperiment_1.32.0 Biobase_2.62.0              MatrixGenerics_1.14.0      
# [16] matrixStats_1.5.0           GenomicRanges_1.54.1        GenomeInfoDb_1.38.8        
# [19] IRanges_2.36.0              S4Vectors_0.40.2            BiocGenerics_0.52.0        
# [22] biomaRt_2.58.2              rstudioapi_0.17.1          
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.3                   bitops_1.0-9                gridExtra_2.3              
# [4] GSEABase_1.64.0             rlang_1.1.6                 compiler_4.3.2             
# [7] RSQLite_2.4.1               DelayedMatrixStats_1.24.0   png_0.1-8                  
# [10] vctrs_0.6.5                 reshape2_1.4.4              stringr_1.5.1              
# [13] pkgconfig_2.0.3             crayon_1.5.3                fastmap_1.2.0              
# [16] dbplyr_2.5.0                XVector_0.42.0              graph_1.80.0               
# [19] bit_4.6.0                   beachmat_2.18.1             zlibbioc_1.48.2            
# [22] cachem_1.1.0                progress_1.2.2              blob_1.2.4                 
# [25] rhdf5filters_1.14.1         DelayedArray_0.32.0         Rhdf5lib_1.24.2            
# [28] BiocParallel_1.36.0         irlba_2.3.5.1               parallel_4.3.2             
# [31] prettyunits_1.2.0           R6_2.6.1                    stringi_1.8.7              
# [34] RColorBrewer_1.1-3          Rcpp_1.1.0                  assertthat_0.2.1           
# [37] Matrix_1.6-5                igraph_2.1.4                tidyselect_1.2.1           
# [40] dichromat_2.0-0.1           abind_1.4-8                 codetools_0.2-19           
# [43] curl_6.3.0                  tibble_3.3.0                plyr_1.8.9                 
# [46] withr_3.0.2                 rARPACK_0.11-0              KEGGREST_1.42.0            
# [49] BiocFileCache_2.10.2        xml2_1.3.8                  Biostrings_2.70.3          
# [52] pillar_1.11.0               filelock_1.0.2              ellipse_0.5.0              
# [55] generics_0.1.3              RCurl_1.98-1.17             hms_1.1.3                  
# [58] sparseMatrixStats_1.14.0    scales_1.4.0                xtable_1.8-4               
# [61] glue_1.8.0                  tools_4.3.2                 ScaledMatrix_1.10.0        
# [64] RSpectra_0.16-1             annotate_1.80.0             locfit_1.5-9.12            
# [67] babelgene_22.9              XML_3.99-0.18               rhdf5_2.46.1               
# [70] grid_4.3.2                  tidyr_1.3.1                 AnnotationDbi_1.64.1       
# [73] colorspace_2.1-1            SingleCellExperiment_1.24.0 GenomeInfoDbData_1.2.11    
# [76] BiocSingular_1.18.0         HDF5Array_1.30.1            rsvd_1.0.5                 
# [79] cli_3.6.5                   rappdirs_0.3.3              S4Arrays_1.6.0             
# [82] corpcor_1.6.10              gtable_0.3.6                digest_0.6.37              
# [85] SparseArray_1.6.2           ggrepel_0.9.6               farver_2.1.2               
# [88] memoise_2.0.1               lifecycle_1.0.4             httr_1.4.7                 
# [91] bit64_4.6.0-1      
