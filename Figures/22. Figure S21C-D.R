####1. LOAD LIBRARIES & SET GLOBAL OPTIONS####
library(GEOquery)
library(rstudioapi)
library(oligo)
library(clariomshumantranscriptcluster.db)
library(limma)
library(decoupleR)
library(progeny)
library(dplyr)
library(tidyr)
library(ggplot2)
options(timeout = 600)

####2. DOWNLOAD MICROARRAY DATA FROM GEO####

#Set directory for download
setwd(dirname(getActiveDocumentContext()$path))
dir <- dirname(getActiveDocumentContext()$path)

#Get data
getGEOSuppFiles("GSE183484")

#De-compress data
tar_path <- file.path(dir, "GSE183484", "GSE183484_RAW.tar")
untar(tar_path, exdir = file.path(dir, "GSE183484_RAW"))

####3. READ MICROARRAY DATA & EXTRACT GENE EXPRESSION####

#Read the CEL files
raw <- read.celfiles(list.celfiles("GSE183484_RAW", full.names=TRUE,listGzipped = T))

#Normalize data with RMA
norm <- rma(raw)

#Extract the expression matrix
expr <- exprs(norm)

#Check processed data
dim(expr)
head(expr[, 1:3])

####4. TRANSFORM MICROARRAY PROBE IDs INTO HGNC SYMBOLS####

#Map probe IDs to gene symbols
annotations <- select(clariomshumantranscriptcluster.db, 
                      keys = rownames(expr), 
                      columns = c("SYMBOL", "GENENAME"), 
                      keytype = "PROBEID")

#Add probe IDs to the expression matrix
expr_df <- as.data.frame(expr)
expr_df$PROBEID <- rownames(expr_df)

#Merge annotations with expression
final_annotated_data <- merge(annotations, expr_df, by = "PROBEID")

#Remove probes without a gene symbol
clean_annot <- final_annotated_data[!is.na(final_annotated_data$SYMBOL), ]

#Remove unnecessary columns
sample_cols <- grep("^GSM", colnames(clean_annot), value = TRUE)
expr_numeric <- as.matrix(clean_annot[, sample_cols])
class(expr_numeric) <- "numeric"

#Collapse multiple probes per gene & check expression matrix
expr_gene <- avereps(expr_numeric, ID = clean_annot$SYMBOL)
head(expr_gene)[, 1:2]
dim(clean_annot)
dim(expr_gene)

####5. FILTER ONLY GROUPS OF INTEREST & PERFORM CONTRAST####

#Samples "GSM5558262","GSM5558264" are normal Ctrl Macrophages and samples "GSM5558263","GSM5558265" are normal C5a treated macrophages
colnames(expr_gene)

#Filter samples of interest & prepare metadata for contrast
expr_gene <- expr_gene[,c(1,3,2,4)]
treatment_group <- factor(c("Ctrl","Ctrl","C5a_Treated","C5a_Treated"))

#Create the design matrix, fit linear model & perform contrast
design <- model.matrix(~ 0 + treatment_group)
colnames(design) <- levels(treatment_group)
fit <- lmFit(expr_gene, design)
cont <- makeContrasts(
  C5a_vs_Control = C5a_Treated - Ctrl,
  levels = design
)
fit <- contrasts.fit(fit, cont)
fit <- eBayes(fit)

#Extract & check the results
C5a_control <- topTable(fit, coef = "C5a_vs_Control", number = Inf)
head(C5a_control)

####6. RUN DecoupleR TO INFER ACTIVATED PATHWAYS IN C5a TREATMENT####

#Extract the limma test-statistics
limma_t <- C5a_control["t"]

#Load PROGENy & rename columns
#progeny <- decoupleR::get_progeny(organism = 'human', 
#                              top = 500)
local_progeny <- as.data.frame(progeny::model_human_full)
net <- local_progeny %>%
  rename(source = pathway, target = gene) %>%
  select(source, target, weight)
#Keep the top 100 genes per pathway
net <- net %>%
  group_by(source) %>%
  slice_max(order_by = abs(weight), n = 100, with_ties = FALSE) %>%
  ungroup()

#Ensure gene names are characters
net$target <- as.character(net$target)

#Infer pathway activities
contrast_acts <- decoupleR::run_mlm(mat  = limma_t, 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor = 'weight', 
                                    minsize = 5)


####7. FIGURE S21C####

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

ggplot(data = contrast_acts[contrast_acts$p_value < 0.05,], 
       mapping = aes(x = stats::reorder(source, score), 
                     y = score)) + 
  geom_bar(mapping = aes(fill = score),
           color = "black",
           stat = "identity") +
  scale_fill_gradient2(low = colors[1], 
                       mid = "whitesmoke", 
                       high = colors[2], 
                       midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size = 10, 
                                   face = "bold"),
        axis.text.y = element_text(size = 10, 
                                   face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

####8. FIGURE S21D####

#Select the pathway to visualize
pathway <- 'NFkB'

#Extract only the genes belonging to the NFkB pathway & make data frame for plot
df <- net %>%
  dplyr::filter(source == pathway) %>%
  dplyr::arrange(target) %>%
  dplyr::mutate(ID = target, 
                color = "3") %>%
  tibble::column_to_rownames('target')

#Keep only genes that are present both in PROGENy and in the differential expression results
inter <- sort(dplyr::intersect(rownames(limma_t), rownames(df)))
df <- df[inter, ]

#Add the test-statistic from the limma analysis
df['t_value'] <- limma_t[inter, ]

#Assign colours according to whether the observed gene regulation
#agrees with the expected direction of regulation in PROGENy
df <- df %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value > 0, '1', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value < 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value > 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value < 0, '1', color))

#Define colours for concordant and discordant genes
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

#FIGURE S21D
ggplot(data = df, 
       mapping = aes(x = weight, 
                     y = t_value, 
                     color = color)) + 
  geom_point(size = 2.5, 
             color = "black") + 
  geom_point(size = 1.5) +
  scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
  geom_label_repel(mapping = aes(label = ID)) + 
  theme_classic() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  ggtitle(pathway)

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
#   [1] ggplot2_3.5.2                           tidyr_1.3.1                            
# [3] dplyr_1.1.4                             progeny_1.24.0                         
# [5] decoupleR_2.8.0                         limma_3.58.1                           
# [7] clariomshumantranscriptcluster.db_8.8.0 org.Hs.eg.db_3.18.0                    
# [9] AnnotationDbi_1.64.1                    oligo_1.66.0                           
# [11] Biostrings_2.70.3                       GenomeInfoDb_1.38.8                    
# [13] XVector_0.42.0                          IRanges_2.36.0                         
# [15] S4Vectors_0.40.2                        oligoClasses_1.64.0                    
# [17] rstudioapi_0.17.1                       GEOquery_2.70.0                        
# [19] Biobase_2.62.0                          BiocGenerics_0.52.0                    
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1            farver_2.1.2                blob_1.2.4                 
# [4] bitops_1.0-9                fastmap_1.2.0               RCurl_1.98-1.17            
# [7] lifecycle_1.0.4             statmod_1.5.0               KEGGREST_1.42.0            
# [10] RSQLite_2.4.1               magrittr_2.0.3              compiler_4.3.2             
# [13] rlang_1.1.6                 tools_4.3.2                 data.table_1.17.4          
# [16] S4Arrays_1.6.0              bit_4.6.0                   DelayedArray_0.32.0        
# [19] plyr_1.8.9                  xml2_1.3.8                  RColorBrewer_1.1-3         
# [22] abind_1.4-8                 withr_3.0.2                 purrr_1.0.4                
# [25] grid_4.3.2                  preprocessCore_1.64.0       affxparser_1.74.0          
# [28] scales_1.4.0                iterators_1.0.14            dichromat_2.0-0.1          
# [31] SummarizedExperiment_1.32.0 cli_3.6.5                   crayon_1.5.3               
# [34] generics_0.1.3              reshape2_1.4.4              httr_1.4.7                 
# [37] tzdb_0.5.0                  DBI_1.2.3                   cachem_1.1.0               
# [40] stringr_1.5.1               zlibbioc_1.48.2             splines_4.3.2              
# [43] BiocManager_1.30.26         matrixStats_1.5.0           vctrs_0.6.5                
# [46] Matrix_1.6-5                hms_1.1.3                   bit64_4.6.0-1              
# [49] ggrepel_0.9.6               foreach_1.5.2               affyio_1.72.0              
# [52] glue_1.8.0                  codetools_0.2-19            stringi_1.8.7              
# [55] gtable_0.3.6                GenomicRanges_1.54.1        tibble_3.3.0               
# [58] pillar_1.11.0               GenomeInfoDbData_1.2.11     R6_2.6.1                   
# [61] ff_4.5.2                    lattice_0.22-7              readr_2.1.4                
# [64] png_0.1-8                   memoise_2.0.1               Rcpp_1.1.0                 
# [67] gridExtra_2.3               SparseArray_1.6.2           MatrixGenerics_1.14.0      
# [70] pkgconfig_2.0.3            

