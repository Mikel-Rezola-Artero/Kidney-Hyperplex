####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(scales)
library(psych)
library(Matrix)
library(qgraph)

####2. LOAD ARTERY DATA####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(getActiveDocumentContext()$path),
                             "/Only arteries/Processed_Artery.rds"))


####3. FIGURE 3G####

#VCAM-1 in Arteries across disease and cell types (FIGURE 3G)
DotPlot(Seurat.obj,split.by = "Disease",
        features = c("VCAM.1.Pos"),cols = c("lightgray","blue","lightblue"),
        scale = FALSE) + 
  labs(title = "VCAM-1 positivity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )


####4. FIGURE 3H####

#Extract data from Seurat object
cluster_counts <-  data.frame("Cells" = Idents(Seurat.obj),
                              "VCAM-1" = Seurat.obj$VCAM.1.Pos,
                              "Artery" = Seurat.obj$Artery,
                              "Disease" = Seurat.obj$Disease)

#Determine which cells are endothelial and VCAM-1 positive
cluster_counts$VCAM.1 <- ifelse(cluster_counts$Cells == "Endothelial" & 
                                  cluster_counts$VCAM.1 == "1",1,0)


#Obtain data VCAM-1 positivity on endothelial cells across Artery and Disease
positive_df <- cluster_counts %>%
  group_by(Disease, Artery) %>%
  summarise(
    "VCAM-1" = mean(VCAM.1),
  ) %>%
  pivot_longer(cols = c(-Artery, -Disease),
               names_to = "Marker", values_to = "Frequency")

#Order factor for plot aesthetics
positive_df$Disease <- factor(positive_df$Disease, levels = c("LN", "Diabetes", "CCE"))

#Check by Kruskal Wallis if marker variability is explained by Disease
compare_means(
  Frequency ~ Disease,
  data = positive_df,
  group.by = c("Marker"),
  method = "kruskal.test",
  p.adjust.method = "bonferroni"
)
#     Marker   .y.           p    p.adj p.format p.signif method        
#      <chr>  <chr>        <dbl>   <dbl>  <chr>    <chr>    <chr>         
#   1 VCAM-1 Frequency 7.58e-12 7.6e-12  7.6e-12   **** Kruskal-Wallis

#Plot data on VCAM-1 positivity on endothelial cells by Artery and Disease
ggplot(positive_df, aes(x = Disease, y = Frequency, fill = Disease)) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.5) +
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.size = 0.2, alpha = 0.7) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("CCE", "LN"), c("CCE", "Diabetes")),
    label = "p.signif",
    tip.length = 0.01
  ) +
  facet_grid( ~ Marker, switch = "y") +
  scale_y_continuous(labels = percent_format(),
                     limits = c(0, 0.15)) +
  scale_fill_manual(values = c("LN" = "red","Diabetes" ="green","CCE" = "magenta")) +
  theme_classic() +
  labs(
    title = "VCAM-1 Endothelial Cells in Arteries across Diseases",
    x = "Disease Group",
    y = "VCAM-1 Pos Endo by Artery (%)",
    fill = "Marker"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


####5. FIGURE 3I####

#Correlations for CCE arteries between C' deposits and VCAM-1

#Filter CCE data and supervised labels
CCE <- Seurat.obj[, Seurat.obj$Disease == "CCE"]
cluster_counts <-  data.frame("Cells" = Idents(CCE),
                              "Artery" = CCE$Artery,
                              "C1q" = CCE$C1q.Pos,
                              "C3c" = CCE$C3c.Pos,
                              "C3d" = CCE$C3d.Pos,
                              "MBL" = CCE$MBL.Pos,
                              "C4d" = CCE$C4d.Pos,
                              "C9" = CCE$C9.Pos,
                              "VCAM-1" = CCE$VCAM.1.Pos)

#Calculate positive frequencies in each cell type within each artery
freq_df <- cluster_counts %>%
  group_by(Cells, Artery) %>%
  summarise(
    C1q = mean(C1q),
    MBL = mean(MBL),
    C3c = mean(C3c),
    C3d = mean(C3d),
    C4d = mean(C4d),
    C9 = mean(C9),
    "VCAM-1" = mean(VCAM.1)
  )
mat <- freq_df[,3:9]

#Calculate spearman correlations for C' and VCAM-1 positivity in each artery per cell type
cor_matrix <- cor(mat,use = "pairwise.complete.obs", method = "spearman")
testRes <- corr.test(mat,
                     use = "pairwise",
                     method="spearman",
                     adjust="bonferroni", 
                     alpha=.05,
                     ci=TRUE,
                     minlength=5,
                     normal=TRUE)

#Calculate make all p-vals adjusted p-vals
testRes <- as.matrix(forceSymmetric(testRes$p,uplo = "U")) > 0.05
for(i in 1:nrow(cor_matrix)){
  for(j in 1:nrow(cor_matrix)){
    if(testRes[i,j] > 0.05){
      cor_matrix[i,j] <- 0.001
    }
  }
}

#Plot only significant correlations with rho > 0.3
colnames(cor_matrix) <- sub("_Pos","",colnames(cor_matrix))
rownames(cor_matrix) <- colnames(cor_matrix)
qgraph(cor_matrix,
       layout="spring",esize = 40,vsize=10,
       color = c("green","royalblue","orange","violet","red",
                 "magenta","yellow","lightyellow","lightblue","magenta","gray","pink"),
       borders = T,threshold =0.3,
       label.cex = 1,
       labels = colnames(cor_matrix),
       graph = "cor",
       edge.width = 0.35,
       node.resolution = 100,
       repulsion = 0.1,
       label.scale=F,
       graph = "association",
       posCol = c("brown"),negCol =c("blue"))

####6. Session Information####

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
#   [1] qgraph_1.9.8       Matrix_1.6-5       psych_2.5.3        scales_1.4.0       ggpubr_0.6.0      
# [6] ggplot2_3.5.2      tidyr_1.3.1        dplyr_1.1.4        Seurat_5.3.0       SeuratObject_5.1.0
# [11] sp_2.2-0           rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.3         spatstat.utils_3.1-4  
# [5] rmarkdown_2.29         farver_2.1.2           vctrs_0.6.5            ROCR_1.0-11           
# [9] spatstat.explore_3.4-3 base64enc_0.1-3        rstatix_0.7.2          htmltools_0.5.8.1     
# [13] broom_1.0.8            Formula_1.2-5          sctransform_0.4.2      parallelly_1.45.0     
# [17] KernSmooth_2.23-26     htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9            
# [21] plotly_4.10.4          zoo_1.8-12             igraph_2.1.4           mime_0.13             
# [25] lifecycle_1.0.4        pkgconfig_2.0.3        R6_2.6.1               fastmap_1.2.0         
# [29] fitdistrplus_1.1-11    future_1.58.0          shiny_1.10.0           fdrtool_1.2.18        
# [33] digest_0.6.37          colorspace_2.1-1       patchwork_1.3.0        tensor_1.5            
# [37] RSpectra_0.16-1        irlba_2.3.5.1          Hmisc_5.2-3            labeling_0.4.3        
# [41] progressr_0.15.1       spatstat.sparse_3.1-0  httr_1.4.7             polyclip_1.10-7       
# [45] abind_1.4-8            compiler_4.3.2         withr_3.0.2            glasso_1.11           
# [49] htmlTable_2.4.3        backports_1.5.0        carData_3.0-5          fastDummies_1.7.5     
# [53] ggsignif_0.6.4         MASS_7.3-60            corpcor_1.6.10         gtools_3.9.5          
# [57] pbivnorm_0.6.0         tools_4.3.2            foreign_0.8-90         lmtest_0.9-40         
# [61] httpuv_1.6.16          future.apply_1.11.0    nnet_7.3-19            goftest_1.2-3         
# [65] quadprog_1.5-8         glue_1.8.0             nlme_3.1-168           promises_1.3.3        
# [69] grid_4.3.2             checkmate_2.3.2        Rtsne_0.17             cluster_2.1.6         
# [73] reshape2_1.4.4         generics_0.1.3         gtable_0.3.6           spatstat.data_3.1-6   
# [77] data.table_1.17.4      utf8_1.2.6             car_3.1-3              spatstat.geom_3.4-1   
# [81] RcppAnnoy_0.0.22       ggrepel_0.9.6          RANN_2.6.1             pillar_1.10.2         
# [85] stringr_1.5.1          spam_2.11-1            RcppHNSW_0.6.0         later_1.4.2           
# [89] splines_4.3.2          lattice_0.22-7         survival_3.8-3         deldir_1.0-9          
# [93] tidyselect_1.2.1       miniUI_0.1.1.1         pbapply_1.7-2          knitr_1.50            
# [97] gridExtra_2.3          scattermore_1.2        stats4_4.3.2           xfun_0.52             
# [101] matrixStats_1.5.0      stringi_1.8.7          lazyeval_0.2.2         evaluate_1.0.3        
# [105] codetools_0.2-19       tibble_3.2.1           cli_3.6.5              uwot_0.2.3            
# [109] rpart_4.1.24           xtable_1.8-4           reticulate_1.42.0      lavaan_0.6-19         
# [113] dichromat_2.0-0.1      Rcpp_1.0.14            globals_0.18.0         spatstat.random_3.4-1 
# [117] png_0.1-8              spatstat.univar_3.1-3  parallel_4.3.2         dotCall64_1.2         
# [121] jpeg_0.1-10            listenv_0.9.0          viridisLite_0.4.2      ggridges_0.5.4        
# [125] purrr_1.0.4            rlang_1.1.6            cowplot_1.1.3          mnormt_2.1.1 