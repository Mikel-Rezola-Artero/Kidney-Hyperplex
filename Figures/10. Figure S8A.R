####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggpubr)
library(scales)
library(ggplot2)

####2. LOAD GLOMERULI DATA####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(getActiveDocumentContext()$path),
                             "/Only glomeruli/Processed_Glomeruli.rds"))


####3. FIGURE S8A####

#Show C deposits on for Endothelial Cells, Mesangial Cells & "PECs/PT-Conv" at Glomeruli level
Endo.Mes.PECs <- Seurat.obj[ ,
                             (Idents(Seurat.obj) == "Endothelial" |
                                Idents(Seurat.obj) == "Mesangial" |
                                Idents(Seurat.obj) == "PECs/PT_VCAM1") ]
cluster_counts <-  data.frame("Cells" = Idents(Endo.Mes.PECs),
                              "Glomeruli" = Endo.Mes.PECs$Glomeruli,
                              "C1q" = Endo.Mes.PECs$C1q.Pos,
                              "C4d" = Endo.Mes.PECs$C4d.Pos,
                              "C3c" = Endo.Mes.PECs$C3c.Pos,
                              "C3d" = Endo.Mes.PECs$C3d.Pos,
                              "C9" = Endo.Mes.PECs$C9.Pos,
                              "Patient" = Endo.Mes.PECs$slide,
                              "Disease" = Endo.Mes.PECs$Disease
)

#Summarize by Glomeruli, Cell Type, Patient and Disease
positive_df <- cluster_counts %>%
  group_by(Cells, Disease, Glomeruli, Patient) %>%
  summarise(
    C3c = mean(C3c),
    C3d = mean(C3d),
    C4d = mean(C4d),
    C9 = mean(C9),
    C1q = mean(C1q),
    .groups = "drop" ) %>% 
  pivot_longer(cols = c(-Cells, -Glomeruli, -Disease, -Patient),
               names_to = "Marker", 
               values_to = "Frequency")
  

#Order factor for plot aesthetics
positive_df$Cells <- factor(positive_df$Cells,levels = c("Endothelial","Mesangial","PECs/PT_VCAM1"))
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
# # A tibble: 15 × 8
#      Cells       Marker    .y.        p      p.adj  p.format p.signif    method        
#      <fct>      <fct>     <chr>     <dbl>    <dbl>    <chr>  <chr>       <chr>         
#   1 Endothelial   C3c    Frequency 2.55e-15 3.8 e-14 2.6e-15  ****     Kruskal-Wallis
# 2 Endothelial   C3d    Frequency 1.71e-23 2.6 e-22 < 2e-16  ****     Kruskal-Wallis
# 3 Endothelial   C4d    Frequency 3.54e-23 5.30e-22 < 2e-16  ****     Kruskal-Wallis
# 4 Endothelial   C9     Frequency 1.32e-20 2   e-19 < 2e-16  ****     Kruskal-Wallis
# 5 Endothelial   C1q    Frequency 3.68e-16 5.50e-15 3.7e-16  ****     Kruskal-Wallis
# 6 Mesangial     C3c    Frequency 2.47e-19 3.7 e-18 < 2e-16  ****     Kruskal-Wallis
# 7 Mesangial     C3d    Frequency 2.82e-20 4.20e-19 < 2e-16  ****     Kruskal-Wallis
# 8 Mesangial     C4d    Frequency 1.72e-21 2.6 e-20 < 2e-16  ****     Kruskal-Wallis
# 9 Mesangial     C9     Frequency 2.38e-18 3.60e-17 < 2e-16  ****     Kruskal-Wallis
# 10 Mesangial     C1q    Frequency 1.47e-16 2.20e-15 < 2e-16  ****     Kruskal-Wallis
# 11 PECs/PT_VCAM1 C3c    Frequency 7.66e-10 1.10e- 8 7.7e-10  ****     Kruskal-Wallis
# 12 PECs/PT_VCAM1 C3d    Frequency 7.37e-12 1.10e-10 7.4e-12  ****     Kruskal-Wallis
# 13 PECs/PT_VCAM1 C4d    Frequency 1.30e-10 1.9 e- 9 1.3e-10  ****     Kruskal-Wallis
# 14 PECs/PT_VCAM1 C9     Frequency 3.69e-11 5.50e-10 3.7e-11  ****     Kruskal-Wallis
# 15 PECs/PT_VCAM1 C1q    Frequency 1.20e-10 1.80e- 9 1.2e-10  ****     Kruskal-Wallis

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
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.2)) +
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
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.2),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
  )


####4. Session Information#####

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
#   [1] scales_1.4.0       tidyr_1.3.1        dplyr_1.1.4        ggpubr_0.6.0       ggplot2_3.5.2     
# [6] Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0           rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.6           
# [5] magrittr_2.0.3         RcppAnnoy_0.0.22       matrixStats_1.5.0      ggridges_0.5.4        
# [9] compiler_4.3.2         spatstat.geom_3.4-1    png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3        fastmap_1.2.0         
# [17] backports_1.5.0        utf8_1.2.6             promises_1.3.3         purrr_1.0.4           
# [21] jsonlite_2.0.0         goftest_1.2-3          later_1.4.2            spatstat.utils_3.1-4  
# [25] broom_1.0.8            irlba_2.3.5.1          parallel_4.3.2         cluster_2.1.6         
# [29] R6_2.6.1               ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3    
# [33] spatstat.data_3.1-6    reticulate_1.42.0      car_3.1-3              parallelly_1.45.0     
# [37] spatstat.univar_3.1-3  lmtest_0.9-40          scattermore_1.2        Rcpp_1.0.14           
# [41] tensor_1.5             future.apply_1.11.0    zoo_1.8-12             sctransform_0.4.2     
# [45] httpuv_1.6.16          Matrix_1.6-5           splines_4.3.2          igraph_2.1.4          
# [49] tidyselect_1.2.1       dichromat_2.0-0.1      abind_1.4-8            codetools_0.2-19      
# [53] spatstat.random_3.4-1  miniUI_0.1.1.1         spatstat.explore_3.4-3 listenv_0.9.0         
# [57] lattice_0.22-7         tibble_3.2.1           plyr_1.8.9             withr_3.0.2           
# [61] shiny_1.10.0           ROCR_1.0-11            Rtsne_0.17             future_1.58.0         
# [65] fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.1-11   
# [69] pillar_1.10.2          carData_3.0-5          KernSmooth_2.23-26     plotly_4.10.4         
# [73] generics_0.1.3         RcppHNSW_0.6.0         globals_0.18.0         xtable_1.8-4          
# [77] glue_1.8.0             lazyeval_0.2.2         tools_4.3.2            data.table_1.17.4     
# [81] RSpectra_0.16-1        ggsignif_0.6.4         RANN_2.6.1             dotCall64_1.2         
# [85] cowplot_1.1.3          grid_4.3.2             nlme_3.1-168           patchwork_1.3.0       
# [89] Formula_1.2-5          cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-1           
# [93] viridisLite_0.4.2      uwot_0.2.3             gtable_0.3.6           rstatix_0.7.2         
# [97] digest_0.6.37          progressr_0.15.1       ggrepel_0.9.6          htmlwidgets_1.6.4     
# [101] farver_2.1.2           htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7            
# [105] mime_0.13              MASS_7.3-60 
