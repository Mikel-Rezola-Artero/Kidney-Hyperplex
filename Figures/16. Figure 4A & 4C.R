####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(scales)

####2. LOAD ARTERY DATA####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(getActiveDocumentContext()$path),
                             "/Only arteries/Processed_Artery.rds"))


####3. FIGURE 4A####

#Plot Immune subcluster frequencies by artery and disease

#Extract data from Seurat object
cluster_counts <-  data.frame("Cells" = Idents(Seurat.obj),
                              "Patient" = Seurat.obj$slide,
                              "CD68" = Seurat.obj$CD68.Pos,
                              "C5aR1" = Seurat.obj$C5aR1.Pos,
                              "CD11b" = Seurat.obj$CD11b.Pos,
                              "CD3" = Seurat.obj$CD3.Pos,
                              "CD20" = Seurat.obj$CD20.Pos,
                              "MPO" = Seurat.obj$MPO.Pos,
                              "Artery" = Seurat.obj$Artery,
                              "Disease" = Seurat.obj$Disease)

#Modify Marker positivity to consider only the Immune cell cluster cells 
cluster_counts$MPO <- ifelse(cluster_counts$Cells == "Immune" & cluster_counts$MPO == "1",1,0)
cluster_counts$CD68 <- ifelse(cluster_counts$Cells == "Immune" & cluster_counts$CD68 == "1",1,0)
cluster_counts$C5aR1 <- ifelse(cluster_counts$Cells == "Immune" & cluster_counts$C5aR1 == "1",1,0)
cluster_counts$CD11b <- ifelse(cluster_counts$Cells == "Immune" & cluster_counts$CD11b == "1",1,0)
cluster_counts$CD3 <- ifelse(cluster_counts$Cells == "Immune" & cluster_counts$CD3 == "1",1,0)
cluster_counts$CD20 <- ifelse(cluster_counts$Cells == "Immune" & cluster_counts$CD20 == "1",1,0)

#Calculate frequency by Artery
positive_df <- cluster_counts %>%
  group_by(Disease, Artery, Patient) %>%
  summarise(
    MPO = mean(MPO),
    CD68 = mean(CD68),
    C5aR1 = mean(C5aR1),
    CD11b = mean(CD11b),
    CD3 = mean(CD3),
    CD20 = mean(CD20),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(-Artery, -Disease, -Patient),
               names_to = "Marker", values_to = "Frequency")

#Order factor for plot aesthetics
positive_df$Disease <- factor(positive_df$Disease, levels = c("LN", "Diabetes", "CCE"))
positive_df$Marker <- factor(positive_df$Marker,levels = c("MPO","CD68","C5aR1","CD11b","CD20","CD3"))


#Check by Kruskal Wallis if all markers' variability is explained by Disease
compare_means(
  Frequency ~ Disease,
  data = positive_df,
  group.by = c("Marker"),
  method = "kruskal.test",
  p.adjust.method = "bonferroni"
)
#  Marker      .y.      p     p.adj   p.format p.signif    method        
#  <fct>      <chr>   <dbl>   <dbl>     <chr>   <chr>       <chr>         
# 1 MPO    Frequency 6.88e- 1 1   e+ 0 0.6875   ns       Kruskal-Wallis
# 2 CD68   Frequency 8.50e- 8 5.10e- 7 8.5e-08  ****     Kruskal-Wallis
# 3 C5aR1  Frequency 2.85e-15 1.7 e-14 2.9e-15  ****     Kruskal-Wallis
# 4 CD11b  Frequency 3.97e- 3 2.4 e- 2 0.0040   **       Kruskal-Wallis
# 5 CD3    Frequency 1.00e- 1 6   e- 1 0.1002   ns       Kruskal-Wallis
# 6 CD20   Frequency 1.60e- 3 9.6 e- 3 0.0016   **       Kruskal-Wallis

#Calculate multiple comparisons using wicox test with Bonferroni Correction
comparisons_df <- compare_means(
  Frequency ~ Disease,
  data = positive_df,
  group.by = c("Marker"),
  ref.group = "CCE",
  method = "wilcox.test",
  p.adjust.method = "bonferroni"
)


#Set y coordinate for p.values
comparisons_df$y.position <- 0.19
comparisons_df$y.position[seq(1, nrow(comparisons_df), by = 2)] <- 0.20

#Plot Comparisons across disease (FIGURE 4A)
ggplot(positive_df, aes(x = Disease, y = Frequency, fill = Marker)) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.5) +
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.size = 0.2, alpha = 0.7) +
  stat_pvalue_manual(
    comparisons_df,
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    hide.ns = F) +
  facet_grid( ~ Marker, switch = "y") +
  scale_y_continuous(labels = percent_format(),
                     limits = c(0, 0.2)) +
  scale_fill_manual(values = c("C5aR1" = "lightblue","CD11b" ="green","MPO" = "magenta",
                               "CD20" = "orange", "CD3" = "red", "CD68" = "violet")) +
  theme_classic() +
  labs(
    title = "Immune Infiltrate in Arteries across Diseases",
    x = "Disease Group",
    y = "Cell Frequency (%)",
    fill = "Marker"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) + NoLegend()


####4. FIGURE 4C####

#Plot Associations between immune cells in CCE and C5aR1 and CD11b

#Obtain CCE Immune cell data
CCE.Imm <- Seurat.obj[, Seurat.obj$Disease == "CCE" & 
                        Idents(Seurat.obj) == "Immune"]

#Extract data from Seurat object
Imune_counts <-  data.frame("MPO" = CCE.Imm$MPO.Pos,
                            "CD20" = CCE.Imm$CD20.Pos,
                            "CD3" = CCE.Imm$CD3.Pos,
                            "CD68" = CCE.Imm$CD68.Pos,
                            "C5aR1" = CCE.Imm$C5aR1.Pos,
                            "CD11b" = CCE.Imm$CD11b.Pos)

#Run fisher test for p.value and O.R. calculation
#Set comparisons
marker_pairs <- list(
  c("C5aR1", "CD3"),c("C5aR1", "CD20"),c("C5aR1", "CD68"),c("C5aR1", "MPO"),
  c("CD11b", "CD3"),c("CD11b", "CD20"),c("CD11b", "CD68"),c("CD11b", "MPO")
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
           Marker2 %in% c("CD3", "CD20", "CD68","MPO"))

#Create bubble plot (FIGURE 4C)
ggplot(plot_df, aes(x = Marker1, y = Marker2)) +
  geom_point(aes(size = OR, fill = p.adj, alpha = log2(OR+1)), 
             shape = 21, color = "black") +
  scale_size_continuous(range = c(3, 12)) +
  scale_fill_gradient(low = "red", high = "white", 
                      name = "Adjusted\np-value") +
  theme_minimal() +
  labs(
    x = "Complement Markers",
    y = "Immune Markers",
    size = "Odds Ratio",
    title = "Complement Markers vs Immune Markers (CCE)"
  ) +
  theme(
    axis.text.x = element_text(angle = 0,size = 10, 
                               hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold",size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


####5. Session Information####

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
#   [1] scales_1.4.0       ggpubr_0.6.0       ggplot2_3.5.2      tidyr_1.3.1       
# [5] dplyr_1.1.4        Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0          
# [9] rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3         
# [4] rlang_1.1.6            magrittr_2.0.3         RcppAnnoy_0.0.22      
# [7] matrixStats_1.5.0      ggridges_0.5.4         compiler_4.3.2        
# [10] spatstat.geom_3.4-1    png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3       
# [16] fastmap_1.2.0          backports_1.5.0        labeling_0.4.3        
# [19] utf8_1.2.6             promises_1.3.3         purrr_1.0.4           
# [22] jsonlite_2.0.0         goftest_1.2-3          later_1.4.2           
# [25] spatstat.utils_3.1-4   broom_1.0.8            irlba_2.3.5.1         
# [28] parallel_4.3.2         cluster_2.1.6          R6_2.6.1              
# [31] ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3    
# [34] spatstat.data_3.1-6    reticulate_1.42.0      car_3.1-3             
# [37] parallelly_1.45.0      spatstat.univar_3.1-3  lmtest_0.9-40         
# [40] scattermore_1.2        Rcpp_1.0.14            tensor_1.5            
# [43] future.apply_1.11.0    zoo_1.8-12             sctransform_0.4.2     
# [46] httpuv_1.6.16          Matrix_1.6-5           splines_4.3.2         
# [49] igraph_2.1.4           tidyselect_1.2.1       dichromat_2.0-0.1     
# [52] abind_1.4-8            codetools_0.2-19       spatstat.random_3.4-1 
# [55] miniUI_0.1.1.1         spatstat.explore_3.4-3 listenv_0.9.0         
# [58] lattice_0.22-7         tibble_3.2.1           plyr_1.8.9            
# [61] withr_3.0.2            shiny_1.10.0           ROCR_1.0-11           
# [64] Rtsne_0.17             future_1.58.0          fastDummies_1.7.5     
# [67] survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.1-11   
# [70] pillar_1.10.2          carData_3.0-5          KernSmooth_2.23-26    
# [73] plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0        
# [76] globals_0.18.0         xtable_1.8-4           glue_1.8.0            
# [79] lazyeval_0.2.2         tools_4.3.2            data.table_1.17.4     
# [82] RSpectra_0.16-1        ggsignif_0.6.4         RANN_2.6.1            
# [85] dotCall64_1.2          cowplot_1.1.3          grid_4.3.2            
# [88] nlme_3.1-168           patchwork_1.3.0        Formula_1.2-5         
# [91] cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-1           
# [94] viridisLite_0.4.2      uwot_0.2.3             gtable_0.3.6          
# [97] rstatix_0.7.2          digest_0.6.37          progressr_0.15.1      
# [100] ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2          
# [103] htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7            
# [106] mime_0.13              MASS_7.3-60         