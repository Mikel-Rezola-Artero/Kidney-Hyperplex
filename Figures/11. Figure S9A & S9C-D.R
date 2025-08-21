####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(scales)

####2. LOAD GLOMERULI DATA####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(getActiveDocumentContext()$path),
                             "/Only glomeruli/Processed_Glomeruli.rds"))


####3. FIGURE S9A####

#Extract data from Seurat object
cluster_counts <-  data.frame("Cells" = Idents(Seurat.obj),
                              "Glomeruli" = Seurat.obj$Glomeruli,
                              "Disease" = Seurat.obj$Disease,
                              "Patient" = Seurat.obj$slide
)

#Calculate immune cell frequency per glomeruli in each slide
Immune_freq <- cluster_counts %>%
  group_by(Disease, Glomeruli, Cells, Patient) %>% 
  summarise(Cell_Count = n(),
            .groups = "drop") %>% 
  group_by(Disease, Glomeruli, Patient) %>% 
  mutate(Frequency = Cell_Count / sum(Cell_Count),
         .groups = "drop") %>% 
  filter(Cells == "Immune")

#Set disease order for plot aesthetics
Immune_freq$Disease <- factor(Immune_freq$Disease, levels = c("LN","Diabetes","CCE"))

#Check by Kruskal Wallis if all markers' variability is explained by Disease
compare_means(
  Frequency ~ Disease,
  data = Immune_freq,
  method = "kruskal.test",
  p.adjust.method = "bonferroni")
# # A tibble: 1 × 6
#   .y.                p       p.adj  p.format p.signif    method        
#   <chr>           <dbl>      <dbl>    <chr>   <chr>       <chr>         
#   1 Frequency 0.000000142 0.00000014 1.4e-07  ****     Kruskal-Wallis

#Plot FIGURE S8A
ggplot(Immune_freq, aes(x = Disease, y = Frequency, fill = Disease)) +
  geom_jitter(size = 0.4) + geom_violin(width = 0.7) + 
  geom_boxplot(width = 0.5,outlier.size = 0.2,alpha = 0.7) +
  scale_y_continuous(labels = percent_format()) +
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


####4. FIGURE S9C####

DotPlot(Seurat.obj[, Idents(Seurat.obj) == "Immune"],
        features = c("CD68.Pos", "CD163.Pos", "CD11b.Pos", "C5aR1.Pos",
                     "MPO.Pos", "CD3.Pos", "CD20.Pos"),
        group.by = "Disease", 
        scale = FALSE) + 
  scale_size(range = c(2, 12)) +  
  scale_color_gradient(low = "white", high = "blue", oob = squish) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  guides(
    size = guide_legend(title = "Percentage Positive"),
    color = guide_colorbar(title = "Average Expression")
  ) #FIGURE S9C


####5. FIGURE S9D####

#Obtain Lupus Immune cell data
Lupus.Imm <- Seurat.obj[ , Seurat.obj$Disease == "LN" & 
                           Idents(Seurat.obj) == "Immune"]

#Extract immune marker data from Seurat object
Imune_counts <-  data.frame("CD20" = Lupus.Imm$CD20.Pos,
                            "CD3" = Lupus.Imm$CD3.Pos,
                            "MPO" = Lupus.Imm$MPO.Pos,
                            "CD68" = Lupus.Imm$CD68.Pos,
                            "C5aR1" = Lupus.Imm$C5aR1.Pos,
                            "CD11b" = Lupus.Imm$CD11b.Pos)

#Run fisher test for p.value and O.R. calculation
#Set comparisons
marker_pairs <- list(
  c("C5aR1", "CD3"),c("C5aR1", "CD20"),c("C5aR1", "MPO"),c("C5aR1", "CD68"),
  c("CD11b", "CD3"),c("CD11b", "CD20"),c("CD11b", "MPO"),c("CD11b", "CD68")
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
           Marker2 %in% c("CD3", "CD20", "MPO", "CD68"))

#Re-order for plot aesthetics
plot_df$Marker2 <- factor(plot_df$Marker2, levels = c("CD3","CD20","MPO","CD68"))

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
  )#FIGURE S9D


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
#   [1] scales_1.4.0       ggpubr_0.6.0       ggplot2_3.5.2      tidyr_1.3.1        dplyr_1.1.4       
# [6] Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0           rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.6           
# [5] magrittr_2.0.3         RcppAnnoy_0.0.22       matrixStats_1.5.0      ggridges_0.5.4        
# [9] compiler_4.3.2         spatstat.geom_3.4-1    png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3        fastmap_1.2.0         
# [17] backports_1.5.0        labeling_0.4.3         utf8_1.2.6             promises_1.3.3        
# [21] purrr_1.0.4            jsonlite_2.0.0         goftest_1.2-3          later_1.4.2           
# [25] spatstat.utils_3.1-4   broom_1.0.8            irlba_2.3.5.1          parallel_4.3.2        
# [29] cluster_2.1.6          R6_2.6.1               ica_1.0-3              stringi_1.8.7         
# [33] RColorBrewer_1.1-3     spatstat.data_3.1-6    reticulate_1.42.0      car_3.1-3             
# [37] parallelly_1.45.0      spatstat.univar_3.1-3  lmtest_0.9-40          scattermore_1.2       
# [41] Rcpp_1.0.14            tensor_1.5             future.apply_1.11.0    zoo_1.8-12            
# [45] sctransform_0.4.2      httpuv_1.6.16          Matrix_1.6-5           splines_4.3.2         
# [49] igraph_2.1.4           tidyselect_1.2.1       dichromat_2.0-0.1      abind_1.4-8           
# [53] codetools_0.2-19       spatstat.random_3.4-1  miniUI_0.1.1.1         spatstat.explore_3.4-3
# [57] listenv_0.9.0          lattice_0.22-7         tibble_3.2.1           plyr_1.8.9            
# [61] withr_3.0.2            shiny_1.10.0           ROCR_1.0-11            Rtsne_0.17            
# [65] future_1.58.0          fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7       
# [69] fitdistrplus_1.1-11    pillar_1.10.2          carData_3.0-5          KernSmooth_2.23-26    
# [73] plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0         globals_0.18.0        
# [77] xtable_1.8-4           glue_1.8.0             lazyeval_0.2.2         tools_4.3.2           
# [81] data.table_1.17.4      RSpectra_0.16-1        ggsignif_0.6.4         RANN_2.6.1            
# [85] dotCall64_1.2          cowplot_1.1.3          grid_4.3.2             nlme_3.1-168          
# [89] patchwork_1.3.0        Formula_1.2-5          cli_3.6.5              spatstat.sparse_3.1-0 
# [93] spam_2.11-1            viridisLite_0.4.2      uwot_0.2.3             gtable_0.3.6          
# [97] rstatix_0.7.2          digest_0.6.37          progressr_0.15.1       ggrepel_0.9.6         
# [101] htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1      lifecycle_1.0.4       
# [105] httr_1.4.7             mime_0.13              MASS_7.3-60  
