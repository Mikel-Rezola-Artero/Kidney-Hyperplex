####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(ggplot2)
library(scales)
library(tidyr)
library(dplyr)
library(ggpubr)
library(glmmTMB)


####2. LOAD ARTERY DATA####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(getActiveDocumentContext()$path),
                             "/Only arteries/Processed_Artery_v2.rds"))


####3. FIGURE S15B####

#Plot annotated clusters as DimPlot (FIGURE S15B)
DimPlot(object = Seurat.obj,repel = T,
        label = T,label.box = T,label.size = 5,
        cols = c("grey","gold","red","violet"),
        pt.size = 0.1,reduction = "umap") + NoLegend()


####4. FIGURE S15C####

#Check distribution of key cell markers as DotPlot (FIGURE S15C)
DotPlot(Seurat.obj,
        c("aSMA",
          "CD31","VCAM1",
          "CD45","CD68","C5aR1","CD11b","CD163","CD3","CD20"
        ),
        scale = T,
        col.min = 0.6,
        col.max = 1.5,
        dot.min = F)+ scale_size(range = c(15, 15)) + guides(size = "none") +
  scale_color_gradient(low = "white",
                       high = "blue",
                       limits = c(0.6, 1.5),
                       oob = squish) + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())

#### 5. PREPARE ARTERY DATA FOR FIGURE S15D ####

# Filter endothelial and vascular smooth muscle cells (vSMC)
Artery.cells <- Seurat.obj[,Idents(Seurat.obj) %in% c("Endothelial","vSMC")]

# Extract data
cluster_counts <- data.frame(
  Cells   = Idents(Artery.cells),
  Artery  = Artery.cells$Artery,
  C1q     = Artery.cells$C1q.Pos,
  MBL     = Artery.cells$MBL.Pos,
  C4d     = Artery.cells$C4d.Pos,
  C3c     = Artery.cells$C3c.Pos,
  C3d     = Artery.cells$C3d.Pos,
  C9      = Artery.cells$C9.Pos,
  Patient = Artery.cells$slide,
  Disease = Artery.cells$Disease
)

# Calculate marker-positive counts within each cell type for each artery
count_df <- cluster_counts %>%
  group_by(Cells, Disease, Artery, Patient) %>%
  summarise(
    # Pos cell counts
    C1q_Pos = sum(C1q, na.rm = TRUE),
    MBL_Pos = sum(MBL, na.rm = TRUE),
    C4d_Pos = sum(C4d, na.rm = TRUE),
    C3c_Pos = sum(C3c, na.rm = TRUE),
    C3d_Pos = sum(C3d, na.rm = TRUE),
    C9_Pos  = sum(C9,  na.rm = TRUE),
    # Total number of specific cell type within specific artery
    Total_Cells = n(),
    .groups = "drop"
  )

# Convert to long format
positive_df <- count_df %>%
  pivot_longer(cols = c(C1q_Pos, MBL_Pos, C4d_Pos, C3c_Pos, C3d_Pos, C9_Pos),
               names_to = "Marker", values_to = "Positive") %>%
  mutate(Marker = sub("_Pos$", "", Marker),
         Negative = Total_Cells - Positive,
         Frequency = Positive / Total_Cells)


# Set factor order for plot
positive_df$Cells <- factor(positive_df$Cells, levels = c("Endothelial","vSMC"))
marker_order <- c("C1q", "MBL", "C4d", "C3c", "C3d", "C9")
positive_df$Marker <- factor(positive_df$Marker, levels = marker_order)

# Set CCE as the reference group for the model testing
positive_df$Disease <- factor(positive_df$Disease, 
                              levels = c("CCE", "LN", "Diabetes"))
positive_df$Patient <- factor(positive_df$Patient)


#### 6. BETA-BINOMIAL MIXED-EFFECTS MODELS FOR FIGURE S15D ####

bb_models <- list()
bb_comparisons <- list()

for (cell_type in levels(positive_df$Cells)) {
  
  for (marker in levels(positive_df$Marker)) {
    
    # Load data for current cell type × marker
    tmp <- positive_df %>% 
      filter(Cells == cell_type,
             Marker == marker,
             !is.na(Positive),
             !is.na(Negative),
             !is.na(Disease),
             !is.na(Patient)) %>%
      droplevels()
    
    #Fit Beta Binomial using counts
    fit <- glmmTMB(cbind(Positive,Negative) ~ Disease + (1 | Patient),
                   family = betabinomial(link = "logit"),
                   data = tmp)
    # Outcome: marker-positive cells | all remaining cells
    # Covariates: 
    #  1) Disease -> CCE as reference for pre-specified comparisons
    #  2) (1 | Patient) -> Patient included as a random effect with random intercept
    
    # Save model name
    model_name <- paste(cell_type, marker, sep = "_")
    
    # Save model for convergence checks
    bb_models[[model_name]] <- fit
    
    # Extract coefficient P values
    comp_df <- data.frame(contrast = c("CCE vs LN", "CCE vs Diabetes"),
                          Cells = cell_type,
                          Marker = marker,
                          p.value = summary(fit)$coefficients$cond[c("DiseaseLN","DiseaseDiabetes"), "Pr(>|z|)"])
    bb_comparisons[[model_name]] <- comp_df
  }
}


# Check convergence
print(
  data.frame(
    Model = names(bb_models),
    pdHess = sapply(
      bb_models,
      function(x) x$sdr$pdHess
    ),
    row.names = NULL
  )
)
#            Model pdHess
#  Endothelial_C1q   TRUE
#  Endothelial_MBL   TRUE
#  Endothelial_C4d   TRUE
#  Endothelial_C3c   TRUE
#  Endothelial_C3d   TRUE
#   Endothelial_C9   TRUE
#         vSMC_C1q   TRUE
#         vSMC_MBL   TRUE
#         vSMC_C4d   TRUE
#         vSMC_C3c   TRUE
#         vSMC_C3d   TRUE
#          vSMC_C9   TRUE

#Apply FDR p value correction
bb_results <- bind_rows(bb_comparisons)
bb_results$p.adj <- p.adjust(bb_results$p.value,method = "BH")

#Create significance labels from adj.p for the plot
bb_results <- bb_results %>%
  mutate(p.signif = case_when(p.adj < 0.0001 ~ "****",
                              p.adj < 0.001  ~ "***",
                              p.adj < 0.01   ~ "**",
                              p.adj < 0.05   ~ "*",
                              TRUE           ~ "ns"))

#### 7. PATIENT-LEVEL VALUES FOR VISUALIZATION FOR FIGURE S15D ####

patient_means <- positive_df %>%
  group_by(Patient, Disease, Cells, Marker) %>%
  summarise(Positive = sum(Positive),
            Total_Cells = sum(Total_Cells),
            Frequency = Positive / Total_Cells,
            .groups = "drop")


#### 8. PREPARE RESULTS FOR PLOTTING FOR FIGURE S15D ####

# Set p value comparisons and y axis locations for plots
comparisons_plot <- bb_results %>% 
  mutate(group1 = "CCE",
         group2 = case_when(contrast == "CCE vs LN" ~ "LN",
                            contrast == "CCE vs Diabetes" ~ "Diabetes"),
         y.position = case_when(contrast == "CCE vs LN"       ~ 1.05,
                                contrast == "CCE vs Diabetes" ~ 1.15))

# Force facet order across dataframes
positive_df$Marker <- factor(positive_df$Marker, levels = marker_order)
patient_means$Marker <- factor(patient_means$Marker, levels = marker_order)
comparisons_plot$Marker <- factor(comparisons_plot$Marker, levels = marker_order)
positive_df$Cells <- factor(positive_df$Cells, levels = c("Endothelial", "vSMC"))
patient_means$Cells <- factor(patient_means$Cells, levels = c("Endothelial", "vSMC"))
comparisons_plot$Cells <- factor(comparisons_plot$Cells, levels = c("Endothelial", "vSMC"))


#### 9. PLOT FIGURE S15D ####

p_artery <- ggplot(positive_df, aes(x = Disease, y = Frequency, fill = Marker)) +
  
  #Individual artery plots
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.3) +
  
  #Patient-level plot
  geom_point(data = patient_means, aes(x = Disease, y = Frequency),
             inherit.aes = FALSE, position = position_jitter(width = 0.08, height = 0),
             shape = 21, size = 2.2, stroke = 0.7, fill = "white") +
  
  #Add Beta-binomial GLMM comparison adjusted p values
  stat_pvalue_manual(comparisons_plot, label = "p.signif",
                     xmin = "group1", xmax = "group2",
                     y.position = "y.position",
                     tip.length = 0.01, bracket.size = 0.5, hide.ns = TRUE) +
  #Facet by cell and marker
  facet_grid(Cells ~ Marker, switch = "y") +
  
  #Y axis
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(-0.01, 1.20),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  
  #Marker colours
  scale_fill_manual(values = c("C1q" = "lightblue",
                               "MBL" = "brown",
                               "C4d" = "violet",
                               "C3c" = "orange",
                               "C3d" = "red",
                               "C9"  = "royalblue")) +
  
  # Plot Labels
  labs(title = "Complement Deposits in Arteries across Diseases",
       x = "Disease Group",
       y = "Complement-Positive Cells (%)",
       fill = "Marker") +
  
  # Theme
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")


p_artery


####10. Session Information#####

sessionInfo()

# R version 4.6.1 (2026-06-24)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 24.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
# 
# locale:
# [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Paris
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] glmmTMB_1.1.14     ggpubr_1.0.0       dplyr_1.2.1        tidyr_1.3.2        scales_1.4.0       ggplot2_4.0.3      Seurat_5.5.1       SeuratObject_5.4.0
# [9] sp_2.2-3           rstudioapi_0.19.0 