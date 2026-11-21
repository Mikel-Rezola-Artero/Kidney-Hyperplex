#### 1. LOAD LIBRARIES ####
library(rstudioapi)
library(tidyr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(glmmTMB)

####2. LOAD GLOMERULI DATA####

#Load Processed data
setwd(dirname(getActiveDocumentContext()$path))
Seurat.obj <- readRDS("../Data2/Processed_Glomeruli_v2.rds")


#### 3. PREPARE DATA FOR FIGURE S11A ####

# Filter endothelial, mesangial and PECs/PT_VCAM1 cells
Endo.Mes.PECs <- Seurat.obj[, Idents(Seurat.obj) %in% c("Endothelial", "Mesangial", "PECs/PT_VCAM1")]

# Extract data
cluster_counts <- data.frame(
  Cells     = Idents(Endo.Mes.PECs),
  Glomeruli = Endo.Mes.PECs$Glomeruli,
  C1q       = Endo.Mes.PECs$C1q.Pos,
  MBL       = Endo.Mes.PECs$MBL.Pos,
  C4d       = Endo.Mes.PECs$C4d.Pos,
  C3c       = Endo.Mes.PECs$C3c.Pos,
  C3d       = Endo.Mes.PECs$C3d.Pos,
  C9        = Endo.Mes.PECs$C9.Pos,
  Patient   = Endo.Mes.PECs$slide,
  Disease   = Endo.Mes.PECs$Disease
)

# Calculate marker-positive counts within each cell type for each glomerulus
count_df <- cluster_counts %>%
  group_by(Cells,Disease,Glomeruli,Patient) %>%
  summarise(
    
    # Pos cell counts
    C1q_Pos = sum(C1q, na.rm = TRUE), 
    MBL_Pos = sum(MBL, na.rm = TRUE),
    C4d_Pos = sum(C4d, na.rm = TRUE),
    C3c_Pos = sum(C3c, na.rm = TRUE),
    C3d_Pos = sum(C3d, na.rm = TRUE),
    C9_Pos  = sum(C9,  na.rm = TRUE),
    
    # Total number of specific cell type within specific glomerulus
    Total_Cells = n(),
    
    .groups = "drop"
  )


# Convert to long format
positive_df <- count_df %>% 
  pivot_longer(cols = c(C1q_Pos, MBL_Pos, C4d_Pos, C3c_Pos, C3d_Pos, C9_Pos),
               names_to = "Marker", values_to = "Positive") %>%
  mutate(Marker = sub("_Pos$", "", Marker), # Clean marker names
         Negative = Total_Cells - Positive, #Rest of specific cells not positive for current marker
         Frequency = Positive / Total_Cells)# Glomerulus-level proportion for visualization

# Set factor order for plots
positive_df$Cells <- factor(positive_df$Cells,
                            levels = c("Endothelial","Mesangial","PECs/PT_VCAM1"))
marker_order <- c("C1q","MBL","C4d","C3c","C3d","C9")
positive_df$Marker <- factor(positive_df$Marker,levels = marker_order)

# Set LN as the reference disease
positive_df$Disease <- factor(positive_df$Disease, levels = c("LN", "Diabetes", "CCE"))

# Make sure patient of origin is a factor
positive_df$Patient <- factor(positive_df$Patient)


#### 4. BETA-BINOMIAL MIXED-EFFECTS MODEL PER CELL TYPE × MARKER ####

bb_models <- list()
bb_comparisons <- list()

for (cell_type in levels(positive_df$Cells)) {
  
  for (marker in levels(positive_df$Marker)) {
    
    # Load data for current cell type × marker
    tmp <- positive_df %>%
      filter(
        Cells == cell_type,
        Marker == marker,
        !is.na(Positive),
        !is.na(Negative),
        !is.na(Disease),
        !is.na(Patient)
      ) %>%
      droplevels()
    
    #Fit Beta Binomial using counts
    fit <- glmmTMB(cbind(Positive,Negative) ~ Disease + (1 | Patient),
                   family = betabinomial(link = "logit"), 
                   data = tmp)
    # Outcome: marker-positive cells | all remaining cells
    # Covariates: 
    #  1) Disease -> LN as reference for pre-specified comparisons
    #  2) (1 | Patient) -> Patient included as a random effect with random intercept
    
    # Save model name
    model_name <- paste(cell_type, marker, sep = "_")
    
    # Save model for convergence checks
    bb_models[[model_name]] <- fit
    
    # Extract coefficient P values
    comp_df <- data.frame(contrast = c("LN vs Diabetes", "LN vs CCE"),
                          Cells = cell_type,
                          Marker = marker,
                          p.value = summary(fit)$coefficients$cond[c("DiseaseDiabetes","DiseaseCCE"),"Pr(>|z|)"])
    
    # Store coefficient contrast results
    bb_comparisons[[model_name]] <- comp_df
  }
}

# Check models' convergences
print(data.frame(Model = names(bb_models),
                 pdHess = sapply(bb_models,
                                 function(x) x$sdr$pdHess),
                 row.names = NULL))
#              Model pdHess
#    Endothelial_C1q   TRUE
#    Endothelial_MBL   TRUE
#    Endothelial_C4d   TRUE
#    Endothelial_C3c   TRUE
#    Endothelial_C3d   TRUE
#     Endothelial_C9   TRUE
#      Mesangial_C1q   TRUE
#      Mesangial_MBL   TRUE
#      Mesangial_C4d   TRUE
#      Mesangial_C3c   TRUE
#      Mesangial_C3d   TRUE
#       Mesangial_C9   TRUE
#  PECs/PT_VCAM1_C1q   TRUE
#  PECs/PT_VCAM1_MBL   TRUE
#  PECs/PT_VCAM1_C4d   TRUE
#  PECs/PT_VCAM1_C3c   TRUE
#  PECs/PT_VCAM1_C3d   TRUE
#   PECs/PT_VCAM1_C9   TRUE

# Combine all raw coefficient P values
bb_results <- bind_rows(bb_comparisons)

# Apply Global BH/FDR correction
bb_results$p.adj <- p.adjust(bb_results$p.value, method = "BH")

# Add significance labels
bb_results <- bb_results %>%
  mutate(p.signif = case_when(p.adj < 0.0001 ~ "****",
                              p.adj < 0.001  ~ "***",
                              p.adj < 0.01   ~ "**",
                              p.adj < 0.05   ~ "*",
                              TRUE           ~ "ns"))


#### 5. PATIENT-LEVEL VALUES FOR VISUALIZATION ####

patient_means <- positive_df %>%
  group_by(Patient, Disease, Cells, Marker) %>%
  summarise(Frequency = sum(Positive) / sum(Total_Cells), .groups = "drop")

#### 6. PREPARE RESULTS FOR PLOTTING ####

# Set p value comparisons and y axis locations for plots
comparisons_plot <- bb_results %>%
  mutate(group1 = "LN",
         group2 = case_when(contrast == "LN vs Diabetes" ~ "Diabetes",
                            contrast == "LN vs CCE"      ~ "CCE"),
         y.position = case_when(contrast == "LN vs Diabetes" ~ 1.05,
                                contrast == "LN vs CCE"      ~ 1.15))

# Force facet order across dataframes
marker_order <- c("C1q", "MBL", "C4d", "C3c", "C3d", "C9")
positive_df$Marker <- factor(positive_df$Marker, levels = marker_order)
patient_means$Marker <- factor(patient_means$Marker, levels = marker_order)
comparisons_plot$Marker <- factor(comparisons_plot$Marker, levels = marker_order)
positive_df$Cells <- factor(positive_df$Cells, 
                            levels = c("Endothelial", "Mesangial", "PECs/PT_VCAM1"))
patient_means$Cells <- factor(patient_means$Cells, levels = levels(positive_df$Cells))
comparisons_plot$Cells <- factor(comparisons_plot$Cells, levels = levels(positive_df$Cells))


#### 7. FIGURE S11A ####

p_S11A <- ggplot(positive_df,
                 aes(x = Disease, y = Frequency, fill = Marker)) +
  # Individual glomerulus plots
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.3) +
  
  # Patient-level values
  geom_point(data = patient_means, aes(x = Disease, y = Frequency),
             inherit.aes = FALSE, 
             position = position_jitter(width = 0.08, height = 0),
             shape = 21, size = 2.2, stroke = 0.7, fill = "white") +
  
  # Add Beta-binomial GLMM-derived globally adjusted P values
  stat_pvalue_manual(comparisons_plot,
                     label = "p.signif",
                     xmin = "group1", xmax = "group2",
                     y.position = "y.position",
                     tip.length = 0.01, bracket.size = 0.5, hide.ns = TRUE) +
  
  # Facet by cell type × complement marker
  facet_grid(Cells ~ Marker, switch = "y") +
  
  # Y axis
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(-0.01, 1.20),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  # Marker colours
  scale_fill_manual(values = c("C1q" = "lightblue",
                               "MBL" = "brown",
                               "C4d" = "violet",
                               "C3c" = "orange",
                               "C3d" = "red",
                               "C9"  = "royalblue")) +
  
  # Labels
  labs(title = "Complement Deposits on Glomeruli across Diseases",
       x = "Disease Group",
       y = "Complement-Positive Cells (%)",
       fill = "Marker") +
  
  # Theme
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))


p_S11A


#### 8. SESSION INFORMATION ####

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
# [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
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
# [1] glmmTMB_1.1.14     scales_1.4.0       ggpubr_1.0.0       ggplot2_4.0.3      dplyr_1.2.1        Seurat_5.5.1      
# [7] SeuratObject_5.4.0 sp_2.2-3           tidyr_1.3.2       rstudioapi_0.17.1  