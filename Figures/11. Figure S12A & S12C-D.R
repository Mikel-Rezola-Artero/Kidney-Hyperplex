####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(scales)
library(glmmTMB)

####2. LOAD GLOMERULI DATA####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(getActiveDocumentContext()$path),
                             "/Only glomeruli/Processed_Glomeruli_v2.rds"))


#### 3. PREPARE DATA FOR BETA BINOMIAL MIXED EFFECTS MODELS & PLOT ####

#Extract counts
cluster_counts <- data.frame(
  Cells = Idents(Seurat.obj),
  Glomeruli = Seurat.obj$Glomeruli,
  Patient   = Seurat.obj$slide,
  Disease   = Seurat.obj$Disease,
  CD68      = Seurat.obj$CD68.Pos,
  CD3       = Seurat.obj$CD3.Pos,
  CD20      = Seurat.obj$CD20.Pos,
  MPO       = Seurat.obj$MPO.Pos
)

#Define marker-positive immune cells
cluster_counts <- cluster_counts %>%
  mutate(
    CD68 = ifelse(Cells == "Immune" & CD68 == "1", 1, 0),
    CD3 = ifelse(Cells == "Immune" & CD3 == "1", 1, 0),
    CD20 = ifelse(Cells == "Immune" & CD20 == "1", 1, 0),
    MPO = ifelse(Cells == "Immune" & MPO == "1", 1, 0))

#Convert to long format (one row per marker)
cell_long <- cluster_counts %>%
  pivot_longer(cols = c(CD68, CD3, CD20, MPO),
               names_to = "Marker",
               values_to = "Marker_Status")

#Calculate number of positive and negative marker immune cells per glomeruli
positive_df <- cell_long %>%
  group_by(Disease, Patient, Glomeruli, Marker) %>%
  summarise(Positive = sum(Marker_Status, na.rm = TRUE), # Number of marker-positive cells
            Total_Cells = n(), # Total number of cells in glomeruli
            .groups = "drop") %>%
  mutate(Negative = Total_Cells - Positive, # Number of marker-negative cells in glomeruli
         Frequency = Positive / Total_Cells) # Proportion of positive immune cells in glomeruli for visualization

#Set LN as the reference disease
positive_df$Disease <- factor(positive_df$Disease,
                              levels = c("LN","Diabetes","CCE"))
#Set immune marker order
marker_order <- c("CD68","MPO","CD3","CD20")
positive_df$Marker <- factor(positive_df$Marker,levels = marker_order)

#Make sure that patient of origin is a factor
positive_df$Patient <- factor(positive_df$Patient)


#### 4. FIT BETA-BINOMIAL MIXED-EFFECTS MODELS ####

bb_models <- list()
bb_comparisons <- list()

for (marker in levels(positive_df$Marker)) {
  
  # Load data for current marker
  tmp <- positive_df %>%
    filter(Marker == marker,
           !is.na(Positive),
           !is.na(Negative),
           !is.na(Disease),
           !is.na(Patient)) %>%
    droplevels()
  
  #Fit Beta Binomial using counts
  fit <- glmmTMB(
    cbind(Positive, Negative) ~ Disease + (1 | Patient),
    family = betabinomial(link = "logit"),
    data = tmp)
  # Outcome: marker-positive immune cells | all remaining cells in glomeruli
  # Covariates: 
  #  1) Disease -> LN as reference for pre-specified comparisons
  #  2) (1 | Patient) -> Patient included as a random effect with random intercept
  
  # Save model for convergence checks
  bb_models[[marker]] <- fit
  
  # Extract coefficient P values
  comp_df <- data.frame(contrast = c("LN vs Diabetes", "LN vs CCE"),
                        Marker = marker,
                        p.value = summary(fit)$coefficients$cond[c("DiseaseDiabetes","DiseaseCCE"),"Pr(>|z|)"])
  
  # Store coefficient contrast results
  bb_comparisons[[marker]] <- comp_df
  
}

#Check models' convergences
print(data.frame(Model = names(bb_models),
                 pdHess = sapply(bb_models,
                                 function(x) x$sdr$pdHess),
                 row.names = NULL))
#   Model pdHess
#    CD68   TRUE
#     CD3   TRUE
#    CD20   TRUE
#     MPO   TRUE

# Combine contrast results per marker
bb_results <- bind_rows(bb_comparisons)

#Adjust p values by BH method
bb_results$p.adj <- p.adjust(bb_results$p.value, method = "BH")

#Adjust significance labels to be displayed on top of the plot
bb_results <- bb_results %>% 
  mutate(p.signif = case_when(p.adj < 0.0001 ~ "****",
                              p.adj < 0.001  ~ "***",
                              p.adj < 0.01   ~ "**",
                              p.adj < 0.05   ~ "*",
                              TRUE           ~ "ns"))


#### 5. PREPARE PATIENT-LEVEL SUMMARIES FOR VISUALIZATION ####

patient_means <- positive_df %>%
  group_by(Patient, Disease, Marker) %>%
  summarise(
    Positive = sum(Positive),
    Total_Cells = sum(Total_Cells),
    Frequency = Positive / Total_Cells,
    n_glomeruli = n(),
    .groups = "drop")


#Set factor order accros patient level dataframes
patient_means$Disease <- factor(patient_means$Disease,
                                levels = c("LN", "Diabetes","CCE"))
patient_means$Marker <- factor(patient_means$Marker, levels = marker_order)


#### 6. PLOT FIGURE S12A ####

# Prepare p value positions for plot
marker_max <- positive_df %>%
  group_by(Marker) %>%
  summarise(max_frequency = max(Frequency, na.rm = TRUE),
            .groups = "drop")
comparisons_plot <- bb_results %>%
  left_join(marker_max, by = "Marker") %>%
  mutate(group1 = "LN",
         group2 = case_when(contrast == "LN vs Diabetes" ~ "Diabetes",
                            contrast == "LN vs CCE" ~ "CCE"),
         
         y.position = case_when(contrast == "LN vs Diabetes" ~ max_frequency * 1.10,
                                contrast == "LN vs CCE" ~ max_frequency * 1.20))

#Force marker order
positive_df$Marker <- factor(positive_df$Marker, levels = marker_order)
patient_means$Marker <- factor(patient_means$Marker, levels = marker_order)
comparisons_plot$Marker <- factor(comparisons_plot$Marker, levels = marker_order)


#Make plot 
p_immune_markers <- ggplot(positive_df,
                           aes(x = Disease, y = Frequency, fill = Marker)) +
  
  # Individual glomeruli plots
  geom_violin(width = 0.7, alpha = 0.3, trim = T) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.3) +
  
  # Patient-level plots
  geom_point(data = patient_means,
             aes(x = Disease, y = Frequency),
             inherit.aes = FALSE,
             position = position_jitter(width = 0.08, height = 0, seed = 123),
             shape = 21, size = 2.2, stroke = 0.7, fill = "white") +
  
  
  # Add Beta-binomial GLMM adjusted P values
  stat_pvalue_manual(comparisons_plot,
                     label = "p.signif",
                     xmin = "group1", xmax = "group2",
                     y.position = "y.position",
                     tip.length = 0.01, bracket.size = 0.5, hide.ns = TRUE) +
  
  # Facet by immune marker
  facet_wrap(~ Marker, ncol = 2, scales = "free_y") +
  
  # Y axis
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02,0.25))) +
  
  # Marker colours
  scale_fill_manual(values = c("CD68" = "violet",
                               "CD3"  = "red",
                               "CD20" = "orange",
                               "MPO"  = "magenta")) +
  
  # Labels
  labs(title = "Immune Cell Populations in Glomeruli across Diseases",
       x = "Disease Group",
       y = "Marker-Positive Immune Cells (% of Glomerular Cells)",
       fill = "Marker") +
  
  # Theme
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,face = "bold"))


p_immune_markers


####7. FIGURE S12C####

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
  ) #FIGURE S12C


####8. FIGURE S12D####

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
  )#FIGURE S12D


####9. Session Information#####

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
# [1] glmmTMB_1.1.14     scales_1.4.0       ggpubr_1.0.0       ggplot2_4.0.3      dplyr_1.2.1        Seurat_5.5.1       SeuratObject_5.4.0 sp_2.2-3          
# [9] tidyr_1.3.2  SeuratObject_5.1.0 sp_2.2-0           rstudioapi_0.17.1 