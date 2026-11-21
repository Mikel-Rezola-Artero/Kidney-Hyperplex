####1. LOAD LIBRARIES####
library(rstudioapi)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(scales)
library(glmmTMB)


####2. LOAD ARTERY DATA####

#Load Processed data
Seurat.obj <- readRDS(paste0(dirname(getActiveDocumentContext()$path),
                             "/Only arteries/Processed_Artery_v2.rds"))


####3. FIGURE 4A####

#Plot Immune cell frequencies by artery and disease

# Extract data
cluster_counts <- data.frame(
  Cells   = Idents(Seurat.obj),
  Patient = Seurat.obj$slide,
  MPO     = Seurat.obj$MPO.Pos,
  CD68    = Seurat.obj$CD68.Pos,
  C5aR1   = Seurat.obj$C5aR1.Pos,
  CD11b   = Seurat.obj$CD11b.Pos,
  CD3     = Seurat.obj$CD3.Pos,
  CD20    = Seurat.obj$CD20.Pos,
  Artery  = Seurat.obj$Artery,
  Disease = Seurat.obj$Disease
)

# Define marker-positive immune cells (MPO, CD68, C5aR1, CD11b, CD20, CD3)
cluster_counts$MPO <- ifelse(
  cluster_counts$Cells == "Immune" & cluster_counts$MPO == "1", 1, 0)
cluster_counts$CD68 <- ifelse(
  cluster_counts$Cells == "Immune" & cluster_counts$CD68 == "1", 1, 0)
cluster_counts$C5aR1 <- ifelse(
  cluster_counts$Cells == "Immune" & cluster_counts$C5aR1 == "1", 1, 0)
cluster_counts$CD11b <- ifelse(
  cluster_counts$Cells == "Immune" & cluster_counts$CD11b == "1", 1, 0)
cluster_counts$CD3 <- ifelse(
  cluster_counts$Cells == "Immune" & cluster_counts$CD3 == "1", 1, 0)
cluster_counts$CD20 <- ifelse(
  cluster_counts$Cells == "Immune" & cluster_counts$CD20 == "1", 1, 0)

# Calculate marker-positive counts per artery
count_df <- cluster_counts %>%
  group_by(Disease, Artery, Patient) %>%
  summarise(
    MPO_Pos   = sum(MPO, na.rm = TRUE), # MPO immune-cell counts
    CD68_Pos  = sum(CD68, na.rm = TRUE), # CD68 immune-cell counts
    C5aR1_Pos = sum(C5aR1, na.rm = TRUE), # C5aR1 immune-cell counts
    CD11b_Pos = sum(CD11b, na.rm = TRUE), # CD11b immune-cell counts
    CD20_Pos  = sum(CD20, na.rm = TRUE), # CD20 immune-cell counts
    CD3_Pos   = sum(CD3, na.rm = TRUE), # CD3 immune-cell counts
    Total_Cells = n(), # Total number of cells in the artery
    .groups = "drop")

# Convert to long format
positive_df <- count_df %>%
  pivot_longer(cols = c(MPO_Pos, CD68_Pos, C5aR1_Pos, CD11b_Pos, CD20_Pos, CD3_Pos),
               names_to = "Marker", values_to = "Positive") %>%
  mutate(Marker = sub("_Pos$","",Marker), # Clean marker names
         Negative = Total_Cells - Positive, # Rest of cells that are not Specific Marker Positive immune cells
         Frequency = Positive / Total_Cells) # Artery-level proportion for visualization

# Set CCE as reference for the beta binomial model
positive_df$Disease <- factor(positive_df$Disease, levels = c("CCE","LN","Diabetes"))

#Order markers for the plot aesthetics
marker_order <- c("MPO", "CD68", "C5aR1", "CD11b", "CD20", "CD3")
positive_df$Marker <- factor(positive_df$Marker, levels = marker_order)

#Make sure Patient of origin is defined as a factor
positive_df$Patient <- factor(positive_df$Patient)

#Beta bnomial mixed-effects model analysis per immue marker
bb_models <- list()
bb_comparisons <- list()


for (marker in levels(positive_df$Marker)) {
  
  # Load data for current immune marker
  tmp <- positive_df %>%
    filter(Marker == marker, 
           !is.na(Positive), !is.na(Negative), !is.na(Disease),!is.na(Patient)
    ) %>% droplevels()
  
  #Fit Beta Binomial using counts
  fit <- glmmTMB(
    cbind(Positive, Negative) ~ Disease + (1 | Patient),
    family = betabinomial(link = "logit"),
    data = tmp)
  # Outcome: Marker+ immune cell counts | all other cell counts
  # Covariates: 
  #  1) Disease -> CCE as reference for pre-specified comparisons
  #  2) (1 | Patient) -> Patient included as a random effect with random intercept
  
  #Save model to check convergence later on
  bb_models[[marker]] <- fit
  
  
  #Get the coeffient contrast p-values from model
  comp_df <- data.frame(contrast = c("CCE vs LN","CCE vs Diabetes"),
                        Marker = marker,
                        p.value = summary(fit)$coefficients$cond[,'Pr(>|z|)'][c("DiseaseLN", "DiseaseDiabetes")])
  
  #Store Coefficient contrast results
  bb_comparisons[[marker]] <- comp_df
}

# Check models' convergences
print(
  data.frame(Marker = names(bb_models),
             pdHess = sapply(bb_models, function(x) x$sdr$pdHess))
)
# Marker pdHess
#    MPO   TRUE
#   CD68   TRUE
#  C5aR1   TRUE
#  CD11b   TRUE
#   CD20   TRUE
#    CD3   TRUE

#FDR p-value adjustment on combined results
bb_results <- bind_rows(bb_comparisons)
bb_results$p.adj <- p.adjust(bb_results$p.value, method = "BH")
bb_results <- bb_results %>%
  mutate(
    p.signif = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001  ~ "***",
      p.adj < 0.01   ~ "**",
      p.adj < 0.05   ~ "*",
      TRUE           ~ "ns"))

#Patient level summaries for plot
patient_means <- positive_df %>%
  group_by(Patient,Disease,Marker) %>%
  summarise(Frequency = sum(Positive) / sum(Total_Cells), .groups = "drop")

#Set disease labels and p.adj y axis positions
comparisons_plot <- bb_results %>% 
  mutate(group1 = "CCE",
         group2 = case_when(contrast == "CCE vs LN"       ~ "LN",
                            contrast == "CCE vs Diabetes" ~ "Diabetes"),
         y.position = case_when(contrast == "CCE vs LN"       ~ 0.21,
                                contrast == "CCE vs Diabetes" ~ 0.23))
# CD68 position p-adj value position
comparisons_plot$y.position[comparisons_plot$Marker == "CD68"] <- c(0.38,0.40)
# C5aR1 position p-adj value position
comparisons_plot$y.position[comparisons_plot$Marker == "C5aR1"] <- c(0.25,0.29)

# Force marker order in all plotting data frames
positive_df$Marker <- factor(positive_df$Marker,levels = marker_order)
patient_means$Marker <- factor(patient_means$Marker,levels = marker_order)
comparisons_plot$Marker <- factor(comparisons_plot$Marker,levels = marker_order)

#Make Figure 4A
p_4A <- ggplot(positive_df, aes(x = Disease, y = Frequency, fill = Marker)) +
  
  # Artery-level plots
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.3) +
  
  # Patient-level plot
  geom_point(data = patient_means, aes(x = Disease, y = Frequency),
             inherit.aes = FALSE,
             position = position_jitter(width = 0.08, height = 0),
             shape = 21, size = 2.2, stroke = 0.7, fill = "white") +
  
  # Beta-binomial GLMM-derived adjusted P values
  stat_pvalue_manual(comparisons_plot,
                     label = "p.signif",
                     xmin = "group1", xmax = "group2",
                     y.position = "y.position",
                     tip.length = 0.01, bracket.size = 0.5, hide.ns = TRUE) +
  #Facet by Marker
  facet_grid(~ Marker, switch = "y") +
  
  #Adjust y axis
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(-0.01,0.40),
                     breaks = seq(0,0.40,by = 0.05)) +
  
  #Set colors
  scale_fill_manual(values = c("C5aR1" = "lightblue",
                               "CD11b" = "green",
                               "MPO"   = "magenta",
                               "CD20"  = "orange",
                               "CD3"   = "red",
                               "CD68"  = "violet")) +
  
  #Set overall theme, add plot title and adjust x and y axis labels 
  theme_classic() +
  labs(title = "Immune Infiltrate in Arteries across Diseases",
       x = "Disease Group",
       y = "Marker+ Immune Cells (% of Arterial Cells)",
       fill = "Marker") +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))
p_4A


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
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] glmmTMB_1.1.14     scales_1.4.0       ggpubr_1.0.0       ggplot2_4.0.3      dplyr_1.2.1        Seurat_5.5.1      
# [7] SeuratObject_5.4.0 sp_2.2-3           tidyr_1.3.2  rstudioapi_0.17.1 