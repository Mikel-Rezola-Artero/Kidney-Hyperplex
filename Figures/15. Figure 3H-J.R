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
library(glmmTMB)


####2. LOAD ARTERY DATA####

#Load Processed data
setwd(dirname(getActiveDocumentContext()$path))
Seurat.obj <- readRDS("../Data2/Processed_Artery_v2.rds")


####3. FIGURE 3H####

#VCAM-1 in Arteries across diseases (FIGURE 3H)
DotPlot(Seurat.obj,
	group.by = "Disease",
        features = c("VCAM.1.Pos"),
	cols = c("lightgray","blue","lightblue"),
        scale = FALSE) + 
  labs(title = "VCAM-1 positivity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )


####4. FIGURE 3I: BETA-BINOMIAL MIXED-EFFECTS MODEL####

# Extract data
cluster_counts <- data.frame(
  Cells   = Idents(Seurat.obj),
  VCAM.1  = Seurat.obj$VCAM.1.Pos,
  Artery  = Seurat.obj$Artery,
  Patient = Seurat.obj$slide,
  Disease = Seurat.obj$Disease
)

# Define endothelial VCAM-1 positivity
cluster_counts$VCAM.1 <- ifelse(
  cluster_counts$Cells == "Endothelial" & cluster_counts$VCAM.1 == "1",
  1, # 1 =  cell is both endothelial AND VCAM-1+
  0 # 0 = every other cell in the artery
)

# Calculate VCAM-1+ endothelial-cell counts and total cell counts (minus VCAM-1+ endothelial-cell counts) per artery
count_df <- cluster_counts %>%
  
  group_by(Disease, Artery, Patient) %>%
  
  summarise(
    VCAM_Pos = sum(VCAM.1, na.rm = TRUE), # Number of VCAM-1+ endothelial cells
    Total_Cells = n(), # Total number of cells in the artery
    .groups = "drop" #remove grouping
  ) %>%
  
  mutate(
    Rest_of_Cells = Total_Cells - VCAM_Pos, # Number of cells that are NOT VCAM-1+ endothelial cells
    Frequency = VCAM_Pos / Total_Cells, # Proportion of all arterial cells that are VCAM-1+ endothelial cells for plot
    Marker = "VCAM-1")

# Set CCE as the reference group
count_df$Disease <- factor(count_df$Disease,
                           levels = c("CCE","LN","Diabetes"))

#Make sure Patient of origin is defined as factor
count_df$Patient <- factor(count_df$Patient)

# Fit beta-binomial mixed-effects model to control for pseudo replication
fit_vcam <- glmmTMB(
  cbind(VCAM_Pos, Rest_of_Cells) ~ Disease + (1 | Patient),
  family = betabinomial(link = "logit"), data = count_df)
# Outcome: VCAM-1+ endothelial cell counts | all other cell counts
# Covariates: 
#  1) Disease -> CCE as reference for pre-specified comparisons
#  2) (1 | Patient) -> Patient included as a random effect with random intercept

# Model summary
summary(fit_vcam)
# Conditional model:
#               Estimate Std.    Error       z      value Pr(>|z|)    
#  (Intercept)      -2.3844     0.2211 -10.786    < 2e-16   ***
#  DiseaseLN        -3.2861     0.4396  -7.475   7.71e-14   ***
#  DiseaseDiabetes  -2.1277     0.3684  -5.776   7.67e-09   ***
fit_vcam$sdr$pdHess # Check model convergence
#TRUE

#Get the coeffient contrast p-values from model
comparisons_vcam <- data.frame(contrast = c("CCE vs LN","CCE vs Diabetes"),
                               p.value = summary(fit_vcam)$coefficients$cond[,'Pr(>|z|)'][c("DiseaseLN", "DiseaseDiabetes")])

#FDR p-value adjustment
comparisons_vcam$p.adj <- p.adjust(comparisons_vcam$p.value, method = "BH")

#Add significance labels for plot
comparisons_vcam <- comparisons_vcam %>%
  mutate(plot_label = case_when(p.adj < 0.0001 ~ "****",
                                p.adj < 0.001  ~ "***",
                                p.adj < 0.01   ~ "**",
                                p.adj < 0.05   ~ "*",
                                TRUE ~ paste0("P = ",format.pval(p.adj, digits = 3,eps = 0.001))))

# Patient-level summaries for visualization
# Calculate per patient: VCAM-1+ endothelial cells / total arterial cells
patient_means <- count_df %>% 
  group_by(Patient,Disease,Marker) %>%
  summarise(Frequency = sum(VCAM_Pos) / sum(Total_Cells),
            .groups = "drop")

# Prepare mixed-model adjusted P values for plotting (coordinate and group)
comparisons_plot <- comparisons_vcam %>%
  mutate(
    group1 = "CCE",
    group2 = c("LN", "Diabetes"),
    y.position = c(0.60, 0.65))

# Plot FIGURE 3I
p_3I <- ggplot(count_df, aes(x = Disease, y = Frequency, fill = Disease)) +
  
  # Individual arteries plots
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.3) +
  
  # Patient-level dotplot
  geom_point(data = patient_means,
             aes(x = Disease, y = Frequency),
             inherit.aes = FALSE, 
             position = position_jitter(width = 0.08, height = 0),
             shape = 21, size = 2.2, stroke = 0.7, fill = "white") +
  
  # Add beta-binomial mixed-model P values
  stat_pvalue_manual(comparisons_plot, label = "plot_label",
                     xmin = "group1", xmax = "group2",
                     y.position = "y.position", 
                     tip.length = 0.01,
                     bracket.size = 0.5,
                     hide.ns = FALSE) +
  
  #Adjust y axis
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(-0.01, 0.68),
                     breaks = c(0, 0.2, 0.4, 0.6)) +
  
  #Set colors
  scale_fill_manual(values = c("LN"       = "red",
                               "Diabetes" = "green",
                               "CCE"      = "magenta")) +
  
  #Set overall theme, add plot title and adjust x and y axis labels 
  theme_classic() +
  labs(title = "VCAM-1+ Endothelial Cells in Arteries across Diseases",
       x = "Disease Group",
       y = "VCAM-1+ Endothelial Cells (% of Arterial Cells)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

#Plot
p_3I

####5. FIGURE 3J####

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
testRes <- as.matrix(forceSymmetric(testRes$p,uplo = "U"))
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
       posCol = c("brown"),negCol =c("blue"))


####6. Session Information####

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
# [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8   
# [6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Paris
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base  
#
# other attached packages:
# [1] glmmTMB_1.1.14     scales_1.4.0       ggpubr_1.0.0       ggplot2_4.0.3      dplyr_1.2.1        Seurat_5.5.1       SeuratObject_5.4.0 sp_2.2-3   
# [8] qgraph_1.9.8       Matrix_1.6-5       psych_2.5.3      tidyr_1.3.1      rstudioapi_0.17.1 