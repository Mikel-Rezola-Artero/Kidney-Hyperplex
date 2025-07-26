####1. LOAD LIBRARIES####
library(readxl)
library(rstudioapi)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(scales)

####2. IMPORT & CLEAN DATA####

#Load data from script directory
setwd(dirname(getActiveDocumentContext()$path))
data <- read_xlsx("LN_Bulgarie_DBfinal20250527.xlsx")

#Are there any duplicated patients?
table(duplicated(data$Abrev.),duplicated(data$Name))
#yes, 5 patients this is due to multiple visits we want visit 1

#Remove second visit duplicated patients
data <- data[-c(2,12,17,19,34),]


####3. FIGURE S9A-D####

#Take glomeruli scoring information and rename columns
df <- data[,c("C3c_parietal_Score",
              "C3c_MesangialScore",
              "C3c_total_score",
              "C5b9_parietal_Score",
              "C5b9_Mesangial_Score",
              "C5b9_totalGlomerular_score",
              "C5aR1_glomeruli (cells/mm2)")]
colnames(df) <- c("C3c.Cap",
                  "C3c.Mes",
                  "C3c.Tot.Glom",
                  "C5b9.Cap",
                  "C5b9.Mes",
                  "C5b9.Tot.Glom",
                  "C5aR1.den")

#Transform scores into factors
df <- df %>%
  mutate(across(.cols = -C5aR1.den, .fns = as.factor))

#Group score levels 
levels(df$C3c.Cap) <- c(0,1,1,1)
levels(df$C3c.Mes) <- c(0,1,1)
levels(df$C3c.Tot.Glom) <- c(0,1,1,1,1)
levels(df$C5b9.Mes) <- c(0,0,1,1)
levels(df$C5b9.Tot.Glom) <- c(0,0,0,1,1,1,1)
levels(df$C5b9.Cap) <- c(0,0,1,1)

#Store clinical measurements of interest 
C3 <-  as.numeric(data$`C3, g/L (N: 0,75-1,65)`)
C4 <- as.numeric(data$`C4, g/L (N: 0,20-0,65)`)
dsDNA <- as.numeric(data$`anti-dsDNA, U/mL (N<20)`)
eGFR <- as.numeric(data$`eGFR, mL/min/1,73 sqm`)

#Add clinical measurements to df
df <- cbind(df,C3,C4,dsDNA,eGFR)

#Modify data frame for graph generation
df_long <- pivot_longer(df, 
                        cols = c(-C5aR1.den,-C3,-C4,-dsDNA,-eGFR), 
                        names_to = "Marker", 
                        values_to = "Score")

#Order factor for plot aesthetics
df_long$Marker <- factor(df_long$Marker,levels = c("C3c.Cap","C3c.Mes","C3c.Tot.Glom",
                                                   "C5b9.Cap","C5b9.Mes","C5b9.Tot.Glom"))

#Calculate multiple comparisons using Wilcox test with Bonferroni Correction
comparisons_df <-compare_means(
  C3 ~ Score,
  data = df_long,
  group.by = c("Marker"),
  method = "wilcox.test",
  p.adjust.method = "BH"
)

#Set y coordinate for p.values
comparisons_df$y.position <- 1.9

#Remove scores with missing values
df_long <- df_long %>% 
  filter(!is.na(Score))

#Plot C3 concentration vs Complement Scores (Low/High)
C3.Plot <- ggplot(df_long, aes(x = factor(Score), y = C3, fill = Marker)) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.5) +
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.size = 0.2, alpha = 0.7) +
  stat_pvalue_manual(
    comparisons_df,
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    hide.ns = F) +
  facet_wrap(~ Marker, nrow = 2) +
  scale_fill_manual(values = c(
    "C3c.Cap" = "lightblue", "C3c.Mes" = "green", "C5b9.Mes" = "orange",
    "C5b9.Tot.Glom" = "red", "C5b9.Cap" = "violet", "C5aR1.den" = "royalblue"
  )) +
  theme_classic() +
  labs(
    title = "Complement Deposits on Glomeruli & C3 Plasma",
    x = "Low = 0 / High = 1",
    y = "C3 Plasma (g/L)",
    fill = "Marker"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) + theme(legend.position = "none")



#Calculate multiple comparisons using Wilcox test with Bonferroni Correction
comparisons_df <-compare_means(
  C4 ~ Score,
  data = df_long,
  group.by = c("Marker"),
  method = "wilcox.test",
  p.adjust.method = "BH"
)

#Set y coordinate for p.values
comparisons_df$y.position <- 0.55

#Remove scores with missing values
df_long <- df_long %>% 
  filter(!is.na(Score))

#Plot C4 concentration vs Complement Scores (Low/High)
C4.Plot <- ggplot(df_long, aes(x = factor(Score), y = C4, fill = Marker)) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.5) +
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.size = 0.2, alpha = 0.7) +
  stat_pvalue_manual(
    comparisons_df,
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    hide.ns = F) +
  facet_wrap(~ Marker, nrow = 2) +
  scale_fill_manual(values = c(
    "C3c.Cap" = "lightblue", "C3c.Mes" = "green", "C5b9.Mes" = "orange",
    "C5b9.Tot.Glom" = "red", "C5b9.Cap" = "violet", "C5aR1.den" = "royalblue"
  )) +
  theme_classic() +
  labs(
    title = "Complement Deposits on Glomeruli & C4 Plasma",
    x = "Low = 0 / High = 1",
    y = "C4 Plasma (g/L)",
    fill = "Marker"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) + theme(legend.position = "none")

#Calculate multiple comparisons using Wilcox test with Bonferroni Correction
comparisons_df <-compare_means(
  dsDNA ~ Score,
  data = df_long,
  group.by = c("Marker"),
  method = "wilcox.test",
  p.adjust.method = "BH"
)

#Set y coordinate for p.values
comparisons_df$y.position <- 204

#Remove scores with missing values
df_long <- df_long %>% 
  filter(!is.na(Score))

#Plot dsDNA concentration vs Complement Scores (Low/High)
dsDNA.Plot <- ggplot(df_long, aes(x = factor(Score), y = dsDNA, fill = Marker)) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.5) +
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.size = 0.2, alpha = 0.7) +
  stat_pvalue_manual(
    comparisons_df,
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    hide.ns = F) +
  facet_wrap(~ Marker, nrow = 2) +
  scale_fill_manual(values = c(
    "C3c.Cap" = "lightblue", "C3c.Mes" = "green", "C5b9.Mes" = "orange",
    "C5b9.Tot.Glom" = "red", "C5b9.Cap" = "violet", "C5aR1.den" = "royalblue"
  )) +
  theme_classic() +
  labs(
    title = "Complement Deposits on Glomeruli & Anti-dsDNA Abs",
    x = "Low = 0 / High = 1",
    y = "Anti-dsDNA Abs (U/ml)",
    fill = "Marker"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) + theme(legend.position = "none")


#Calculate multiple comparisons using Wilcox test with Bonferroni Correction
comparisons_df <-compare_means(
  eGFR ~ Score,
  data = df_long,
  group.by = c("Marker"),
  method = "wilcox.test",
  p.adjust.method = "BH"
)

#Set y coordinate for p.values
comparisons_df$y.position <- 146

#Remove scores with missing values
df_long <- df_long %>% 
  filter(!is.na(Score))

#Plot dsDNA concentration vs Complement Scores (Low/High)
eGFR.Plot <- ggplot(df_long, aes(x = factor(Score), y = eGFR, fill = Marker)) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.5) +
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.size = 0.2, alpha = 0.7) +
  stat_pvalue_manual(
    comparisons_df,
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    hide.ns = F) +
  facet_wrap(~ Marker, nrow = 2) +
  scale_fill_manual(values = c(
    "C3c.Cap" = "lightblue", "C3c.Mes" = "green", "C5b9.Mes" = "orange",
    "C5b9.Tot.Glom" = "red", "C5b9.Cap" = "violet", "C5aR1.den" = "royalblue"
  )) +
  theme_classic() +
  labs(
    title = "Complement Deposits on Glomeruli & eGFR",
    x = "Low = 0 / High = 1",
    y = "estimated Glomerular Filtration Rate",
    fill = "Marker"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) + theme(legend.position = "none")

#PLOT THE DATA (FIGURE S9A-D)
(C3.Plot | C4.Plot) / (dsDNA.Plot | eGFR.Plot)


####4. FIGURE S9E-G####

# Variables to correlate with C5aR1 infiltrate
variables <- c("C3", "C4", "dsDNA")

# Create correlation plots
custom_breaks <- c(0, 25, 50, 75, 125, 250, 500, 1500, 2000)
x_labels <- c("C3 Plasma", "C4 Plasma","Anti-dsDNA Abs")

plots <- lapply(seq_along(variables), function(i) {
  ggscatter(df, x = variables[i], y = "C5aR1.den",
            add = "reg.line"#"loess"
            , conf.int = TRUE,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = x_labels[i], ylab = "C5aR1.den") +
    scale_x_continuous(trans = pseudo_log_trans(base = 10)) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10),
                       breaks = custom_breaks)+
    theme(
      axis.text.x = element_text(size = 8),  
      axis.text.y = element_text(size = 8)   
    )
})

#Arrange plots (FIGURE S9E-G)
ggarrange(plotlist = plots, ncol = 3, nrow = 1)


####5. Session Information#####

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
#   [1] scales_1.4.0      patchwork_1.3.0   ggpubr_0.6.0      ggplot2_3.5.2     tidyr_1.3.1       dplyr_1.1.4      
# [7] rstudioapi_0.17.1 readxl_1.4.5     
# 
# loaded via a namespace (and not attached):
#   [1] Matrix_1.6-5       gtable_0.3.6       compiler_4.3.2     ggsignif_0.6.4     tidyselect_1.2.1  
# [6] dichromat_2.0-0.1  splines_4.3.2      lattice_0.22-7     R6_2.6.1           labeling_0.4.3    
# [11] generics_0.1.3     Formula_1.2-5      backports_1.5.0    tibble_3.2.1       car_3.1-3         
# [16] pillar_1.10.2      RColorBrewer_1.1-3 rlang_1.1.6        broom_1.0.8        cli_3.6.5         
# [21] mgcv_1.9-3         withr_3.0.2        magrittr_2.0.3     grid_4.3.2         cowplot_1.1.3     
# [26] nlme_3.1-168       lifecycle_1.0.4    vctrs_0.6.5        rstatix_0.7.2      glue_1.8.0        
# [31] farver_2.1.2       cellranger_1.1.0   abind_1.4-8        carData_3.0-5      purrr_1.0.4       
# [36] tools_4.3.2        pkgconfig_2.0.3   
# > sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3;  LAPACK version 3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Paris
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods  
# [7] base     
# 
# other attached packages:
#   [1] scales_1.4.0      patchwork_1.3.0   ggpubr_0.6.0     
# [4] ggplot2_3.5.2     tidyr_1.3.1       dplyr_1.1.4      
# [7] rstudioapi_0.17.1 readxl_1.4.5     
# 
# loaded via a namespace (and not attached):
#   [1] Matrix_1.6-5       gtable_0.3.6       compiler_4.3.2    
# [4] ggsignif_0.6.4     tidyselect_1.2.1   dichromat_2.0-0.1 
# [7] splines_4.3.2      lattice_0.22-7     R6_2.6.1          
# [10] labeling_0.4.3     generics_0.1.3     Formula_1.2-5     
# [13] backports_1.5.0    tibble_3.2.1       car_3.1-3         
# [16] pillar_1.10.2      RColorBrewer_1.1-3 rlang_1.1.6       
# [19] broom_1.0.8        cli_3.6.5          mgcv_1.9-3        
# [22] withr_3.0.2        magrittr_2.0.3     grid_4.3.2        
# [25] cowplot_1.1.3      nlme_3.1-168       lifecycle_1.0.4   
# [28] vctrs_0.6.5        rstatix_0.7.2      glue_1.8.0        
# [31] farver_2.1.2       cellranger_1.1.0   abind_1.4-8       
# [34] carData_3.0-5      purrr_1.0.4        tools_4.3.2       
# [37] pkgconfig_2.0.3  
