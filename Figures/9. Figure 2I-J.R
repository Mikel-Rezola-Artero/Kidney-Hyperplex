####1. LOAD LIBRARIES####
library(readxl)
library(rstudioapi)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(pROC)
library(RColorBrewer)
library(ggrepel)

####2. IMPORT & CLEAN DATA####

#Load data from script directory
setwd(dirname(getActiveDocumentContext()$path))
data <- read_xlsx("LN_Bulgarie_DBfinal20250527.xlsx")

#Are there any duplicated patients?
table(duplicated(data$Abrev.),duplicated(data$Name))
#yes, 5 patients this is due to multiple visits we want visit 1

#Remove second visit duplicated patients
data <- data[-c(2,12,17,19,34),]


####3. FIGURE 2I####

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

#Modify data frame for graph generation
df_long <- pivot_longer(df, 
                        cols = -C5aR1.den, 
                        names_to = "Marker", 
                        values_to = "Score")

#Order factor for plot aesthetics
df_long$Marker <- factor(df_long$Marker,
                         levels = c("C3c.Cap","C3c.Mes","C3c.Tot.Glom",
                                    "C5b9.Cap","C5b9.Mes","C5b9.Tot.Glom")
                         )

#Calculate multiple comparisons using Wilcox test with Bonferroni Correction
comparisons_df <- compare_means(
  C5aR1.den ~ Score,
  data = df_long,
  group.by = c("Marker"),
  method = "wilcox.test",
  p.adjust.method = "BH"
)

#Set y coordinate for adjusted p.values
comparisons_df$y.position <- 1725

#Remove rows with missing values for scores (avoid "NA levels")
df_long <- df_long %>% 
  filter(!is.na(Score))

#Plot C5aR1 density in glomeruli vs Complement Scores (Low/High)
ggplot(df_long, aes(x = factor(Score), y = C5aR1.den, fill = Marker)) +
  geom_jitter(size = 0.5, width = 0.2, alpha = 0.5) +
  geom_violin(width = 0.7, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.size = 0.2, alpha = 0.7) +
  stat_pvalue_manual(
    comparisons_df,
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    hide.ns = TRUE) +
  facet_wrap(~ Marker, nrow = 2) +
  scale_fill_manual(values = c(
    "C3c.Cap" = "lightblue", "C3c.Mes" = "green", "C5b9.Mes" = "orange",
    "C5b9.Tot.Glom" = "red", "C5b9.Par" = "violet", "C5aR1.den" = "royalblue"
  )) +
  theme_classic() +
  labs(
    title = "Complement Deposits on Glomeruli & C5aR1 Infiltrate",
    x = "Low = 0 / High = 1",
    y = "C5aR1 (cells/mm²)",
    fill = "Marker"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) + theme(legend.position = "none")


####4. FIGURE 2J####

#Store clinical measurements of interest 
C3 <-  as.numeric(data$`C3, g/L (N: 0,75-1,65)`)
C4 <- as.numeric(data$`C4, g/L (N: 0,20-0,65)`)
dsDNA <- as.numeric(data$`anti-dsDNA, U/mL (N<20)`)

#Store clinical endpoint (proteinuria, g/day)
PU <-  factor(ifelse(data$`PU, g/L` > 0.5 , 1,0))

#Add clinical measurements and endpoint to data
data <- cbind(df,C3,C4,dsDNA,PU)

#Transform C5aR1 density with UQ cut-off
data$C5aR1.den <- ifelse(df$C5aR1.den > quantile(df$C5aR1.den,probs = 0.75,
                                                 na.rm = T), 1, 0)

#Set clinical end point & predictor names
y <- "PU" #Proteinuria g/24h
data$y <- data[[y]] #Assign end point data
predictors <- colnames(df) #Assign Predictors: Clinical Measurements + Complement Scores


#Fit univariate logistic regression models & extract ROC data & coef p-values
roc_data <- lapply(predictors, function(var) {
  
  #Fit Logistic Regression
  model <- glm(as.formula(paste("PU ~", var)), data = data, family = binomial)
  #Obtain predicted probabilities of y = 1
  probs <- predict(model, type = "response")
  #Extract ROC data
  roc_obj <- roc(data$y[as.numeric(names(probs))], probs)
  
  #Get coefficient p-value for the predictor variable (from summary)
  p_val <- summary(model)$coefficients[2, 4]
  
  #Calculate and Store False Positive & True Positive Rates / Store p-values
  data.frame(
    fpr = rev(1 - roc_obj$specificities),
    tpr = rev(roc_obj$sensitivities),
    predictor = var,
    auc = as.numeric(auc(roc_obj)),
    p.value = p_val
  )
}) %>% bind_rows()

#Create significance column for the models
roc_data <- roc_data %>%
  mutate(sig = ifelse(p.value < 0.05, "yes", "no"))

#Select significant labels for plotting (at location tpr ~0.75)
label_points <- roc_data %>%
  group_by(predictor) %>%
  arrange(abs(tpr - 0.75)) %>%
  slice(1) %>%
  filter(sig == "yes") %>% #label only significant
  ungroup() %>%
  mutate(
    predictor_label = predictor, #since significant predictors have no "(ns)"
    display_color = predictor_label
  )

#Prepare color palettes
sig_predictors <- unique(roc_data$predictor[roc_data$sig == "yes"])
nonsig_predictors <- unique(roc_data$predictor[roc_data$sig == "no"])

roc_data <- roc_data %>%
  mutate(
    predictor_label = ifelse(sig == "yes", predictor, paste0(predictor, " (ns)")),
    display_color = predictor_label,
    display_linetype = predictor_label  #linetype keys must match
  )

#Colors for significant predictors
if (length(sig_predictors) >= 3 && length(sig_predictors) <= 9) {
  sig_colors <- setNames(brewer.pal(length(sig_predictors), "Set1"), sig_predictors)
} else {
  sig_colors <- setNames(hue_pal()(length(sig_predictors)), sig_predictors)
}

#Gray shades for non-significant
nonsig_colors <- setNames(scales::grey_pal()(length(nonsig_predictors)),
                          paste0(nonsig_predictors, " (ns)"))

color_palette <- c(sig_colors, nonsig_colors)

#Set line types to help distinguishing between predictors
linetype_palette <- c(
  setNames(rep("solid", length(sig_predictors)), sig_predictors),
  setNames(
    rep_len(c("dashed", "dotted", "dotdash", "twodash"), length(nonsig_predictors)),
    paste0(nonsig_predictors, " (ns)")
  )
)

#Plot the ROC curves
p1 <- ggplot(roc_data, aes(x = fpr, y = tpr, group = predictor_label)) +
  geom_line(aes(color = display_color, linetype = display_linetype), linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_label_repel(
    data = label_points,
    aes(label = paste0(predictor, "\nAUC = ", round(auc, 3))),
    size = 4,
    show.legend = FALSE,
    segment.color = "gray50",
    segment.size = 0.5,
    box.padding = 0.5,
    point.padding = 0.25,
    nudge_y = 0.03
  ) +
  scale_color_manual(values = color_palette, name = "Predictor") +
  scale_linetype_manual(values = linetype_palette, name = "Predictor") +
  labs(
    title = "ROC Curves for Univariate Logistic Models",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.box = "vertical",
    legend.position = "right"
  )

p1

#Test if combining with infiltrates increases predicting power
#Create combined labels
data$C5b9.Cap_C5aR1.den <- interaction(data$C5b9.Cap,data$C5aR1.den)
levels(data$C5b9.Cap_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")
data$C5b9.Tot.Glom_C5aR1.den <- interaction(data$C5b9.Tot.Glom,data$C5aR1.den)
levels(data$C5b9.Tot.Glom_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")
data$C5b9.Mes_C5aR1.den <- interaction(data$C5b9.Mes,data$C5aR1.den)
levels(data$C5b9.Mes_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")
data$C3c.Cap_C5aR1.den <- interaction(data$C3c.Cap,data$C5aR1.den)
levels(data$C3c.Cap_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")
data$C3c.Mes_C5aR1.den <- interaction(data$C3c.Mes,data$C5aR1.den)
levels(data$C3c.Mes_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")
data$C3c.Tot.Glom_C5aR1.den <- interaction(data$C3c.Tot.Glom,data$C5aR1.den)
levels(data$C3c.Tot.Glom_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")

predictors <- c("C3c.Cap_C5aR1.den",
                "C3c.Mes_C5aR1.den",
                "C3c.Tot.Glom_C5aR1.den",
                "C5b9.Mes_C5aR1.den",
                "C5b9.Tot.Glom_C5aR1.den",
                "C5b9.Cap_C5aR1.den",
                "C3","C4","dsDNA")  # Your variables



#Fit univariate logistic regression models & extract ROC data & coef p-values
roc_data <- lapply(predictors, function(var) {
  
  #Fit Logistic Regression
  model <- glm(as.formula(paste("y ~", var)), data = data, family = binomial)
  #Obtain predicted probabilities of y = 1
  probs <- predict(model, type = "response")
  #Extract ROC data
  roc_obj <- roc(data$y[as.numeric(names(probs))], probs)
  
  # Get coefficient p-value for the predictor variable (from summary)
  p_val <- summary(model)$coefficients[2, 4]
  
  #Calculate and Store False Positive & True Positive Rates / Store p-values
  data.frame(
    fpr = rev(1 - roc_obj$specificities),
    tpr = rev(roc_obj$sensitivities),
    predictor = var,
    auc = as.numeric(auc(roc_obj)),
    p.value = p_val
  )
}) %>% bind_rows()

#Create significance column for the models
roc_data <- roc_data %>%
  mutate(sig = ifelse(p.value < 0.05, "yes", "no"))

#Select significant labels for plotting (at location tpr ~0.75)
label_points <- roc_data %>%
  group_by(predictor) %>%
  arrange(abs(tpr - 0.75)) %>%
  slice(1) %>%
  filter(sig == "yes") %>% # label only significant
  ungroup() %>%
  mutate(
    predictor_label = predictor, # since significant predictors have no "(ns)"
    display_color = predictor_label
  )

#Prepare color palettes
sig_predictors <- unique(roc_data$predictor[roc_data$sig == "yes"])
nonsig_predictors <- unique(roc_data$predictor[roc_data$sig == "no"])

roc_data <- roc_data %>%
  mutate(
    predictor_label = ifelse(sig == "yes", predictor, paste0(predictor, " (ns)")),
    display_color = predictor_label,
    display_linetype = predictor_label  # linetype keys must match
  )

# Colors for significant predictors
if (length(sig_predictors) >= 3 && length(sig_predictors) <= 9) {
  sig_colors <- setNames(brewer.pal(length(sig_predictors), "Set1"), sig_predictors)
} else {
  sig_colors <- setNames(hue_pal()(length(sig_predictors)), sig_predictors)
}

# Gray shades for non-significant
nonsig_colors <- setNames(scales::grey_pal()(length(nonsig_predictors)),
                          paste0(nonsig_predictors, " (ns)"))

color_palette <- c(sig_colors, nonsig_colors)

#Set line types to help distinguishing between predictors
linetype_palette <- c(
  setNames(rep("solid", length(sig_predictors)), sig_predictors),
  setNames(
    rep_len(c("dashed", "dotted", "dotdash", "twodash"), length(nonsig_predictors)),
    paste0(nonsig_predictors, " (ns)")
  )
)

#Plot the ROC curves
p2 <- ggplot(roc_data, aes(x = fpr, y = tpr, group = predictor_label)) +
  geom_line(aes(color = display_color, linetype = display_linetype), linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_label_repel(
    data = label_points,
    aes(label = paste0(predictor, "\nAUC = ", round(auc, 3))),
    size = 4,
    show.legend = FALSE,
    segment.color = "gray50",
    segment.size = 0.5,
    box.padding = 0.5,
    point.padding = 0.25,
    nudge_y = 0.03
  ) +
  scale_color_manual(values = color_palette, name = "Predictor") +
  scale_linetype_manual(values = linetype_palette, name = "Predictor") +
  labs(
    title = "ROC Curves for Univariate Logistic Models",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.box = "vertical",
    legend.position = "right"
  )

p2

#FIGURE 2J
p1 | p2



####5. EXTRA OR and p-values for Logistic Regressions####

# Run Univariate Logistic Regressions in non-combined predictors
dependent <- "PU"  # Replace with your binary outcome
explanatory <- colnames(df)  # Your variables

# Prepare storage
results <- data.frame(
  Predictor = character(),
  Estimate = numeric(),
  OR = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(explanatory)) {
  pred <- explanatory[i]
  
  # Fit univariate logistic regression
  model <- glm(as.formula(paste("PU ~", pred)), data = data, family = binomial)
  summ <- summary(model)
  
  # Extract coefficient for the predictor (usually 2nd row)
  coef_est <- summ$coefficients[2, "Estimate"]
  p_val <- summ$coefficients[2, "Pr(>|z|)"]
  or <- exp(coef_est)
  
  # Append results
  results <- rbind(results, data.frame(
    Predictor = pred,
    Estimate = coef_est,
    OR = or,
    P_value = p_val,
    stringsAsFactors = FALSE
  ))
}

# Round numeric columns nicely
results$Estimate <- round(results$Estimate, 3)
results$OR <- round(results$OR, 3)
results$P_value <- signif(results$P_value, 3)

print(results)#Confirmation that only C5b-9 Capillary, Mesangial and Total are significant
#        Predictor Estimate  OR   P_value
# 1       C3c.Cap    0.474 1.607  0.5120
# 2       C3c.Mes   -0.788 0.455  0.3820
# 3  C3c.Tot.Glom   -0.223 0.800  0.7560
# 4      C5b9.Cap    1.572 4.815  0.0462
# 5      C5b9.Mes    1.572 4.815  0.0462
# 6 C5b9.Tot.Glom    1.573 4.821  0.0373
# 7     C5aR1.den    0.722 2.059  0.4200

#Test if combining with infiltrates increases predicting power
#Create combined labels
data$C5b9.Cap_C5aR1.den <- interaction(data$C5b9.Cap,data$C5aR1.den)
levels(data$C5b9.Cap_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")
data$C5b9.Tot.Glom_C5aR1.den <- interaction(data$C5b9.Tot.Glom,data$C5aR1.den)
levels(data$C5b9.Tot.Glom_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")
data$C5b9.Mes_C5aR1.den <- interaction(data$C5b9.Mes,data$C5aR1.den)
levels(data$C5b9.Mes_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")
data$C3c.Cap_C5aR1.den <- interaction(data$C3c.Cap,data$C5aR1.den)
levels(data$C3c.Cap_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")
data$C3c.Mes_C5aR1.den <- interaction(data$C3c.Mes,data$C5aR1.den)
levels(data$C3c.Mes_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")

explanatory <- c("C3c.Mes_C5aR1.den",
                 "C3c.Cap_C5aR1.den","C5b9.Mes_C5aR1.den",
                 "C5b9.Tot.Glom_C5aR1.den","C5b9.Cap_C5aR1.den")  # Your variables


results <- data.frame(
  Predictor = character(),
  Level = character(),
  Estimate = numeric(),
  OR = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (pred in explanatory) {
  model <- glm(as.formula(paste("PU ~", pred)), data = data, family = binomial)
  summ <- summary(model)
  
  # Extract all coefficients except intercept
  coefs <- summ$coefficients[-1, , drop = FALSE]   # drop intercept row
  
  # Extract levels from row names like 'predInt' and 'predHigh'
  levels <- gsub(paste0("^", pred), "", rownames(coefs))
  
  # Append each level info
  for (j in seq_len(nrow(coefs))) {
    coef_est <- coefs[j, "Estimate"]
    p_val <- coefs[j, "Pr(>|z|)"]
    or <- exp(coef_est)
    
    results <- rbind(results, data.frame(
      Predictor = pred,
      Level = levels[j],
      Estimate = round(coef_est, 3),
      OR = round(or, 3),
      P_value = signif(p_val, 3),
      stringsAsFactors = FALSE
    ))
  }
}

print(results)#Confirmation that only C5b-9 Capillary, Mesangial and Total + C5aR1 are significant
#               Predictor    Level    Estimate  OR  P_value
# 1       C3c.Mes_C5aR1.den Int/High    0.000 1.000  1.0000
# 2       C3c.Cap_C5aR1.den Int/High    0.405 1.500  0.6000
# 3      C5b9.Mes_C5aR1.den Int/High    1.584 4.875  0.0488
# 4 C5b9.Tot.Glom_C5aR1.den Int/High    2.015 7.500  0.0144
# 5      C5b9.Cap_C5aR1.den Int/High    1.584 4.875  0.0488


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
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
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
#   [1] ggrepel_0.9.6      RColorBrewer_1.1-3 pROC_1.18.5        ggpubr_0.6.0       ggplot2_3.5.2     
# [6] tidyr_1.3.1        dplyr_1.1.4        readxl_1.4.5       rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6      compiler_4.3.2    ggsignif_0.6.4    Rcpp_1.0.14       tidyselect_1.2.1  dichromat_2.0-0.1
# [7] scales_1.4.0      R6_2.6.1          plyr_1.8.9        labeling_0.4.3    generics_0.1.3    Formula_1.2-5    
# [13] backports_1.5.0   tibble_3.2.1      car_3.1-3         pillar_1.10.2     rlang_1.1.6       broom_1.0.8      
# [19] cli_3.6.5         withr_3.0.2       magrittr_2.0.3    grid_4.3.2        lifecycle_1.0.4   vctrs_0.6.5      
# [25] rstatix_0.7.2     glue_1.8.0        farver_2.1.2      cellranger_1.1.0  abind_1.4-8       carData_3.0-5    
# [31] purrr_1.0.4       tools_4.3.2       pkgconfig_2.0.3  
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
#   [1] ggrepel_0.9.6      RColorBrewer_1.1-3 pROC_1.18.5       
# [4] ggpubr_0.6.0       ggplot2_3.5.2      tidyr_1.3.1       
# [7] dplyr_1.1.4        readxl_1.4.5       rstudioapi_0.17.1 
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6      compiler_4.3.2    ggsignif_0.6.4   
# [4] Rcpp_1.0.14       tidyselect_1.2.1  dichromat_2.0-0.1
# [7] scales_1.4.0      R6_2.6.1          plyr_1.8.9       
# [10] labeling_0.4.3    generics_0.1.3    Formula_1.2-5    
# [13] backports_1.5.0   tibble_3.2.1      car_3.1-3        
# [16] pillar_1.10.2     rlang_1.1.6       broom_1.0.8      
# [19] cli_3.6.5         withr_3.0.2       magrittr_2.0.3   
# [22] grid_4.3.2        lifecycle_1.0.4   vctrs_0.6.5      
# [25] rstatix_0.7.2     glue_1.8.0        farver_2.1.2     
# [28] cellranger_1.1.0  abind_1.4-8       carData_3.0-5    
# [31] purrr_1.0.4       tools_4.3.2       pkgconfig_2.0.3 