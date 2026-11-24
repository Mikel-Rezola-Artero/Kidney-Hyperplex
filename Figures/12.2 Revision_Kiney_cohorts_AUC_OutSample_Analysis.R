####1. LOAD LIBRARIES & SET WORKING DIRECTORY####
library(readxl)
library(rstudioapi)
library(dplyr)
library(pROC)
library(ggplot2)
library(ggrepel)

# Set working directory
setwd(dirname(getActiveDocumentContext()$path))

####2. Logistic Regression Fit on Exploratory Cohort####

#Load data from script directory
data <- read_xlsx("LN_Bulgarie_DBfinal20250527.xlsx")

#Are there any duplicated patients? (anonymised names)
table(duplicated(data$Abrev.),duplicated(data$Name))
#yes, 5 patients this is due to multiple visits we want visit 1

#Remove second visit duplicated patients
data <- data[-c(2,12,17,19,34),]

#Rename columns
C_markers <- c("C3c.Cap", "C3c.Mes", "C3c.Tot.Glom",
               "C5b9.Cap", "C5b9.Mes", "C5b9.Tot.Glom",
               "C5aR1.den")
wanted <- c("C3c_parietal_Score", "C3c_MesangialScore", "C3c_total_score",
            "C5b9_parietal_Score", "C5b9_Mesangial_Score", "C5b9_totalGlomerular_score",
            "C5aR1_glomeruli (cells/mm2)")
colnames(data)[match(wanted, colnames(data))] <- C_markers

#Transform C5aR1 density with UQ cut-off
data$C5aR1.den <- ifelse(data$C5aR1.den > quantile(data$C5aR1.den,probs = 0.75,
                                                   na.rm = T), 1, 0)
#Transform scores into factors
for (i in C_markers) {
  data[,colnames(data) == i] <- as.factor(unlist(data[,which(colnames(data) == i)]))
}

#Group score levels 
levels(data$C3c.Cap) <- c(0,1,1,1)
levels(data$C3c.Mes) <- c(0,1,1)
levels(data$C3c.Tot.Glom) <- c(0,1,1,1,1)
levels(data$C5b9.Mes) <- c(0,0,1,1)
levels(data$C5b9.Tot.Glom) <- c(0,0,0,1,1,1,1)
levels(data$C5b9.Cap) <- c(0,0,1,1)

#Add clinical endpoint (proteinuria, g/day) to data
data$PU <- factor(ifelse(data$`PU, g/L` > 0.5 , 1,0))

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

# Run Univariate Logistic Regressions in non-combined predictors
explanatory <- c(C_markers,paste0(C_markers[-7],"_",C_markers[7]))

# Prepare storage
results <- data.frame(
  Predictor = character(),
  Estimate = numeric(),
  OR = numeric(),
  CI_low = numeric(),
  CI_high = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(explanatory)) {
  pred <- explanatory[i]
  
  # Fit univariate logistic regression
  model <- glm(
    as.formula(paste("PU ~", pred)),
    data = data,
    family = binomial
  )
  
  summ <- summary(model)
  
  # Extract coefficient for predictor
  coef_est <- summ$coefficients[2, "Estimate"]
  p_val <- summ$coefficients[2, "Pr(>|z|)"]
  
  # Odds ratio
  or <- exp(coef_est)
  
  # 95% CI for OR from fitted model
  ci <- exp(confint(model)[2, ])
  
  # Append results
  results <- rbind(
    results,
    data.frame(
      Predictor = pred,
      Estimate = coef_est,
      OR = or,
      CI_low = ci[1],
      CI_high = ci[2],
      P_value = p_val,
      stringsAsFactors = FALSE
    )
  )
}

# Round numeric columns
results$Estimate <- round(results$Estimate, 3)
results$OR <- round(results$OR, 3)
results$CI_low <- round(results$CI_low, 3)
results$CI_high <- round(results$CI_high, 3)
results$P_value <- signif(results$P_value, 3)

print(results)
# Predictor            Estimate     OR  CI_low  CI_high  P_value
# C3c.Cap                 0.474  1.607   0.377    6.713   0.5120
# C3c.Mes                -0.788  0.455   0.072    2.829   0.3820
# C3c.Tot.Glom           -0.223  0.800   0.180    3.176   0.7560
# C5b9.Cap                1.572  4.815   1.113   26.325   0.0462
# C5b9.Mes                1.572  4.815   1.113   26.325   0.0462
# C5b9.Tot.Glom           1.573  4.821   1.158   23.394   0.0373
# C5aR1.den               0.722  2.059   0.400   15.666   0.4200
# C5b9.Cap_C5aR1.den      1.584  4.875   1.090   27.367   0.0488
# C5b9.Tot.Glom_C5aR1.den 2.015  7.500   1.629   43.911   0.0144
# C5b9.Mes_C5aR1.den      1.584  4.875   1.090   27.367   0.0488
# C3c.Cap_C5aR1.den       0.405  1.500   0.311    6.846   0.6000
# C3c.Mes_C5aR1.den       0.000  1.000   0.245    4.247   1.0000
# C3c.Tot.Glom_C5aR1.den -0.511  0.600   0.078    3.201   0.5730

####3. Use Exploration Cohort Logistic Regression Fit on Validation Cohort####

#Load data from script directory
data2 <- read_xlsx("../Data2/Matrice_C9C5aR_ROC_validation.xlsx")

#Rename columns
colnames(data2) <- c("Patients",
                     "C5b9.Tot.Glom",
                     "C5b9.Cap",
                     "C5b9.Mes",
                     "C5aR1.den",
                     "PU")

#Transform scores into factors
for (i in colnames(data2)) {
  data2[,colnames(data2) == i] <- as.factor(unlist(data2[,which(colnames(data2) == i)]))
}

#Create combined labels
data2$C5b9.Cap_C5aR1.den <- interaction(data2$C5b9.Cap,data2$C5aR1.den)
levels(data2$C5b9.Cap_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")
data2$C5b9.Tot.Glom_C5aR1.den <- interaction(data2$C5b9.Tot.Glom,data2$C5aR1.den)
levels(data2$C5b9.Tot.Glom_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")
data2$C5b9.Mes_C5aR1.den <- interaction(data2$C5b9.Mes,data2$C5aR1.den)
levels(data2$C5b9.Mes_C5aR1.den) <- c("Low","Int/High","Int/High","Int/High")


#Create LR models to calculate ROC curves on test data (validation cohort)
# Model 1: Combined C5b9.Tot.Glom + C5aR1.den
fit_combined <- glm(PU ~ C5b9.Tot.Glom_C5aR1.den,
                    data = data, family = binomial)
# Model 2: C5aR1.den alone
fit_c5ar1 <- glm(PU ~ C5aR1.den,
                 data = data, family = binomial)
# Model 3: C5b9.Tot.Glom alone
fit_c9 <- glm(PU ~ C5b9.Tot.Glom,
              data = data, family = binomial)

#Apply trained models to test data (validation cohort)
data2$pred_combined <- predict(fit_combined,
                               newdata = data2,
                               type = "response") # Model 1
data2$pred_c5ar1 <- predict(fit_c5ar1,
                            newdata = data2,
                            type = "response") # Model 2
data2$pred_c9 <- predict(fit_c9,
                         newdata = data2,
                         type = "response") # Model 3

#Perform ROC analysis in validation cohort
roc_combined_data2 <- roc(response = data2$PU,
                          predictor = data2$pred_combined,
                          ci = TRUE,
                          boot.n = 1000) # Model 1 
roc_c5ar1_data2 <- roc(response = data2$PU,
                       predictor = data2$pred_c5ar1,
                       ci = TRUE,
                       boot.n = 1000) # Model 2
roc_c9_data2 <- roc(response = data2$PU,
                    predictor = data2$pred_c9,
                    ci = TRUE,
                    boot.n = 1000) # Model 3

#DeLong test combined vs C5aR1 or C9 on test data
p.adjust(c(roc.test(roc_combined_data2, roc_c5ar1_data2, method = "delong")$p.val,
           roc.test(roc_combined_data2, roc_c9_data2, method = "delong")$p.val),
         method = "BH")

# #Apply trained models to train data (exploratory cohort)
# data$pred_combined <- predict(fit_combined,newdata = data,type = "response") # Model 1
# data$pred_c5ar1 <- predict(fit_c5ar1,newdata = data,type = "response") # Model 2
# data$pred_c9 <- predict(fit_c9,newdata = data,type = "response") # Model 3
# #Perform ROC analysis in exploratory cohort
# roc_combined_data <- roc(response = data$PU,predictor = data$pred_combined,ci = TRUE,boot.n = 1000) # Model 1 
# roc_c5ar1_data <- roc(response = data$PU,predictor = data$pred_c5ar1,ci = TRUE,boot.n = 1000) # Model 2
# roc_c9_data <- roc(response = data$PU,predictor = data$pred_c9,ci = TRUE,boot.n = 1000) # Model 3
# #DeLong test combined vs C5aR1 or C9 on train data
# p.adjust(c(roc.test(roc_combined_data, roc_c5ar1_data, method = "delong")$p.val,
#            roc.test(roc_combined_data, roc_c9_data, method = "delong")$p.val),
#          method = "BH")

#Plot ROC curves for external validation

# Convert ROC objects to data frames
roc_combined_df <- data.frame(
  fpr = rev(1 - roc_combined_data2$specificities),
  tpr = rev(roc_combined_data2$sensitivities),
  predictor = "C5b9.Tot.Glom + C5aR1.den",
  auc = as.numeric(auc(roc_combined_data2))
)

roc_c5ar1_df <- data.frame(
  fpr = rev(1 - roc_c5ar1_data2$specificities),
  tpr = rev(roc_c5ar1_data2$sensitivities),
  predictor = "C5aR1.den",
  auc = as.numeric(auc(roc_c5ar1_data2))
)

roc_c9_df <- data.frame(
  fpr = rev(1 - roc_c9_data2$specificities),
  tpr = rev(roc_c9_data2$sensitivities),
  predictor = "C5b9.Tot.Glom",
  auc = as.numeric(auc(roc_c9_data2))
)

# Combine
roc_data_validation <- bind_rows(
  roc_combined_df,
  roc_c5ar1_df,
  roc_c9_df
)

# Choose one point on each curve for direct labels
label_points <- roc_data_validation %>%
  group_by(predictor) %>%
  arrange(abs(tpr - 0.75)) %>%
  slice(1) %>%
  ungroup()

# Colours and line types
color_palette <- c(
  "C5b9.Tot.Glom + C5aR1.den" = "#C8202A",
  "C5aR1.den"                  = "#3A7DC9",
  "C5b9.Tot.Glom"              = "#2B2B2B"
)

linetype_palette <- c(
  "C5b9.Tot.Glom + C5aR1.den" = "solid",
  "C5aR1.den"                  = "solid",
  "C5b9.Tot.Glom"              = "dashed"
)

# Plot
ggplot(
  roc_data_validation,
  aes(x = fpr, y = tpr, group = predictor)
) +
  
  # ROC curves
  geom_line(
    aes(color = predictor, linetype = predictor),
    linewidth = 1.2
  ) +
  
  # Chance line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.6
  ) +
  
  # Direct labels + AUC
  geom_label_repel(
    data = label_points,
    aes(
      label = paste0(
        predictor,
        "\nAUC = ",
        sprintf("%.3f", auc)
      ),
      color = predictor
    ),
    size = 4,
    show.legend = FALSE,
    segment.color = "gray50",
    segment.size = 0.5,
    box.padding = 0.5,
    point.padding = 0.25,
    nudge_y = 0.03,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = color_palette,
    name = "Predictor"
  ) +
  
  scale_linetype_manual(
    values = linetype_palette,
    name = "Predictor"
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  labs(
    title = "ROC Curves: Validation Cohort",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  
  theme_minimal(base_size = 20) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.box = "vertical",
    legend.position = "none"
  )
