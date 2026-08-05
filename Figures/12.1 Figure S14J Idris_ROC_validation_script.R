# ============================================================
# ROC curve analysis — Validation cohort
# Markers: C5b9.Tot.Glom, C5aR1.den, and their combination
# Outcome: proteinuria > 0.5 g/g (binary: 1 = yes, 0 = no)
# ============================================================

# --- Packages -----------------------------------------------
library(readxl)
library(pROC)
library(ggplot2)

# --- Data ---------------------------------------------------
df <- read_excel("matrice_C9C5aR.xlsx", sheet = "real data")

# Clean column names if needed
colnames(df) <- trimws(colnames(df))

# Variables:
#   C5b9.Tot.Glom : binary (0 = low, 1 = high)
#   C5aR1.den     : binary (0 = low, 1 = high; cut-off = 500.9012)
#   proteinuria   : binary outcome (0 = <0.5 g/g, 1 = ≥0.5 g/g)

# Remove rows with missing values
df <- df[complete.cases(df[, c("C5b9.Tot.Glom", "C5aR1.den", "proteinuria")]), ]
cat("n =", nrow(df), "\n")
cat("Outcome: ", sum(df$proteinuria == 1), "positive /",
    sum(df$proteinuria == 0), "negative\n\n")

# --- Logistic regression models -----------------------------

# Model 1: Combined C5b9.Tot.Glom + C5aR1.den
fit_combined <- glm(proteinuria ~ C5b9.Tot.Glom + C5aR1.den,
                    data = df, family = binomial)

# Model 2: C5aR1.den alone
fit_c5ar1 <- glm(proteinuria ~ C5aR1.den,
                 data = df, family = binomial)

# Model 3: C5b9.Tot.Glom alone
fit_c5b9 <- glm(proteinuria ~ C5b9.Tot.Glom,
                data = df, family = binomial)

# --- ROC curves ---------------------------------------------
roc_combined <- roc(df$proteinuria,
                    predict(fit_combined, type = "response"),
                    ci = TRUE, boot.n = 1000)

roc_c5ar1    <- roc(df$proteinuria,
                    predict(fit_c5ar1, type = "response"),
                    ci = TRUE, boot.n = 1000)

roc_c5b9     <- roc(df$proteinuria,
                    predict(fit_c5b9, type = "response"),
                    ci = TRUE, boot.n = 1000)

# --- Print AUC + 95% CI ------------------------------------
cat("AUC — C5b9.Tot.Glom + C5aR1.den: ",
    round(auc(roc_combined), 3),
    " [", round(ci(roc_combined)[1], 3), "–",
    round(ci(roc_combined)[3], 3), "]\n")

cat("AUC — C5aR1.den alone:            ",
    round(auc(roc_c5ar1), 3),
    " [", round(ci(roc_c5ar1)[1], 3), "–",
    round(ci(roc_c5ar1)[3], 3), "]\n")

cat("AUC — C5b9.Tot.Glom alone:        ",
    round(auc(roc_c5b9), 3),
    " [", round(ci(roc_c5b9)[1], 3), "–",
    round(ci(roc_c5b9)[3], 3), "]\n\n")

# --- DeLong tests -------------------------------------------
cat("DeLong test — Combined vs C5b9.Tot.Glom alone:\n")
print(roc.test(roc_combined, roc_c5b9, method = "delong"))

cat("\nDeLong test — Combined vs C5aR1.den alone:\n")
print(roc.test(roc_combined, roc_c5ar1, method = "delong"))

# --- Plot ---------------------------------------------------
# Convert ROC objects to data frames for ggplot2
roc_to_df <- function(roc_obj, label) {
  data.frame(
    FPR   = 1 - roc_obj$specificities,
    TPR   = roc_obj$sensitivities,
    Model = label
  )
}

plot_df <- rbind(
  roc_to_df(roc_combined, paste0("C5b9.Tot.Glom + C5aR1.den\nAUC = ",
            round(auc(roc_combined), 3),
            " [", round(ci(roc_combined)[1], 2), "–",
            round(ci(roc_combined)[3], 2), "]")),
  roc_to_df(roc_c5ar1,    paste0("C5aR1.den\nAUC = ",
            round(auc(roc_c5ar1), 3),
            " [", round(ci(roc_c5ar1)[1], 2), "–",
            round(ci(roc_c5ar1)[3], 2), "]")),
  roc_to_df(roc_c5b9,     paste0("C5b9.Tot.Glom\nAUC = ",
            round(auc(roc_c5b9), 3),
            " [", round(ci(roc_c5b9)[1], 2), "–",
            round(ci(roc_c5b9)[3], 2), "]"))
)

cols <- c("#C8202A", "#3A7DC9", "#2B2B2B")
names(cols) <- unique(plot_df$Model)

p <- ggplot(plot_df, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey70", linewidth = 0.5) +
  scale_color_manual(values = cols) +
  scale_linetype_manual(values = c("solid", "solid", "dashed")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  labs(x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       color = NULL, linetype = NULL) +
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.60, 0.20),
        legend.text     = element_text(size = 7),
        legend.key.size = unit(0.8, "lines"))

ggsave("ROC_validation_cohort.png", plot = p,
       width = 4.13, height = 5.10, units = "cm",
       dpi = 400, bg = "white")

cat("Figure saved as ROC_validation_cohort.png\n")
