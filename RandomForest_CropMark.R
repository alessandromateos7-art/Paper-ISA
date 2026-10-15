# =============================================================================
# Random Forest crop-mark discrimination using Sentinel-2 vegetation indices
# Study area: Lempeș region, north of Hărman, Romania
# Author: Alessandro Mateos Sierra
# Version: 1.0 (original methodological workflow, documented for publication)
# =============================================================================
# This script evaluates all 255 combinations of 1–4 predictors drawn from nine
# vegetation indices. A fixed, non-stratified 80/20 random split is used for
# model comparison, as described in the associated manuscript.
#
# IMPORTANT: Validation data are used for both comparison and model selection.
# The resulting metrics are exploratory model-selection performance, NOT an
# independent assessment of generalisation. The split is not spatially grouped.
# =============================================================================

# 1. PACKAGES ------------------------------------------------------------------
library(terra)
library(randomForest)
library(pROC)
library(caret)
library(viridis)

# 2. INPUT AND OUTPUT PATHS -----------------------------------------------------
# Replace these example relative paths with the locations of your own files.
# Required inputs: a nine-band GeoTIFF and two reference point shapefiles.
indices_path <- "data/Indice_Lempes_Mask2.tif"
crop_path <- "data/Crop.shp"
background_path <- "data/Fondo.shp"
output_path <- "results/RandomForest_crop_mark_score_scenario2.tif"

# 3. LOAD DATA -----------------------------------------------------------------
Indices <- rast(indices_path)
stopifnot(nlyr(Indices) == 9)
names(Indices) <- c(
  "NDVI", "GNDVI", "MSR", "SAVI", "CLgreen", "RVI", "GRVI", "NRVI", "MTVI"
)
crop <- vect(crop_path)
fondo <- vect(background_path)

# 4. EXTRACT VALUES AND PREPARE REFERENCE DATA ---------------------------------
# terra::extract() returns an ID column. It is intentionally retained to
# preserve the original script's data structure and random sampling procedure.
val_crop <- extract(Indices, crop)
val_fondo <- extract(Indices, fondo)
val_crop$clase <- "crop"
val_fondo$clase <- "fondo"
df <- rbind(val_crop, val_fondo)
df <- na.omit(df)
df$clase <- as.factor(df$clase)
# Factor levels are explicitly specified to make the class interpretation
# stable across environments. These match the original alphabetical ordering.
df$clase <- factor(df$clase, levels = c("crop", "fondo"))

# 5. FIXED TRAINING / VALIDATION PARTITION --------------------------------------
# Keep the original random split (without stratification) to reproduce the
# manuscript workflow. Do not replace with createDataPartition() here.
set.seed(123)
train_idx <- sample(1:nrow(df), 0.8 * nrow(df))
train <- df[train_idx, ]
test <- df[-train_idx, ]

# 6. ALL COMBINATIONS OF ONE TO FOUR INDICES -----------------------------------
vars_all <- names(Indices)
comb_list <- list()
max_vars <- 4
for (i in 1:max_vars) {
  comb_list <- c(comb_list, combn(vars_all, i, simplify = FALSE))
}
stopifnot(length(comb_list) == 255)

# 7. FIT AND EVALUATE RANDOM FOREST MODELS -------------------------------------
resultados <- data.frame()
for (vars in comb_list) {
  cat("\n=============================\n")
  cat("Variables:", paste(vars, collapse = ", "), "\n")

  model_formula <- as.formula(
    paste("clase ~", paste(vars, collapse = " + "))
  )

  # Reset the seed before EACH model, preserving the original fitting strategy.
  # randomForest defaults: bootstrap sampling with replacement, nodesize = 1,
  # mtry = floor(sqrt(number of predictors)) for classification.
  set.seed(123)
  model <- randomForest(model_formula, data = train, ntree = 500)

  # The score is the fraction of tree votes for 'crop', not a calibrated
  # probability of archaeological crop-mark occurrence.
  prob <- predict(model, test, type = "prob")[, "crop"]
  roc_obj <- roc(
    test$clase,
    prob,
    levels = c("fondo", "crop"),
    direction = "<",
    quiet = TRUE
  )
  auc_val <- as.numeric(auc(roc_obj))

  pred_class <- predict(model, test, type = "response")
  cm <- confusionMatrix(pred_class, test$clase, positive = "crop")
  acc <- as.numeric(cm$overall["Accuracy"])
  kappa <- as.numeric(cm$overall["Kappa"])
  sens <- as.numeric(cm$byClass["Sensitivity"])
  spec <- as.numeric(cm$byClass["Specificity"])
  prec <- as.numeric(cm$byClass["Pos Pred Value"])
  f1 <- as.numeric(cm$byClass["F1"])
  bal_acc <- as.numeric(cm$byClass["Balanced Accuracy"])

  resultados <- rbind(
    resultados,
    data.frame(
      variables = paste(vars, collapse = "+"),
      n_vars = length(vars),
      AUC = auc_val,
      Accuracy = acc,
      Kappa = kappa,
      Sensitivity = sens,
      Specificity = spec,
      Precision = prec,
      F1 = f1,
      Balanced_Accuracy = bal_acc
    )
  )
}

# 8. RANK MODELS BY VALIDATION AUC ----------------------------------------------
resultados <- resultados[order(-resultados$AUC), ]
print(head(resultados, 20))

# 9. SELECT SCENARIO 2 FOR FINAL MAPPING ---------------------------------------
# The second-ranked AUC model is intentionally selected rather than the first.
# In the reported top-five comparison, Scenario 2 (GNDVI + RVI + MTVI) had the
# highest Accuracy, Kappa and Specificity, with the same Sensitivity as the
# other four scenarios. Scenario 1 had a marginally higher AUC.
# IMPORTANT: This assumes the second row remains GNDVI+RVI+MTVI with the same
# input data, factor levels, package behaviour and random partition.
best <- resultados[2, ]
best_vars <- unlist(strsplit(as.character(best$variables), "\\+"))
if (!identical(best_vars, c("GNDVI", "RVI", "MTVI"))) {
  warning("Second-ranked model differs from the manuscript's Scenario 2; verify inputs and package versions.")
}
cat("\n=============================\n")
cat("SELECTED MODEL:", paste(best_vars, collapse = ", "), "\n")
cat("AUC:", best$AUC, "\n")

# 10. REFIT SELECTED MODEL ON ALL REFERENCE OBSERVATIONS ------------------------
formula_best <- as.formula(
  paste("clase ~", paste(best_vars, collapse = " + "))
)
set.seed(123)
model_best <- randomForest(formula_best, data = df, ntree = 500)

# 11. MAP RANDOM FOREST CROP-MARK SCORE -----------------------------------------
pred_raster <- predict(Indices[[best_vars]], model_best, type = "prob")
prob_crop <- pred_raster[["crop"]]
names(prob_crop) <- "RF_crop_mark_score"
plot(prob_crop, main = paste(best_vars, collapse = " + "), col = viridis(100))

# 12. EXPORT --------------------------------------------------------------------
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
writeRaster(prob_crop, output_path, overwrite = TRUE)
# Optional: save the complete model-comparison table for auditability.
# write.csv(resultados, "results/RandomForest_all_255_models.csv", row.names = FALSE)
