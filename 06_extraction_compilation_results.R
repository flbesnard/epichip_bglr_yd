library(data.table)
library(ggplot2)

## Load trait list
Liste_car <- fread("/g2b/fbesnard/GITHUB/HAPCAR/Choix_caracteres_ibl")
Liste_car <- Liste_car[Race == RACE]

RACE=66
## Create output directory structure
output_dir <- paste0("/travail/", Sys.getenv("USER"), "/bglr2_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_dir, "/plots"), showWarnings = FALSE)
dir.create(paste0(output_dir, "/tables"), showWarnings = FALSE)

## Initialize results table with all columns
Results <- data.frame(
  Caractere = character(),
  Var_G_model1 = numeric(),
  Var_G_model2 = numeric(),
  Var_Epi_model1 = numeric(),
  Var_Residuelle_model1 = numeric(),
  Var_Residuelle_model2 = numeric(),
  H2_model1 = numeric(),
  H2_model2 = numeric(),
  Cor_Model1bis = numeric(),
  Cor_Model2bis = numeric(),
  RMSE_Model1bis = numeric(),
  RMSE_Model2bis = numeric(),
  Bias_Model1bis = numeric(),
  Bias_Model2bis = numeric(),
  stringsAsFactors = FALSE
)

## Loop through each trait
for (car in Liste_car$id_caractere_ctig) {
  
  cat("\n========================================\n")
  cat("Processing trait:", car, "\n")
  cat("========================================\n")
  
  ## Create trait-specific directory
  trait_dir <- paste0(output_dir, "/", car)
  dir.create(trait_dir, showWarnings = FALSE)
  
  ## -------------------------------------------------------------------
  ## ANALYSIS 1: Variance components comparison (Model 1 vs Model 2)
  ## -------------------------------------------------------------------
  
  # Load models
  load(paste0("/travail/", Sys.getenv("USER"), "/bglr2/", car, "/model1.RData"))
  load(paste0("/travail/", Sys.getenv("USER"), "/bglr2/", car, "/model2.RData"))
  
  # Extract variance components
  G_model1 <- model1$ETA$gebv$varU
  G_model2 <- model2$ETA$gebv2$varU
  Epi_model1 <- model1$ETA$met$varU
  ResVar_model1 <- model1$varE
  ResVar_model2 <- model2$varE
  
  # Calculate heritabilities
  H2_model1 <- (G_model1 + Epi_model1) / (G_model1 + Epi_model1 + ResVar_model1)
  H2_model2 <- G_model2 / (G_model2 + ResVar_model2)
  
  ## -------------------------------------------------------------------
  ## ANALYSIS 2: Top methylation marks (Model 3)
  ## -------------------------------------------------------------------
  
  load(paste0("/travail/", Sys.getenv("USER"), "/bglr2/", car, "/model3.RData"))
  
  # Extract methylation weights
  met_weights <- model3$ETA$metcpg$b
  
  # Select top 1000 marks
  top_met_indices <- order(abs(met_weights), decreasing = TRUE)[1:1000]
  top_met_weights <- met_weights[top_met_indices]
  
  # Save top methylation marks
  top_met_df <- data.frame(
    Mark_Index = top_met_indices, 
    Weight = top_met_weights,
    Abs_Weight = abs(top_met_weights)
  )
  write.table(
    top_met_df, 
    file = paste0(trait_dir, "/top_1000_met_marks.txt"), 
    row.names = FALSE, 
    sep = "\t", 
    quote = FALSE
  )
  
  ## -------------------------------------------------------------------
  ## ANALYSIS 3: Predictive ability (Model 1bis vs Model 2bis)
  ## -------------------------------------------------------------------
  
  # Load validation models
  load(paste0("/travail/", Sys.getenv("USER"), "/bglr2/", car, "/model1bis.RData"))
  load(paste0("/travail/", Sys.getenv("USER"), "/bglr2/", car, "/model2bis.RData"))
  
  # Prepare validation data
  merged_data <- data.table(
    y = model1$y,
    y_masked = model1bis$y,
    yhat1bis = model1bis$yHat,
    yhat2bis = model2bis$yHat
  )
  
  # Keep only validation animals (have y but not y_masked)
  Validation_population <- merged_data[!is.na(y) & is.na(y_masked)]
  
  # Calculate accuracy metrics
  corr_model1bis <- cor(Validation_population$y, Validation_population$yhat1bis, use = "complete.obs")
  corr_model2bis <- cor(Validation_population$y, Validation_population$yhat2bis, use = "complete.obs")
  
  rmse1 <- sqrt(mean((Validation_population$y - Validation_population$yhat1bis)^2, na.rm = TRUE))
  rmse2 <- sqrt(mean((Validation_population$y - Validation_population$yhat2bis)^2, na.rm = TRUE))
  
  bias1 <- coef(lm(y ~ yhat1bis, data = Validation_population))[2]
  bias2 <- coef(lm(y ~ yhat2bis, data = Validation_population))[2]
  
  # Calculate ranking metrics
  Validation_population[, Rank_Model1bis := frank(-yhat1bis, ties.method = "average")]
  Validation_population[, Rank_Model2bis := frank(-yhat2bis, ties.method = "average")]
  Validation_population[, Rank_Shift := Rank_Model1bis - Rank_Model2bis]
  Validation_population[, Delta_EBV := yhat1bis - yhat2bis]
  
  ## -------------------------------------------------------------------
  ## SAVE PLOTS
  ## -------------------------------------------------------------------
  
  # Plot 1: GEBV comparison (Model 1bis vs Model 2bis)
  p1 <- ggplot(Validation_population, aes(x = yhat1bis, y = yhat2bis)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_smooth(method = "lm", color = "darkblue", se = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = paste("GEBV Comparison:", car),
      subtitle = paste0("Cor = ", round(cor(Validation_population$yhat1bis, Validation_population$yhat2bis), 3)),
      x = "GEBV Model 1bis (G + Epi)",
      y = "GEBV Model 2bis (G only)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(
    filename = paste0(trait_dir, "/plot1_GEBV_comparison.png"),
    plot = p1,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # Plot 2: Distribution of ΔEBV
  p2 <- ggplot(Validation_population, aes(x = Delta_EBV)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(
      title = paste("Distribution of ΔEBV:", car),
      subtitle = paste0("Mean = ", round(mean(Validation_population$Delta_EBV), 3), 
                       " | SD = ", round(sd(Validation_population$Delta_EBV), 3)),
      x = "ΔEBV (G + Epi - G only)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(
    filename = paste0(trait_dir, "/plot2_Delta_EBV_distribution.png"),
    plot = p2,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # Plot 3: Rank shift distribution
  p3 <- ggplot(Validation_population, aes(x = Rank_Shift)) +
    geom_histogram(bins = 30, fill = "coral", color = "black", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(
      title = paste("Rank Shift Distribution:", car),
      subtitle = paste0("Mean = ", round(mean(Validation_population$Rank_Shift), 2), 
                       " | Max = ", max(abs(Validation_population$Rank_Shift))),
      x = "Rank Shift (Model 1bis - Model 2bis)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(
    filename = paste0(trait_dir, "/plot3_Rank_shift_distribution.png"),
    plot = p3,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # Plot 4: Prediction accuracy - Model 1bis
  p4 <- ggplot(Validation_population, aes(x = yhat1bis, y = y)) +
    geom_point(color = "darkgreen", alpha = 0.6) +
    geom_smooth(method = "lm", color = "darkgreen", se = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = paste("Prediction Accuracy - Model 1bis (G + Epi):", car),
      subtitle = paste0("Cor = ", round(corr_model1bis, 3), " | RMSE = ", round(rmse1, 3), 
                       " | Bias = ", round(bias1, 3)),
      x = "Predicted GEBV",
      y = "Observed Y"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(
    filename = paste0(trait_dir, "/plot4_Accuracy_Model1bis.png"),
    plot = p4,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # Plot 5: Prediction accuracy - Model 2bis
  p5 <- ggplot(Validation_population, aes(x = yhat2bis, y = y)) +
    geom_point(color = "darkorange", alpha = 0.6) +
    geom_smooth(method = "lm", color = "darkorange", se = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = paste("Prediction Accuracy - Model 2bis (G only):", car),
      subtitle = paste0("Cor = ", round(corr_model2bis, 3), " | RMSE = ", round(rmse2, 3), 
                       " | Bias = ", round(bias2, 3)),
      x = "Predicted GEBV",
      y = "Observed Y"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(
    filename = paste0(trait_dir, "/plot5_Accuracy_Model2bis.png"),
    plot = p5,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ## -------------------------------------------------------------------
  ## SAVE VALIDATION DATA
  ## -------------------------------------------------------------------
  
  write.table(
    Validation_population,
    file = paste0(trait_dir, "/validation_data.txt"),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
  
  ## -------------------------------------------------------------------
  ## STORE RESULTS IN GLOBAL TABLE
  ## -------------------------------------------------------------------
  
  Results <- rbind(
    Results,
    data.frame(
      Caractere = car,
      Var_G_model1 = G_model1,
      Var_G_model2 = G_model2,
      Var_Epi_model1 = Epi_model1,
      Var_Residuelle_model1 = ResVar_model1,
      Var_Residuelle_model2 = ResVar_model2,
      H2_model1 = H2_model1,
      H2_model2 = H2_model2,
      Cor_Model1bis = corr_model1bis,
      Cor_Model2bis = corr_model2bis,
      RMSE_Model1bis = rmse1,
      RMSE_Model2bis = rmse2,
      Bias_Model1bis = bias1,
      Bias_Model2bis = bias2,
      stringsAsFactors = FALSE
    )
  )
  
  cat("Completed trait:", car, "\n")
}

## -------------------------------------------------------------------
## SAVE FINAL RESULTS TABLE
## -------------------------------------------------------------------

# Save as text file
write.table(
  Results,
  file = paste0(output_dir, "/tables/Complete_Results_Summary.txt"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# Save as CSV
write.csv(
  Results,
  file = paste0(output_dir, "/tables/Complete_Results_Summary.csv"),
  row.names = FALSE
)

# Save as RData
save(Results, file = paste0(output_dir, "/tables/Complete_Results_Summary.RData"))

## -------------------------------------------------------------------
## CREATE SUMMARY PLOTS ACROSS ALL TRAITS
## -------------------------------------------------------------------

# Plot: H2 comparison
p_h2 <- ggplot(Results, aes(x = H2_model2, y = H2_model1)) +
  geom_point(size = 3, color = "darkblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(label = Caractere), vjust = -0.5, size = 2.5) +
  labs(
    title = "Heritability Comparison: G vs G+Epi",
    x = "H² Model 2 (G only)",
    y = "H² Model 1 (G + Epi)"
  ) +
  theme_minimal()

ggsave(
  filename = paste0(output_dir, "/plots/Summary_H2_comparison.png"),
  plot = p_h2,
  width = 10,
  height = 8,
  dpi = 300
)

# Plot: Correlation comparison
p_cor <- ggplot(Results, aes(x = Cor_Model2bis, y = Cor_Model1bis)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(label = Caractere), vjust = -0.5, size = 2.5) +
  labs(
    title = "Prediction Accuracy Comparison",
    x = "Correlation Model 2bis (G only)",
    y = "Correlation Model 1bis (G + Epi)"
  ) +
  theme_minimal()

ggsave(
  filename = paste0(output_dir, "/plots/Summary_Correlation_comparison.png"),
  plot = p_cor,
  width = 10,
  height = 8,
  dpi = 300
)

# Plot: Epigenetic variance contribution
p_epi <- ggplot(Results, aes(x = reorder(Caractere, Var_Epi_model1), y = Var_Epi_model1)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Epigenetic Variance Contribution by Trait",
    x = "Trait",
    y = "Epigenetic Variance"
  ) +
  theme_minimal()

ggsave(
  filename = paste0(output_dir, "/plots/Summary_Epigenetic_variance.png"),
  plot = p_epi,
  width = 10,
  height = 8,
  dpi = 300
)

cat("\n========================================\n")
cat("Analysis completed successfully!\n")
cat("Results saved in:", output_dir, "\n")
cat("========================================\n")