######################
### SI - Half a century of per- and polyfluoroalkyl substances in Northern Gannet eggs: Impact of regulations
### Step 6 - Temporal trends of PFAS in NOGA of Bonaventure Island, QC using GAMs
### Created by Anaïs Fournier on 05-14-2025
#####################

# Load libraries
library(ggplot2)
library(mgcv)
library(dplyr)
library(tidyr)
library(gratia)
library(MuMIn)
library(tibble)
library(purrr)
library(patchwork)
library(car)

# Import datasets
dat <- read.csv2("/Users/anaisfournier/Desktop/Stats/Step 3&4 - Correlation of Moisture, Lipid, d15N, d13CS and PFAS Concentrations.csv")
dat[, 3:7] <- lapply(dat[, 3:7], as.numeric)
data <- read.csv2("/Users/anaisfournier/Desktop/Stats/Step 5 - Lipid-Corrected Concentrations for d15N Corection.csv")
data[, 3:7] <- lapply(data[, 3:7], as.numeric)
data_corr <- read.csv2("/Users/anaisfournier/Desktop/Stats/Step 6 - Lipid+d15N-Corrected Concentrations for GAMs.csv")
data_corr[, 3:7] <- lapply(data_corr[, 3:7], as.numeric)


### Comparison of top Generalized Additive Models (GAM) for each compound ###

# List of compounds
compounds <- c("PFOA", "PFNA", "PFDA", "PFUdA", "PFDoA", "PFTrDA", "PFTeDA",
               "PFHxS", "PFOS", "PFDS", "PFCA", "PFSA", "PFAS")

# List of covariates
vars_final <- c("Year", "Moisture", "Lipid", "d15N", "d13CS")

# Function to generate 3 choices for each variable: not included, linear, or smooth
options_per_var <- function(var) {
  c("", var, paste0("s(", var, ")"))
}

# Generate all possible combinations
all_combos <- expand.grid(map(vars_final, options_per_var), stringsAsFactors = FALSE)

# Remove combinations that contain both "var" and "s(var)"
valid_combos <- all_combos[
  apply(all_combos, 1, function(row) {
    all(sapply(vars_final, function(var) {
      !(var %in% row & paste0("s(", var, ")") %in% row)
    }))
  }),
]

# Remove null model if all terms are ""
valid_combos <- valid_combos[rowSums(valid_combos == "") < length(vars_final), ]

# Create model formulas
model_formulas_data <- apply(valid_combos, 1, function(row) {
  terms <- row[row != ""]
  paste(terms, collapse = " + ")
})

# Add intercept-only model
model_formulas_data <- c("1", model_formulas_data)

# Total number of models
length(model_formulas_data)


analyze_compound_data <- function(compound, dat) {
  # Select relevant variables and remove rows with missing values
  all_vars <- unique(c("Year", "Moisture", "Lipid", "d15N", "d13CS", compound))
  data_clean <- dat %>%
    select(all_of(gsub("s\\(|\\)", "", all_vars))) %>%  # Remove 's()' for selection
    na.omit()
  
  models <- list()
  
  # Fit GAMs for each model formula
  for (formula in model_formulas_data) {
    formula_text <- paste0(compound, " ~ ", formula)
    model <- tryCatch({
      gam(as.formula(formula_text), data = data_clean, method = "REML")
    }, error = function(e) NULL)
    
    if (!is.null(model)) {
      models[[formula]] <- model
    }
  }
  
  # Skip compound if no valid model could be fitted
  if (length(models) == 0) {
    message(paste("No valid models could be fitted for compound:", compound))
    return(NULL)
  }
  
  # Generate AICc model selection table
  aictab <- model.sel(models, rank = "AICc")
  model_names <- rownames(as.data.frame(aictab))
  deltaAICc <- aictab$delta
  
  # Extract null model details
  null_model_included <- "1" %in% model_names
  null_model_position <- which(model_names == "1")
  null_model_AICc <- if (length(null_model_position) > 0) aictab$AICc[null_model_position] else NA
  null_model_weight <- if (length(null_model_position) > 0) aictab$weight[null_model_position] else NA
  null_model_in_top <- null_model_position %in% which(deltaAICc < 2)
  
  # Select models within ΔAICc < 2
  top_models <- which(deltaAICc < 2)
  if (length(top_models) == 0) return(NULL)
  
  selected <- tibble(
    Model = model_names[top_models],
    AICc = aictab$AICc[top_models],
    deltaAICc = deltaAICc[top_models],
    Weight = aictab$weight[top_models],
    edf = sapply(model_names[top_models], function(name) sum(summary(models[[name]])$edf)),
    R2 = sapply(model_names[top_models], function(name) summary(models[[name]])$r.sq)
  )
  
  # Compile final result
  result <- tibble(
    Compound = compound,
    Models = paste(selected$Model, collapse = "; "),
    AICc = paste(round(selected$AICc, 2), collapse = "; "),
    deltaAICc = paste(round(selected$deltaAICc, 2), collapse = "; "),
    Weight = paste(round(selected$Weight, 3), collapse = "; "),
    edf = paste(round(selected$edf, 2), collapse = "; "),
    R2 = paste(round(selected$R2, 3), collapse = "; "),
    Null_model_AICc = round(null_model_AICc, 2),
    Null_model_weight = format(null_model_weight, scientific = TRUE),
    Null_in_top = null_model_in_top
  )
  
  return(result)
}

# Run for data (lipid-corrected only)
results <- lapply(compounds, function(comp) analyze_compound_data(comp, dat))
models_GAMs <- bind_rows(Filter(Negate(is.null), results))

view(models_GAMs)
write.csv2(models_GAMs, "/Users/anaisfournier/Desktop/Stats/Models_GAMs.csv", row.names = FALSE)

# Create an empty list to store averaged models for each compound
models_avg_list <- list()

# Loop over each target compound
for (compound in compounds) {
  
  # Select all relevant covariates including the compound of interest
  all_vars <- unique(c("Year", "Moisture", "Lipid", "d15N", "d13CS", compound))
  
  # Subset and clean the dataset by removing rows with missing values
  data_clean <- dat %>%
    select(all_of(all_vars)) %>%
    na.omit()
  
  models <- list()
  
  # Fit generalized additive models (GAMs) using each candidate formula
  for (formula in model_formulas_data) {
    formula_text <- paste0(compound, " ~ ", formula)
    model <- tryCatch({
      gam(as.formula(formula_text), data = data_clean, method = "REML")
    }, error = function(e) NULL)  # Handle models that fail to converge
    
    if (!is.null(model)) {
      models[[formula]] <- model
    }
  }
  
  # Skip to next compound if no models were successfully fitted
  if (length(models) == 0) next
  
  # Create an AICc model selection table
  aictab <- model.sel(models, rank = "AICc")
  df_tab <- as.data.frame(aictab)
  df_tab$Model <- rownames(df_tab)
  
  # Sort models by weight and retain those contributing up to 95% cumulative weight
  df_tab <- df_tab %>%
    arrange(desc(weight)) %>%
    mutate(cum_weight = cumsum(weight)) %>%
    filter(cum_weight <= 0.95 | row_number() == 1)
  
  selected_model_names <- df_tab$Model
  
  # Retrieve selected models for model averaging
  selected_models <- models[selected_model_names]
  
  # Perform model averaging using selected models (within 95% cumulative AICc weight)
  model_avg <- model.avg(selected_models)
  
  # Store averaged model in the list
  models_avg_list[[compound]] <- model_avg
}

# Extract variable importance (sum of weights) from each averaged model
importance_data <- lapply(names(models_avg_list), function(compound) {
  model_avg <- models_avg_list[[compound]]
  imp <- sw(model_avg)  # sw() returns sum of weights for each predictor
  
  tibble(
    Compound = compound,
    Variable = names(imp),
    Importance = as.numeric(imp)
  )
})

# Combine all results into a single dataframe
importance_df <- bind_rows(importance_data)

# Convert raw importance values into percentages per compound
importance_df <- importance_df %>%
  mutate(ImportancePercent = 100 * Importance)

# Reconvert the 'Compound' column into a factor with custom order for plotting
importance_df$Compound <- factor(importance_df$Compound,
                                 levels = c("PFOA", "PFNA", "PFDA", "PFUdA", "PFDoA",
                                            "PFTrDA", "PFTeDA", "PFHxS", "PFOS", "PFDS",
                                            "PFCA", "PFSA", "PFAS"))

# Define custom facet labels with parsed mathematical expressions for grouping
compound_labels <- c(
  "PFOA" = "PFOA  (C8)",
  "PFNA" = "PFNA  (C9)",
  "PFDA" = "PFDA  (C10)",
  "PFUdA" = "PFUdA  (C11)",
  "PFDoA" = "PFDoA  (C12)",
  "PFTrDA" = "PFTrDA  (C13)",
  "PFTeDA" = "PFTeDA  (C14)",
  "PFHxS" = "PFHxS  (C6)",
  "PFOS"  = "PFOS  (C8)",
  "PFDS"  = "PFDS  (C10)",
  "PFCA"  = "Sigma[13]*PFCAs",
  "PFSA"  = "Sigma[4]*PFSAs",
  "PFAS"  = "Sigma[17]*PFASs"
)

 variable_labels <- c(
  "d15N"      = expression(delta^15 * "N"),
  "d13CS"     = expression(delta^13 * "C"),
  "s(d15N)"   = expression("s(" * delta^15 * "N" * ")"),
  "s(d13CS)"  = expression("s(" * delta^13 * "C" * ")")
)

# Create a barplot of variable importance per compound with custom facet labels
ggplot(importance_df, aes(x = ImportancePercent, y = Variable, fill = Variable)) +
  geom_col(show.legend = FALSE) +  
  scale_y_discrete(labels = variable_labels) + 
  scale_x_continuous(limits = c(0, 100)) +
  labs(x = "Relative importance (%)", y = "Variable") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank()
  ) +
  # Apply parsed expressions as custom labels for each compound facet
  facet_wrap(
    ~ Compound,
    scales = "free_y",
    labeller = labeller(Compound = as_labeller(compound_labels, label_parsed))  
  )


variable_means <- aggregate(ImportancePercent ~ Variable, data = importance_df, FUN = mean)
ordered_vars <- variable_means[order(-variable_means$ImportancePercent), "Variable"]
importance_df$Variable <- factor(importance_df$Variable, levels = ordered_vars)

# Create a heatmap
ggplot(importance_df, aes(x = Compound, y = Variable, fill = ImportancePercent)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(ImportancePercent, 0)), size = 3, color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Relative importance (%)") +
  scale_x_discrete(labels = function(x) parse(text = compound_labels[x])) +  
  scale_y_discrete(labels = variable_labels, limits = rev(levels(importance_df$Variable))) +
  labs(x = "PFAS", y = "Variable") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid = element_blank(),
    legend.position = "right"
  )



### Temporal trends ###

data$PFOS_PFHxS <- data$PFOS + data$PFHxS

data_corr$PFOS_PFHxS <- data_corr$PFOS + data_corr$PFHxS

# Definition of plot titles for each compound
plot_titles <- list(
  "PFOA"   = "PFOA", 
  "PFNA"   = "PFNA", 
  "PFDA"   = "PFDA", 
  "PFUdA"  = "PFUdA", 
  "PFDoA"  = "PFDoA", 
  "PFTrDA" = "PFTrDA", 
  "PFTeDA" = "PFTeDA",
  "PFHxS"  = "PFHxS", 
  "PFOS"   = "PFOS",
  "PFDS"   = "PFDS",
  "PFOS_PFHxS" = "PFOS + PFHxS",
  "PFAS"   = expression(Σ[17] * "PFASs"),
  "PFCA"   = expression(Σ[13] * "PFCAs"),
  "PFSA"   = expression(Σ[4]  * "PFSAs")
)


### PFOS
compounds <- "PFOS"

for (compound in compounds) {
  cat("\n Analysis for:", compound, "\n")
  
  # Remove missing values
  data_clean <- na.omit(data[, c("Year", compound)])
  
  # Fit a GAM model
  gam_model <- gam(as.formula(paste(compound, "~ s(Year, k=10)")), data = data_clean, method = "REML")
  
  # Compute derivatives and 95% confidence intervals
  deriv_values <- derivatives(gam_model, term = "s(Year)", type = "central", level = 0.95) %>%
    select(Year, .derivative, .lower_ci, .upper_ci) %>%
    mutate(
      trend = case_when(
        .lower_ci > 0 ~ "Significantly increasing",
        .upper_ci < 0 ~ "Significantly decreasing",
        TRUE ~ "Neutral"
      )
    )
  
  # Detect local maximum: transition from increasing to decreasing trend
  critical_points <- deriv_values %>%
    mutate(sign_change = lag(.derivative) > 0 & .derivative < 0) %>%
    filter(sign_change)
  
  if (nrow(critical_points) > 0) {
    
    # Select the year with the highest value of the compound
    max_year_data <- data_clean %>%
      filter(Year %in% critical_points$Year) %>%
      arrange(desc(!!sym(compound))) %>%
      slice(1)
    
    year_max <- max_year_data$Year
  
  }
    
  # Estimate confidence interval around the GAM maximum using predict()
  pred_data <- data.frame(Year = seq(min(data_clean$Year), max(data_clean$Year), length.out = 500))
  pred <- predict(gam_model, newdata = pred_data, se.fit = TRUE, type = "response")
  pred_data$fit <- pred$fit
  pred_data$se  <- pred$se.fit
  
  # Year corresponding to the maximum predicted GAM value
  year_max_gam <- pred_data$Year[which.max(pred_data$fit)]
    
  # Confidence interval for the GAM maximum
  ci_range <- pred_data %>%
      filter(fit >= max(fit) - 1.96 * se) %>%
      summarize(min_year = min(Year), max_year = max(Year))
    
  cat("\nLipid corrected GAM: detected maximum for", compound, "in year:", round(year_max_gam, 2), "\n")
  cat("95% confidence interval for the year of maximum: [", round(ci_range$min_year, 2), ", ", round(ci_range$max_year, 2), "]\n")
  
  
  # Get predicted concentration and 95% CI at the year of the maximum
  new_data <- data.frame(Year = year_max_gam)
  pred <- predict(gam_model, newdata = new_data, se.fit = TRUE)
  
  fit_value <- pred$fit           
  se_value <- pred$se.fit 
  lower_ci <- fit_value - 1.96 * se_value
  upper_ci <- fit_value + 1.96 * se_value
  
  ci_range$fit <- fit_value
  year_gam_max <- new_data$Year
  ci_range$year_gam_max <- year_gam_max
  
  cat("Lipid corrected GAM: concentration at year of maximum (", round(year_max_gam, 2), "): ", round(fit_value, 2), "ng/g ww\n")
  cat("95% confidence interval: [", round(lower_ci, 2), ", ", round(upper_ci, 2), "]\n")
  
  
  # Segments distinguishing increasing and decreasing trends
  deriv_segments <- deriv_values %>%
    arrange(Year) %>%
    mutate(
      is_neutral = (trend == "Neutral"),
      trend_change = (trend != lag(trend, default = first(trend))),
      segment_id = cumsum(trend_change | is_neutral)
    ) %>%
    filter(trend != "Neutral")
  
  deriv_segments$fit <- predict(gam_model, newdata = data.frame(Year = deriv_segments$Year))
  
  # Remove missing values in corrected dataset
  data_clean2 <- na.omit(data_corr[, c("Year", compound)])
  
  # Fit corrected GAM model
  gam_model2 <- gam(as.formula(paste(compound, "~ s(Year, k=10)")), data = data_clean2, method = "REML")
  
  # Compute derivatives and 95% CI
  deriv_values2 <- derivatives(gam_model2, term = "s(Year)", type = "central", level = 0.95) %>%
    select(Year, .derivative, .lower_ci, .upper_ci) %>%
    mutate(
      trend = case_when(
        .lower_ci > 0 ~ "Significantly increasing",
        .upper_ci < 0 ~ "Significantly decreasing",
        TRUE ~ "Neutral"
      )
    )
  
  # Detect local maximum 
  critical_points2 <- deriv_values2 %>%
    mutate(sign_change = lag(.derivative) > 0 & .derivative < 0) %>%
    filter(sign_change)
  
  if (nrow(critical_points2) > 0) {
    
    max_year_data2 <- data_clean2 %>%
      filter(Year %in% critical_points2$Year) %>%
      arrange(desc(!!sym(compound))) %>%
      slice(1)
    
    year_max2 <- max_year_data2$Year
    
  }
    
  pred_data2 <- data.frame(Year = seq(min(data_clean2$Year), max(data_clean2$Year), length.out = 500))
  pred_data2$fit <- predict(gam_model2, newdata = pred_data2, se.fit = TRUE)$fit
  pred_data2$se <- predict(gam_model2, newdata = pred_data2, se.fit = TRUE)$se.fit
    
  year_max_gam2 <- pred_data2$Year[which(pred_data2$fit == max(pred_data2$fit, na.rm = TRUE))]
    
  ci_range2 <- pred_data2 %>%
      filter(fit >= max(fit) - 1.96 * se) %>%
      summarize(min_year2 = min(Year), max_year2 = max(Year))
    
  cat("\nLipid-d15N corrected GAM: detected maximum for", compound, "in year:", round(year_max_gam2, 2), "\n")
  cat("95% CI for year of maximum: [", round(ci_range2$min_year2, 2), ", ", round(ci_range2$max_year2, 2), "]\n")
  
  # Predicted value and CI at the maximum year
  new_data2 <- data.frame(Year = year_max_gam2)
  pred2 <- predict(gam_model2, newdata = new_data2, se.fit = TRUE)
  
  fit_value2 <- pred2$fit           
  se_value2 <- pred2$se.fit 
  lower_ci2 <- fit_value2 - 1.96 * se_value2
  upper_ci2 <- fit_value2 + 1.96 * se_value2
  
  ci_range2$fit <- fit_value2
  year_gam_max2 <- new_data2$Year
  ci_range2$year_gam_max2 <- year_gam_max2
  
  cat("Lipid-d15N corrected GAM: concentration at year of maximum (", round(year_max_gam2, 2), "): ", round(fit_value2, 2), "ng/g ww\n")
  cat("95% confidence interval: [", round(lower_ci2, 2), ", ", round(upper_ci2, 2), "]\n")
  
  # Compute start and end of each significant trend period
  trend_periods <- deriv_segments %>%
    group_by(segment_id, trend) %>%
    summarise(
      start_year = min(Year),
      end_year = max(Year),
      .groups = "drop"
    ) %>%
    mutate(
      start_value = map_dbl(start_year, ~ predict(gam_model, newdata = data.frame(Year = .x))),
      end_value   = map_dbl(end_year, ~ predict(gam_model, newdata = data.frame(Year = .x)))
    )
  
  
  cat("\nSignificant trend periods (Lipid corrected model) for", compound, ":\n")
  print(trend_periods)
  
  # Same for corrected model
  deriv_segments2 <- deriv_values2 %>%
    arrange(Year) %>%
    mutate(
      is_neutral = (trend == "Neutral"),
      trend_change = (trend != lag(trend, default = first(trend))),
      segment_id = cumsum(trend_change | is_neutral)
    ) %>%
    filter(trend != "Neutral")
  
  deriv_segments2$fit <- predict(gam_model2, newdata = data.frame(Year = deriv_segments2$Year))
  
  trend_periods2 <- deriv_segments2 %>%
    group_by(segment_id, trend) %>%
    summarise(
      start_year = min(Year),
      end_year = max(Year),
      .groups = "drop"
    ) %>%
    mutate(
      start_value = map_dbl(start_year, ~ predict(gam_model2, newdata = data.frame(Year = .x))),
      end_value   = map_dbl(end_year, ~ predict(gam_model2, newdata = data.frame(Year = .x)))
    )
  
  cat("\nSignificant trend periods (Lipid-d15N corrected model) for", compound, ":\n")
  print(trend_periods2)
  
  # Plot 1: First derivative of the GAM without δ15N correction
  p0 <- ggplot(deriv_values, aes(x = Year, y = .derivative)) +
    geom_line(color = "#4682B4", linewidth = 0.8) +  # Main derivative curve
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#4682B4", alpha = 0.2) +  # Confidence interval
    geom_line(data = deriv_segments, aes(y = .derivative, group = segment_id, color = trend), linewidth = 1) +  # Significant trends
    scale_color_manual(values = c("Significantly increasing" = "#E63946", 
                                  "Significantly decreasing" = "#2A9D8F")) +  # Custom colors for trend segments
    labs(title = plot_titles[[compound]], x = "Year", y = expression("f ' (x)"), color = "") +
    theme_minimal(base_size = 14) +  
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  print(p0)
  
  # Plot 2: First derivative of the GAM with δ15N correction
  p1 <- ggplot(deriv_values2, aes(x = Year, y = .derivative)) +
    geom_line(color = "#C5A8D1", linewidth = 0.8) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#C5A8D1", alpha = 0.3) +
    geom_line(data = deriv_segments2, aes(y = .derivative, group = segment_id, color = trend), linewidth = 1) +
    scale_color_manual(values = c("Significantly increasing" = "#E63946", 
                                  "Significantly decreasing" = "#2A9D8F")) +
    labs(title = plot_titles[[compound]], x = "Year", y = expression("f ' (x)"), color = "") +
    theme_minimal(base_size = 14) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  print(p1)
  
  
  # Plot 3: Raw and smoothed concentrations over time with derivative trend segments and critical points
  p2 <- ggplot(data_clean, aes(x = Year, y = .data[[compound]])) +
    geom_point(alpha = 0.5, color = "#4682B4") +
    geom_smooth(method = "gam", formula = y ~ s(x), 
                aes(fill = "Without"), color = "#4682B4", alpha = 0.2, level = 0.95) +  
    geom_line(data = deriv_segments,   
              aes(x = Year,   
                  y = predict(gam_model, newdata = deriv_segments),   
                  color = trend,   
                  group = segment_id),   
              linewidth = 1.5)
  
  # Add points and smoothing for δ15N-corrected data
  p2 <- p2 +
    geom_point(data = data_clean2, alpha = 0.5, color = "#C5A8D1") +
    geom_smooth(data = data_clean2, method = "gam", formula = y ~ s(x), 
                aes(fill = "With"), color = "#C5A8D1", alpha = 0.3, level = 0.95) +  
    geom_line(data = deriv_segments2,   
              aes(x = Year,   
                  y = predict(gam_model2, newdata = deriv_segments2),   
                  color = trend,   
                  group = segment_id),   
              linewidth = 1.5) 
  
  # Add GAM1 and GAM2 change points with horizontal confidence intervals
  p2 <- p2 + 
    geom_point(
      aes(x = year_max_gam, 
          y = predict(gam_model, newdata = data.frame(Year = year_max_gam)), 
          color = "GAM1 change point"), 
      size = 3, shape = 21, fill = "#3A6A94"
    ) +
    geom_errorbarh(
      aes(xmin = ci_range$min_year, xmax = ci_range$max_year, 
          y = predict(gam_model, newdata = data.frame(Year = year_max_gam))), 
      height = 3, size = 0.6, color = "#3A6A94"
    ) +
    geom_point(
      aes(x = year_max_gam2, 
          y = predict(gam_model2, newdata = data.frame(Year = year_max_gam2)), 
          color = "GAM2 change point"), 
      size = 3, shape = 21, fill = "#5A3D77"
    ) +
    geom_errorbarh(
      aes(xmin = ci_range2$min_year2, xmax = ci_range2$max_year2, 
          y = predict(gam_model2, newdata = data.frame(Year = year_max_gam2))), 
      height = 3, size = 0.6, color = "#5A3D77"
    )
  
  # Add regulatory and historical periods
  p2 <- p2 + 
    geom_errorbarh(aes(xmin = 2002, xmax = 2024, 
                       y = 25, 
                       color = "Regulatory period"), 
                   height = 3, size = 0.6)
  
  # Phase-out PFOS (2000-2003)
  p2 <- p2 + 
    geom_errorbarh(aes(xmin = 2000, xmax = 2003, 
                       y = 40, 
                       color = "Phase-out of PFOS"), 
                   height = 3, size = 0.6)
  
  # Add LOAEL threshold for PFOS
  p2 <- p2 +
    geom_hline(aes(yintercept = 92.4, color = "LOAEL PFOS"), linetype = "dashed", linewidth = 0.8)
  
  # Customize colors and legends
  p2 <- p2 + 
    scale_color_manual(
      name = "", 
      values = c(
        "Significantly increasing" = "#E63946",  
        "Significantly decreasing" = "#2A9D8F",
        "Regulatory period" = "#66C2B1",
        "Phase-out of PFOS" = "darkred",
        "LOAEL PFOS" = "#333333"
      ),
      labels = c(
        "LOAEL PFOS" = "LOAEL for PFOS individual exposure",
        "Phase-out of PFOS" = "Phase-out of PFOS"
      )
    ) +
    scale_fill_manual(
      name = "", 
      values = c("Without" = "#4682B4", "With" = "#C5A8D1"),
      labels = c(
        "Without" = expression("Without " * delta^15 * "N correction"),
        "With" = expression("With " * delta^15 * "N correction")
      )
    )
  
  # Improve theme and axis aesthetics
  p2 <- p2 +
    theme_minimal(base_size = 14) +  
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    labs(
      title = plot_titles[[compound]],
      x = "Year",
      y = "Lipid-standardized concentration (ng/g ww)"
    )
  
  # Print plot
  print(p2)
  
}


### PFOS + PFHxS
compounds <- "PFOS_PFHxS"

for (compound in compounds) {
  cat("\n Analysis for:", compound, "\n")
  
  # Remove missing values
  data_clean <- na.omit(data[, c("Year", compound)])
  
  # Fit a GAM model
  gam_model <- gam(as.formula(paste(compound, "~ s(Year, k=10)")), data = data_clean, method = "REML")
  
  # Compute derivatives and 95% confidence intervals
  deriv_values <- derivatives(gam_model, term = "s(Year)", type = "central", level = 0.95) %>%
    select(Year, .derivative, .lower_ci, .upper_ci) %>%
    mutate(
      trend = case_when(
        .lower_ci > 0 ~ "Significantly increasing",
        .upper_ci < 0 ~ "Significantly decreasing",
        TRUE ~ "Neutral"
      )
    )
  
  # Detect local maximum: transition from increasing to decreasing trend
  critical_points <- deriv_values %>%
    mutate(sign_change = lag(.derivative) > 0 & .derivative < 0) %>%
    filter(sign_change)
  
  if (nrow(critical_points) > 0) {
    
    # Select the year with the highest value of the compound
    max_year_data <- data_clean %>%
      filter(Year %in% critical_points$Year) %>%
      arrange(desc(!!sym(compound))) %>%
      slice(1)
    
    year_max <- max_year_data$Year
    
    # Estimate confidence interval around the GAM maximum using predict()
    pred_data <- data.frame(Year = seq(min(data_clean$Year), max(data_clean$Year), length.out = 500))
    pred_data$fit <- predict(gam_model, newdata = pred_data, se.fit = TRUE)$fit
    pred_data$se <- predict(gam_model, newdata = pred_data, se.fit = TRUE)$se.fit
    
    # Year corresponding to the maximum predicted GAM value
    year_max_gam <- pred_data$Year[which.max(pred_data$fit)]
    
    # Confidence interval for the GAM maximum
    ci_range <- pred_data %>%
      filter(fit >= max(fit) - 1.96 * se) %>%
      summarize(min_year = min(Year), max_year = max(Year))
    
    cat("\nLipid corrected GAM: detected maximum for", compound, "in year:", round(year_max_gam, 2), "\n")
    cat("95% confidence interval for the year of maximum: [", round(ci_range$min_year, 2), ", ", round(ci_range$max_year, 2), "]\n")
  }
  
  # Get predicted concentration and 95% CI at the year of the maximum
  new_data <- data.frame(Year = year_max_gam)
  pred <- predict(gam_model, newdata = new_data, se.fit = TRUE)
  
  fit_value <- pred$fit
  se_value <- pred$se.fit
  lower_ci <- fit_value - 1.96 * se_value
  upper_ci <- fit_value + 1.96 * se_value
  
  cat("Lipid corrected GAM: concentration at year of maximum (", round(year_max_gam, 2), "): ", round(fit_value, 2), "ng/g ww\n")
  cat("95% confidence interval: [", round(lower_ci, 2), ", ", round(upper_ci, 2), "]\n")
  
  
  # Segments distinguishing increasing and decreasing trends
  deriv_segments <- deriv_values %>%
    arrange(Year) %>%
    mutate(
      is_neutral = (trend == "Neutral"),
      trend_change = (trend != lag(trend, default = first(trend))),
      segment_id = cumsum(trend_change | is_neutral)
    ) %>%
    filter(trend != "Neutral")
  
  # Remove missing values in corrected dataset
  data_clean2 <- na.omit(data_corr[, c("Year", compound)])
  
  # Fit corrected GAM model
  gam_model2 <- gam(as.formula(paste(compound, "~ s(Year, k=10)")), data = data_clean2, method = "REML")
  
  # Compute derivatives and 95% CI
  deriv_values2 <- derivatives(gam_model2, term = "s(Year)", type = "central", level = 0.95) %>%
    select(Year, .derivative, .lower_ci, .upper_ci) %>%
    mutate(
      trend = case_when(
        .lower_ci > 0 ~ "Significantly increasing",
        .upper_ci < 0 ~ "Significantly decreasing",
        TRUE ~ "Neutral"
      )
    )
  
  # Detect local maximum 
  critical_points2 <- deriv_values2 %>%
    mutate(sign_change = lag(.derivative) > 0 & .derivative < 0) %>%
    filter(sign_change)
  
  if (nrow(critical_points2) > 0) {
    
    max_year_data2 <- data_clean2 %>%
      filter(Year %in% critical_points2$Year) %>%
      arrange(desc(!!sym(compound))) %>%
      slice(1)
    
    year_max2 <- max_year_data2$Year
    
    pred_data2 <- data.frame(Year = seq(min(data_clean2$Year), max(data_clean2$Year), length.out = 500))
    pred_data2$fit <- predict(gam_model2, newdata = pred_data2, se.fit = TRUE)$fit
    pred_data2$se <- predict(gam_model2, newdata = pred_data2, se.fit = TRUE)$se.fit
    
    year_max_gam2 <- pred_data2$Year[which(pred_data2$fit == max(pred_data2$fit, na.rm = TRUE))]
    
    ci_range2 <- pred_data2 %>%
      filter(fit >= max(fit) - 1.96 * se) %>%
      summarize(min_year2 = min(Year), max_year2 = max(Year))
    
    cat("\nLipid-d15N corrected GAM: detected maximum for", compound, "in year:", round(year_max_gam2, 2), "\n")
    cat("95% CI for year of maximum: [", round(ci_range2$min_year2, 2), ", ", round(ci_range2$max_year2, 2), "]\n")
  }
  
  # Predicted value and CI at the maximum year
  new_data2 <- data.frame(Year = year_max_gam2)
  pred2 <- predict(gam_model2, newdata = new_data2, se.fit = TRUE)
  
  fit_value2 <- pred2$fit
  se_value2 <- pred2$se.fit
  lower_ci2 <- fit_value2 - 1.96 * se_value2
  upper_ci2 <- fit_value2 + 1.96 * se_value2
  
  cat("Lipid-d15N corrected GAM: concentration at year of maximum (", round(year_max_gam2, 2), "): ", round(fit_value2, 2), "ng/g ww\n")
  cat("95% confidence interval: [", round(lower_ci2, 2), ", ", round(upper_ci2, 2), "]\n")
  
  # Compute start and end of each significant trend period
  trend_periods <- deriv_segments %>%
    group_by(segment_id, trend) %>%
    summarise(
      start_year = min(Year),
      end_year = max(Year),
      start_value = predict(gam_model, newdata = data.frame(Year = min(Year))),
      end_value = predict(gam_model, newdata = data.frame(Year = max(Year)))
    ) %>%
    ungroup()
  
  cat("\nSignificant trend periods (Lipid corrected model) for", compound, ":\n")
  print(trend_periods)
  
  # Same for corrected model
  deriv_segments2 <- deriv_values2 %>%
    arrange(Year) %>%
    mutate(
      is_neutral = (trend == "Neutral"),
      trend_change = (trend != lag(trend, default = first(trend))),
      segment_id = cumsum(trend_change | is_neutral)
    ) %>%
    filter(trend != "Neutral")
  
  trend_periods2 <- deriv_segments2 %>%
    group_by(segment_id, trend) %>%
    summarise(
      start_year = min(Year),
      end_year = max(Year),
      start_value = predict(gam_model2, newdata = data.frame(Year = min(Year))),
      end_value = predict(gam_model2, newdata = data.frame(Year = max(Year)))
    ) %>%
    ungroup()
  
  cat("\nSignificant trend periods (Lipid-d15N corrected model) for", compound, ":\n")
  print(trend_periods2)
  
  # Plot 1: First derivative of the GAM without δ15N correction
  p3 <- ggplot(deriv_values, aes(x = Year, y = .derivative)) +
    geom_line(color = "#4682B4", linewidth = 0.8) +  # Main derivative curve
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#4682B4", alpha = 0.2) +  # Confidence interval
    geom_line(data = deriv_segments, aes(y = .derivative, group = segment_id, color = trend), linewidth = 1) +  # Significant trends
    scale_color_manual(values = c("Significantly increasing" = "#E63946", 
                                  "Significantly decreasing" = "#2A9D8F")) +  # Custom colors for trend segments
    labs(title = plot_titles[[compound]], x = "Year", y = expression("f ' (x)"), color = "") +
    theme_minimal(base_size = 14) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")  
    )
  
  print(p3)
  
  # Plot 2: First derivative of the GAM with δ15N correction
  p4 <- ggplot(deriv_values2, aes(x = Year, y = .derivative)) +
    geom_line(color = "#C5A8D1", linewidth = 0.8) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#C5A8D1", alpha = 0.3) +
    geom_line(data = deriv_segments2, aes(y = .derivative, group = segment_id, color = trend), linewidth = 1) +
    scale_color_manual(values = c("Significantly increasing" = "#E63946", 
                                  "Significantly decreasing" = "#2A9D8F")) +
    labs(title = plot_titles[[compound]], x = "Year", y = expression("f ' (x)"), color = "") +
    theme_minimal(base_size = 14) +  
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")  
    )
  
  print(p4)
  
  # Plot 3: Raw and smoothed concentrations over time with derivative trend segments and critical points
  p5 <- ggplot(data_clean, aes(x = Year, y = .data[[compound]])) +
    geom_point(alpha = 0.5, color = "#4682B4") +
    geom_smooth(method = "gam", formula = y ~ s(x), 
                aes(fill = "Without"), color = "#4682B4", alpha = 0.2, level = 0.95) +  
    geom_line(data = deriv_segments,   
              aes(x = Year,   
                  y = predict(gam_model, newdata = deriv_segments),   
                  color = trend,   
                  group = segment_id),   
              linewidth = 1.5)
  
  # Add points and smoothing for δ15N-corrected data
  p5 <- p5 +
    geom_point(data = data_clean2, alpha = 0.5, color = "#C5A8D1") +
    geom_smooth(data = data_clean2, method = "gam", formula = y ~ s(x), 
                aes(fill = "With"), color = "#C5A8D1", alpha = 0.3, level = 0.95) +  
    geom_line(data = deriv_segments2,   
              aes(x = Year,   
                  y = predict(gam_model2, newdata = deriv_segments2),   
                  color = trend,   
                  group = segment_id),   
              linewidth = 1.5) 
  
  # Add GAM1 and GAM2 change points with horizontal confidence intervals
  p5 <- p5 + 
    geom_point(
      aes(x = year_max_gam, 
          y = predict(gam_model, newdata = data.frame(Year = year_max_gam)), 
          color = "GAM1 change point"), 
      size = 2, shape = 21, fill = "#3A6A94"
    ) +
    geom_errorbarh(
      aes(xmin = ci_range$min_year, xmax = ci_range$max_year, 
          y = predict(gam_model, newdata = data.frame(Year = year_max_gam))), 
      height = 2, size = 0.6, color = "#3A6A94"
    ) +
    geom_point(
      aes(x = year_max_gam2, 
          y = predict(gam_model2, newdata = data.frame(Year = year_max_gam2)), 
          color = "GAM2 change point"), 
      size = 2, shape = 21, fill = "#5A3D77"
    ) +
    geom_errorbarh(
      aes(xmin = ci_range2$min_year2, xmax = ci_range2$max_year2, 
          y = predict(gam_model2, newdata = data.frame(Year = year_max_gam2))), 
      height = 2, size = 0.6, color = "#5A3D77"
    )
  
    # Add historical and regulatory context
  p5 <- p5 +
    geom_errorbarh(aes(xmin = 2002, xmax = 2024, y = 20, color = "Regulatory period"), height = 2, size = 0.6) +
    geom_errorbarh(aes(xmin = 2000, xmax = 2003, y = 35, color = "Phase-out of PFOS"), height = 2, size = 0.6)
    
    # Add LOAEL thresholds
  p5 <- p5 +
    geom_hline(aes(yintercept = 81.3, linetype = "LOAEL PFOS:PFHxS"), color = "#333333", linewidth = 1.2)
  
    # Customize colors and labels
  p5 <- p5 +
    scale_color_manual(
      name = "", 
      values = c(
        "Significantly increasing" = "#E63946",
        "Significantly decreasing" = "#2A9D8F",
        "Regulatory period" = "#66C2B1",
        "Phase-out of PFOS" = "darkred"
        ),
      labels = c(
        "Significantly increasing" = "Significantly increasing",
        "Significantly decreasing" = "Significantly decreasing",
        "Regulatory period" = "Regulatory period",
        "Phase-out of PFOS" = "Phase-out of PFOS"
        )
    ) +
  
    scale_linetype_manual(
      name = "",
      values = c(
        "LOAEL PFOS:PFHxS" = "dotted"
      ),
      labels = c(
        "LOAEL PFOS:PFHxS" = "LOAEL for PFOS:PFHxS mixture"
      )
    ) +
    
    scale_fill_manual(
      name = "",
      values = c("Without" = "#4682B4", "With" = "#C5A8D1"),
      labels = c("Without" = expression("Without " * delta^15 * "N correction"),
                 "With" = expression("With " * delta^15 * "N correction"))
    ) +
    
    # Aesthetic theme and axis/plot title
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    labs(
      title = plot_titles[[compound]],
      x = "Year",
      y = "Lipid-standardized concentration (ng/g ww)"
    )
  
  # Display plot
  print(p5)
  
}



### PFAS
compounds <- "PFAS"

for (compound in compounds) {
  cat("\n Analysis for:", compound, "\n")
  
  # Remove missing values
  data_clean <- na.omit(data[, c("Year", compound)])
  
  # Fit a GAM model
  gam_model <- gam(as.formula(paste(compound, "~ s(Year, k=10)")), data = data_clean, method = "REML")
  
  # Compute derivatives and 95% confidence intervals
  deriv_values <- derivatives(gam_model, term = "s(Year)", type = "central", level = 0.95) %>%
    select(Year, .derivative, .lower_ci, .upper_ci) %>%
    mutate(
      trend = case_when(
        .lower_ci > 0 ~ "Significantly increasing",
        .upper_ci < 0 ~ "Significantly decreasing",
        TRUE ~ "Neutral"
      )
    )
  
  # Detect local maximum: transition from increasing to decreasing trend
  critical_points <- deriv_values %>%
    mutate(sign_change = lag(.derivative) > 0 & .derivative < 0) %>%
    filter(sign_change)
  
  if (nrow(critical_points) > 0) {
    
    # Select the year with the highest value of the compound
    max_year_data <- data_clean %>%
      filter(Year %in% critical_points$Year) %>%
      arrange(desc(!!sym(compound))) %>%
      slice(1)
    
    year_max <- max_year_data$Year
    
    # Estimate confidence interval around the GAM maximum using predict()
    pred_data <- data.frame(Year = seq(min(data_clean$Year), max(data_clean$Year), length.out = 500))
    pred_data$fit <- predict(gam_model, newdata = pred_data, se.fit = TRUE)$fit
    pred_data$se <- predict(gam_model, newdata = pred_data, se.fit = TRUE)$se.fit
    
    # Year corresponding to the maximum predicted GAM value
    year_max_gam <- pred_data$Year[which.max(pred_data$fit)]
    
    # Confidence interval for the GAM maximum
    ci_range <- pred_data %>%
      filter(fit >= max(fit) - 1.96 * se) %>%
      summarize(min_year = min(Year), max_year = max(Year))
    
    cat("\nLipid corrected GAM: detected maximum for", compound, "in year:", round(year_max_gam, 2), "\n")
    cat("95% confidence interval for the year of maximum: [", round(ci_range$min_year, 2), ", ", round(ci_range$max_year, 2), "]\n")
  }
  
  # Get predicted concentration and 95% CI at the year of the maximum
  new_data <- data.frame(Year = year_max_gam)
  pred <- predict(gam_model, newdata = new_data, se.fit = TRUE)
  
  fit_value <- pred$fit
  se_value <- pred$se.fit
  lower_ci <- fit_value - 1.96 * se_value
  upper_ci <- fit_value + 1.96 * se_value
  
  cat("Lipid corrected GAM: concentration at year of maximum (", round(year_max_gam, 2), "): ", round(fit_value, 2), "ng/g ww\n")
  cat("95% confidence interval: [", round(lower_ci, 2), ", ", round(upper_ci, 2), "]\n")
  
  
  # Segments distinguishing increasing and decreasing trends
  deriv_segments <- deriv_values %>%
    arrange(Year) %>%
    mutate(
      is_neutral = (trend == "Neutral"),
      trend_change = (trend != lag(trend, default = first(trend))),
      segment_id = cumsum(trend_change | is_neutral)
    ) %>%
    filter(trend != "Neutral")
  
  # Remove missing values in corrected dataset
  data_clean2 <- na.omit(data_corr[, c("Year", compound)])
  
  # Fit corrected GAM model
  gam_model2 <- gam(as.formula(paste(compound, "~ s(Year, k=10)")), data = data_clean2, method = "REML")
  
  # Compute derivatives and 95% CI
  deriv_values2 <- derivatives(gam_model2, term = "s(Year)", type = "central", level = 0.95) %>%
    select(Year, .derivative, .lower_ci, .upper_ci) %>%
    mutate(
      trend = case_when(
        .lower_ci > 0 ~ "Significantly increasing",
        .upper_ci < 0 ~ "Significantly decreasing",
        TRUE ~ "Neutral"
      )
    )
  
  # Detect local maximum 
  critical_points2 <- deriv_values2 %>%
    mutate(sign_change = lag(.derivative) > 0 & .derivative < 0) %>%
    filter(sign_change)
  
  if (nrow(critical_points2) > 0) {
    
    max_year_data2 <- data_clean2 %>%
      filter(Year %in% critical_points2$Year) %>%
      arrange(desc(!!sym(compound))) %>%
      slice(1)
    
    year_max2 <- max_year_data2$Year
    
    pred_data2 <- data.frame(Year = seq(min(data_clean2$Year), max(data_clean2$Year), length.out = 500))
    pred_data2$fit <- predict(gam_model2, newdata = pred_data2, se.fit = TRUE)$fit
    pred_data2$se <- predict(gam_model2, newdata = pred_data2, se.fit = TRUE)$se.fit
    
    year_max_gam2 <- pred_data2$Year[which(pred_data2$fit == max(pred_data2$fit, na.rm = TRUE))]
    
    ci_range2 <- pred_data2 %>%
      filter(fit >= max(fit) - 1.96 * se) %>%
      summarize(min_year2 = min(Year), max_year2 = max(Year))
    
    cat("\nLipid-d15N corrected GAM: detected maximum for", compound, "in year:", round(year_max_gam2, 2), "\n")
    cat("95% CI for year of maximum: [", round(ci_range2$min_year2, 2), ", ", round(ci_range2$max_year2, 2), "]\n")
  }
  
  # Predicted value and CI at the maximum year
  new_data2 <- data.frame(Year = year_max_gam2)
  pred2 <- predict(gam_model2, newdata = new_data2, se.fit = TRUE)
  
  fit_value2 <- pred2$fit
  se_value2 <- pred2$se.fit
  lower_ci2 <- fit_value2 - 1.96 * se_value2
  upper_ci2 <- fit_value2 + 1.96 * se_value2
  
  cat("Lipid-d15N corrected GAM: concentration at year of maximum (", round(year_max_gam2, 2), "): ", round(fit_value2, 2), "ng/g ww\n")
  cat("95% confidence interval: [", round(lower_ci2, 2), ", ", round(upper_ci2, 2), "]\n")
  
  # Compute start and end of each significant trend period
  trend_periods <- deriv_segments %>%
    group_by(segment_id, trend) %>%
    summarise(
      start_year = min(Year),
      end_year = max(Year),
      start_value = predict(gam_model, newdata = data.frame(Year = min(Year))),
      end_value = predict(gam_model, newdata = data.frame(Year = max(Year)))
    ) %>%
    ungroup()
  
  cat("\nSignificant trend periods (Lipid corrected model) for", compound, ":\n")
  print(trend_periods)
  
  # Same for corrected model
  deriv_segments2 <- deriv_values2 %>%
    arrange(Year) %>%
    mutate(
      is_neutral = (trend == "Neutral"),
      trend_change = (trend != lag(trend, default = first(trend))),
      segment_id = cumsum(trend_change | is_neutral)
    ) %>%
    filter(trend != "Neutral")
  
  trend_periods2 <- deriv_segments2 %>%
    group_by(segment_id, trend) %>%
    summarise(
      start_year = min(Year),
      end_year = max(Year),
      start_value = predict(gam_model2, newdata = data.frame(Year = min(Year))),
      end_value = predict(gam_model2, newdata = data.frame(Year = max(Year)))
    ) %>%
    ungroup()
  
  cat("\nSignificant trend periods (Lipid-d15N corrected model) for", compound, ":\n")
  print(trend_periods2)
  
  # Plot 1: First derivative of the GAM without δ15N correction
  p3 <- ggplot(deriv_values, aes(x = Year, y = .derivative)) +
    geom_line(color = "#4682B4", linewidth = 0.8) +  # Main derivative curve
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#4682B4", alpha = 0.2) +  # Confidence interval
    geom_line(data = deriv_segments, aes(y = .derivative, group = segment_id, color = trend), linewidth = 1) +  # Significant trends
    scale_color_manual(values = c("Significantly increasing" = "#E63946", 
                                  "Significantly decreasing" = "#2A9D8F")) +  # Custom colors for trend segments
    labs(title = plot_titles[[compound]], x = "Year", y = expression("f ' (x)"), color = "") +
    theme_minimal(base_size = 14) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")  
    )
  
  print(p3)
  
  # Plot 2: First derivative of the GAM with δ15N correction
  p4 <- ggplot(deriv_values2, aes(x = Year, y = .derivative)) +
    geom_line(color = "#C5A8D1", linewidth = 0.8) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#C5A8D1", alpha = 0.3) +
    geom_line(data = deriv_segments2, aes(y = .derivative, group = segment_id, color = trend), linewidth = 1) +
    scale_color_manual(values = c("Significantly increasing" = "#E63946", 
                                  "Significantly decreasing" = "#2A9D8F")) +
    labs(title = plot_titles[[compound]], x = "Year", y = expression("f ' (x)"), color = "") +
    theme_minimal(base_size = 14) +  
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")  
    )
  
  print(p4)
  
  # Plot 3: Raw and smoothed concentrations over time with derivative trend segments and critical points
  p5 <- ggplot(data_clean, aes(x = Year, y = .data[[compound]])) +
    geom_point(alpha = 0.5, color = "#4682B4") +
    geom_smooth(method = "gam", formula = y ~ s(x), 
                aes(fill = "Without"), color = "#4682B4", alpha = 0.2, level = 0.95) +  
    geom_line(data = deriv_segments,   
              aes(x = Year,   
                  y = predict(gam_model, newdata = deriv_segments),   
                  color = trend,   
                  group = segment_id),   
              linewidth = 1.5)
  
  # Add points and smoothing for δ15N-corrected data
  p5 <- p5 +
    geom_point(data = data_clean2, alpha = 0.5, color = "#C5A8D1") +
    geom_smooth(data = data_clean2, method = "gam", formula = y ~ s(x), 
                aes(fill = "With"), color = "#C5A8D1", alpha = 0.3, level = 0.95) +  
    geom_line(data = deriv_segments2,   
              aes(x = Year,   
                  y = predict(gam_model2, newdata = deriv_segments2),   
                  color = trend,   
                  group = segment_id),   
              linewidth = 1.5) 
  
  # Add GAM1 and GAM2 change points with horizontal confidence intervals
  p5 <- p5 + 
    geom_point(
      aes(x = year_max_gam, 
          y = predict(gam_model, newdata = data.frame(Year = year_max_gam)), 
          color = "GAM1 change point"), 
      size = 3, shape = 21, fill = "#3A6A94"
    ) +
    geom_errorbarh(
      aes(xmin = ci_range$min_year, xmax = ci_range$max_year, 
          y = predict(gam_model, newdata = data.frame(Year = year_max_gam))), 
      height = 3, size = 0.6, color = "#3A6A94"
    ) +
    geom_point(
      aes(x = year_max_gam2, 
          y = predict(gam_model2, newdata = data.frame(Year = year_max_gam2)), 
          color = "GAM2 change point"), 
      size = 3, shape = 21, fill = "#5A3D77"
    ) +
    geom_errorbarh(
      aes(xmin = ci_range2$min_year2, xmax = ci_range2$max_year2, 
          y = predict(gam_model2, newdata = data.frame(Year = year_max_gam2))), 
      height = 3, size = 0.6, color = "#5A3D77"
    )
  
  # Add historical and regulatory context
  p5 <- p5 +
    geom_errorbarh(aes(xmin = 2002, xmax = 2024, y = 50, color = "Regulatory period"), height = 3, size = 0.6) +
    geom_errorbarh(aes(xmin = 2010, xmax = 2015, y = 57, color = "Phase-out of PFOA and LC-PFCAs"), height = 3, size = 0.6) +
    geom_errorbarh(aes(xmin = 2000, xmax = 2003, y = 65, color = "Phase-out of PFOS"), height = 3, size = 0.6)
  
  
  # Customize colors and labels
  p5 <- p5 +
    scale_color_manual(
      name = "", 
      values = c(
        "Significantly increasing" = "#E63946",
        "Significantly decreasing" = "#2A9D8F",
        "Regulatory period" = "#66C2B1",
        "Phase-out of PFOS" = "darkred",
        "Phase-out of PFOA and LC-PFCAs" = "darkorange"
      ),
      labels = c(
        "Significantly increasing" = "Significantly increasing",
        "Significantly decreasing" = "Significantly decreasing",
        "Regulatory period" = "Regulatory period",
        "Phase-out of PFOS" = "Phase-out of PFOS",
        "Phase-out of PFOA and LC-PFCAs" = "Phase-out of PFOA and LC-PFCAs"
      )
    ) +
    
    scale_fill_manual(
      name = "",
      values = c("Without" = "#4682B4", "With" = "#C5A8D1"),
      labels = c("Without" = expression("Without " * delta^15 * "N correction"),
                 "With" = expression("With " * delta^15 * "N correction"))
    ) +
    
    # Aesthetic theme and axis/plot title
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    labs(
      title = plot_titles[[compound]],
      x = "Year",
      y = "Lipid-standardized concentration (ng/g ww)"
    )
  
  # Display plot
  print(p5)
  
}


### PFSA
compounds <- "PFSA"

for (compound in compounds) {
  cat("\n Analysis for:", compound, "\n")
  
  # Remove missing values
  data_clean <- na.omit(data[, c("Year", compound)])
  
  # Fit a GAM model
  gam_model <- gam(as.formula(paste(compound, "~ s(Year, k=10)")), data = data_clean, method = "REML")
  
  # Compute derivatives and 95% confidence intervals
  deriv_values <- derivatives(gam_model, term = "s(Year)", type = "central", level = 0.95) %>%
    select(Year, .derivative, .lower_ci, .upper_ci) %>%
    mutate(
      trend = case_when(
        .lower_ci > 0 ~ "Significantly increasing",
        .upper_ci < 0 ~ "Significantly decreasing",
        TRUE ~ "Neutral"
      )
    )
  
  # Detect local maximum: transition from increasing to decreasing trend
  critical_points <- deriv_values %>%
    mutate(sign_change = lag(.derivative) > 0 & .derivative < 0) %>%
    filter(sign_change)
  
  if (nrow(critical_points) > 0) {
    
    # Select the year with the highest value of the compound
    max_year_data <- data_clean %>%
      filter(Year %in% critical_points$Year) %>%
      arrange(desc(!!sym(compound))) %>%
      slice(1)
    
    year_max <- max_year_data$Year
    
    # Estimate confidence interval around the GAM maximum using predict()
    pred_data <- data.frame(Year = seq(min(data_clean$Year), max(data_clean$Year), length.out = 500))
    pred_data$fit <- predict(gam_model, newdata = pred_data, se.fit = TRUE)$fit
    pred_data$se <- predict(gam_model, newdata = pred_data, se.fit = TRUE)$se.fit
    
    # Year corresponding to the maximum predicted GAM value
    year_max_gam <- pred_data$Year[which.max(pred_data$fit)]
    
    # Confidence interval for the GAM maximum
    ci_range <- pred_data %>%
      filter(fit >= max(fit) - 1.96 * se) %>%
      summarize(min_year = min(Year), max_year = max(Year))
    
    cat("\nLipid corrected GAM: detected maximum for", compound, "in year:", round(year_max_gam, 2), "\n")
    cat("95% confidence interval for the year of maximum: [", round(ci_range$min_year, 2), ", ", round(ci_range$max_year, 2), "]\n")
  }
  
  # Get predicted concentration and 95% CI at the year of the maximum
  new_data <- data.frame(Year = year_max_gam)
  pred <- predict(gam_model, newdata = new_data, se.fit = TRUE)
  
  fit_value <- pred$fit
  se_value <- pred$se.fit
  lower_ci <- fit_value - 1.96 * se_value
  upper_ci <- fit_value + 1.96 * se_value
  
  cat("Lipid corrected GAM: concentration at year of maximum (", round(year_max_gam, 2), "): ", round(fit_value, 2), "ng/g ww\n")
  cat("95% confidence interval: [", round(lower_ci, 2), ", ", round(upper_ci, 2), "]\n")
  
  
  # Segments distinguishing increasing and decreasing trends
  deriv_segments <- deriv_values %>%
    arrange(Year) %>%
    mutate(
      is_neutral = (trend == "Neutral"),
      trend_change = (trend != lag(trend, default = first(trend))),
      segment_id = cumsum(trend_change | is_neutral)
    ) %>%
    filter(trend != "Neutral")
  
  # Remove missing values in corrected dataset
  data_clean2 <- na.omit(data_corr[, c("Year", compound)])
  
  # Fit corrected GAM model
  gam_model2 <- gam(as.formula(paste(compound, "~ s(Year, k=10)")), data = data_clean2, method = "REML")
  
  # Compute derivatives and 95% CI
  deriv_values2 <- derivatives(gam_model2, term = "s(Year)", type = "central", level = 0.95) %>%
    select(Year, .derivative, .lower_ci, .upper_ci) %>%
    mutate(
      trend = case_when(
        .lower_ci > 0 ~ "Significantly increasing",
        .upper_ci < 0 ~ "Significantly decreasing",
        TRUE ~ "Neutral"
      )
    )
  
  # Detect local maximum 
  critical_points2 <- deriv_values2 %>%
    mutate(sign_change = lag(.derivative) > 0 & .derivative < 0) %>%
    filter(sign_change)
  
  if (nrow(critical_points2) > 0) {
    
    max_year_data2 <- data_clean2 %>%
      filter(Year %in% critical_points2$Year) %>%
      arrange(desc(!!sym(compound))) %>%
      slice(1)
    
    year_max2 <- max_year_data2$Year
    
    pred_data2 <- data.frame(Year = seq(min(data_clean2$Year), max(data_clean2$Year), length.out = 500))
    pred_data2$fit <- predict(gam_model2, newdata = pred_data2, se.fit = TRUE)$fit
    pred_data2$se <- predict(gam_model2, newdata = pred_data2, se.fit = TRUE)$se.fit
    
    year_max_gam2 <- pred_data2$Year[which(pred_data2$fit == max(pred_data2$fit, na.rm = TRUE))]
    
    ci_range2 <- pred_data2 %>%
      filter(fit >= max(fit) - 1.96 * se) %>%
      summarize(min_year2 = min(Year), max_year2 = max(Year))
    
    cat("\nLipid-d15N corrected GAM: detected maximum for", compound, "in year:", round(year_max_gam2, 2), "\n")
    cat("95% CI for year of maximum: [", round(ci_range2$min_year2, 2), ", ", round(ci_range2$max_year2, 2), "]\n")
  }
  
  # Predicted value and CI at the maximum year
  new_data2 <- data.frame(Year = year_max_gam2)
  pred2 <- predict(gam_model2, newdata = new_data2, se.fit = TRUE)
  
  fit_value2 <- pred2$fit
  se_value2 <- pred2$se.fit
  lower_ci2 <- fit_value2 - 1.96 * se_value2
  upper_ci2 <- fit_value2 + 1.96 * se_value2
  
  cat("Lipid-d15N corrected GAM: concentration at year of maximum (", round(year_max_gam2, 2), "): ", round(fit_value2, 2), "ng/g ww\n")
  cat("95% confidence interval: [", round(lower_ci2, 2), ", ", round(upper_ci2, 2), "]\n")
  
  # Compute start and end of each significant trend period
  trend_periods <- deriv_segments %>%
    group_by(segment_id, trend) %>%
    summarise(
      start_year = min(Year),
      end_year = max(Year),
      start_value = predict(gam_model, newdata = data.frame(Year = min(Year))),
      end_value = predict(gam_model, newdata = data.frame(Year = max(Year)))
    ) %>%
    ungroup()
  
  cat("\nSignificant trend periods (Lipid corrected model) for", compound, ":\n")
  print(trend_periods)
  
  # Same for corrected model
  deriv_segments2 <- deriv_values2 %>%
    arrange(Year) %>%
    mutate(
      is_neutral = (trend == "Neutral"),
      trend_change = (trend != lag(trend, default = first(trend))),
      segment_id = cumsum(trend_change | is_neutral)
    ) %>%
    filter(trend != "Neutral")
  
  trend_periods2 <- deriv_segments2 %>%
    group_by(segment_id, trend) %>%
    summarise(
      start_year = min(Year),
      end_year = max(Year),
      start_value = predict(gam_model2, newdata = data.frame(Year = min(Year))),
      end_value = predict(gam_model2, newdata = data.frame(Year = max(Year)))
    ) %>%
    ungroup()
  
  cat("\nSignificant trend periods (Lipid-d15N corrected model) for", compound, ":\n")
  print(trend_periods2)
  
  # Plot 1: First derivative of the GAM without δ15N correction
  p6 <- ggplot(deriv_values, aes(x = Year, y = .derivative)) +
    geom_line(color = "#4682B4", linewidth = 0.8) +  # Main derivative curve
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#4682B4", alpha = 0.2) +  # Confidence interval
    geom_line(data = deriv_segments, aes(y = .derivative, group = segment_id, color = trend), linewidth = 1) +  # Significant trends
    scale_color_manual(values = c("Significantly increasing" = "#E63946", 
                                  "Significantly decreasing" = "#2A9D8F")) +  # Custom colors for trend segments
    labs(title = plot_titles[[compound]], x = "Year", y = expression("f ' (x)"), color = "") +
    theme_minimal(base_size = 14) +  
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  
  
  print(p6)
  
  # Plot 2: First derivative of the GAM with δ15N correction
  p7 <- ggplot(deriv_values2, aes(x = Year, y = .derivative)) +
    geom_line(color = "#C5A8D1", linewidth = 0.8) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#C5A8D1", alpha = 0.3) +
    geom_line(data = deriv_segments2, aes(y = .derivative, group = segment_id, color = trend), linewidth = 1) +
    scale_color_manual(values = c("Significantly increasing" = "#E63946", 
                                  "Significantly decreasing" = "#2A9D8F")) +
    labs(title = plot_titles[[compound]], x = "Year", y = expression("f ' (x)"), color = "") +
    theme_minimal(base_size = 14) +  
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  
  
  print(p7)
  
  # Generate predicted values for trend segments (without and with δ15N correction)
  deriv_segments$predicted <- predict(gam_model, newdata = deriv_segments)
  deriv_segments2$predicted <- predict(gam_model2, newdata = deriv_segments2)
  
  # Plot 3: Raw and smoothed concentrations over time, including trend segments and critical points
  p8 <- ggplot(data_clean, aes(x = Year, y = .data[[compound]])) +
    geom_point(alpha = 0.5, color = "#4682B4") +
    geom_smooth(method = "gam", formula = y ~ s(x), 
                aes(fill = "Without"), color = "#4682B4", alpha = 0.2, level = 0.95) +  
    geom_line(data = deriv_segments,   
              aes(x = Year,   
                  y = predict(gam_model, newdata = deriv_segments),   
                  color = trend,   
                  group = segment_id),   
              linewidth = 1.5)
  
  # Add points and smoothing for δ15N-corrected data
  p8 <- p8 +
    geom_point(data = data_clean2, alpha = 0.5, color = "#C5A8D1") +
    geom_smooth(data = data_clean2, method = "gam", formula = y ~ s(x), 
                aes(fill = "With"), color = "#C5A8D1", alpha = 0.3, level = 0.95) +  
    geom_line(data = deriv_segments2,   
              aes(x = Year,   
                  y = predict(gam_model2, newdata = deriv_segments2),   
                  color = trend,   
                  group = segment_id),   
              linewidth = 1.5) 
  
  # Add GAM1 and GAM2 change points with horizontal confidence intervals
  p8 <- p8 +
    geom_point(
      aes(x = year_max_gam, 
          y = predict(gam_model, newdata = data.frame(Year = year_max_gam)), 
          color = "GAM1 change point"), 
      size = 3, shape = 21, fill = "#3A6A94"
    ) +
    geom_errorbarh(
      aes(xmin = ci_range$min_year, xmax = ci_range$max_year, 
          y = predict(gam_model, newdata = data.frame(Year = year_max_gam))), 
      height = 2, size = 0.6, color = "#3A6A94"
    ) +
    geom_point(
      aes(x = year_max_gam2, 
          y = predict(gam_model2, newdata = data.frame(Year = year_max_gam2)), 
          color = "GAM2 change point"), 
      size = 3, shape = 21, fill = "#5A3D77"
    ) +
    geom_errorbarh(
      aes(xmin = ci_range2$min_year2, xmax = ci_range2$max_year2, 
          y = predict(gam_model2, newdata = data.frame(Year = year_max_gam2))), 
      height = 2, size = 0.6, color = "#5A3D77"
    )
  
  # Add historical and regulatory context
  p8 <- p8 +
    geom_errorbarh(aes(xmin = 2002, xmax = 2024, y = 35, color = "Regulatory period"),
                   height = 2, size = 0.6) +  # Regulatory actions timeline
    geom_errorbarh(aes(xmin = 2000, xmax = 2003, y = 45, color = "Phase-out of PFOS"),
                   height = 2, size = 0.6)  # PFOS phase-out period
  
  
  # Customize colors and labels
  p8 <- p8 +
    scale_color_manual(
      name = "", 
      values = c(
        "Significantly increasing" = "#E63946",
        "Significantly decreasing" = "#2A9D8F",
        "Regulatory period" = "#66C2B1",
        "Phase-out of PFOS" = "darkred"
      ),
      labels = c(
        "Significantly increasing" = "Significantly increasing",
        "Significantly decreasing" = "Significantly decreasing",
        "Regulatory period" = "Regulatory period",
        "Phase-out of PFOS" = "Phase-out of PFOS"
      )
    ) +
    
    
    scale_fill_manual(
      name = "",
      values = c("Without" = "#4682B4", "With" = "#C5A8D1"),
      labels = c("Without" = expression("Without " * delta^15 * "N correction"),
                 "With" = expression("With " * delta^15 * "N correction"))
    ) +
    
    # Aesthetic theme and axis/plot title
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    labs(
      title = plot_titles[[compound]],
      x = "Year",
      y = "Lipid-standardized concentration (ng/g ww)"
    )
  
  # Display plot
  print(p8)
  
}


### PFCA
compounds <- "PFCA"

for (compound in compounds) {
  cat("\n Analysis for:", compound, "\n")
  
  # Remove missing values
  data_clean <- na.omit(data[, c("Year", compound)])
  
  # Fit a GAM model
  gam_model <- gam(as.formula(paste(compound, "~ s(Year, k=10)")), data = data_clean, method = "REML")
  
  # Compute derivatives and 95% confidence intervals
  deriv_values <- derivatives(gam_model, term = "s(Year)", type = "central", level = 0.95) %>%
    select(Year, .derivative, .lower_ci, .upper_ci) %>%
    mutate(
      trend = case_when(
        .lower_ci > 0 ~ "Significantly increasing",
        .upper_ci < 0 ~ "Significantly decreasing",
        TRUE ~ "Neutral"
      )
    )
  
  # Detect local maximum: transition from increasing to decreasing trend
  critical_points <- deriv_values %>%
    mutate(sign_change = lag(.derivative) > 0 & .derivative < 0) %>%
    filter(sign_change)
  
  if (nrow(critical_points) > 0) {
    
    # Select the year with the highest value of the compound
    max_year_data <- data_clean %>%
      filter(Year %in% critical_points$Year) %>%
      arrange(desc(!!sym(compound))) %>%
      slice(1)
    
    year_max <- max_year_data$Year
    
    # Estimate confidence interval around the GAM maximum using predict()
    pred_data <- data.frame(Year = seq(min(data_clean$Year), max(data_clean$Year), length.out = 500))
    pred_data$fit <- predict(gam_model, newdata = pred_data, se.fit = TRUE)$fit
    pred_data$se <- predict(gam_model, newdata = pred_data, se.fit = TRUE)$se.fit
    
    # Year corresponding to the maximum predicted GAM value
    year_max_gam <- pred_data$Year[which.max(pred_data$fit)]
    
    # Confidence interval for the GAM maximum
    ci_range <- pred_data %>%
      filter(fit >= max(fit) - 1.96 * se) %>%
      summarize(min_year = min(Year), max_year = max(Year))
    
    cat("\nLipid corrected GAM: detected maximum for", compound, "in year:", round(year_max_gam, 2), "\n")
    cat("95% confidence interval for the year of maximum: [", round(ci_range$min_year, 2), ", ", round(ci_range$max_year, 2), "]\n")
  }
  
  # Get predicted concentration and 95% CI at the year of the maximum
  new_data <- data.frame(Year = year_max_gam)
  pred <- predict(gam_model, newdata = new_data, se.fit = TRUE)
  
  fit_value <- pred$fit
  se_value <- pred$se.fit
  lower_ci <- fit_value - 1.96 * se_value
  upper_ci <- fit_value + 1.96 * se_value
  
  cat("Lipid corrected GAM: concentration at year of maximum (", round(year_max_gam, 2), "): ", round(fit_value, 2), "ng/g ww\n")
  cat("95% confidence interval: [", round(lower_ci, 2), ", ", round(upper_ci, 2), "]\n")
  
  # Segments distinguishing increasing and decreasing trends
  deriv_segments <- deriv_values %>%
    arrange(Year) %>%
    mutate(
      is_neutral = (trend == "Neutral"),
      trend_change = (trend != lag(trend, default = first(trend))),
      segment_id = cumsum(trend_change | is_neutral)
    ) %>%
    filter(trend != "Neutral")
  
  # Remove missing values in corrected dataset
  data_clean2 <- na.omit(data_corr[, c("Year", compound)])
  
  # Fit corrected GAM model
  gam_model2 <- gam(as.formula(paste(compound, "~ s(Year, k=10)")), data = data_clean2, method = "REML")
  
  # Compute derivatives and 95% CI
  deriv_values2 <- derivatives(gam_model2, term = "s(Year)", type = "central", level = 0.95) %>%
    select(Year, .derivative, .lower_ci, .upper_ci) %>%
    mutate(
      trend = case_when(
        .lower_ci > 0 ~ "Significantly increasing",
        .upper_ci < 0 ~ "Significantly decreasing",
        TRUE ~ "Neutral"
      )
    )
  
  # Detect local maximum 
  critical_points2 <- deriv_values2 %>%
    mutate(sign_change = lag(.derivative) > 0 & .derivative < 0) %>%
    filter(sign_change)
  
  if (nrow(critical_points2) > 0) {
    
    max_year_data2 <- data_clean2 %>%
      filter(Year %in% critical_points2$Year) %>%
      arrange(desc(!!sym(compound))) %>%
      slice(1)
    
    year_max2 <- max_year_data2$Year
    
    pred_data2 <- data.frame(Year = seq(min(data_clean2$Year), max(data_clean2$Year), length.out = 500))
    pred_data2$fit <- predict(gam_model2, newdata = pred_data2, se.fit = TRUE)$fit
    pred_data2$se <- predict(gam_model2, newdata = pred_data2, se.fit = TRUE)$se.fit
    
    year_max_gam2 <- pred_data2$Year[which(pred_data2$fit == max(pred_data2$fit, na.rm = TRUE))]
    
    ci_range2 <- pred_data2 %>%
      filter(fit >= max(fit) - 1.96 * se) %>%
      summarize(min_year2 = min(Year), max_year2 = max(Year))
    
    cat("\nLipid-d15N corrected GAM: detected maximum for", compound, "in year:", round(year_max_gam2, 2), "\n")
    cat("95% CI for year of maximum: [", round(ci_range2$min_year2, 2), ", ", round(ci_range2$max_year2, 2), "]\n")
  }
  
  # Predicted value and CI at the maximum year
  new_data2 <- data.frame(Year = year_max_gam2)
  pred2 <- predict(gam_model2, newdata = new_data2, se.fit = TRUE)
  
  fit_value2 <- pred2$fit
  se_value2 <- pred2$se.fit
  lower_ci2 <- fit_value2 - 1.96 * se_value2
  upper_ci2 <- fit_value2 + 1.96 * se_value2
  
  cat("Lipid-d15N corrected GAM: concentration at year of maximum (", round(year_max_gam2, 2), "): ", round(fit_value2, 2), "ng/g ww\n")
  cat("95% confidence interval: [", round(lower_ci2, 2), ", ", round(upper_ci2, 2), "]\n")
  
  # Compute start and end of each significant trend period
  trend_periods <- deriv_segments %>%
    group_by(segment_id, trend) %>%
    summarise(
      start_year = min(Year),
      end_year = max(Year),
      start_value = predict(gam_model, newdata = data.frame(Year = min(Year))),
      end_value = predict(gam_model, newdata = data.frame(Year = max(Year)))
    ) %>%
    ungroup()
  
  cat("\nSignificant trend periods (Lipid corrected model) for", compound, ":\n")
  print(trend_periods)
  
  # Same for corrected model
  deriv_segments2 <- deriv_values2 %>%
    arrange(Year) %>%
    mutate(
      is_neutral = (trend == "Neutral"),
      trend_change = (trend != lag(trend, default = first(trend))),
      segment_id = cumsum(trend_change | is_neutral)
    ) %>%
    filter(trend != "Neutral")
  
  trend_periods2 <- deriv_segments2 %>%
    group_by(segment_id, trend) %>%
    summarise(
      start_year = min(Year),
      end_year = max(Year),
      start_value = predict(gam_model2, newdata = data.frame(Year = min(Year))),
      end_value = predict(gam_model2, newdata = data.frame(Year = max(Year)))
    ) %>%
    ungroup()
  
  cat("\nSignificant trend periods (Lipid-d15N corrected model) for", compound, ":\n")
  print(trend_periods2)
  
  # Plot 1: First derivative of the GAM without δ15N correction
  p9 <- ggplot(deriv_values, aes(x = Year, y = .derivative)) +
    geom_line(color = "#4682B4", linewidth = 0.8) +  # Main derivative curve
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#4682B4", alpha = 0.2) +  # Confidence interval
    geom_line(data = deriv_segments, aes(y = .derivative, group = segment_id, color = trend), linewidth = 1) +  # Significant trends
    scale_color_manual(values = c("Significantly increasing" = "#E63946", 
                                  "Significantly decreasing" = "#2A9D8F")) +  # Custom colors for trend segments
    labs(title = plot_titles[[compound]], x = "Year", y = expression("f ' (x)"), color = "") +
    theme_minimal(base_size = 14) +  
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  
  
  print(p9)
  
  # Plot 2: First derivative of the GAM with δ15N correction
  p10 <- ggplot(deriv_values2, aes(x = Year, y = .derivative)) +
    geom_line(color = "#C5A8D1", linewidth = 0.8) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#C5A8D1", alpha = 0.3) +
    geom_line(data = deriv_segments2, aes(y = .derivative, group = segment_id, color = trend), linewidth = 1) +
    scale_color_manual(values = c("Significantly increasing" = "#E63946", 
                                  "Significantly decreasing" = "#2A9D8F")) +
    labs(title = plot_titles[[compound]], x = "Year", y = expression("f ' (x)"), color = "") +
    theme_minimal(base_size = 14) +  
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  
  
  print(p10)
  
  deriv_segments$y_pred <- predict(gam_model, newdata = deriv_segments)
  
  # Precompute predicted values for GAM trend segments
  deriv_segments$predicted <- predict(gam_model, newdata = deriv_segments)
  deriv_segments2$predicted <- predict(gam_model2, newdata = deriv_segments2)
  
  # Precompute predicted concentrations at GAM change points
  pred_gam1 <- predict(gam_model, newdata = data.frame(Year = year_max_gam))
  pred_gam2 <- predict(gam_model2, newdata = data.frame(Year = year_max_gam2))
  
  # Plot 3: Raw and smoothed concentrations over time with derivative trend segments and change points
  p11 <- ggplot(data_clean, aes(x = Year, y = .data[[compound]])) +
    geom_point(alpha = 0.5, color = "#4682B4") +
    geom_smooth(method = "gam", formula = y ~ s(x), 
                aes(fill = "Without"), color = "#4682B4", alpha = 0.2, level = 0.95) +  
    geom_line(data = deriv_segments,   
              aes(x = Year,   
                  y = predict(gam_model, newdata = deriv_segments),   
                  color = trend,   
                  group = segment_id),   
              linewidth = 1.5)
  
  # Add points and smoothing for δ15N-corrected data
  p11 <- p11 +
    geom_point(data = data_clean2, alpha = 0.5, color = "#C5A8D1") +
    geom_smooth(data = data_clean2, method = "gam", formula = y ~ s(x), 
                aes(fill = "With"), color = "#C5A8D1", alpha = 0.3, level = 0.95) +  
    geom_line(data = deriv_segments2,   
              aes(x = Year,   
                  y = predict(gam_model2, newdata = deriv_segments2),   
                  color = trend,   
                  group = segment_id),   
              linewidth = 1.5) 
  
  # Add GAM1 and GAM2 change points with horizontal confidence intervals
  p11 <- p11 +
    geom_point(
      aes(x = year_max_gam, 
          y = predict(gam_model, newdata = data.frame(Year = year_max_gam)), 
          color = "GAM1 change point"), 
      size = 3, shape = 21, fill = "#3A6A94"
    ) +
    geom_errorbarh(
      aes(xmin = ci_range$min_year, xmax = ci_range$max_year, 
          y = predict(gam_model, newdata = data.frame(Year = year_max_gam))), 
      height = 1, size = 0.6, color = "#3A6A94"
    ) +
    geom_point(
      aes(x = year_max_gam2, 
          y = predict(gam_model2, newdata = data.frame(Year = year_max_gam2)), 
          color = "GAM2 change point"), 
      size = 3, shape = 21, fill = "#5A3D77"
    ) +
    geom_errorbarh(
      aes(xmin = ci_range2$min_year2, xmax = ci_range2$max_year2, 
          y = predict(gam_model2, newdata = data.frame(Year = year_max_gam2))), 
      height = 1, size = 0.6, color = "#5A3D77"
    )
  
  # Add historical and regulatory context
  p11 <- p11 +
    geom_errorbarh(aes(xmin = 2002, xmax = 2024, y = 10, color = "Regulatory period"), 
                   height = 1, size = 0.6) +
    geom_errorbarh(aes(xmin = 2010, xmax = 2015, y = 20, color = "Phase-out of PFOA and LC-PFCAs"), 
                   height = 1, size = 0.6)
  
  # Customize colors for trend segments and regulatory annotations
  p11 <- p11 +
    scale_color_manual(
      name = "", 
      values = c(
        "Significantly increasing" = "#E63946",  
        "Significantly decreasing" = "#2A9D8F",
        "Phase-out of PFOS" = "darkred",
        "Phase-out of PFOA and LC-PFCAs" = "darkorange",
        "Regulatory period" = "#66C2B1"))
  
  # Customize fill for smoothed GAM curves
  p11 <- p11 +
    scale_fill_manual(
      name = "", 
      values = c("Without" = "#4682B4", "With" = "#C5A8D1"),
      labels = c(
        "Without" = expression("Without " * delta^15 * "N correction"),
        "With" = expression("With " * delta^15 * "N correction")))
  
  # Final plot customization
  p11 <- p11 +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
    labs(
      title = plot_titles[[compound]],
      x = "Year",
      y = "Lipid-standardized concentration (ng/g ww)")
  
  # Display plot
  print(p11)
  
}


### Relative contributions ###

# Create a table indicating the peak year of concentration for each compound
peak_years_table <- data.frame(
  Compound = c("PFOS", "PFHxS", "PFDS", "PFSA",
               "PFOA", "PFNA", "PFUdA", "PFDoA", "PFTrDA", "PFTeDA", "PFDA", 
               "PFCA", "PFAS"),
  Peak_Year = c(1994, 1999, 1994, 1994,
                1999, 2006, 2007, 2007, 2005, 2004, 2024, 
                2006, 1996)
)


create_proportion_table <- function(data, compounds) {
  compounds <- setdiff(compounds, "PFAS")  
  data %>%
    select(Year, all_of(compounds), PFAS) %>%
    pivot_longer(
      cols = all_of(compounds),
      names_to = "Compound",
      values_to = "Concentration"
    ) %>%
    mutate(Proportion = (Concentration / PFAS) * 100) %>%
    select(Year, Compound, Proportion) %>%
    arrange(Year, Compound) %>%
    replace_na(list(Proportion = 0))
}

# Compute the relative contribution of selected PFAS compounds
df_proportion <- create_proportion_table(data, c("PFOA", "PFNA", "PFDA", "PFUdA", "PFDoA", "PFTrDA", "PFTeDA", "PFHxS", "PFOS", "PFDS"))  

# Recode compound names to include carbon chain length or groupings
df_proportion$Compound <- dplyr::recode(df_proportion$Compound,
                                        "PFOA" = "PFOA (C8)",
                                        "PFNA" = "PFNA (C9)",
                                        "PFDA" = "PFDA (C10)",
                                        "PFUdA" = "PFUdA (C11)",
                                        "PFDoA" = "PFDoA (C12)",
                                        "PFTrDA" = "PFTrDA (C13)",
                                        "PFTeDA" = "PFTeDA (C14)",
                                        "PFHxS" = "PFHxS (C6)",
                                        "PFOS" = "PFOS (C8)",
                                        "PFDS" = "PFDS (C10)",
                                        "PFAS" = "Σ[17] PFASs",
                                        "PFCA" = "Σ[13] PFCAs",
                                        "PFSA" = "Σ[4] PFSAs",
                                        .default = df_proportion$Compound
)


# Plot for compounds with proportions < 20%
p_bottom <- ggplot(df_proportion[df_proportion$Proportion < 20, ], 
                   aes(x = Year, y = Proportion, color = Compound)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  labs(x = NULL, y = NULL, color = "") +
  scale_color_manual(
    values = c(
      "PFOA (C8)" = "#1f77b4",
      "PFNA (C9)" = "#ff7f0e",
      "PFDA (C10)" = "#2ca02c",
      "PFUdA (C11)" = "#d62728",
      "PFDoA (C12)" = "#9467bd",
      "PFTrDA (C13)" = "#8c564b",
      "PFTeDA (C14)" = "#e377c2",
      "PFHxS (C6)" = "#7f7f7f",
      "PFDS (C10)" = "#17becf"
    ),
    breaks = c("PFOA (C8)", "PFNA (C9)", "PFDA (C10)", "PFUdA (C11)", 
               "PFDoA (C12)", "PFTrDA (C13)", "PFTeDA (C14)", 
               "PFHxS (C6)", "PFDS (C10)")
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0, 20))

# Plot for compounds with proportions > 40%
p_top <- ggplot(df_proportion[df_proportion$Proportion > 40, ], 
                aes(x = Year, y = Proportion, color = Compound)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  labs(x = NULL, y = NULL, color = "") +
  scale_color_manual(
    values = c("PFOS (C8)" = "#cab2d6"),
    breaks = c("PFOS (C8)")
  ) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "right"
  ) +
  scale_y_continuous(limits = c(40, 80))

# Combine top and bottom plots using patchwork
p_final <- p_top / p_bottom + 
  plot_layout(heights = c(1, 2)) +
  labs(x = "Year", y = "Relative contribution (%)", color = "") +
  scale_color_manual(
    values = c(
      "PFOA (C8)" = "#1f77b4",
      "PFNA (C9)" = "#ff7f0e",
      "PFDA (C10)" = "#2ca02c",
      "PFUdA (C11)" = "#d62728",
      "PFDoA (C12)" = "#9467bd",
      "PFTrDA (C13)" = "#8c564b",
      "PFTeDA (C14)" = "#e377c2",
      "PFHxS (C6)" = "#7f7f7f",
      "PFOS (C8)" = "#cab2d6",
      "PFDS (C10)" = "#17becf"
    ),
    breaks = c("PFOA (C8)", "PFNA (C9)", "PFDA (C10)", "PFUdA (C11)", 
               "PFDoA (C12)", "PFTrDA (C13)", "PFTeDA (C14)", 
               "PFHxS (C6)", "PFOS (C8)", "PFDS (C10)")
  ) +
  theme(
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 1.2),
    legend.position = "right",
    legend.justification = c(0.5, 1.2),
    legend.box.margin = margin(l = 10)
  )

# Display the final plot
p_final


# Recode compound names in the peak years table
peak_years_table$Compound <- dplyr::recode(peak_years_table$Compound,
                                           "PFOA" = "PFOA (C8)",
                                           "PFNA" = "PFNA (C9)",
                                           "PFUdA" = "PFUdA (C11)",
                                           "PFDoA" = "PFDoA (C12)",
                                           "PFTrDA" = "PFTrDA (C13)",
                                           "PFTeDA" = "PFTeDA (C14)",
                                           "PFHxS" = "PFHxS (C6)",
                                           "PFOS" = "PFOS (C8)",
                                           "PFDS" = "PFDS (C10)",
                                           "PFAS" = "Σ[17] PFASs",
                                           "PFCA" = "Σ[13] PFCAs",
                                           "PFSA" = "Σ[4] PFSAs",
                                           .default = peak_years_table$Compound
)


compare_before_after_with_tests <- function(df_proportion, peak_years) {
  # Filter relevant compounds
  useful_compounds <- unique(df_proportion$Compound)
  peak_years <- peak_years %>% filter(Compound %in% useful_compounds)
  
  results <- data.frame()
  
  for (i in 1:nrow(peak_years)) {
    compound <- peak_years$Compound[i]
    peak <- peak_years$Peak_Year[i]
    
    subset <- df_proportion %>% filter(Compound == compound)
    
    before <- subset %>% filter(Year < peak) %>% pull(Proportion)
    after <- subset %>% filter(Year >= peak) %>% pull(Proportion)
    
    if (length(before) >= 3 & length(after) >= 3) {
      test_df <- data.frame(Proportion = c(before, after),
                            Group = rep(c("before", "after"), c(length(before), length(after))))
      
      # Normality tests
      norm_before <- shapiro.test(before)
      norm_after <- shapiro.test(after)
      normality_ok <- norm_before$p.value > 0.05 & norm_after$p.value > 0.05
      
      # Homogeneity of variances
      levene <- car::leveneTest(Proportion ~ Group, data = test_df)
      
      # Selection of appropriate statistical test
      if (normality_ok) {
        if (levene$`Pr(>F)`[1] > 0.05) {
          test <- t.test(before, after, var.equal = TRUE)
          test_used <- "t-test"
        } else {
          test <- t.test(before, after, var.equal = FALSE)
          test_used <- "Welch t-test"
        }
      } else {
        test <- wilcox.test(before, after)
        test_used <- "Wilcoxon"
      }
      
      # Append results
      results <- rbind(results, data.frame(
        Compound = compound,
        Peak_Year = peak,
        Normality_Before_p = round(norm_before$p.value, 4),
        Normality_After_p = round(norm_after$p.value, 4),
        Variance_Equal_p = round(levene$`Pr(>F)`[1], 4),
        Test = test_used,
        p_value = round(test$p.value, 4),
        Median_Before = round(median(before, na.rm = TRUE), 2),
        Median_After = round(median(after, na.rm = TRUE), 2)
      ))
    }
  }
  
  return(results)
}

# Apply the comparison function and export results
stat_results <- compare_before_after_with_tests(df_proportion, peak_years_table)
view(stat_results)
write.csv2(stat_results, "/Users/anaisfournier/Desktop/Stats/stat_results.csv", row.names = FALSE)


### Ecological half lives ###

# Define the peak year for each compound
peak_years_table <- data.frame(
  Compound = c("PFOS", "PFHxS", "PFDS", "PFSA",
               "PFOA", "PFNA", "PFUdA", "PFDoA", "PFTrDA", "PFTeDA", "PFDA", 
               "PFCA", "PFAS"),
  Peak_Year = c(1994, 1999, 1994, 1994,
                1999, 2006, 2007, 2007, 2005, 2004, 2024, 
                2007, 1996)
)

# Initialize result table
half_lives <- data.frame(Compound = character(), Half_Life = numeric(), stringsAsFactors = FALSE)

# Select compound columns (columns 8 onward)
compound_cols <- colnames(data)[8:ncol(data)]

# Estimate half-life for each compound
for (compound in compound_cols) {
  peak_year <- peak_years_table$Peak_Year[peak_years_table$Compound == compound]
  
  if (length(peak_year) == 0) {
    warning(paste("Peak year not found for compound:", compound))
    next
  }
  
  sub_data <- data %>%
    select(Year, Concentration = all_of(compound)) %>%
    filter(!is.na(Concentration), Concentration > 0, Year >= peak_year, Year <= 2024)
  
  if (nrow(sub_data) >= 5) {
    lm_model <- lm(log(Concentration) ~ Year, data = sub_data)
    slope <- coef(lm_model)[2]
    
    pred_peak <- exp(predict(lm_model, newdata = data.frame(Year = peak_year)))
    pred_2024 <- exp(predict(lm_model, newdata = data.frame(Year = 2024)))
    
    r_est <- -log(pred_2024 / pred_peak) / (2024 - peak_year)
    t_half <- ifelse(r_est > 0, log(2) / r_est, NA)
  } else {
    t_half <- NA
  }
  
  half_lives <- rbind(half_lives, data.frame(Compound = compound, Half_Life = t_half))
}

# Output results
print(half_lives)
write.csv2(half_lives, "/Users/anaisfournier/Desktop/Stats/half_lives_lipid.csv", row.names = FALSE)


peak_years_table <- data.frame(
  Compound = c("PFOS", "PFHxS", "PFDS", "PFSA",
               "PFOA", "PFNA", "PFUdA", "PFDoA", "PFTrDA", "PFTeDA", "PFDA", 
               "PFCA", "PFAS"),
  Peak_Year = c(1997, 1999, 1997, 1996,
                2000, 2006, 2006, 2006, 2005, 2004, 2009, 
                2006, 1998)
)

half_lives <- data.frame(Compound = character(), Half_Life = numeric(), stringsAsFactors = FALSE)

compound_cols <- colnames(data_corr)[8:ncol(data_corr)]

for (compound in compound_cols) {
  peak_year <- peak_years_table$Peak_Year[peak_years_table$Compound == compound]
  
  if (length(peak_year) == 0) {
    warning(paste("Peak year not found for compound:", compound))
    next
  }
  
  sub_data <- data_corr %>%
    select(Year, Concentration = all_of(compound)) %>%
    filter(!is.na(Concentration), Concentration > 0, Year >= peak_year, Year <= 2024)
  
  if (nrow(sub_data) >= 5) {
    lm_model <- lm(log(Concentration) ~ Year, data = sub_data)
    slope <- coef(lm_model)[2]
    
    pred_peak <- exp(predict(lm_model, newdata = data.frame(Year = peak_year)))
    pred_2024 <- exp(predict(lm_model, newdata = data.frame(Year = 2024)))
    
    r_est <- -log(pred_2024 / pred_peak) / (2024 - peak_year)
    t_half <- ifelse(r_est > 0, log(2) / r_est, NA)
  } else {
    t_half <- NA
  }
  
  half_lives <- rbind(half_lives, data.frame(Compound = compound, Half_Life = t_half))
}

# Output results
print(half_lives)
write.csv2(half_lives, "/Users/anaisfournier/Desktop/Stats/half_lives_lipid_d15N.csv", row.names = FALSE)


