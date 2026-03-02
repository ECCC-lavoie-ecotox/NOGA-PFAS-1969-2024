###############################
### SI - Half a century of per- and polyfluoroalkyl substances in Northern Gannet eggs: Impact of regulations
### Step 4 - Lipid correction
### Created by Anaïs Fournier on 05-14-2025
###############################

# Import dataset 
data <- read.csv2("/Users/anaisfournier/Desktop/Stats/Step 3&4 - Correlation of Moisture, Lipid, d15N, d13CS and PFAS Concentrations.csv")
data[, 3:7] <- lapply(data[, 3:7], as.numeric)

## Method used in the study :  normalizing concentrations as if all samples had the same lipid content
# More details are provided in the Supporting Information file (page S4)
# Lipid content maximum value (L_max)
L_max <- 5.54

# Select columns 8 to 20 for the analysis
selected_columns <- c(8:20)

# Initialize a data frame to store the results
results <- data.frame(Compound = character(), Intercept = numeric(), Slope = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Initialize a data frame to store corrected concentrations
stand <- data.frame(FID = data$FID, 
                             Location = data$Location, 
                             Year = data$Year, 
                             Moisture = data$Moisture,
                             Lipid = data$Lipid, 
                             d15N = data$d15N, 
                             d13CS = data$d13CS)

# Loop through each selected column
for (col in selected_columns) {
  compound_name <- colnames(data)[col]  # Retrieve the compound name
  
  # Remove only the rows where the compound is NA
  subset_data <- data[!is.na(data[[compound_name]]), c(compound_name, "Lipid")]
  
  # Fit a linear model
  model <- lm(subset_data[[compound_name]] ~ subset_data$Lipid)
  
  # Extract the intercept, slope, and p-value
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  pvalue <- summary(model)$coefficients[2, 4]  
  
  # Calculate corrected concentrations using the formula
  corrected_concentrations <- subset_data[[compound_name]] + slope * (L_max - subset_data$Lipid)
  
  # Ensure corrected_concentrations matches the full data set by aligning on Sample_ID
  corrected_concentrations_full <- rep(NA, nrow(data))  # Initialize a vector of NAs for all rows
  corrected_concentrations_full[!is.na(data[[compound_name]])] <- corrected_concentrations  # Fill in the corrected values for rows with data
  
  # Add the corrected concentrations to the corrected_data frame
  stand[[compound_name]] <- corrected_concentrations_full
  
  # Add the results to the results data frame
  results <- rbind(results, data.frame(Compound = compound_name, Intercept = intercept, Slope = slope, p_value = pvalue))
}

# Reset row names
rownames(results) <- NULL

# To see the results
view(results)
view(stand)


## The typical method used in the literature : to divide by the lipid content
# Create a new data frame named 'div' with selected columns and lipid-corrected concentrations
div <- data.frame(
  FID = data$FID, 
  Location = data$Location, 
  Year = data$Year, 
  Moisture = data$Moisture,
  Lipid = data$Lipid, 
  d15N = data$d15N, 
  d13CS = data$d13CS
)

# Select columns 8 to 20 (which contain the concentrations)
selected_columns <- 8:20

# Loop through each of the selected columns to calculate the lipid-corrected concentrations
for (col in selected_columns) {
  compound_name <- colnames(data)[col]
  
  # Create a new column with the corrected concentration by dividing by Lipid content
  div[[paste(compound_name)]] <- data[[compound_name]] / (data$Lipid/100)
}

# To see the new data frame 'div' with lipid-corrected concentrations
view(div)


# Comparison of the two methods (with PFOS because it is the mojor compound)
# Scatter plot comparing the two correction methods
plot(div$PFOS ~ stand$PFOS,
     xlab = "Method 2", 
     ylab = "Method 1")

# Fit a linear model to assess the relationship between the two methods
mod <- lm(div$PFOS ~ stand$PFOS) 
# Add the regression line to the plot
abline(mod)

# Calculate and display the coefficient R²
R2 <- summary(mod)$r.squared
text(x = min(stand$PFOS), y = max(div$PFOS), 
     labels = paste("R² =", round(R2, 3)), pos = 4)

write.csv2(stand, "/Users/anaisfournier/Desktop/Stats/Step 5 - Lipid-Corrected Concentrations for d15N Corection.csv", row.names = FALSE)
