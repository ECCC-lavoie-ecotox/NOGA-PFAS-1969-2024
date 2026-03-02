###############################
### SI - Half a century of per- and polyfluoroalkyl substances in Northern Gannet eggs: Impact of regulations
### Step 3 - Analysis of temporal trends and linear relationships between PFOS concentrations and Lipid, Moisture, δ15N or δ13CS
### Created by Anaïs Fournier on 05-14-2025
###############################

# Load libraries
library(ggplot2)

# Import dataset 
dat <- read.csv2("/Users/anaisfournier/Desktop/Stats/Step 3&4 - Correlation of Moisture, Lipid, d15N, d13CS and PFAS Concentrations.csv")
dat[, 3:7] <- lapply(dat[, 3:7], as.numeric)

# Plot of Lipid content over time
ggplot(dat, aes(x = Year, y = Lipid)) +
  geom_point(color = "black") +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black", fill = "black", alpha = 0.2) +
  labs(title = "Lipid Content",
       x = "Year",
       y = "Lipid Content (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# Plot of Moisture content over time
ggplot(dat, aes(x = Year, y = Moisture)) +
  geom_point(color = "black") +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black", fill = "black", alpha = 0.2) +
  labs(title = "Moisture Content",
       x = "Year",
       y = "Moisture Content (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# Plot of δ15N values over time
ggplot(dat, aes(x = Year, y = d15N)) +
  geom_point(color = "black") +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black", fill = "black", alpha = 0.2) +
  labs(title = expression(delta^15 * "N"),
       x = "Year",
       y = expression(delta^15 * "N (‰)")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# Plot of δ13CS values over time
ggplot(dat, aes(x = Year, y = d13CS)) +
  geom_point(color = "black") +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black", fill = "black", alpha = 0.2) +
  labs(title = expression(delta^13 * "CS"),
       x = "Year",
       y = expression(delta^13 * "CS (‰)")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )


# Linear regression between PFOS concentrations (because it is the major compound) and the 4 parameters
lm_lipid <- lm(PFOS ~ Lipid, data = dat)
lm_moisture <- lm(PFOS ~ Moisture, data = dat)
lm_d15N <- lm(PFOS ~ d15N, data = dat)
lm_d13CS <- lm(PFOS ~ d13CS, data = dat)

# Summary of the regression models
summary(lm_lipid)
summary(lm_moisture)
summary(lm_d15N)
summary(lm_d13CS)

# Extracting coefficients and other information from the model summaries
coef_lipid <- summary(lm_lipid)$coefficients
coef_moisture <- summary(lm_moisture)$coefficients
coef_d15N <- summary(lm_d15N)$coefficients
coef_d13CS <- summary(lm_d13CS)$coefficients

# Creating a table to display the results
results <- data.frame(
  Variable = c("Lipid", "Moisture", "δ15N", "δ13CS"),
  Intercept = c(coef_lipid[1,1], coef_moisture[1,1], coef_d15N[1,1], coef_d13CS[1,1]),
  Slope = c(coef_lipid[2,1], coef_moisture[2,1], coef_d15N[2,1], coef_d13CS[2,1]),
  P_value = c(coef_lipid[2,4], coef_moisture[2,4], coef_d15N[2,4], coef_d13CS[2,4])
)

# Displaying the results
print(results)

# Scatter plots with linear regression between PFOS concentrations and variables
# PFOS vs Lipid
ggplot(dat, aes(x = Lipid, y = PFOS)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "black", alpha = 0.2) +
  labs(x = "Lipid Content (%)",
       y = "PFOS Concentration (ng/g ww)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# PFOS vs Moisture
ggplot(dat, aes(x = Moisture, y = PFOS)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "black", alpha = 0.2) +
  labs(x = "Moisture Content (%)",
       y = "PFOS Concentration (ng/g ww)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# PFOS vs δ15N
ggplot(dat, aes(x = d15N, y = PFOS)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "black", alpha = 0.2) +
  labs(x = expression(delta^15 * "N (‰)"),
       y = "PFOS Concentration (ng/g ww)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# PFOS vs δ13CS
ggplot(dat, aes(x = d13CS, y = PFOS)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "black", alpha = 0.2) +
  labs(x = expression(delta^13 * "CS (‰)"),
       y = "PFOS Concentration (ng/g ww)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

