# Load required libraries
rm()
library(raster)
library(sf)
library(dplyr)
library(caret)
library(ranger)
library(ggplot2)
library(lime)
library(reshape2)
library(stats)
library(corrplot)
library(spdep)
library(ape)
library(gstat)
library(sp)
library(viridis)
library(RColorBrewer)
library(ggrepel)

# Load datasets (120000 samples)
load("/gpfs1/work/acostave/120000/shannon_samples.RData")
load("/gpfs1/work/acostave/120000/utm32.RData")
load("/gpfs1/work/acostave/120000/tsm_matrix.RData")

####

load("/gpfs1/work/acostave/2020/traindata.RData")
load("/gpfs1/work/acostave/2020/testdata.RData")
load("/gpfs1/work/acostave/2020/forestmodel.RData")

sampled_coords_df <- data.frame(coordinates(sample_coordinates_utm32))

# Load raster layers
  slope <- raster("/gpfs1/data/satellite/dgm_5m/mosaic/slope_res20_horn.tif")
  aspect <- raster("/gpfs1/data/satellite/dgm_5m/mosaic/aspect_res20_horn.tif")
  soil <- raster("/gpfs1/work/acostave/soilinformation_1000_20mutm32.tif")
  entropy <- raster("/gpfs1/work/acostave/Entropy2/mosaic_entropy.tif")
  forestcover <- raster("/gpfs1/work/acostave/Forest_Count/forestcount.tif")
  fca <- raster("/gpfs1/data/satellite/sentinel/germany/forest-condition/v0007-0005/2020/tif-Germany_cog/Forest_condition_v0007-0005_Germany_2020_full_bmean_R20m_cog.tif")
  elevation <- raster("/gpfs1/work/acostave/elevation_cleaned.tif")
  twi <- raster("/gpfs1/schlecker/home/acostave/twi_res20_cog.tif")
  #precip <- raster("/gpfs1/work/acostave/precipitation/precipitation_2017.tif")
  dominant_sp <- raster("/gpfs1/work/acostave/Species_Dominant/species_dominant_mosaic.tif")
  #gr <- raster("/gpfs1/data/rs-validation/validation_forest_condition/growing_areas/raster/wg_raster_20m.tif")
  
print ("rasters imported")
# Extract values from raster layers
  slope_ <- extract(slope, sample_coordinates_utm32)
  aspect_ <- extract(aspect, sample_coordinates_utm32)
  soil_ <- extract(soil, sample_coordinates_utm32)
  entropy_ <- extract(entropy, sample_coordinates_utm32)
  forestcover_ <- extract(forestcover, sample_coordinates_utm32)
  elevation_ <- extract(elevation, sample_coordinates_utm32)
  fca_ <- extract(fca, sample_coordinates_utm32)
  twi_ <- extract(twi, sample_coordinates_utm32)
  #precip_ <- extract(precip, sample_coordinates_utm32)
  dominant_sp_ <- extract(dominant_sp, sample_coordinates_utm32)
  #growingreg_ <- extract(gr, sample_coordinates_utm32)
  
# Create dataframe
  shannon_41x41 <- shannon_matrix[,"41x41"] / log(11)
  shannon_df <- as.data.frame(shannon_41x41)
  soil_df <- as.data.frame(soil_)
  elevation_df <- as.data.frame(elevation_)
  slope_df <- as.data.frame(slope_)
  aspect_df <- as.data.frame(aspect_)
  entropy_df <- as.data.frame(entropy_)
  forestcover_df <- as.data.frame(forestcover_)
  fca_df <- as.data.frame(fca_)
  fca_df2 <- fca_df / 10000
  twi_df <- as.data.frame(twi_)
  dominant_spdf<- as.data.frame(dominant_sp_)
  #growingreg_df <- as.data.frame(growingreg_)
 
#Combine data
tablemodel <- data.frame(shannon_df, soil_df, elevation_df, slope_df, aspect_df, entropy_df, twi_df, dominant_spdf, forestcover_df, fca_df2)
tablemodel$soil_ <- as.factor(tablemodel$soil_)
tablemodel$elevation_ <- as.factor(tablemodel$elevation_)
tablemodel$dominant_sp_ <- as.factor(tablemodel$dominant_sp_)
#tablemodel$growingreg_ <- as.factor(tablemodel$growingreg_)

str(tablemodel)

# Add coordinates
table_model_coord <- data.frame(x = sampled_coords_df[, 1], y = sampled_coords_df[, 2], tablemodel)
write.csv(table_model_coord, "tablemodel_scaled.csv", row.names = FALSE)


# Remove missing values
fully_cleaned_data <- na.omit(table_model_coord)
save(fully_cleaned_data, file = "/gpfs1/work/acostave/FINALRESULTS/fully_cleaned_data2020.RData")
str(fully_cleaned_data)
# Correlation matrix
categorical_vars <- c("elevation_", "dominant_sp_", "soil_")
columns_to_scale <- setdiff(names(fully_cleaned_data), c("x", "y", "fca_", categorical_vars))
numeric_columns <- setdiff(columns_to_scale, categorical_vars)
correlation_matrix <- cor(fully_cleaned_data[, numeric_columns], use = "complete.obs")

png("/gpfs1/work/acostave/correlation_matrix2020.png", width = 1000, height = 800, res = 150)

corrplot(correlation_matrix, 
         method = "color", col = viridis(200), addCoef.col = "black", number.cex = 0.6,                
         number.font = 1, diag = TRUE, tl.cex = 0.8, cl.cex = 0.8, tl.col = "black",  
         mar = c(1, 1, 1, 1))    

# Guardar la imagen de la matriz de correlación
print(correlation_matrix)

dev.off()

fully_cleaned_data[columns_to_scale] <- as.data.frame(lapply(
  fully_cleaned_data[columns_to_scale], 
  function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
))

# Split dataset
set.seed(123)
trainIndex <- sample(1:nrow(fully_cleaned_data), size = 0.8 * nrow(fully_cleaned_data))
trainData <- fully_cleaned_data[trainIndex, ]
testData <- fully_cleaned_data[-trainIndex, ]


output_dir <- "/work/acostave/FINALRESULTS"
trainData_file <- file.path(output_dir, "traindata2020.RData")
dir_path_td <- dirname(trainData_file)
if (!dir.exists(dir_path_td)) {
  dir.create(dir_path_td, recursive = TRUE)
  print(paste("Created directory for traindata2020.RData:", dir_path_td))
}
print("traindatacreated")
save(trainData, file=trainData_file)


testData_file <- file.path(output_dir, "testdata2020.RData")
dir_path_ted <- dirname(testData_file)
if (!dir.exists(dir_path_ted)) {
  dir.create(dir_path_ted, recursive = TRUE)
  print(paste("Created directory for testdata2020.RData:", dir_path_ted))
}
print("testdatacreated")
save(testData, file=testData_file)
str(fully_cleaned_data)
print(dim(fully_cleaned_data))
print(dim(trainData))
print(dim(testData))
print(names(trainData))

###Recursive Feature Elimination (RFE) RFE to select important features

# columns_rfe <- columns_to_scale
# control <- rfeControl(functions = rfFuncs, method = "cv", number = 3)
# rfe_results <- rfe(
# x = trainData[, columns_rfe, drop = FALSE],  # Features for RFE
# y = as.vector(trainData$fca_),                              # Target variable
# sizes = c(1:length(columns_rfe)),            # Subset sizes
# rfeControl = control
# )
# # 
# # # Print RFE results
# print(rfe_results)
# 
# 
# ## Extract selected features
# selected_features <- predictors(rfe_results)
# cat("Selected Features:\n")
# print(selected_features)
# # 
# # 
# final_features <- c(selected_features, "soil_", "elevation_", "dominant_sp_" )
# trainData <- trainData[, c(final_features, "fca_"), drop = FALSE]
# testData <- testData[, c(final_features, "fca_"), drop = FALSE]

tuneGrid <- expand.grid(
  splitrule = c("variance", "extratrees"),
  min.node.size = c(1, 5),
  mtry = c(2,3)
)

# Train Random Forest model
random_forest <- train(
  fca_ ~ . -x -y,
  data = trainData,
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5),
  importance = "permutation",
  num.trees = 700
)

print(random_forest)

save(random_forest, file = "/gpfs1/work/acostave/FINALRESULTS/forestmodel2020.RData")

final_model <- random_forest$finalModel

# Extract residuals
coordinates(testData) <- ~x + y
predictions <- predict(random_forest, newdata = testData)
testData$residuals <- testData$fca_ - predictions

# Performance metrics
var_explained <- final_model$r.squared
print(paste("Percentage of variance explained:", var_explained * 100, "%"))

oob_error <- final_model$prediction.error
print(paste("Out-of-Bag Error:", oob_error))

ss_res <- sum((testData$fca_ - predictions)^2)
ss_tot <- sum((testData$fca_ - mean(testData$fca_))^2)
r_squared <- 1 - (ss_res / ss_tot)
rmse <- sqrt(mean((predictions - testData$fca_)^2))
mae <- mean(abs(predictions - testData$fca_))
mse <- mean((predictions - testData$fca_)^2)

cat("Testing - R2 Score:", r_squared, "\n")
cat("Testing - MAE:", mae, "\n")
cat("Testing - MSE:", mse, "\n")
cat("Testing - RMSE:", rmse, "\n")

train_predictions <- predict(random_forest, newdata = trainData)

ss_res_train <- sum((trainData$fca_ - train_predictions)^2)
ss_tot_train <- sum((trainData$fca_ - mean(trainData$fca_))^2)
r_squared_train <- 1 - (ss_res_train / ss_tot_train)
rmse_train <- sqrt(mean((train_predictions - trainData$fca_)^2))

# Métricas adicionales en entrenamiento
mae_train <- mean(abs(train_predictions - trainData$fca_))
mse_train <- mean((train_predictions - trainData$fca_)^2)

# Mostrar métricas del conjunto de entrenamiento
cat("Training - R2 Score:", r_squared_train, "\n")
cat("Training - MAE:", mae_train, "\n")
cat("Training - MSE:", mse_train, "\n")
cat("Training - RMSE:", rmse_train, "\n")

# Semivariogram analysis

semivariogram_residuals <- variogram(residuals ~ 1, data = testData)
semivario_df_residuals <- as.data.frame(semivariogram_residuals)
semivariogram_fca <- variogram(fca_ ~ 1, data = testData)
semivario_df_fca <- as.data.frame(semivariogram_fca)

p_residuals <- ggplot(semivario_df_residuals, aes(x = dist, y = gamma)) +
  geom_point(color = "#4B0082", size = 2, alpha = 1) +  # Puntos oscuros y grandes
  labs(title = "Semivariogram of Residuals", 
       x = "Distance (meters)", 
       y = "Semivariance") +
  scale_y_continuous(limits = c(0, max(semivario_df_residuals$gamma) * 1.2)) +  # Ajuste del eje Y
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    panel.grid.major = element_line(linetype = "solid", color = "gray80", linewidth = 0.4),
    panel.grid.minor = element_line(linetype = "solid", color = "gray90", linewidth = 0.2),
    legend.position = "bottom"  
  ) +
  guides(color = guide_legend(title = "Residuals"))

p_fca <- ggplot(semivario_df_fca, aes(x = dist, y = gamma)) +
  geom_point(color = "#4B0082", size = 2, alpha = 1) +
  labs(title = "Semivariogram of fca", 
       x = "Distance (meters)", 
       y = "Semivariance") +
  scale_y_continuous(limits = c(0, max(semivario_df_fca$gamma) * 1.2)) +  # Ajuste del eje Y
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    panel.grid.major = element_line(linetype = "solid", color = "gray80", linewidth = 0.4),
    panel.grid.minor = element_line(linetype = "solid", color = "gray90", linewidth = 0.2),
    legend.position = "bottom"
  )+
  guides(color = guide_legend(title = "FCA Values"))


ggsave("semivariogram_residuals2020.png", plot = p_residuals, width = 8, height = 6, dpi = 300)
print(p_residuals)
ggsave("semivariogram_fca2020.png", plot = p_fca, width = 8, height = 6, dpi = 300)
print(p_fca)


# Moran's I test
nb <- dnearneigh(coordinates(testData), 0, max(dist(coordinates(testData)))/10)
lw <- nb2listw(nb, style = "W")
morans_test <- moran.test(testData$residuals, lw)
print(morans_test)

# Spatial Residuals map
library(ggplot2)
library(viridis)
library(sf)

# Ensure residuals column exists and remove NAs
if (!"residuals" %in% names(testData_df)) stop("Column 'residuals' not found in dataset!")
testData_df <- testData_df[!is.na(testData_df$residuals), ]

# Define breakpoints for the color scale
residual_breaks <- seq(floor(min(testData_df$residuals, na.rm = TRUE)), 
                       ceiling(max(testData_df$residuals, na.rm = TRUE)), 
                       length.out = 5)  # Adjust number of intervals if needed

# Create Spatial Residuals Map with improved title formatting
plot2 <- ggplot(testData_df, aes(x = x, y = y, color = residuals)) +
  geom_point(size = 0.3, alpha = 0.7) +
  scale_color_viridis(option = "D", name = "Residuals", 
                      breaks = residual_breaks, 
                      labels = function(x) sprintf("%.1f", x)) +  # Format labels to 1 decimal
  labs(title = "Spatial Distribution of Residuals", x = "Coordinate X", y = "Coordinate Y") +
  theme_minimal() +
  theme(
    legend.position = "right", 
    panel.grid.major = element_line(color = "gray80"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Title centered, bold, bigger
  ) +
  coord_fixed()  # Ensures the aspect ratio is maintained properly

# Save the plot
ggsave("/gpfs1/work/acostave/FINALRESULTS/ggplot_Residuals_Map_Annotated.png", plot2, width = 8, height = 8, dpi = 300)


####BLINDFRIENDLY PALETTE#########
color1 <- "#1b9e77"  # Verde accesible
color2 <- "#d95f02"  # Naranja accesible
color3 <- "#7570b3"  # Azul accesible

#  Predicted vs Actual Values
dataTest <- data.frame(Actual = testData$fca_, Predicted = predictions)

pred_actual_plot <- ggplot(dataTest, aes(x = Actual, y = Predicted)) +
  geom_point(color = color3, alpha = 0.9, size = 1) +  # Puntos más visibles
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = color2, linewidth = 1.2) +  # Línea ideal
  labs(title = "Predicted vs Actual Values", x = "Actual Values", y = "Predicted Values") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold"))

ggsave("predicted_vs_actual2020.png", plot = pred_actual_plot, width = 8, height = 6, dpi = 300)

# Residuals vs Fitted Values
residuals <- testData$fca_ - predictions
residuals_plot <- ggplot(data.frame(Fitted = predictions, Residuals = residuals), aes(x = Fitted, y = Residuals)) +
  geom_point(color = color1, alpha = 0.9, size = 1) +  # Puntos más oscuros
  geom_hline(yintercept = 0, linetype = "dashed", color = color2, linewidth = 1.2) +  # Línea en 0
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold"))

ggsave("residuals_vs_fitted2020.png", plot = residuals_plot, width = 8, height = 6, dpi = 300)

#   Histogram of Residuals
hist_residuals_plot <- ggplot(data.frame(Residuals = residuals), aes(x = Residuals)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = color3, color = "black", alpha = 0.8) +
  geom_density(color = color1, linewidth = 1) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = color2, linewidth = 1.2) +  # Línea en 0
  labs(title = "Histogram and Density of Residuals", x = "Residuals", y = "Density") +
  theme_bw(base_size = 14)+
  scale_fill_manual(name = "Legend", values = c("Residuals" = "#7570b3")) +
  theme(legend.position = "bottom")


ggsave("histogram_residuals2020.png", plot = hist_residuals_plot, width = 8, height = 6, dpi = 300)

#   Histogram of fca
hist_fca_plot <- ggplot(data.frame(FCA = testData$fca_), aes(x = FCA)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = color3, color = "black", alpha = 0.8) +
  geom_density(color = color1, linewidth = 1) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = color2, linewidth = 1.2) +  # Línea en 0
  labs(title = "Histogram and Density of FCA", x = "FCA", y = "Density") +
  theme_bw(base_size = 14)+
  scale_fill_manual(name = "Legend", values = c("FCA Values" = "#7570b3")) +
  theme(legend.position = "bottom")


ggsave("hist_fca_plot2020.png", hist_fca_plot, width = 8, height = 6, dpi = 300)

#  Normal Q-Q Plot of Residuals
qq_plot_residuals <- ggplot(data.frame(sample = residuals), aes(sample = sample)) +
  stat_qq(color = color3, size = 2) +
  stat_qq_line(color = color2, linewidth = 1.2) +  # Línea de referencia
  labs(title = "Normal Q-Q Plot of Residuals") +
  theme_bw(base_size = 14)+
  scale_color_manual(name = "Legend", values = c("Residuals" = "#7570b3")) +
  theme(legend.position = "bottom")

ggsave("qq_plot_residuals2020.png", plot = qq_plot_residuals, width = 8, height = 6, dpi = 300)

#  Normal Q-Q plot of FCA
qq_plot_fca <- ggplot(data.frame(sample = testData$fca_), aes(sample = sample)) +
  stat_qq(color = color1, size = 2) +
  stat_qq_line(color = color2, linewidth = 1.2) +  # Línea de referencia
  labs(title = "Normal Q-Q Plot of FCA") +
  theme_bw(base_size = 14)+
  scale_color_manual(name = "Legend", values = c("FCA Values" = "#1b9e77")) +
  theme(legend.position = "bottom")


ggsave("qq_plot_fca2020.png", plot = qq_plot_fca, width = 8, height = 6, dpi = 300)



# Variable importance
importance_values <- final_model$variable.importance

soil_vars <- grep("^soil_", names(importance_values), value = TRUE)
elev_vars <- grep("^elevation_", names(importance_values), value = TRUE)
dominant_vars <- grep("^dominant_sp_", names(importance_values), value = TRUE)
#growingreg_vars <- grep("^dominant_sp_", names(importance_values), value = TRUE)

soil_importance <- sum(importance_values[soil_vars])
elevation_importance <- sum(importance_values[elev_vars])
dominantsp_importance <- sum(importance_values[dominant_vars])
#growingreg_importance <- sum(importance_values[growingreg_vars])

num_soil <- length(soil_vars)    
num_elevation <- length(elev_vars)  
num_dominant_sp <- length(dominant_vars)
#num_growingreg <- length(growingreg_vars)

soil_importance <- soil_importance / num_soil
elevation_importance <- elevation_importance / num_elevation
dominantsp_importance <- dominantsp_importance / num_dominant_sp
#growingreg_importance <- growingreg_importance / num_growingreg


continuous_importance <- importance_values[c("entropy_", "slope_", "aspect_", "twi_", "shannon_41x41", "forestcover_")]

variable_importance <- data.frame(
  Variable = c("entropy_", "slope_", "aspect_", "twi_", "shannon_41x41", 
               "dominant_sp_", "soil_", "elevation_", "forestcover_"),
  Importance = c(
    continuous_importance["entropy_"],
    continuous_importance["slope_"],
    continuous_importance["aspect_"],
    continuous_importance["twi_"],
    continuous_importance["shannon_41x41"],
    continuous_importance["forestcover_"],
    dominantsp_importance,
    soil_importance,
    elevation_importance
  )
)

variable_importance$Importance <- variable_importance$Importance / sum(variable_importance$Importance) * 100


cat("Variable Importance:\n")
print(variable_importance)


variable_importance <- variable_importance[order(-variable_importance$Importance), ]

importance_plot <- ggplot(variable_importance, aes(x = reorder(Variable, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = 'identity', color = "black", width = 0.7) +  
  coord_flip() +  
  scale_fill_viridis(option = "D", direction = -1) +  
  geom_text(aes(label = round(Importance, 1)), hjust = -0.2, size = 4) +  # Mostrar valores en cada barra
  labs(
    title = "Variable Importance in Random Forest",
    x = "Variables",
    y = "Importance (%)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  
  )

ggsave("variable_importance2020.png", plot = importance_plot, width = 8, height = 6, dpi = 300)
print(importance_plot)
write.csv(variable_importance, "variable_importance.csv", row.names = FALSE)



cat("Modeling process completed successfully.\n")
