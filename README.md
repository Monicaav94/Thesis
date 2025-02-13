# Thesis
# 🌳 Forest Condition Anomaly (FCA) Prediction - Thesis Repository

## 📖 Overview
This repository contains the code and results for my **Master’s thesis**, where I applied machine learning techniques to predict **Forest Condition Anomalies (FCA)** using environmental variables. The study focuses on analyzing the factors influencing FCA across German forests.

## 📂 Repository Structure
- **📁 Codes/**: R script used for data processing, model training, and evaluation.
- **📁 Data/**: Processed datasets used in the study.
- **📁 plots/**:  Visualizations such as histograms, residual plots, and variable importance.
- **📁 results/**: Random Forest model results.

## 🖥️ Methodology
- **Machine Learning Model:** Random Forest (Ranger package)
- **Environmental Variables:** Soil properties, elevation, biodiversity indices, forest fragmentation (entropy), slope, aspect, and Shannon diversity index.
- **Evaluation Metrics:** R², RMSE, residual analysis, and spatial autocorrelation (Moran’s I).

## 📊 Key Findings
- The Random Forest model achieved an **R² ≈ 0.08**, indicating that FCA is influenced by additional environmental factors beyond those included in the model.
- Residuals showed **weak spatial autocorrelation**, suggesting that FCA variations are primarily driven by **localized environmental stressors**.
- A **linear regression model** was also tested, but it showed **even lower predictive power (R² ≈ 0.03)**.
- A **classification approach** based on FCA thresholds (Lange et al., 2024) achieved **41% accuracy**, reinforcing the complexity of FCA prediction.
- Integration of **Climate data, management Data, specific tree attributes and evaluating FCA gradientes over time**, will improve the prediction of FCA values.

## 🔧 How to Run the Code
To run the analysis, ensure you have the required R libraries installed:
```r
install.packages(c("raster", "sf", "dplyr", "caret", "ranger", "ggplot2", "sp", "gstat", "corrplot", "spdep", "ape"))

Since the original raster files are not included in this repository due to their size, the .RData files contain pre-processed data used in the study:

    fully_cleaned_data.RData: Contains the cleaned dataset used for training and testing the models (fully_cleaned_data).
    traindata.RData: The training subset (trainData) used to fit the Random Forest model.
    testdata.RData: The test subset (testData) used to evaluate model performance.
    forestmodel.RData: Stores the trained Random Forest model (randomForest) for further analysis.
