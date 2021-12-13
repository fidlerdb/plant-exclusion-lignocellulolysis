rm(list = ls())

library(randomForest)

Hydroxybenz_RF <- readRDS("Metabolomics/Results/4-HydroxybenzoicAcid_RandomForest.rds")
Benz_RF <- readRDS("Metabolomics/Results/BenzoicAcid_RandomForest.rds")
Fucose_RF <- readRDS("Metabolomics/Results/Fucose_RandomForest.rds")
Galactose_RF <- readRDS("Metabolomics/Results/Galactose_RandomForest.rds")
Glucose_RF <- readRDS("Metabolomics/Results/Glucose_RandomForest.rds")
Vanillin_RF <- readRDS("Metabolomics/Results/VanillicAcid_RandomForest.rds")
Xylose_RF <- readRDS("Metabolomics/Results/Xylose_RandomForest.rds")

######## Get Random Forest Results ########

# Cellulose breakdown
# Get the results
Glucose_RF$rf_best                    # r^2 = 13.49%
Glucose_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Glucose_RF$mtry_Plot                  # Visualization of this
varImpPlot(Glucose_RF$rf_best)        # Vizualization of importance of results

getTree(Glucose_RF$rf_best, 1, labelVar = TRUE)
getTree(Glucose_RF$rf_best, 2, labelVar = TRUE)
getTree(Glucose_RF$rf_best, 3, labelVar = TRUE)
getTree(Glucose_RF$rf_best, 4, labelVar = TRUE)
getTree(Glucose_RF$rf_best, 5, labelVar = TRUE)

### Glucose
# probably not the best predictor of cellulose degradation (especially in control plots)
# Might be useful. 35 predictors in best model R^2  = 0.39

# Hemicellulose breakdown
### Xylose
# Get the results
Xylose_RF$rf_best                    # r^2 = 62.73%, 42 predictors
Xylose_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Xylose_RF$mtry_Plot                  # Visualization of this
varImpPlot(Xylose_RF$rf_best)        # Vizualization of importance of results

# It appears that the Acidithiobacillales, Acidobacteriales, Pseudomonadalesand  
#    are the most important orders for xylose production (hemicellulose breakdown)
#    Used all 42 predictors for the best model. Suggests that maybe 

### fucose
# Get the results
Fucose_RF$rf_best                    # r^2 = 54.2%, 41 predictors
Fucose_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Fucose_RF$mtry_Plot                  # Visualization of this
varImpPlot(Fucose_RF$rf_best)        # Vizualization of importance of results

### 3,6-anhydro-D-galactose
# Get the results
Galactose_RF$rf_best                    # r^2 = -30.69%, 16 predictors
Galactose_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Galactose_RF$mtry_Plot                  # Visualization of this
varImpPlot(Galactose_RF$rf_best)        # Vizualization of importance of results

#### Lignin Breakdown

### Vanillic Acid
# Get the results
Vanillin_RF$rf_best                    # r^2 = -25%, 26 predictors
Vanillin_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Vanillin_RF$mtry_Plot                  # Visualization of this
varImpPlot(Vanillin_RF$rf_best)        # Vizualization of importance of results

### 4-Hydroxybenzoic Acid
# Get the results
Hydroxybenz_RF$rf_best                    # r^2 = 12%, 42 predictors
Hydroxybenz_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Hydroxybenz_RF$mtry_Plot                  # Visualization of this
varImpPlot(Hydroxybenz_RF$rf_best)        # Vizualization of importance of results

### Benzoic Acid
# Get the results
Benz_RF$rf_best                    # r^2 = -36%, 2 predictors
Benz_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Benz_RF$mtry_Plot                  # Visualization of this
varImpPlot(Benz_RF$rf_best)        # Vizualization of importance of results
