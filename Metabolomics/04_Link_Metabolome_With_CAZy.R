rm(list = ls())

# This script investigates the relative abundances of 
# breakdown products of cellulose hemicellulose and liugnin

# Load packages
library(data.table)
library(ggplot2)
library(cowplot)
library(randomForest)

#### Data Import and Cleaning #### 

# Load in the data
mf <- data.frame(fread("Metabolomics/Data/Metabolome_Blackout_Raw.csv"))
# of <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_Orders_CAZymes.rds")
of <- data.frame(fread("CAZy/Cleaned_Data/Henfaes_CAZyFamilyRichness.csv"))


head(mf)
head(of)
head(sf[1:5])

mf <- mf[, c(1,7:20)] # Keep only the most useful columns

# Lignocellulose breakdown products
BreakdownProducts <- c('glucose'
                       , 'xylose', 'fucose', '3,6-anhydro-D-galactose'
                       , 'vanillic acid', '4-hydroxybenzoic acid', 'benzoic acid')

# Keep only those products
mf <- mf[mf$BinBase.name %in% BreakdownProducts,]

head(mf)

head(mf)
AbundData <- data.frame(t(mf[2:length(mf)]))
names(AbundData) <- mf[,1]
AbundData$Sample <- row.names(AbundData)


ofm <- merge(of, AbundData, by = "Sample")
head(ofm[1:5])
head(ofm[,300:311])
#names(ofm)[44] <- "Bacteroidetes_Order_II"
# Save this as a file

#write.csv(ofm, "Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.csv", row.names = FALSE)

saveRDS(ofm, file = "Metabolomics/Cleaned Data/Henfaes_CAZy_Metabolome.rds")

#### Analysis of correlations between orders and metabolites ####

rm(list = ls())
source("Functions/CrossValidate_rf.R") # Load a function for cross-validated random forests
ofm <- readRDS("Metabolomics/Cleaned Data/Henfaes_CAZy_Metabolome.rds")
#ofm <- data.frame(fread("Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.csv"))


#### Find Orders which were not rare              #### 
#### (median abundance >10%ile of all abundances) ####

head(ofm[1:10])


# Include <- lapply(ofm[4:05], FUN = median) > quantile(unlist(ofm[4:305]), 0.1)

# Names of all orders
CAZymes <- names(ofm[4:304])

# Which of these are above the 10% threshold for further analysis
CAZymesToInclude <- CAZymes#[CAZymes %in% names(Include[Include == TRUE])]

# Keep these
CAZyFrame <- ofm[, names(ofm) %in% c(CAZymesToInclude)]
UnvaryingPredictors <- which(apply(CAZyFrame, 2, function(x){min(x) - max(x)}) == 0)
CAZyFrame <- CAZyFrame[, -UnvaryingPredictors]

## Make all the useful data frames
# Cellulose
df_glucose <- cbind(CAZyFrame, ofm$glucose)
  #ofm[, names(ofm) %in% c(CAZymesToInclude, 'glucose')]

# Hemicellulose
df_xylose <- cbind(CAZyFrame, ofm$xylose)# ofm[, names(ofm) %in% c(CAZymesToInclude, 'xylose')]
df_fucose <- cbind(CAZyFrame, ofm$fucose)# ofm[, names(ofm) %in% c(CAZymesToInclude, 'fucose')]
df_galactose <- cbind(CAZyFrame, ofm$galactose)#<- ofm[, names(ofm) %in% c(CAZymesToInclude, '3,6-anhydro-D-galactose')]

# Lignin
df_vanillin <- cbind(CAZyFrame, ofm$`vanillic acid`)# ofm[, names(ofm) %in% c(CAZymesToInclude, 'vanillic acid')]
df_hydroxybenz <- cbind(CAZyFrame, ofm$`4-hydroxybenzoic acid`)# ofm[, names(ofm) %in% c(CAZymesToInclude, '4-hydroxybenzoic acid')]
df_benz <- cbind(CAZyFrame, ofm$`benzoic acid`)# ofm[, names(ofm) %in% c(CAZymesToInclude, 'benzoic acid')]

# Get results for which orders were most important for lignocellulose degradation
#     along with a measure of the variability in the predictability of those 
#     measurements

#### Use parallel processing for this as it is fairly intensive ####

library(doParallel)
cl <- makePSOCKcluster(4)
registerDoParallel(cl)

#### Random Forests for cellulose ####

Glucose_RF <- CrossValidate_rf(y = "glucose", data = df_glucose, k = 3, repeats = 50)

# Save the output so it can be used/viewed again
saveRDS(Glucose_RF, file = "Metabolomics/Results/Glucose_CAZymes_RandomForest.rds")

#### Random Forests for hemicellulose ####

Xylose_RF <- CrossValidate_rf(y = "xylose", data = df_xylose, k = 3, repeats = 50)
Fucose_RF <- CrossValidate_rf(y = "fucose", data = df_fucose, k = 3, repeats = 50)
Galactose_RF <- CrossValidate_rf(y = "`3,6-anhydro-D-galactose`", data = df_galactose, k = 3, repeats = 50)

# Save the output so it can be used/viewed again
saveRDS(Xylose_RF, file = "Metabolomics/Results/Xylose_CAZymes_RandomForest.rds")
saveRDS(Fucose_RF, file = "Metabolomics/Results/Fucose_CAZymes_RandomForest.rds")
saveRDS(Galactose_RF, file = "Metabolomics/Results/Galactose_CAZymes_RandomForest.rds")


#### Random Forests for Lignin #### 

Vanillin_RF <- CrossValidate_rf(y = "`vanillic acid`", data = df_vanillin, k = 3, repeats = 50)
Hydroxybenz_RF <- CrossValidate_rf(y = "`4-hydroxybenzoic acid`", data = df_hydroxybenz, k = 3, repeats = 50)
Benz_RF <- CrossValidate_rf(y = "`benzoic acid`", data = df_benz, k = 3, repeats = 50)

# Save the output so it can be used/viewed again
saveRDS(Vanillin_RF, file = "Metabolomics/Results/VanillicAcid_CAZymes_RandomForest.rds")
saveRDS(Hydroxybenz_RF, file = "Metabolomics/Results/4-HydroxybenzoicAcid_CAZymes_RandomForest.rds")
saveRDS(Benz_RF, file = "Metabolomics/Results/BenzoicAcid_CAZymes_RandomForest.rds")


######## Get Random Forest Results ########

# Cellulose breakdown
# Get the results
Glucose_RF$rf_best                    # r^2 = -2%, 2/302 predictors
Glucose_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Glucose_RF$mtry_Plot                  # Visualization of this
varImpPlot(Glucose_RF$rf_best)        # Vizualization of importance of results

### Glucose
# probably not the best predictor of cellulose degradation (especially in control plots)
# Might be useful. 35 predictors in best model R^2  = 0.39

# Hemicellulose breakdown
### Xylose
# Get the results
Xylose_RF$rf_best                    # r^2 = 59%, 292/301 predictors
Xylose_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Xylose_RF$mtry_Plot                  # Visualization of this
varImpPlot(Xylose_RF$rf_best)        # Vizualization of importance of results

# It appears that the Acidithiobacillales, Acidobacteriales, Pseudomonadalesand  
#    are the most important orders for xylose production (hemicellulose breakdown)
#    Used all 42 predictors for the best model. Suggests that maybe 

### fucose
# Get the results
Fucose_RF$rf_best                    # r^2 = 48%, 247/301 predictors
Fucose_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Fucose_RF$mtry_Plot                  # Visualization of this
varImpPlot(Fucose_RF$rf_best)        # Vizualization of importance of results

### 3,6-anhydro-D-galactose
# Get the results
Galactose_RF$rf_best                    # r^2 = 14%
Galactose_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Galactose_RF$mtry_Plot                  # Visualization of this
varImpPlot(Galactose_RF$rf_best)        # Vizualization of importance of results

#### Lignin Breakdown

### Vanillic Acid
# Get the results
Vanillin_RF$rf_best                    # r^2 = 14%
Vanillin_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Vanillin_RF$mtry_Plot                  # Visualization of this
varImpPlot(Vanillin_RF$rf_best)        # Vizualization of importance of results

### 4-Hydroxybenzoic Acid
# Get the results
Hydroxybenz_RF$rf_best                    # r^2 = 14%
Hydroxybenz_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Hydroxybenz_RF$mtry_Plot                  # Visualization of this
varImpPlot(Hydroxybenz_RF$rf_best)        # Vizualization of importance of results

### Benzoic Acid
# Get the results
Benz_RF$rf_best                    # r^2 = 14%
Benz_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Benz_RF$mtry_Plot                  # Visualization of this
varImpPlot(Benz_RF$rf_best)        # Vizualization of importance of results