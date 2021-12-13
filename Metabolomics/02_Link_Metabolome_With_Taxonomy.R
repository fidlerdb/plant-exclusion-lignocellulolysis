rm(list = ls())

# This script investigates the relationships between relative
# abundances of lignocellulose breakdown products and the 
# relative abundances of different microbial orders and species,
# total microbial abundance, species richness, and Shannon's H metrics 

# Load packages
library(data.table)
library(ggplot2)
library(ggforce)
library(cowplot)
#library(randomForest)

#### Data Import and Cleaning #### 

# Load in the data
mf <- data.frame(fread("Metabolomics/Data/Metabolome_Blackout_Raw.csv"))
of <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_Orders_CAZymes.rds")
sf <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_Species_CAZymes.rds")

head(mf) # Metabolites
head(of) # Orders
head(sf) # Species
length(sf)
head(sf[1:5])
min(sf[3:length(sf)])
hist(rowSums(sf[3:length(sf)]))

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
names(ofm)[44] <- "Bacteroidetes_Order_II"
# Save this as a file

#write.csv(ofm, "Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.csv", row.names = FALSE)

saveRDS(ofm, file = "Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.rds")

# Species variation
sfm <- merge(sf, AbundData, by = "Sample")
saveRDS(sfm, file = "Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome_Species.rds")

#### Analysis of correlations between orders and metabolites ####

rm(list = ls())
#source("Functions/CrossValidate_rf.R") # Load a function for cross-validated random forests
source("Functions/range01.R")
ofm <- readRDS("Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.rds")
#ofm <- data.frame(fread("Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.csv"))

names(ofm)
# Scale all response and predictor vaviables between 0 and 1
# Metabolites
ofm[,49:55] <- apply(ofm[,49:55], 2,range01)

# orders
length(3:48)
any(ofm[,3:48] == 0)
ofm[,3:48] <- apply(ofm[,3:48], 2, log2)
ofm[,3:48] <- apply(ofm[,3:48], 2,range01)
hist(ofm[51], breaks = 200)

###### Now fit individual models to search for correlations as simpler is better. #####

library(pwr)
mod <- lm(log10(ofm[,49]+1) ~ 
            ofm[,3]
          , data = ofm)
summary(mod)
# Quick power analysis
pwr.f2.test(u = 1, v = 12, f2 = 0.02, sig.level = 0.05)
pwr.f2.test(u = 1, v = 12, f2 = 0.15, sig.level = 0.05)
pwr.f2.test(u = 1, v = 12, f2 = 0.35, sig.level = 0.05)

# Create a function which gets the results for each metabolite
GetMetabResults <- function(ORDER){
  MetabResults <- lapply(seq_along(ofm[,49:55])
                         , FUN = function(METABOLITE){
                           # Fit a model for each metabolite
                           mod <- lm(log10(unlist(ofm[,49:55][,METABOLITE]+1)) ~ 
                                       ORDER
                                     , data = ofm)
                           # Store the summary
                           Result <- summary(mod)
                           
                           # Extract metabolite name, p value, and r value
                           data.frame(Metabolite = names(ofm[,49:55][METABOLITE])
                                      , p = anova(mod)[5][1,]
                                      , r = Result$coefficients[2,1])
                         })
  rbindlist(MetabResults)
}

MetaboliteResults <- rbindlist(lapply(ofm[, 3:48]
                                      , FUN = function(ORDER){
                                        GetMetabResults(ORDER)}
                                      ), idcol = "Order"
                               )
head(MetaboliteResults)
MetaboliteResults$p_adj <- double(length = nrow(MetaboliteResults))

for(i in MetaboliteResults$Metabolite){
  MetaboliteResults[MetaboliteResults$Metabolite == i,]$p_adj <- 
    p.adjust(MetaboliteResults[MetaboliteResults$Metabolite == i,]$p
             , method = 'fdr')
}



# Metabolite-wise p adjustment
# MetaboliteResults$p_adj <- apply(MetaboliteResults, 1
#                                  , FUN = function(x){
#                                    p.adjust(x[3], method = "fdr"
#                                             , n = length(unique(MetaboliteResults$Order)))
#                                    })

MetaboliteResults[MetaboliteResults$p_adj < 0.05,] # No significant correlations
siggy <- MetaboliteResults[MetaboliteResults$p < 0.05,] # 14 significant correlations
siggy[order(siggy$Metabolite),]


#####
hist(MetaboliteResults$p_adj, breaks = 100)
nrow(MetaboliteResults)
BlackoutPalette <- c('#c2a5cf'
                     ,'#7b3294'
                     ,'#a6dba0'
                     ,'#008837')

plot(xylose ~ Acidithiobacillales, data = ofm
     , col = BlackoutPalette[ofm$Treatment])
plot(fucose~ Acidithiobacillales, data = ofm
     , col = BlackoutPalette[ofm$Treatment])


########
#### Analysis of correlations between species and metabolites ####

rm(list = ls())
#source("Functions/CrossValidate_rf.R") # Load a function for cross-validated random forests
source("Functions/range01.R")
sfm <- readRDS("Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome_Species.rds")
#ofm <- data.frame(fread("Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.csv"))

names(sfm[1:10])
names(sfm[,(length(sfm) - 10) : length(sfm)])
names(sfm[1708]) # Last species

# 7 metabolites at the end. 1706 species

# Scale all response and predictor vaviables between 0 and 1
# metabolites
sfm[,1709:length(sfm)] <- apply(sfm[,1709:length(sfm)], 2,range01)
# species
# log2+1 transform, then put them all on the same relative abundance scale
# for easy visualization
sfm[,3:1708] <- apply(sfm[,3:1708]+1, 2, log2)
sfm[,3:1708] <- apply(sfm[,3:1708], 2, range01)

any(is.na(sfm[,3:1708]))
hist(sfm[51], breaks = 200)

# Now fit individual models to search for correlations as simpler is better.

# Create a function which gets the results for each metabolite
GetMetabResults <- function(SPECIES){
  MetabResults <- lapply(seq_along(sfm[,1709:length(sfm)])
                         , FUN = function(METABOLITE){
                           # Fit a model for each metabolite
                           mod <- lm(unlist(sfm[,1709:length(sfm)][,METABOLITE]) ~ 
                                       SPECIES
                                     , data = sfm)
                           # Store the summary
                           Result <- summary(mod)
                           
                           # Extract metabolite name, p value, and r value
                           data.frame(Metabolite = names(sfm[,1709:length(sfm)][METABOLITE])
                                      , p = anova(mod)[5][1,]
                                      , r = Result$coefficients[2,1])
                         })
  rbindlist(MetabResults)
}

MetaboliteResults <- rbindlist(lapply(sfm[, 3:1708]
                                      , FUN = function(SPECIES){
                                        GetMetabResults(SPECIES)}
                                      ), idcol = "Species"
                               )

# Metabolite-wise p adjustment
MetaboliteResults$p_adj <- double(length = nrow(MetaboliteResults))
for(i in MetaboliteResults$Metabolite){
  MetaboliteResults[MetaboliteResults$Metabolite == i,]$p_adj <- 
    p.adjust(MetaboliteResults[MetaboliteResults$Metabolite == i,]$p
             , method = 'fdr')
}

MetaboliteResults[MetaboliteResults$p_adj < 0.05,] # No significant correlations
MetaboliteResults[MetaboliteResults$p < 0.05,] # No significant correlations

hist(MetaboliteResults$p_adj, breaks = 100)
nrow(MetaboliteResults)
1706*7

########

# Is there a relationship between total microbial abundance and 
# any breakdown product?

AllReadsResults <- lapply(seq_along(sfm[,1709:length(sfm)])
       , FUN = function(METABOLITE){
         # Fit a model for each metabolite
         mod <- lm(unlist(sfm[,1709:length(sfm)][,METABOLITE]) ~ 
                     rowSums(sfm[3:1708])
                   )
         # Store the summary
         Result <- summary(mod)
         
         # Extract metabolite name, p value, and r value
         data.frame(Metabolite = names(sfm[,1709:length(sfm)][METABOLITE])
                    , p = anova(mod)[5][1,]
                    , r = Result$coefficients[2,1])
       })

rbindlist(AllReadsResults) # No.

# Is there a relationship between microbial species richness and 
# any breakdown product?

vegan::specnumber(sfm[3:1708])
RichnessResults <- lapply(seq_along(sfm[,1709:length(sfm)])
                          , FUN = function(METABOLITE){
                            # Fit a model for each metabolite
                            mod <- lm(unlist(sfm[,1709:length(sfm)][,METABOLITE]) ~ 
                                        vegan::specnumber(sfm[3:1708])
                            )
                            # Store the summary
                            Result <- summary(mod)
                            
                            # Extract metabolite name, p value, and r value
                            data.frame(Metabolite = names(sfm[,1709:length(sfm)][METABOLITE])
                                       , p = anova(mod)[5][1,]
                                       , r = Result$coefficients[2,1])
                          })

rbindlist(RichnessResults) # No.

# Is there a relationship between microbial species richness and 
# any breakdown product?

DiversityResults <- lapply(seq_along(sfm[,1709:length(sfm)])
                          , FUN = function(METABOLITE){
                            # Fit a model for each metabolite
                            mod <- lm(unlist(sfm[,1709:length(sfm)][,METABOLITE]) ~ 
                                        vegan::diversity(sfm[3:1708], index = 'shannon')
                            )
                            # Store the summary
                            Result <- summary(mod)
                            
                            # Extract metabolite name, p value, and r value
                            data.frame(Metabolite = names(sfm[,1709:length(sfm)][METABOLITE])
                                       , p = anova(mod)[5][1,]
                                       , r = Result$coefficients[2,1])
                          })

rbindlist(DiversityResults) # No.







######## OLD #####


# 
# 
# ofm[3:48]
# par(mfrow = c(2,2));plot(lm(as.matrix(ofm[,49:55]) ~ Corynebacteriales, data = ofm));par(mfrow = c(1,1))
# anova(lm(as.matrix(ofm[,49:55]) ~ Corynebacteriales, data = ofm))
# 
# 
# 
# pca_hem <- prcomp(ofm[names(ofm) %in% 
#          c('xylose', 'fucose', '3,6-anhydro-D-galactose')]
#          , scale. = TRUE)
# 
# pca_lig <- prcomp(ofm[names(ofm) %in% 
#                           c('vanillic acid', '4-hydroxybenzoic acid', 'benzoic acid')]
#                   , scale. = TRUE)
# 
# plot(pca_hem)
# plot(pca_lig)
# 
# plot(pca_hem$x, col = ofm$Treatment, pch = 16)
# plot(pca_lig$x, col = ofm$Treatment, pch = 16)
# 
# ofm$Hem_deg <- pca_hem$x[,1]
# ofm$Lig_deg <- pca_lig$x[,1]
# 
# plot(Hem_deg ~ log10(Streptomycetales), data = ofm)
# plot(Lig_deg ~ log10(Streptomycetales), data = ofm)
# 
# names(ofm)
# ofm$Shannon_H <- vegan::diversity(ofm[3:8], index = 'shannon')
# #ofm$Shannon_H_Species <- vegan::diversity(sf[3:length(sf)], index = 'shannon')
# 
# # 
# 
# ofm
# anova(lm(as.matrix(ofm[,49:55]) ~ Corynebacteriales, data = ofm))
# 
# 
# plot(glucose ~ Shannon_H, data = ofm)
# plot(Hem_deg ~ Shannon_H, data = ofm);abline(lm(Hem_deg ~ Shannon_H, data = ofm))
# plot(Lig_deg ~ Shannon_H, data = ofm)
# 
# # plot(glucose ~ Shannon_H_Species, data = ofm)
# # plot(Hem_deg ~ Shannon_H_Species, data = ofm);abline(lm(Hem_deg ~ Shannon_H_Species, data = ofm))
# # plot(Lig_deg ~ Shannon_H_Species, data = ofm)
# 
# plot(Shannon_H ~ Treatment, data = ofm)
# 
# # Significant negative correlation between hemicellulose degradation product
# # abundance and diversity
# drop1(lm(Hem_deg ~ Shannon_H, data = ofm), test = 'F')
# drop1(lm(xylose ~ Shannon_H, data = ofm), test = 'F')
# drop1(lm(fucose ~ Shannon_H, data = ofm), test = 'F') 
# drop1(lm(`3,6-anhydro-D-galactose` ~ Shannon_H, data = ofm), test = 'F')
# # Does not apply to 3,6-anhydro-D-galactose
# 
# names(ofm)
# 
# #### Which orders are correlated? ####
# CorrMatrix <- Hmisc::rcorr(as.matrix(ofm[3:48]))
# CorrMatrix
# 
# CorrMatrix$r[CorrMatrix$P < 0.05]
# 
# # Which have an r value greater than 5?
# apply(CorrMatrix$r, 2, FUN = function(x){
#   pos <- x[x > 0.7]
#   neg <- x[x < -0.7]
#   c(pos,neg)
# })

# Many of the orders are highly correlated. Weak a-priori hypotheses so 
# go with it and see what comes out as important later on--then try to 
# understand if this is sensible later on (other lines of evidence. e.g. 
# CAZyme richness/composition in these orders, geochip, cultivation work)


head(ofm)

#### Find Orders which were not rare              #### 
#### (median abundance >10%ile of all abundances) ####

Include <- lapply(ofm[3:48], FUN = median) > quantile(unlist(ofm[3:48]), 0.1)

# Names of all orders
Orders <- names(ofm[3:48])

# Which of these are above the 10% threshold for further analysis
OrdersToInclude <- Orders[Orders %in% names(Include[Include == TRUE])]

## Make all the useful data frames
# Cellulose
df_glucose <- ofm[, names(ofm) %in% c(OrdersToInclude, 'glucose')]

# Hemicellulose
df_xylose <- ofm[, names(ofm) %in% c(OrdersToInclude, 'xylose')]
df_fucose <- ofm[, names(ofm) %in% c(OrdersToInclude, 'fucose')]
df_galactose <- ofm[, names(ofm) %in% c(OrdersToInclude, '3,6-anhydro-D-galactose')]

# Lignin
df_vanillin <- ofm[, names(ofm) %in% c(OrdersToInclude, 'vanillic acid')]
df_hydroxybenz <- ofm[, names(ofm) %in% c(OrdersToInclude, '4-hydroxybenzoic acid')]
df_benz <- ofm[, names(ofm) %in% c(OrdersToInclude, 'benzoic acid')]

# Get results for which orders were most important for lignocellulose degradation
#     along with a measure of the variability in the predictability of those 
#     measurements


#### Random Forests for cellulose ####
Glucose_RF <- CrossValidate_rf(y = "glucose", data = df_glucose, k = 3, repeats = 50)

# Save the output so it can be used/viewed again
saveRDS(Glucose_RF, file = "Metabolomics/Results/Glucose_RandomForest.rds")

#### Random Forests for hemicellulose ####

Xylose_RF <- CrossValidate_rf(y = "xylose", data = df_xylose, k = 3, repeats = 50)
Fucose_RF <- CrossValidate_rf(y = "fucose", data = df_fucose, k = 3, repeats = 50)
Galactose_RF <- CrossValidate_rf(y = "`3,6-anhydro-D-galactose`", data = df_galactose, k = 3, repeats = 50)

# Save the output so it can be used/viewed again
saveRDS(Xylose_RF, file = "Metabolomics/Results/Xylose_RandomForest.rds")
saveRDS(Fucose_RF, file = "Metabolomics/Results/Fucose_RandomForest.rds")
saveRDS(Galactose_RF, file = "Metabolomics/Results/Galactose_RandomForest.rds")


#### Random Forests for Lignin #### 

Vanillin_RF <- CrossValidate_rf(y = "`vanillic acid`", data = df_vanillin, k = 3, repeats = 50)
Hydroxybenz_RF <- CrossValidate_rf(y = "`4-hydroxybenzoic acid`", data = df_hydroxybenz, k = 3, repeats = 50)
Benz_RF <- CrossValidate_rf(y = "`benzoic acid`", data = df_benz, k = 3, repeats = 50)

# Save the output so it can be used/viewed again
saveRDS(Vanillin_RF, file = "Metabolomics/Results/VanillicAcid_RandomForest.rds")
saveRDS(Hydroxybenz_RF, file = "Metabolomics/Results/4-HydroxybenzoicAcid_RandomForest.rds")
saveRDS(Benz_RF, file = "Metabolomics/Results/BenzoicAcid_RandomForest.rds")


######## Get Random Forest Results ########

# Cellulose breakdown
# Get the results
Glucose_RF$rf_best                    # r^2 = 14%
Glucose_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Glucose_RF$mtry_Plot                  # Visualization of this
varImpPlot(Glucose_RF$rf_best)        # Vizualization of importance of results

### Glucose
# probably not the best predictor of cellulose degradation (especially in control plots)
# Might be useful. 35 predictors in best model R^2  = 0.39

# Hemicellulose breakdown
### Xylose
# Get the results
Xylose_RF$rf_best                    # r^2 = 14%
Xylose_RF$CrossValidation_Results    # Effect of number of predictors on RMSE etc.
Xylose_RF$mtry_Plot                  # Visualization of this
varImpPlot(Xylose_RF$rf_best)        # Vizualization of importance of results

# It appears that the Acidithiobacillales, Acidobacteriales, Pseudomonadalesand  
#    are the most important orders for xylose production (hemicellulose breakdown)
#    Used all 42 predictors for the best model. Suggests that maybe 

### fucose
# Get the results
Fucose_RF$rf_best                    # r^2 = 14%
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

###########################################



# cowplot::plot_grid(Glucose_Importance, Xylose_Importance
#                    , Fucose_Importance, Galactose_Importance
#                    , Vanillin_Importance, Hydroxybenz_Importance
#                    , Benz_Importance
#                    , ncol = 2)

##########################

# mtry = sqrt(p) is worse than mtry = p/3
# https://arxiv.org/pdf/0811.3619.pdf



# Does species richness predict abundance of any of the 
# degradation products?

# ofm$Richness <- rowSums(apply(sf[3:length(sf)]
#               , 2
#               ,  function(x)ifelse(x == 0, yes = 0, no = 1)
# ))
# 
# plot(glucose ~ Richness, data = ofm)
# plot(Hem_deg ~ Richness, data = ofm)
# plot(Lig_deg ~ Richness, data = ofm)
# 
# drop1(lm(Lig_deg ~ Richness, data = ofm), test = 'F')
# drop1(lm(glucose ~ Richness, data = ofm), test = 'F')
# drop1(lm(Hem_deg ~ Richness, data = ofm), test = 'F')
# 
# # What about specific compounds
# drop1(lm(xylose ~ Richness, data = ofm), test = 'F')
# drop1(lm(fucose ~ Richness, data = ofm), test = 'F') 
# drop1(lm(`3,6-anhydro-D-galactose` ~ Richness, data = ofm), test = 'F')
# # Not for hemicellulose...
# drop1(lm(`vanillic acid` ~ Richness, data = ofm), test = 'F')
# drop1(lm(`4-hydroxybenzoic acid` ~ Richness, data = ofm), test = 'F') 
# drop1(lm(`benzoic acid` ~ Richness, data = ofm), test = 'F')
# # And not for lignin
