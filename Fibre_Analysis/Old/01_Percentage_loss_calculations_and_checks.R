
# Clean workspace 
rm( list = ls())
gc(full = TRUE)

library(data.table)

#### Import data ####
df <- fread("../Data/Fibre Analysis/Ankom_Raw_Data.csv")

df

#### Calculate probable actual mass of c7 sample ####
# added into bag, based on difference of others in the control group

# View the distribution of differences for all samples and the control group
par(mfrow = c(2,1))
with(df, hist(na.omit(Mass_added_minus_actual)))
with(df, abline(v = mean(na.omit(Mass_added_minus_actual)), col = "red"))
with(df, abline(v = median(na.omit(Mass_added_minus_actual)), col = "blue"))

with(df[grep("c", Sample)], hist(na.omit(Mass_added_minus_actual)))
with(df[grep("c", Sample)], abline(v = mean(na.omit(Mass_added_minus_actual)), col = "red"))
with(df[grep("c", Sample)], abline(v = median(na.omit(Mass_added_minus_actual)), col = "blue"))
par(mfrow = c(1,1))

# Ratio of sample mass to sample mass added
hist(df[grep("c", Sample), na.omit(Sample_Mass / Sample_Mass_Added)])
with(df[grep("c", Sample)], abline(v = mean(na.omit(Sample_Mass / Sample_Mass_Added)), col = "red"))
with(df[grep("c", Sample)], abline(v = median(na.omit(Sample_Mass / Sample_Mass_Added)), col = "blue"))

# Calculate the mean proportion of sample mass lost between sample weighing and sample+bag weighing
df[grep("c", Sample), mean(na.omit(Sample_Mass / Sample_Mass_Added))]

Mass_added_correction <- df[grep("c", Sample), median(na.omit(Sample_Mass / Sample_Mass_Added))]
df[Sample == "c7", Sample_Mass := Sample_Mass_Added * Mass_added_correction]

#### Calculate blank bag correction factor ####

bf <- df[grep("blank", Sample),]

bf

# What happened to bag mass between procedures?
plot(value ~ variable, data = melt(bf[,c("Sample", "Bag_Mass", "NDF", "ADF")], id = "Sample"))

NDF_blank_Correction <- mean(bf$NDF / bf$Bag_Mass)
ADF_blank_Correction <- mean(bf$ADF / bf$NDF) # Use NDF as t0 as the analysis is sequential

#### Calculate % NDF mass loss ####

# 100 * (dried NDF weight - (Bag tare weight * ADF blank correction)) / Sample mass
df[, Percentage_NDF := 100 * (NDF - (Bag_Mass * NDF_blank_Correction)) / Sample_Mass]
df[, Proportion_NDF := (NDF - (Bag_Mass * NDF_blank_Correction)) / Sample_Mass]


#### Calculate % ADF mass loss ####

# New initial bag mass for the samples post-NDF
with(df, Bag_Mass * NDF_blank_Correction)

df[, Post_NDF_Bag_Mass := Bag_Mass * NDF_blank_Correction]

# What happened to bag mass after NDF
plot(value ~ variable, data = melt(df[,c("Sample", "Bag_Mass", "Post_NDF_Bag_Mass")], id = "Sample"))

# 100 * (dried ADF weight - (Bag tare weight * ADF blank correction)) / Sample mass (input after NDF?)
# If bag mass after ADF then it's not % loss from total.
df[, Percentage_ADF := 100 * (ADF - (Post_NDF_Bag_Mass * ADF_blank_Correction)) / Sample_Mass]
df[, Proportion_ADF := (ADF - (Post_NDF_Bag_Mass * ADF_blank_Correction)) / NDF]


#### Check Sample comparisons ####

af <- df[grep("alfalfa", Sample),]

af$Percentage_NDF # -- 35.5 - 36.9
af$Percentage_ADF # -- 30.7 - 31.3

mean(af$Percentage_NDF) # -- 35.5 - 36.9
mean(af$Percentage_ADF) # -- 30.7 - 31.3

between(x = mean(af$Percentage_NDF), lower = 36.2 - 0.7, upper = 36.2 + 0.7)
between(x = mean(af$Percentage_ADF), lower = 31.0 - 0.3, upper = 31.0 + 0.3)

par(mfrow = c(2,1))
hist(af$Percentage_NDF, xlim = c(0,100)
     , main = "NDF % mass loss"
     , xlab = "% mass loss"
     , ylab = "Number of alfalfa samples")
abline(v = c(30.7, 31.3), col = "red")
hist(af$Percentage_ADF, xlim = c(0,100)
     , main = "ADF % mass loss"
     , xlab = "% mass loss"
     , ylab = "Number of alfalfa samples")
abline(v = c(30.7, 31.3), col = "red")
par(mfrow = c(1,1))

## Hmmmmm.
# TODO -- understand why the % mass loss from the alfalfa does not reflect what it should. Samples fall around
# the expected values

#### Check results ####

# Create a subset with a treatment column
        sf <- df[!grep("alfalfa|blank", Sample),]
sf$Treatment <- as.factor(paste(sub("[0-9]", "", sf$Sample)
                                , ifelse(sub("[a-z]", "", sf$Sample) > 3
                                         , yes = "New", no = "Old"
                                         )
                                , sep= "-")
                          )
sf$Treatment <- relevel(sf$Treatment, ref = "b-Old")

# Write this to an output file
fwrite(sf, "Fibre_Analysis/Output_Data/Fibre_Analysis_Data_Processed.csv")

# View the data distribution
hist(sf$NDF)
hist(sf$ADF)

# From here on I started analysis. Stick to what it says on the tin!

# #### Model % NDF mass loss--redo this after calculation of hemicellulose content ####
# 
# # Plot the % mass loss from NDF treatment
# plot(Percentage_NDF ~ Treatment, data = sf, ylab = "Percentage Mass Loss (NDF)", ylim = c(0,100))
# points(Percentage_NDF ~ Treatment, data = sf, pch = 16)
# 
# # Nonparametric tests suggests no difference--low powered
# kruskal.test(x = sf$Percentage_NDF, g = sf$Treatment)
# 
# library(FSA)
# dunnTest(Percentage_NDF ~ Treatment
#          , data = sf
#          , method = "bh") 
# 
# # According to: 
# # https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13234
# # Fitting a classical binomial model does not deal with proportional data 
# # from continuous measurements well. Try fitting a beta regression.
# # 
# # Also note that logit, loglog, cauchit, probit, and cloglog models were tried
# # These had identical model AIC and p values
# library(betareg)
# library(lmtest)
# NDF_mod_beta_logit <- betareg(Proportion_NDF ~ Treatment, data = sf, link = "logit")
# 
# # Check model assumptions
# par(mfrow = c(2,2));plot(NDF_mod_beta_logit);par(mfrow = c(1,1))
# 
# # Model output
# summary(NDF_mod_beta_logit) # Significant differences at the factor level -- partial wald tests
# lrtest(NDF_mod_beta_logit) # Overall model performs poorly
# waldtest(NDF_mod_beta_logit) # But wald test says it's not so bad?
# NDF_mod_beta_logit$pseudo.r.squared
# 
# #### Get model predictions ####
# 
# Predictions_NDF_m <- data.frame(NDF_predicted = predict(NDF_mod_beta_logit 
#                                         , newdata = data.frame(Treatment = levels(sf$Treatment))
#                                         , type = "response"))
# Predictions_NDF_q <- data.frame(predict(NDF_mod_beta_logit
#                             , newdata = data.frame(Treatment = levels(sf$Treatment))
#                             , type = "quantile", at = c(0.025, 0.5, 0.975)))
# Predictions_NDF_q$Treatment <- factor(levels(sf$Treatment))
# Predictions_NDF_m$Treatment <- factor(levels(sf$Treatment))
# 
# library(ggplot2)
# colstokeep = c("Sample", "Treatment", "Percentage_NDF")
# 
# sf[,..colstokeep]
# Predictions_NDF_q$Percentage_NDF <- 1
# 
# ggplot(sf[,..colstokeep], aes(x = Treatment, y = Percentage_NDF)) + 
#   geom_point(position = position_jitter(0.15)) +
#   geom_errorbar(data = Predictions_NDF_q, aes(x = Treatment, ymin = q_0.025*100, ymax = q_0.975*100))
# 
# ####
# 
# #### ADf analysis--redo when cellulose is calculated ####
# 
# # Plot the % mass loss from ADF treatment
# plot(Percentage_ADF ~ Treatment, data = sf, ylab = "Percentage Mass Loss (ADF)", ylim = c(0,100))
# points(Percentage_ADF ~ Treatment, data = sf, pch = 16)
# 
# # Nonparametric tests suggests no difference--low powered
# kruskal.test(x = sf$Percentage_ADF, g = sf$Treatment)
# 
# library(FSA)
# dunnTest(Percentage_ADF ~ Treatment
#          , data = sf
#          , method = "bh") 
# 
# # According to: 
# # https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13234
# # Fitting a classical binomial model does not deal with proportional data 
# # from continuous measurements well. Try fitting a beta regression.
# # 
# # Also note that logit, loglog, cauchit, probit, and cloglog models were tried
# # These had identical model AIC and p values
# library(betareg)
# library(lmtest)
# ADF_mod_beta_logit <- betareg(Proportion_ADF ~ Treatment, data = sf, link = "logit")
# 
# # Check model assumptions
# par(mfrow = c(2,2));plot(ADF_mod_beta_logit);par(mfrow = c(1,1))
# 
# # Model output
# summary(ADF_mod_beta_logit) # Significant differences at the factor level -- partial wald tests
# lrtest(ADF_mod_beta_logit) # Overall model performs poorly
# waldtest(ADF_mod_beta_logit) # But wald test says it's not so bad?
# ADF_mod_beta_logit$pseudo.r.squared
# 
# #### Get model predictions ####
# 
# Predictions_ADF_m <- data.frame(ADF_predicted = predict(ADF_mod_beta_logit 
#                                                         , newdata = data.frame(Treatment = levels(sf$Treatment))
#                                                         , type = "response"))
# Predictions_ADF_q <- data.frame(predict(ADF_mod_beta_logit
#                                         , newdata = data.frame(Treatment = levels(sf$Treatment))
#                                         , type = "quantile", at = c(0.025, 0.5, 0.975)))
# Predictions_ADF_q$Treatment <- factor(levels(sf$Treatment))
# Predictions_ADF_m$Treatment <- factor(levels(sf$Treatment))
# 
# library(ggplot2)
# colstokeep = c("Sample", "Treatment", "Percentage_ADF")
# 
# sf[,..colstokeep]
# Predictions_ADF_q$Percentage_ADF <- 1
# 
# ggplot(sf[,..colstokeep], aes(x = Treatment, y = Percentage_ADF)) + 
#   geom_point(position = position_jitter(0.15)) +
#   geom_errorbar(data = Predictions_ADF_q, aes(x = Treatment, ymin = q_0.025*100, ymax = q_0.975*100))
# 
# 
# 
# 
# 
# 
# 
