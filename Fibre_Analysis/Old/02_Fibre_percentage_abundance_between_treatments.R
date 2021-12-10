#### Data import ####
rm(list = ls())

sf <- fread("Fibre_Analysis/Output_Data/Fibre_Analysis_Data_Processed.csv")
sf$Treatment <- factor(sf$Treatment)

#### Model % NDF mass loss--redo this after calculation of hemicellulose content ####

# Plot the % mass loss from NDF treatment
plot(Percentage_NDF ~ Treatment, data = sf, ylab = "Percentage Mass Loss (NDF)", ylim = c(0,100))
points(Percentage_NDF ~ Treatment, data = sf, pch = 16)

# Nonparametric tests suggests no difference--low powered
kruskal.test(x = sf$Percentage_NDF, g = sf$Treatment)

library(FSA)
dunnTest(Percentage_NDF ~ Treatment
         , data = sf
         , method = "bh") 

# According to: 
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13234
# Fitting a classical binomial model does not deal with proportional data 
# from continuous measurements well. Try fitting a beta regression.
# 
# Also note that logit, loglog, cauchit, probit, and cloglog models were tried
# These had identical model AIC and p values
library(betareg)
library(lmtest)
NDF_mod_beta_logit <- betareg(Proportion_NDF ~ Treatment, data = sf, link = "logit")

# Check model assumptions
par(mfrow = c(2,2));plot(NDF_mod_beta_logit);par(mfrow = c(1,1))

# Model output
summary(NDF_mod_beta_logit) # Significant differences at the factor level -- partial wald tests
lrtest(NDF_mod_beta_logit) # Overall model performs poorly
waldtest(NDF_mod_beta_logit) # But wald test says it's not so bad?
NDF_mod_beta_logit$pseudo.r.squared

#### Get model predictions ####

Predictions_NDF_m <- data.frame(NDF_predicted = predict(NDF_mod_beta_logit 
                                                        , newdata = data.frame(Treatment = levels(sf$Treatment))
                                                        , type = "response"))
Predictions_NDF_q <- data.frame(predict(NDF_mod_beta_logit
                                        , newdata = data.frame(Treatment = levels(sf$Treatment))
                                        , type = "quantile", at = c(0.025, 0.5, 0.975)))
Predictions_NDF_q$Treatment <- factor(levels(sf$Treatment))
Predictions_NDF_m$Treatment <- factor(levels(sf$Treatment))

library(ggplot2)
colstokeep = c("Sample", "Treatment", "Percentage_NDF")

sf[,..colstokeep]
Predictions_NDF_q$Percentage_NDF <- 1

ggplot(sf[,..colstokeep], aes(x = Treatment, y = Percentage_NDF)) + 
  geom_point(position = position_jitter(0.15)) +
  geom_errorbar(data = Predictions_NDF_q, aes(x = Treatment, ymin = q_0.025*100, ymax = q_0.975*100), width = 0.3) + 
  geom_point(data = Predictions_NDF_m, aes(x = Treatment, y = NDF_predicted*100), size = 5) +
  ylab("Percentage mass loss (NDF)")


####

#### ADf analysis--redo when cellulose is calculated ####

# Plot the % mass loss from ADF treatment
plot(Percentage_ADF ~ Treatment, data = sf, ylab = "Percentage Mass Loss (ADF)", ylim = c(0,100))
points(Percentage_ADF ~ Treatment, data = sf, pch = 16)

# Nonparametric tests suggests no difference--low powered
kruskal.test(x = sf$Percentage_ADF, g = sf$Treatment)

library(FSA)
dunnTest(Percentage_ADF ~ Treatment
         , data = sf
         , method = "bh") 

# According to: 
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13234
# Fitting a classical binomial model does not deal with proportional data 
# from continuous measurements well. Try fitting a beta regression.
# 
# Also note that logit, loglog, cauchit, probit, and cloglog models were tried
# These had identical model AIC and p values
library(betareg)
library(lmtest)
ADF_mod_beta_logit <- betareg(Proportion_ADF ~ Treatment|Treatment, data = sf, link = "logit")

# Check model assumptions
par(mfrow = c(2,2));plot(ADF_mod_beta_logit);par(mfrow = c(1,1))

# Model output
summary(ADF_mod_beta_logit) # Significant differences at the factor level -- partial wald tests
lrtest(ADF_mod_beta_logit) # Overall model performs poorly
waldtest(ADF_mod_beta_logit) # But wald test says it's not so bad?
ADF_mod_beta_logit$pseudo.r.squared

#### Get model predictions ####

Predictions_ADF_m <- data.frame(ADF_predicted = predict(ADF_mod_beta_logit 
                                                        , newdata = data.frame(Treatment = levels(sf$Treatment))
                                                        , type = "response"))
Predictions_ADF_q <- data.frame(predict(ADF_mod_beta_logit
                                        , newdata = data.frame(Treatment = levels(sf$Treatment))
                                        , type = "quantile", at = c(0.025, 0.5, 0.975)))
Predictions_ADF_q$Treatment <- factor(levels(sf$Treatment))
Predictions_ADF_m$Treatment <- factor(levels(sf$Treatment))

library(ggplot2)
colstokeep = c("Sample", "Treatment", "Percentage_ADF")

sf[,..colstokeep]
Predictions_ADF_q$Percentage_ADF <- 1

ggplot(sf[,..colstokeep], aes(x = Treatment, y = Percentage_ADF)) + 
  geom_point(position = position_jitter(0.15)) +
  geom_errorbar(data = Predictions_ADF_q, aes(x = Treatment, ymin = q_0.025*100, ymax = q_0.975*100), width = 0.3) +
  geom_point(data = Predictions_ADF_m, aes(x = Treatment, y = ADF_predicted*100), size = 5)








