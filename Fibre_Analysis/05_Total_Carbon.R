## This script assesses the abundance of different lignocellulosic ##
## elements between treatments ##

library(data.table)
library(ggplot2)

#### Data import ####
rm(list = ls())

cf <- fread("../Data/Fibre Analysis/Total_Carbon_Nitrogen.csv")

cf

#### Merge with Ankom data ####

sf <- fread("Fibre_Analysis/Output_Data/Fibre_Analysis_Data_Processed_2021-04-29.csv")

# Sort out the treatment variable
sf$Treatment <- factor(sf$Treatment)
sf$Sample
levels(sf$Treatment) <- c("1-Year\nBare", "10-Year\nBare", "1-Year\nGrassland", "10-Year\nGrassland")
sf$Treatment <- relevel(sf$Treatment, ref = "10-Year\nBare")

# Colour palette for the treatment variable
BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

# Function for 95% confidence intervals
mean_ci <- function(x){
  mean_se(x, mult = 1.96)
}

# Check how the data look
sf

# Scale by input mass (without ash content)
sf[, Proportion_Hemicellulose := Hemicellulose/(Sample_Mass - Ash)]
sf[, Proportion_Cellulose := Cellulose/(Sample_Mass - Ash)]
sf[, Proportion_Lignin := Lignin/(Sample_Mass - Ash)]

####

cf <- cf[, c("Sample", "Percentage_Carbon")]
af <- cf[sf, on = "Sample"]

plot(Percentage_Carbon*Proportion_Cellulose ~ Treatment, data = af)
plot(Percentage_Carbon*Proportion_Hemicellulose ~ Treatment, data = af)
plot(Percentage_Carbon*Proportion_Lignin ~ Treatment, data = af)

cidata <- af[, mean_ci(Percentage_Carbon), by = Treatment]
cidata[,Percentage_Carbon := 1]

ggplot(af, aes(x = Treatment, y = Percentage_Carbon, fill = Treatment)) + 
  stat_summary(fun = "mean", geom = 'bar') +
  geom_errorbar(data = cidata, aes(ymin = ymin, ymax = ymax), width = 0.3, colour = "black") + 
  #stat_summary(fun.data = "mean_ci", size = 1.1) +
  geom_point(colour = "black", position = position_jitter(0.1)) +
  scale_fill_manual(values = BlackoutPalette
                    , breaks = levels(af$Treatment)) +
  theme_classic() +
  theme(legend.position = "none") +
  ylim(0,5) +
  ylab("Total carbon (% mass)")

kruskal.test(af$Percentage_Carbon, af$Treatment)
library(FSA)
dunnTest(af$Percentage_Carbon, af$Treatment)

af$Treatment <- relevel(af$Treatment, ref = "10-Year\nBare")

mod1 <- lm(Percentage_Carbon ~ Treatment, data = af)
par(mfrow=c(2,2));plot(mod1);par(mfrow=c(1,1))
car::qqPlot(residuals(mod1, type = 'pearson'))
anova(mod1)
summary(mod1)

af$Treatment <- relevel(af$Treatment, ref = "1-Year\nBare")
mod1 <- update(mod1)
summary(mod1)

af$Treatment <- relevel(af$Treatment, ref = "1-Year\nGrassland")
mod1 <- update(mod1)
summary(mod1)

# 10-Year Bare - 1-Year Bare              0.09331 .            
# 10-Year Bare - 1-Year Grassland         0.00696 ** 
# 10-Year Bare - 10-Year Grassland        0.00439 **
# 1-Year Bare - 1-Year Grassland          0.1297  
# 1-Year Bare - 10-Year Grassland         0.0667 .
# 1-Year Grassland - 10-Year Grassland    0.60851

# 10-Year Bare        a
# 1-Year Bare         a b
# 1-Year Grassland    - b
# 10-Year Grassland   - b

paste(af[, mean(Percentage_Carbon), by = Treatment]$Treatment
      , round(af[, mean(Percentage_Carbon), by = Treatment]$V1,2)
      , "+-"
      , round(af[, sd(Percentage_Carbon), by = Treatment]$V1,2))


