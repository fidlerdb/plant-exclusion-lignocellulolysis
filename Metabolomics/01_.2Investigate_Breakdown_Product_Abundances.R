rm(list = ls())

# This script investigates the relative abundances of 
# breakdown products of cellulose hemicellulose and lignin

# Load packages
library(data.table)
library(ggplot2)
library(cowplot)
library(multcomp)
source("Functions/get_Model_Letters.R")
mean_ci <- function(x){
  mean_se(x, mult = 1.96)
}
source("Functions/capStr.R")
source("Functions/ggcorrplot.R")
library(dplyr)
library(magrittr)
library(lme4)

#### Data Import and Cleaning #### 

# Load in the data
mf <- data.frame(fread("Metabolomics/Data/Metabolome_Blackout_Raw.csv"))
mf <- mf[, c(1,7:20)] # Keep only the most useful columns

# Put it in an easy format to work with
mfl <- melt(mf, id = "BinBase.name")

# Rename columns sensibly
names(mfl)[2:3] <- c("Sample", "Abundance")

# Create a treatment variable
mfl$Treatment <- NA
mfl$Treatment[grep("b1|b2|b3", mfl$Sample)] <- "b-Old"
mfl$Treatment[grep("c1|c2|c3", mfl$Sample)] <- "c-Old"
mfl$Treatment[grep("b4|b5|b6|b7", mfl$Sample)] <- "b-New"
mfl$Treatment[grep("c4|c5|c6|c7", mfl$Sample)] <- "c-New"
mfl$Treatment <- factor(mfl$Treatment)
head(mfl)

# AllMetabolites <- unique(mfl$BinBase.name)
# AllMetabolites[grep("acet", AllMetabolites)] # Gone through all
                                             # breakdown products 
                                             # shown in paper

# Create a function to normalize abundance data so that 
# "lignin breakdown products" can all be weighted equally
normalize <- function(x){
  nd <- (x-min(x))/(max(x)-min(x))
  return(nd)
}

#### Set up data for plotting and modelling ####

# Colour palette
BlackoutPalette <- c('#c2a5cf'
                     ,'#7b3294'
                     ,'#a6dba0'
                     ,'#008837')

#### Cellulose breakdown ####

plot(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == 'glucose',])
points(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == 'glucose',]
       , col = 'red', pch = 16, cex = 1.5)

pf_cellulose <- mfl[mfl$BinBase.name == 'glucose',]
pf_cellulose$NormAbund <- normalize(pf_cellulose$Abundance)

#### Hemicellulose breakdown ####

plot(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == 'xylose',])
points(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == 'xylose',]
       , col = 'red', pch = 16, cex = 1.5)
plot(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == 'fucose',])
points(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == 'fucose',]
       , col = 'red', pch = 16, cex = 1.5)
plot(log2(Abundance+1) ~ Treatment, data = mfl[mfl$BinBase.name == '3,6-anhydro-D-galactose',])
points(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == '3,6-anhydro-D-galactose',]
       , col = 'red', pch = 16, cex = 1.5)
# plot(log2(Abundance+1) ~ Treatment, data = mfl[mfl$BinBase.name == 'galactinol',])
# points(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == 'galactinol',]
#        , col = 'red', pch = 16, cex = 1.5)

# Make a ggplot2-able dataframe
pf_hemicellulose <- mfl[mfl$BinBase.name %in% c('xylose', 'fucose', '3,6-anhydro-D-galactose'),]
pf_hemicellulose$Abundance[pf_hemicellulose$BinBase.name == 'xylose']

pf_hemicellulose$NormAbund <- NA
pf_hemicellulose$NormAbund[pf_hemicellulose$BinBase.name == "xylose"] <- 
  normalize(pf_hemicellulose$Abundance[pf_hemicellulose$BinBase.name == "xylose"])
pf_hemicellulose$NormAbund[pf_hemicellulose$BinBase.name == "fucose"] <- 
  normalize(pf_hemicellulose$Abundance[pf_hemicellulose$BinBase.name == "fucose"])
pf_hemicellulose$NormAbund[pf_hemicellulose$BinBase.name == '3,6-anhydro-D-galactose'] <- 
  normalize(pf_hemicellulose$Abundance[pf_hemicellulose$BinBase.name == '3,6-anhydro-D-galactose'])
plot(NormAbund ~ Treatment, data = pf_hemicellulose)


#### Lignin

# https://biotechnologyforbiofuels.biomedcentral.com/articles/10.1186/s13068-017-0735-y

plot(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == 'vanillic acid',])
points(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == 'vanillic acid',]
       , col = 'red', pch = 16, cex = 1.5)
plot(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == '4-hydroxybenzoic acid',])
points(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == '4-hydroxybenzoic acid',]
       , col = 'red', pch = 16, cex = 1.5)
plot(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == 'benzoic acid',])
points(log2(Abundance) ~ Treatment, data = mfl[mfl$BinBase.name == 'benzoic acid',]
       , col = 'red', pch = 16, cex = 1.5)

# Make a ggplot2-able dataframe
pf_lignin <- mfl[mfl$BinBase.name %in% c('vanillic acid', '4-hydroxybenzoic acid', 'benzoic acid'),]

pf_lignin$NormAbund <- NA
pf_lignin$NormAbund[pf_lignin$BinBase.name == 'vanillic acid'] <- 
  normalize(pf_lignin$Abundance[pf_lignin$BinBase.name == 'vanillic acid'])
pf_lignin$NormAbund[pf_lignin$BinBase.name == '4-hydroxybenzoic acid'] <- 
  normalize(pf_lignin$Abundance[pf_lignin$BinBase.name == '4-hydroxybenzoic acid'])
pf_lignin$NormAbund[pf_lignin$BinBase.name == 'benzoic acid'] <- 
  normalize(pf_lignin$Abundance[pf_lignin$BinBase.name == 'benzoic acid'])
plot(NormAbund ~ Treatment, data = pf_lignin)

pf_lignin$Treatment
levels(pf_lignin$Treatment)

pf_lignin$Treatment



#### Model the different abundances ####

#### Cellulose
head(pf_cellulose) # Only glucose
cell_mod <- glm(NormAbund ~ Treatment, data = pf_cellulose
                , family = quasibinomial(link = 'logit'))
deviance(cell_mod)/df.residual(cell_mod)
par(mfrow=c(2,2));plot(cell_mod);par(mfrow=c(1,1))
drop1(cell_mod, test = 'Chisq')
summary(cell_mod)
#pf_cellulose$Treatment <- relevel(pf_cellulose$Treatment, ref = 'b-Old')
#cell_mod <- update(cell_mod)
#summary(cell_mod)

cell_letters <- get_Model_Letters(data = pf_cellulose, model = cell_mod)
cell_letters <- data.frame(Treatment = names(cell_letters)
                          , Letters = cell_letters
                          , NormAbund = 0.9)
cell_letters$Treatment <- as.character(cell_letters$Treatment)
cell_letters$Treatment[cell_letters$Treatment == "b-Old"] <- "10-Year\nBare"
cell_letters$Treatment[cell_letters$Treatment == "b-New"] <- "1-Year\nBare"
cell_letters$Treatment[cell_letters$Treatment == "c-Old"] <- "10-Year\nGrassland"
cell_letters$Treatment[cell_letters$Treatment == "c-New"] <- "1-Year\nGrassland"
cell_letters$Treatment <- factor(cell_letters$Treatment, levels = c("1-Year\nBare","10-Year\nBare"
                                                                  ,"1-Year\nGrassland","10-Year\nGrassland"))
#### Hemicellulose
head(pf_hemicellulose)

hem_mod_int <- glm(NormAbund ~ Treatment*BinBase.name, data = pf_hemicellulose
    , family = quasibinomial(link = 'logit'))
deviance(hem_mod_int)/df.residual(hem_mod_int)
par(mfrow=c(2,2));plot(hem_mod_int);par(mfrow=c(1,1))
drop1(hem_mod_int, test = 'Chisq')
hem_mod_sep <- update(hem_mod_int
                  , formula = "NormAbund ~ Treatment + BinBase.name"
               )
drop1(hem_mod_sep, test = 'Chisq')
hem_mod <- update(hem_mod_sep
                  , formula = "NormAbund ~ Treatment")
drop1(hem_mod, test = 'Chisq')
summary(hem_mod)

#####
library(lme4)
hem_mod_int <- glmer(pf_hemicellulose$NormAbund
                           ~ Treatment#+BinBase.name 
                     + (1|Sample)
                     , data = pf_hemicellulose
                   , family = binomial()
                   , control = glmerControl(optCtrl=list(maxfun=2.5e4))
)

hem_mod_int <- lmer(Abundance
                     ~ Treatment + BinBase.name 
                    + Treatment : BinBase.name 
                     + (1|Sample)
                     , data = pf_hemicellulose
                     #, family = binomial()
                     #, control = glmerControl(optCtrl=list(maxfun=2.5e4))
)
plot(hem_mod_int)

library(sjPlot)
summary(hem_mod_int)
drop1(hem_mod_int, test = 'Chisq')

plot_model(hem_mod_int)

# Singular fit. Model too complex, try refitting with the random effect as a fixed effect 

hem_mod_int <- glm(pf_hemicellulose$NormAbund
                   ~ Treatment + Sample#+BinBase.name 
                   #+ (1|Sample)
                   , data = pf_hemicellulose
                   , family = quasibinomial()
                   #, control = glmerControl(optCtrl=list(maxfun=2.5e4)
)

summary(hem_mod_int)
drop1(hem_mod_int, test = 'Chisq') # Model still too complex. But no effect of sample.

hem_mod_int <- glm(pf_hemicellulose$NormAbund
                     ~ Treatment #+ Sample#+BinBase.name 
                     #+ (1|Sample)
                     , data = pf_hemicellulose
                     , family = quasibinomial()
                     #, control = glmerControl(optCtrl=list(maxfun=2.5e4)
                   )

summary(hem_mod_int)
drop1(hem_mod_int, test = 'Chisq')

####

### Does a weighted model say anything different? ###
hem_weights <- melt(tapply(pf_hemicellulose$Abundance, pf_hemicellulose$BinBase.name, mean), value.name = "Weight")
names(hem_weights)[1] <- "BinBase.name"
pf_hemicellulose_wt <- merge(pf_hemicellulose, hem_weights, by = "BinBase.name")

hem_mod_wt <- glm(NormAbund ~ Treatment, data = pf_hemicellulose_wt
                  , family = quasibinomial(link = 'logit')
                  , weights = Weight)
drop1(hem_mod_wt, test = "F")
summary(hem_mod_wt)
par(mfrow = c(2,2));plot(hem_mod_wt);par(mfrow = c(1,1))

# Same outcome

# Get Tukey-esque comparisons from the model
hem_letters <- get_Model_Letters(data = pf_hemicellulose, model = hem_mod)
hem_letters <- data.frame(Treatment = names(hem_letters)
                           , Letters = hem_letters
                           , NormAbund = 0.9)
hem_letters$Treatment <- as.character(hem_letters$Treatment)
hem_letters$Treatment[hem_letters$Treatment == "b-Old"] <- "10-Year\nBare"
hem_letters$Treatment[hem_letters$Treatment == "b-New"] <- "1-Year\nBare"
hem_letters$Treatment[hem_letters$Treatment == "c-Old"] <- "10-Year\nGrassland"
hem_letters$Treatment[hem_letters$Treatment == "c-New"] <- "1-Year\nGrassland"
hem_letters$Treatment <- factor(hem_letters$Treatment, levels = c("1-Year\nBare","10-Year\nBare"
                                                                    ,"1-Year\nGrassland","10-Year\nGrassland"))

pf_hemicellulose$Treatment <- relevel(pf_hemicellulose$Treatment
                                      , ref = "1-Year\nBare")
hem_mod <- update(hem_mod)
summary(hem_mod)

#### Lignin
head(pf_lignin)
#pf_lignin$Treatment <- relevel(pf_lignin$Treatment, ref = 'b-New')

lig_mod_int <- glm(NormAbund ~ Treatment*BinBase.name, data = pf_lignin
                   , family = quasibinomial(link = 'logit'))
deviance(lig_mod_int)/df.residual(lig_mod_int)
car::qqPlot(residuals(lig_mod_int, type = 'deviance'))


lig_mod_int <- glmer(Abundance
                    ~ Treatment + BinBase.name 
                    + Treatment : BinBase.name 
                    + (1|Sample)
                    , data = pf_lignin
                    , family = Gamma(link = 'log')
                    #, control = glmerControl(optCtrl=list(maxfun=2.5e4))
)
plot(lig_mod_int)
summary(lig_mod_int)
drop1(lig_mod_int, test = 'Chisq')

plot_model(lig_mod_int)

par(mfrow=c(2,2));plot(lig_mod_int);par(mfrow=c(1,1))
drop1(lig_mod_int, test = 'Chisq')
lig_mod_sep <- update(lig_mod_int
                      , formula = "NormAbund ~ Treatment + BinBase.name"
)
drop1(lig_mod_sep, test = 'Chisq')
lig_mod <- update(lig_mod_int
                  , formula = "NormAbund ~ Treatment")
deviance(lig_mod)/df.residual(lig_mod)
car::qqPlot(residuals(lig_mod, type = 'deviance'))
par(mfrow=c(2,2));plot(lig_mod);par(mfrow=c(1,1))

drop1(lig_mod, test = 'Chisq')
summary(lig_mod)
# pf_lignin$Treatment <- relevel(pf_lignin$Treatment, ref = 'b-Old')
# lig_mod <- update(lig_mod)
# summary(lig_mod)
# pf_lignin$Treatment <- relevel(pf_lignin$Treatment, ref = 'c-Old')
# lig_mod <- update(lig_mod)
# summary(lig_mod)

lig_letters <- get_Model_Letters(data = pf_lignin, model = lig_mod)
lig_letters <- data.frame(Treatment = names(lig_letters)
                           , Letters = lig_letters
                           , NormAbund = 0.9)
lig_letters$Treatment <- as.character(lig_letters$Treatment)
lig_letters$Treatment[lig_letters$Treatment == "b-Old"] <- "10-Year\nBare"
lig_letters$Treatment[lig_letters$Treatment == "b-New"] <- "1-Year\nBare"
lig_letters$Treatment[lig_letters$Treatment == "c-Old"] <- "10-Year\nGrassland"
lig_letters$Treatment[lig_letters$Treatment == "c-New"] <- "1-Year\nGrassland"
lig_letters$Treatment <- factor(lig_letters$Treatment, levels = c("1-Year\nBare","10-Year\nBare"
                                                                  ,"1-Year\nGrassland","10-Year\nGrassland"))
### Does a weighted model say anything different? ###
lig_weights <- melt(tapply(pf_lignin$Abundance, pf_lignin$BinBase.name, mean), value.name = "Weight")
names(lig_weights)[1] <- "BinBase.name"
pf_lignin_wt <- merge(pf_lignin, lig_weights, by = "BinBase.name")

lig_mod_wt <- glm(NormAbund ~ Treatment, data = pf_lignin_wt
                  , family = quasibinomial(link = 'logit')
                  , weights = Weight)

par(mfrow = c(2,2));plot(lig_mod_wt);par(mfrow = c(1,1))

drop1(lig_mod_wt, test = "F")
summary(lig_mod_wt)

pf_lignin_wt$Treatment <- relevel(pf_lignin_wt$Treatment, ref = 'c-New')
lig_mod_wt <- update(lig_mod_wt)
summary(lig_mod_wt)
pf_lignin_wt$Treatment <- relevel(pf_lignin_wt$Treatment, ref = 'c-Old')
lig_mod_wt <- update(lig_mod_wt)
summary(lig_mod_wt)

# -- no, same outcome

#### Ensure that the plotting data frames are formatted correctly for plotting ####

#levels(pf_cellulose$Treatment) <- levels(mfl$Treatment)
#levels(pf_hemicellulose$Treatment) <- levels(mfl$Treatment)
#levels(pf_lignin$Treatment) <- levels(mfl$Treatment)


#### R hates me for some reason--getting all of the hemicellulose levels muddled, ####
##  so redefine all of the pf_*s  now and check they are correct...
par(mfrow = c(3,1))

pf_cellulose <- mfl[mfl$BinBase.name == 'glucose',]
pf_cellulose$NormAbund <- normalize(pf_cellulose$Abundance)
plot(NormAbund ~ Treatment, data = pf_cellulose)

pf_cellulose$Treatment <- as.character(pf_cellulose$Treatment)
pf_cellulose$Treatment[pf_cellulose$Treatment == "b-Old"] <- "10-Year\nBare"
pf_cellulose$Treatment[pf_cellulose$Treatment == "b-New"] <- "1-Year\nBare"
pf_cellulose$Treatment[pf_cellulose$Treatment == "c-Old"] <- "10-Year\nGrassland"
pf_cellulose$Treatment[pf_cellulose$Treatment == "c-New"] <- "1-Year\nGrassland"
pf_cellulose$Treatment <- factor(pf_cellulose$Treatment, levels = c("1-Year\nBare","10-Year\nBare"
                                                ,"1-Year\nGrassland","10-Year\nGrassland"))


pf_hemicellulose <- mfl[mfl$BinBase.name %in% c('xylose', 'fucose', '3,6-anhydro-D-galactose'),]
pf_hemicellulose$Abundance[pf_hemicellulose$BinBase.name == 'xylose']

pf_hemicellulose$NormAbund <- NA
pf_hemicellulose$NormAbund[pf_hemicellulose$BinBase.name == "xylose"] <- 
  normalize(pf_hemicellulose$Abundance[pf_hemicellulose$BinBase.name == "xylose"])
pf_hemicellulose$NormAbund[pf_hemicellulose$BinBase.name == "fucose"] <- 
  normalize(pf_hemicellulose$Abundance[pf_hemicellulose$BinBase.name == "fucose"])
pf_hemicellulose$NormAbund[pf_hemicellulose$BinBase.name == '3,6-anhydro-D-galactose'] <- 
  normalize(pf_hemicellulose$Abundance[pf_hemicellulose$BinBase.name == '3,6-anhydro-D-galactose'])
plot(NormAbund ~ Treatment, data = pf_hemicellulose)

pf_hemicellulose$Treatment <- as.character(pf_hemicellulose$Treatment)
pf_hemicellulose$Treatment[pf_hemicellulose$Treatment == "b-Old"] <- "10-Year\nBare"
pf_hemicellulose$Treatment[pf_hemicellulose$Treatment == "b-New"] <- "1-Year\nBare"
pf_hemicellulose$Treatment[pf_hemicellulose$Treatment == "c-Old"] <- "10-Year\nGrassland"
pf_hemicellulose$Treatment[pf_hemicellulose$Treatment == "c-New"] <- "1-Year\nGrassland"
pf_hemicellulose$Treatment <- factor(pf_hemicellulose$Treatment, levels = c("1-Year\nBare","10-Year\nBare"
                                                                    ,"1-Year\nGrassland","10-Year\nGrassland"))
# To make the legend pretty
pf_hemicellulose$BinBase.name <- factor(pf_hemicellulose$BinBase.name
                                        , levels = c("fucose", "xylose", "3,6-anhydro-D-galactose")) 

pf_lignin <- mfl[mfl$BinBase.name %in% c('vanillic acid', '4-hydroxybenzoic acid', 'benzoic acid'),]

pf_lignin$NormAbund <- NA
pf_lignin$NormAbund[pf_lignin$BinBase.name == 'vanillic acid'] <- 
  normalize(pf_lignin$Abundance[pf_lignin$BinBase.name == 'vanillic acid'])
pf_lignin$NormAbund[pf_lignin$BinBase.name == '4-hydroxybenzoic acid'] <- 
  normalize(pf_lignin$Abundance[pf_lignin$BinBase.name == '4-hydroxybenzoic acid'])
pf_lignin$NormAbund[pf_lignin$BinBase.name == 'benzoic acid'] <- 
  normalize(pf_lignin$Abundance[pf_lignin$BinBase.name == 'benzoic acid'])
plot(NormAbund ~ Treatment, data = pf_lignin)

pf_lignin$Treatment <- as.character(pf_lignin$Treatment)
pf_lignin$Treatment[pf_lignin$Treatment == "b-Old"] <- "10-Year\nBare"
pf_lignin$Treatment[pf_lignin$Treatment == "b-New"] <- "1-Year\nBare"
pf_lignin$Treatment[pf_lignin$Treatment == "c-Old"] <- "10-Year\nGrassland"
pf_lignin$Treatment[pf_lignin$Treatment == "c-New"] <- "1-Year\nGrassland"
pf_lignin$Treatment <- factor(pf_lignin$Treatment, levels = c("1-Year\nBare","10-Year\nBare"
                                                                    ,"1-Year\nGrassland","10-Year\nGrassland"))

#### Plots ####

pc <- ggplot(pf_cellulose
             , aes(x = Treatment, y = NormAbund)
             ) + 
  geom_point(aes(shape = BinBase.name)
             , alpha = 0.5) + 
  stat_summary(fun.data = 'mean_ci'
               , aes(colour = Treatment)
               , size = 1.1)+
  guides(shape = guide_legend(title = ""
                              , order = 2)
         , colour = 'none') + 
  scale_color_manual(values = BlackoutPalette
                     , breaks = levels(pf_cellulose$Treatment)) + 
  ylab("Abundance") +
  ggtitle("Cellulose") +
  geom_text(data = cell_letters, aes(label = Letters
                                     , x = Treatment
                                     , y = NormAbund
                                     , fontface = 'bold'), colour = 'black') + 
  theme(legend.position = "bottom")

ph <- ggplot(pf_hemicellulose
             , aes(x = Treatment, y = NormAbund)
             ) + 
  geom_point(aes(shape = BinBase.name)
             , position = position_jitter(0.05)
             , alpha = 0.5) + 
  stat_summary(fun.data = 'mean_ci'
               , aes(colour = Treatment)
               , size = 1.1) +
  guides(shape = guide_legend(title = "", nrow = 2, byrow = TRUE)
         , colour = 'none') + 
  scale_color_manual(values = BlackoutPalette
                     , breaks = levels(pf_hemicellulose$Treatment)) + 
  ylab("Abundance") +
  ggtitle("Hemicellulose") +
  geom_text(data = hem_letters, aes(label = Letters
                                    , x = Treatment
                                    , y = NormAbund
                                    , fontface = 'bold'), colour = 'black') + 
  theme(legend.position = "bottom")

pl <- ggplot(pf_lignin
             , aes(x = Treatment, y = NormAbund)
             ) + 
  geom_point(aes(shape = BinBase.name)
             , position = position_jitter(0.05)
             , alpha = 0.5
             ) + 
  stat_summary(fun.data = 'mean_ci'
               , aes(colour = Treatment
               )
               , size = 1.1
               ) +
  guides(shape = guide_legend(title = "", nrow = 2, byrow = TRUE)
         , colour = 'none') + 
  scale_color_manual(values = BlackoutPalette
                     , breaks = levels(pf_lignin$Treatment)) +
  ylab("Abundance") +
  ggtitle("Lignin") +
  geom_text(data = lig_letters, aes(label = Letters
                                    , x = Treatment
                                    , y = NormAbund
                                    , fontface = 'bold'), colour = 'black') + 
  theme(legend.position = "bottom") 



#### put all plots together ####
plot_grid(pc, ph, pl
          , align = "hv"
          # labels = c("Cellulose", "Hemicellulose", "Lignin")
          , ncol = 1
)

# Check the letters make sense and the orders are real because sometimes they're nnot right.......

# Now save because you have the correct results (once in a blue moon)
# ggsave("Figures/Metabolomics/2020-06-21_Metabolite_Abundance.pdf"
#        , width = 7.06
#        , height = 9
#        , units = "in")

# Save plot data to a list for combination with fibre analysis work

saveRDS(list(pf_cellulose, pf_hemicellulose, pf_lignin
     , cell_letters, hem_letters, lig_letters), file = "Fibre_Analysis/Output_Data/Metabolite_Data_For_Figure.rds")


# Not useful here but useful for future color replacement things :)
# library(colorspace)
# LabelPal <- c("#066C0D", "#066C0D", "#066C0D"
#               ,"#6C1706", "#6C1706"
#               , "#0D066C" 
#               , "#6C1706")
# LabelPal <- replace(LabelPal
#                     , list = c("#066C0D", "#6C1706", "#0D066C")
#                     , values = lighten(c("#066C0D", "#6C1706", "#0D066C") , c(1000,0,0))
# )

#### Now add the correlation plot in with these ####

# Check how the metabolites co-occur

# Load in the data
mf <- data.frame(fread("Metabolomics/Data/Metabolome_Blackout_Raw.csv"))
mf <- mf[, c(1,7:20)] # Keep only the most useful columns and rows
mf <- mf[mf$BinBase.name %in% c("glucose"
                                , "fucose", "xylose", "3,6-anhydro-D-galactose"
                                , "vanillic acid", "4-hydroxybenzoic acid", "benzoic acid"),]
# Clean metabolite names
mf[mf$BinBase.name %in% c("3,6-anhydro-D-galactose"
                          , "4-hydroxybenzoic acid"),]$BinBase.name <- c("HBA", "galactose")
mf$BinBase.name <- capStr(mf$BinBase.name)

# Make this into a named matrix with the focus on metabolites
mfs <- as.matrix(mf[,2:length(mf)]) %>% `dimnames<-` (list(mf$BinBase.name, colnames(mf[2:length(mf)])))
mfs <- t(mfs)
mfs

# Correlation matrix with p values
res2 <- Hmisc::rcorr(mfs)



# Create the correlation plot

p_corr2 <- ggcorrplot(res2$r
                      ,type = "lower"
                      #, diag = TRUE
                      , method = "circle"  #"ellipse"
                      , colourscale.direction = "horizontal"
                      , p.mat = res2$P
                      , show.diag = FALSE
                      #, y.label.angle = 10
                      , y.label.hjust = 0.2 #"left"
                      )

p_corr2 <- p_corr2 + theme(legend.position = "bottom") 

# Create the final plot
#p <- plot_grid(ph, pl, NULL, pc, p_corr2, NULL
          #, align = "hv"
          # labels = c("Cellulose", "Hemicellulose", "Lignin")
#          , ncol = 3
#          , labels = c("a)", "b)", "","c)", "d)", "")
 #         , rel_widths = c(1,1,0.2)
#)

#p

save_plot(filename = "Figures/Metabolomics/2020-10-13_Metabolite_Correlation.pdf"
          , plot = p_corr2, base_height = 5, base_width = 5
          )







