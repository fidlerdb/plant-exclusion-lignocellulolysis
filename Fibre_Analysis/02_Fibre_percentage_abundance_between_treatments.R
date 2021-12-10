## This script assesses the abundance of different lignocellulosic ##
## elements between treatments ##

library(data.table)
library(ggplot2)

#### Data import ####
rm(list = ls())

sf <- fread("Fibre_Analysis/Output_Data/Fibre_Analysis_Data_Processed_2021-04-29.csv")

# Sort out the treatment variable
sf$Treatment <- factor(sf$Treatment)
levels(sf$Treatment) <- c("1-Year Bare", "10-Year Bare", "1-Year Grassland", "10-Year Grassland")
sf$Treatment <- relevel(sf$Treatment, ref = "10-Year Bare")

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

sf$Treatment

library(pwr)
pwr.anova.test(k = 4, f = 0.5, sig.level = 0.05, n = 4) # 26% chance of finding a large effect overall
pwr.anova.test(k = 2, f = 0.5, sig.level = 0.05, n = 3) # 15-23% chance of finding a large effect

library(MultNonParam)
kwpower(rep(3:4,2), rep(0.5,4), "normal")

#### Plot hemicellulose content ####

ph <- ggplot(sf, aes(x = Treatment, y = Proportion_Hemicellulose * 100
               , fill = Treatment)) + 
  stat_summary(fun = "median", size = 1.1, geom = "bar") + 
  geom_point(position = position_jitter(0.15), colour = "black") +
  ylab("Hemicellulose (% of biomass )") + 
  scale_fill_manual(values = BlackoutPalette
                     , breaks = levels(sf$Treatment)) +
  theme_classic() +
  theme(legend.position = "none")

ph

## Were the differences between groups?

kruskal.test(sf$Proportion_Hemicellulose, sf$Treatment)
# No differences 

#dunnTest(sf$Proportion_Hemicellulose, sf$Treatment, method = "bh")
dunn.test::dunn.test(sf$Proportion_Hemicellulose, sf$Treatment)
# No differences -- using unadjusted p-values due to tiny (0.05) power

# 10-Year Bare - 1-Year Bare              0.6726038           
# 10-Year Bare - 1-Year Grassland         0.1514264
# 10-Year Bare - 10-Year Grassland        0.4349672
# 1-Year Bare - 1-Year Grassland          0.6726038
# 1-Year Bare - 10-Year Grassland         0.8347166
# 1-Year Grassland - 10-Year Grassland    0.5485818

# 10-Year Bare        a
# 1-Year Bare         a
# 1-Year Grassland    a
# 10-Year Grassland   a

# Multiple comparison letters
hem_letters <- data.table(Treatment = c("10-Year Bare" , "1-Year Bare" , "1-Year Grassland", "10-Year Grassland")
                          , Letters = c("a", "a", "a", "a")
                          , Proportion_Hemicellulose = 50#sf[, max(Proportion_Cellulose) + 
                          #  max(Proportion_Cellulose)*0.1, by = Treatment][,2]*100
)

ph
ph2 <- ph + geom_text(data = hem_letters, aes(x = Treatment, y = Proportion_Hemicellulose, label = Letters), colour = "black") +
  ylim(0,100)

ph2
#### Cellulose ####

pc <- ggplot(sf, aes(x = Treatment, y = Proportion_Cellulose * 100
               , fill = Treatment)) + 
  stat_summary(fun = "median", size = 1.1, geom = 'bar') +   
  geom_point(position = position_jitter(0.15), colour = "black") +
  ylab("Cellulose (% of biomass)") + 
  scale_fill_manual(values = BlackoutPalette
                     , breaks = levels(sf$Treatment)) +
  theme_classic() +
  theme(legend.position = "none")

pc

## Were there differences between groups?

kruskal.test(sf$Proportion_Cellulose, sf$Treatment) 
# Possibility of differences between groups

#dunnTest(sf$Proportion_Cellulose, sf$Treatment, method = "bh") 
dunn.test::dunn.test(sf$Proportion_Cellulose, sf$Treatment, list = TRUE)
# 10-Year Bare - 1-Year Bare              0.03463124 *            
# 10-Year Bare - 1-Year Grassland         0.04191148 * 
# 10-Year Bare - 10-Year Grassland        0.01917248 *
# 1-Year Bare - 1-Year Grassland          0.93264664
# 1-Year Bare - 10-Year Grassland         0.69562694
# 1-Year Grassland - 10-Year Grassland    0.63872909

sf[, median(Proportion_Cellulose), by = Treatment]
sf[, quantile(Proportion_Cellulose, 0.25), by = Treatment]
sf[, quantile(Proportion_Cellulose, 0.75), by = Treatment]


# Multiple comparison letters
cel_letters <- data.table(Treatment = c("10-Year Bare" , "1-Year Bare" , "1-Year Grassland", "10-Year Grassland")
                          , Letters = c("a", "b", "b", "b")
                          , Proportion_Cellulose = 50#sf[, max(Proportion_Cellulose) + 
                                                      #  max(Proportion_Cellulose)*0.1, by = Treatment][,2]*100
                          )

pc
pc2 <- pc + geom_text(data = cel_letters, aes(x = Treatment, y = Proportion_Cellulose, label = Letters), colour = "black") +
  ylim(0,100)

pc2
#### Lignin ####

pl <- ggplot(sf, aes(x = Treatment, y = Proportion_Lignin * 100
               , fill = Treatment)) + 
  stat_summary(fun = "median", size = 1.1, geom = "bar") + 
  geom_point(position = position_jitter(0.15), colour = "black") +
  ylab("Lignin (% of biomass)") + 
  scale_fill_manual(values = BlackoutPalette
                     , breaks = levels(sf$Treatment)) +
  theme_classic() +
  theme(legend.position = "none")

pl

kruskal.test(sf$Proportion_Lignin, sf$Treatment)
# No overall treatment effect -- post hoc testing means nothing
dunnTest(sf$Proportion_Lignin, sf$Treatment, method = "bh")

# 10-Year Bare - 1-Year Bare              0.13710053            
# 10-Year Bare - 1-Year Grassland         0.63872909 
# 10-Year Bare - 10-Year Grassland        0.43496716
# 1-Year Bare - 1-Year Grassland          0.03461056 * 
# 1-Year Bare - 10-Year Grassland         0.51436849
# 1-Year Grassland - 10-Year Grassland    0.19219904

# a = no difference from 10 year bare
# b = different from 10 year bare, not different from 1 year bare

# 10-Year Bare        ab   
# 1-Year Bare         b     
# 1-Year Grassland    a  
# 10-Year Grassland   ab  


# Multiple comparison letters
lig_letters <- data.table(Treatment = c("10-Year Bare" , "1-Year Bare" , "1-Year Grassland", "10-Year Grassland")
                          , Letters = c("a", "a", "a", "a")
                          , Proportion_Lignin = 25#sf[, max(Proportion_Cellulose) + 
                          #  max(Proportion_Cellulose)*0.1, by = Treatment][,2]*100
)

pl
pl2 <- pl + geom_text(data = lig_letters, aes(x = Treatment, y = 12, label = Letters), colour = "black") +
  ylim(0,15)

pl2

#### Plot all polymer % abundances together
plot_grid(pc2, ph2, pl2, align = "v", ncol =  1)

# P values for the analysis obtained. Save the objects which make the plots to an rds to combine with the metabolite abundance plots

saveRDS(list(sf, cel_letters, hem_letters, lig_letters), file = "Fibre_Analysis/Output_Data/Fibre_Data_For_Figure_2021-05-05.rds")





