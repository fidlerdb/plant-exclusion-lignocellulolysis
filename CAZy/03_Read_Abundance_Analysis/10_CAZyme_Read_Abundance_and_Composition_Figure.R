rm(list = ls())

# Create a figure showing total CAZyme abundance and CAZy community composition
library(cowplot)

##### Total CAZyme Abundance ####

# Data import  and checking 
df <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZyFamilyReadsPerKbp.csv")
rft <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZymeResponse_p_FPR.csv")

library(roperators)
#library(vegan)
library(data.table)
library(ggplot2)
#library(Hmisc)
library(doBy)
library(multcomp)

#### Data Cleaning and checking ####
# Normalise per Gbp
df[, 4:length(df)] <- df[, 4:length(df)] / (df$Reads/1e9) 

# Create a total CAZyme column
tf <- df
tf$CAZymes <- rowSums(df[,4:length(df)])
tf$Treatment <- factor(tf$Treatment)
levels(tf$Treatment) <- c("1-Year\nBare"
                       , "10-Year\nBare"
                       , "1-Year\nGrassland"
                       , "10-Year\nGrassland")
tf$Treatment <- relevel(tf$Treatment, ref = "10-Year\nBare")

# Create a model to test if there are differences
mod1 <- lm((CAZymes) ~ Treatment, data = tf)

# Create a sensible 
mean_ci <- function(x){
  mean_se(x, mult = 1.96)
}

# Save the multiple comparison output for plotting
TukeyLabels <- cld(glht(mod1, linfct = mcp(Treatment = "Tukey")), level = 0.05)
TukeyLabels <- data.frame(Treatment = names(TukeyLabels$mcletters$Letters)
                          , Letters = TukeyLabels$mcletters$Letters
                          , CAZymes = rbindlist(tapply(tf$CAZymes, tf$Treatment, mean_ci))$ymax
)

BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

#### CAZyme abundance ####

p_abundance <- ggplot(tf, aes(x = Treatment, y = CAZymes
                    , colour = Treatment)) + 
  geom_point(colour = 'grey50', position = position_jitter(0.15)) +
  stat_summary(fun.data = mean_ci, size = 1.1) +
  theme_classic() + 
  scale_colour_manual(breaks = levels(tf$Treatment)
                      , values = BlackoutPalette) +
  scale_x_discrete(breaks=levels(tf$Treatment)
                   #, labels=c("1-Year\nFallow", "10-Year\nFallow"
                              #, "1-Year\nGrassland", "10-Year\nGrassland")
) +
  guides(colour = 'none') +
  ylab("CAZyme abundance") +
  geom_text(data = TukeyLabels, aes(label = Letters
                                    , x = Treatment
                                    , y = CAZymes + 30), colour = 'black') +
  theme(axis.text = element_text(size = rel(1)))

p_abundance  

ggsave("Final_figure_making/To_Put_in_Inkscape/Figure_2b_CAZyme_abundance.pdf"
       , plot = p_abundance
       , width = 5, height = 4)

#### Community Composition ####

# Data import  and checking 
df <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZyFamilyReadsPerKbp.csv")

library(vegan)
library(scales)
# Data Cleaning and checking

# Normalise per Gbp
df[, 4:length(df)] <- df[, 4:length(df)] / (df$Reads/1e9) 
# CAZyme community matrix
cf <- df[,4:length(df)]

# Treatment variable
Treatment <- factor(df$Treatment)
levels(Treatment) <- c("1-Year Bare"
                          , "10-Year Bare"
                          , "1-Year Grassland"
                          , "10-Year Grassland")
Treatment <- relevel(Treatment, ref = "10-Year Bare")

# Make an NMDS (no multivariate normality assumptions)
nmds1 <- metaMDS(cf, wascores = TRUE)

# Create useful vectors for plotting
BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')
colours <- BlackoutPalette[Treatment]

# how correlated species are. scaled sum of the correlations to each axis
source("Functions/range01.R")
magnitude <- function(x){x <- sqrt(x^2); return(x)}
speciesAlpha <- range01(magnitude(nmds1$species[,1]) + magnitude(nmds1$species[,2]))

#### Plot of NMDS ####
dev.off()
pdf(file = "Final_figure_making/To_Put_in_Inkscape/Figure_2c_CAZyme_Community.pdf"
    , width = 5, height = 4)

# Allow external legend
par(xpd = TRUE)

# Base plot
plot(nmds1$points, type = 'n'
     , xlim = c(-0.3, 0.3)
     , ylim = c(-0.3, 0.3)
     ,bty = 'L')

# CAZyme labels
orditorp(nmds1, display="species"
         , col = alpha("black", speciesAlpha)
)

# Samples CAZy composition values
points(nmds1$points
       , pch = 16
       , col = colours)

# Add convex hulls
for(i in levels(Treatment)){
  ordihull(nmds1$points[grep(i, Treatment),]
           , groups = Treatment[Treatment == i]
           , draw = "polygon"
           , col = colours[grep(i, Treatment)]
           , border = FALSE
           , alpha = 0.7
  )
}

# Information about the performance of the NMDS model
text(x = 0.25, y = -0.28
     , labels = paste0("Stress = ",round(nmds1$stress, 3))
)

# Add legend
legend(#"topleft"
       x = -0.325, y = 0.55
       , legend = levels(Treatment)[1:2]
       , col = BlackoutPalette
       , pch = 16
       , bty = 'n')
legend(#"topleft"
  x = 0, y = 0.55
  , legend = levels(Treatment)[3:4]
  , col = BlackoutPalette[3:4]
  , pch = 16
  , bty = 'n')

dev.off()
#Plot_NMDS <- recordPlot()

#### Put the two plots together into a single figure ####

# Turns out this is nigh on impossible to get the sizes right. Just combine them manually, then add a), b) labels
