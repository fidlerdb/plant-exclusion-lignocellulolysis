rm(list = ls())

#### Data import  and checking ####
df <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZyFamilyReadsPerKbp.csv")

library(vegan)
library(scales)
#### Data Cleaning and checking ####

# Normalise per Gbp
df[, 4:length(df)] <- df[, 4:length(df)] / (df$Reads/1e9) 

# CAZyme community matrix
cf <- df[,4:length(df)]

# Treatment variable
Treatment <- df$Treatment
levels(df$Treatment) <- c("1-Year Fallow"
                          , "10-Year Fallow"
                          , "1-Year Grassland"
                          , "10-Year Grassland")
cf[,1:5]
# Read abundance based diversity scores

#### Did treatment affect CAZyme diversity? ####

# Calculate diversity
cd <- cbind(df$Treatment, cf)
cd$Shannon <- diversity(cf)
cd$Richness <- specnumber(cf)
source("Functions/plotGlm.R")

# Shannon diversity (CAZyme homogeneity)
plot(exp(Shannon) ~ Treatment, data = cd) # Plot values

div_mod <- lm(exp(Shannon) ~ Treatment, data = cd) # Create model

plotGlm(div_mod) # Check assumptions
car::qqPlot(residuals(div_mod,type = "pearson"))

# Test hypotesis
drop1(div_mod, test = 'F') # No treatment effect on shannon's diversity
summary(div_mod)
kruskal.test(Shannon ~ Treatment, data = cd) # check the nonparametric version gives the same output

# CAZy family richness
plot(Richness ~ factor(Treatment), data = cd) # Plot values

rich_mod <- glm(Richness ~ Treatment, data = cd, family = poisson()) # Create model
rich_mod <- glm(Richness ~ Treatment, data = cd, family = quasipoisson()) # Create model

plotGlm(rich_mod) # Check assumptions

car::qqPlot(residuals(rich_mod,type = "deviance"))
deviance(rich_mod)/df.residual(rich_mod) # underdispersion -- conservative model, use quasipoisson

# Test hypotesis
drop1(rich_mod, test = 'Chi') # Strong treatment effect on CAZy family richness
summary(rich_mod)
kruskal.test(Richness ~ Treatment, data = cd) # Rank transformation suggests there is

library(ggplot2)
ggplot(cd, aes(x = Treatment, y = Richness)) + geom_point()

cd2 <- cd
library(data.table)
setDT(cd2)

cd2[, list(Richness.mean = mean(Richness)
           , Richness.sd = sd(Richness)), by = factor(Treatment)]

#### CAZyme community composition ####
# Was there a significant difference between the CAZyme communities in the different treatments?
adonis(cf ~ Treatment)

# Pairwise differencves between treatments
source("Functions/pairwise_adonis.R")
pairwise.adonis(cf, Treatment, p.adjust.m = "fdr")
# Check these weird results with an NMDS -- doesn't work, explains weird results
nmds2 <- metaMDS(cf[Treatment %in% c("b-Old", "c-Old"),]
                 , wascores = TRUE)
library(emmeans)

dat2 <- cbind(cf, Treatment)
dat2 <- dat2#[Treatment %in% c("b-Old", "c-Old"),]
mlm1 <- lm(log10(as.matrix(dat2[names(dat2) %in% names(cf)])+1)
           ~ Treatment, data = dat2)

library(MVN)
mvn(residuals(mlm1), mvnTest = "mardia", multivariatePlot = "qq")


library(micompr)
library(biotools)
cmp <- cmpoutput("CAZymes", 3, cf, Treatment)
cmp
plot(cmp)

cmp <- cmpoutput("CAZymes", 3, cf, Treatment)

comparisons <- data.frame(comp1 = c("b-Old"
                                    , "b-Old"
                                    , "b-Old"
                                    , "b-New"
                                    , "b-New"
                                    , "c-Old")
                          , comp2 = c("b-New"
                                      ,"c-Old"
                                      ,"c-New"
                                      ,"c-New"
                                      ,"c-Old"
                                      ,"c-New")
                          )
library(data.table)

min(table(Treatment))

pairwise_multivariate_comparisons <- rbindlist(apply(comparisons, 1, FUN = function(x){
  Trts <- x[1:2]
  cf2 <- cf[Treatment %in% Trts,]
  Treat <- factor(Treatment[Treatment %in% Trts])
  NPCs <- min(table(Treat))
  cmp <- cmpoutput(paste(Trts, sep = "_"), NPCs, cf2, Treat)
  cmp <- summary(cmp)
  results <- with(cmp, data.frame(Treatment_1 = Trts[1]
                                  , Treatment_2 = Trts[2]
                                  , Test = "MANOVA"
                                  , n_PCs = num.pcs
                                  , p = manova.pvals))
  
  return(results)
}))

pairwise_multivariate_comparisons
p.adjust(pairwise_multivariate_comparisons$p, method = "fdr")

# Really complicated and there's no consensus in the p-values leave it with this dataset.


car::Manova(mlm1)
mlm1# Get the least square means
mlm1.lsm <- emmeans(mlm1, "Treatment")
# Multiple comparisons with fdr p value adjustment
test(contrast(mlm1.lsm, "pairwise"), side = "=",  adjust = "bonferroni")
test(contrast(mlm1.lsm, "pairwise"), side = "=",  adjust = "holm")
test(contrast(mlm1.lsm, "pairwise"), side = "=",  adjust = "fdr")
test(contrast(mlm1.lsm, "pairwise"), side = "=",  adjust = "none")

# Make an NMDS (no multivariate normality assumptions)
nmds1 <- metaMDS(cf, wascores = TRUE) # Low stress so should be a good fit

# Create useful vectors for plotting
BlackoutPalette <- c('#c2a5cf'
                     ,'#7b3294'
                     ,'#a6dba0'
                     ,'#008837')
colours <- BlackoutPalette[Treatment]

# how correlated species are. scaled sum of the correlations to each axis
source("Functions/range01.R")
magnitude <- function(x){x <- sqrt(x^2); return(x)}
speciesAlpha <- range01(magnitude(nmds1$species[,1]) + magnitude(nmds1$species[,2]))

#### Plot of NMDS ####

# Create a plot and save it
pdf("Figures/CAZy/CAZyme_Community_Composition.pdf", width = 6, height = 4.84)

# Allow external legend
par(xpd = TRUE)

# Base plot
plot(nmds1$points, type = 'n'
     , xlim = c(-0.3, 0.3)
     , ylim = c(-0.3, 0.3)
     , bty = 'L')

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
legend(x = -0.325, y = 0.45, legend = levels(df$Treatment)
       , col = BlackoutPalette
       , pch = 16)

# Save the plot

dev.off()


#### Importance of variables along PC1 ####

CAZymeImportance <- na.omit(nmds1$species)
MI <- nrow(CAZymeImportance)

# NMDS Axis 1
sort(CAZymeImportance[,1])[1:10] # Most negative correlations
sort(CAZymeImportance[,1])[MI:(MI-10)] # Most positive correlations

# NMDS Axis 2
sort(CAZymeImportance[,2])[1:10] # Most negative correlations
sort(CAZymeImportance[,2])[MI:(MI-10)] # Most positive correlations

#pwcomp <- pairwise.adonis(cf, Treatment, p.adjust.m = "fdr", perm = 1e8)
#summary(pwcomp)

################
# Make an NMDS (no multivariate normality assumptions)
nmds2 <- metaMDS(cf[Treatment %in% c("b-Old", "c-Old"),], wascores = TRUE)

# Create useful vectors for plotting
BlackoutPalette <- c('#c2a5cf'
                     ,'#7b3294'
                     ,'#a6dba0'
                     ,'#008837')
colours <- BlackoutPalette[Treatment]

# how correlated species are. scaled sum of the correlations to each axis
source("Functions/range01.R")
magnitude <- function(x){x <- sqrt(x^2); return(x)}
speciesAlpha <- range01(magnitude(nmds1$species[,1]) + magnitude(nmds1$species[,2]))

#### Plot of NMDS ####

# Create a plot and save it
pdf("Figures/CAZy/CAZyme_Community_Composition.pdf", width = 6, height = 4.84)

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
legend(x = -0.325, y = 0.45, legend = levels(df$Treatment)
       , col = BlackoutPalette
       , pch = 16)

# Save the plot

dev.off()

##### Choosing the optimal k value ####

# nmds_k2 <- metaMDS(cf, wascores = TRUE, k = 2)
# nmds_k3 <- metaMDS(cf, wascores = TRUE, k = 3)
# nmds_k4 <- metaMDS(cf, wascores = TRUE, k = 4)
# nmds_k5 <- metaMDS(cf, wascores = TRUE, k = 5)
# nmds_k6 <- metaMDS(cf, wascores = TRUE, k = 6)
# 
# 
# stressframe <- data.frame(Model = c("nmds_k2", "nmds_k3", "nmds_k4", "nmds_k5", "nmds_k6")
#            , Stress = c(nmds_k2$stress, nmds_k3$stress
#                         , nmds_k4$stress, nmds_k5$stress, nmds_k6$stress)
#            )
# plot(stressframe)
# 
# 
# # Allow external legend
# par(xpd = TRUE)
# 
# # Base plot
# plot(nmds_k3$points, type = 'n'
#      , xlim = c(-0.3, 0.3)
#      , ylim = c(-0.3, 0.3)
#      ,bty = 'L')
# 
# # CAZyme labels
# orditorp(nmds_k3, display="species"
#          , col = alpha("black", speciesAlpha)
# )
# 
# # Samples CAZy composition values
# points(nmds_k3$points
#        , pch = 16
#        , col = colours)
# 
# # Add convex hulls
# for(i in levels(Treatment)){
#   ordihull(nmds_k3$points[grep(i, Treatment),]
#            , groups = Treatment[Treatment == i]
#            , draw = "polygon"
#            , col = colours[grep(i, Treatment)]
#            , border = FALSE
#            , alpha = 0.7
#   )
# }
# 
# # Information about the performance of the NMDS model
# text(x = 0.25, y = -0.28
#      , labels = paste0("Stress = ",round(nmds_k3$stress, 3))
# )
# 
# # Add legend
# legend(x = -0.325, y = 0.45, legend = levels(df$Treatment)
#        , col = BlackoutPalette
#        , pch = 16)
# 
# 
# sort(nmds_k3$species[,1])
