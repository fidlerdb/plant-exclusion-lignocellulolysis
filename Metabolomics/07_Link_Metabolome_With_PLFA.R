# Discover relationships between PLFA data and metabolite data

#### Data Import and Cleaning #### 
library(data.table)
library(ggplot2)

rm(list = ls())

# Load in the data
mf <- data.frame(fread("Metabolomics/Data/Metabolome_Blackout_Raw.csv"))
pf <- fread("Metabolomics/Data/PLFA_For_Analysis.csv")

head(mf) # Metabolites
head(pf) # PLFA

mf <- mf[, c(1,7:20)] # Keep only the most useful columns

# Lignocellulose breakdown products
BreakdownProducts <- c('glucose'
                       , 'xylose', 'fucose', '3,6-anhydro-D-galactose'
                       , 'vanillic acid', '4-hydroxybenzoic acid', 'benzoic acid')

# Keep only those products
mf <- mf[mf$BinBase.name %in% BreakdownProducts,]

head(mf)

# Create the mergable-metabolite data frame
AbundData <- data.frame(t(mf[2:length(mf)]))
names(AbundData) <- mf[,1]
AbundData$Sample <- row.names(AbundData)

# Tidy the PLFA data
head(pf)
pf$Percent_All_Fungi <- pf$Percent_AM_Fungi + pf$Percent_Fungi
pf$Percent_All_Bacteria <- pf$Percent_Gram_Negative + pf$Percent_Gram_Positive +
  pf$Percent_Actinomycetes + pf$Percent_Anaerobe
pf$FB_Ratio <- pf$Percent_All_Fungi / pf$Percent_All_Bacteria

pm <- merge(pf, AbundData, by = 'Sample')

# Save this as a file
#write.csv(ofm, "Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.csv", row.names = FALSE)
#saveRDS(ofm, file = "Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.rds")


names(pm)
pm <- data.frame(pm)

#### What governs metabolite abundance in soils? ####

pmb <- pm

# Changing model family improves AIC, switch to a Gamma model
lapply(seq_along(pmb[,14:length(pmb)])
       , FUN = function(METABOLITE){
         mod1 <- glm(unlist(pmb[,14:length(pmb)][,METABOLITE]) ~ FB_Ratio + Total_PLFA
                     , data = pmb
                     , family = Gamma())
         mod2 <- glm(unlist(pmb[,14:length(pmb)][,METABOLITE]) ~ FB_Ratio + Total_PLFA
                     , data = pmb
                     , family = Gamma(link = 'log'))
         mod3 <- lm(unlist(pmb[,14:length(pmb)][,METABOLITE]) ~ FB_Ratio + Total_PLFA
                    , data = pmb)
         # mod4 <- glm(unlist(pmb[,14:length(pmb)][,METABOLITE]) ~ FB_Ratio + Total_PLFA
         #            , data = pmb
         #            , family = binomial())
         AIC(mod1,mod2, mod3)
       })

# Refit the models with a gamma family
MetabResults <- lapply(seq_along(pmb[,14:length(pmb)])
                       , FUN = function(METABOLITE){
                         # Fit a model for each metabolite
                         mod <- glm(unlist(pmb[,14:length(pmb)][,METABOLITE]) ~ 
                                     FB_Ratio + Total_PLFA
                                   , data = pmb
                                   , family = Gamma)
                         
                         Result <- summary(mod)
                         
                         # Extract metabolite name, p value, and r value
                         data.frame(Metabolite = names(pmb[,14:length(pmb)][METABOLITE])
                                    , Predictor = dimnames(drop1(mod))[[1]][2:3]
                                    , p = drop1(mod, test = "F")[5][2:3,]
                                    , r = Result$coefficients[2:3,1])
                       })
MetabResults <- rbindlist(MetabResults)
MetabResults 

MetabResults[MetabResults$Predictor == "Total_PLFA",] # None are significant. Refit without these
MetabResults[MetabResults$Predictor == "FB_Ratio",]
p.adjust(MetabResults[MetabResults$Predictor == "FB_Ratio",]$p, method = "fdr")
# just to check

#### Results with the gaussian model, N.B. no change in outcome
# MetabResults <- lapply(seq_along(pmb[,14:length(pmb)])
#                        , FUN = function(METABOLITE){
#                          # Fit a model for each metabolite
#                          mod <- glm(unlist(pmb[,14:length(pmb)][,METABOLITE]) ~ 
#                                      FB_Ratio + Total_PLFA
#                                    , data = pmb
#                                    , family = Gamma)
#                          
#                          Result <- summary(mod)
#                          
#                          # Extract metabolite name, p value, and r value
#                          data.frame(Metabolite = names(pmb[,14:length(pmb)][METABOLITE])
#                                     , Predictor = dimnames(anova(mod))[[1]][1:2]
#                                     , p = anova(mod)[5][1:2,]
#                                     , r = Result$coefficients[2:3,1])
#                        })
# MetabResults <- rbindlist(MetabResults)
# MetabResults 
####

MetabResults[MetabResults$Predictor == "Total_PLFA",] # None are significant. Refit without these
MetabResults[MetabResults$Predictor == "FB_Ratio",]

# Total microbial abundance was not a predictor of lignocellulose breakdown product abundance
# Fungal:Bacterial ratio was a predictor of lignocellulose breakdown product abundance 
#     in a number of cases

# Gamma models perform significantly better
lapply(seq_along(pmb[,14:length(pmb)])
       , FUN = function(METABOLITE){
         mod1 <- glm(unlist(pmb[,14:length(pmb)][,METABOLITE]) ~ FB_Ratio
                     , data = pmb
                     , family = Gamma())
         mod2 <- lm(unlist(pmb[,14:length(pmb)][,METABOLITE]) ~ FB_Ratio
                    , data = pmb)
         data.frame(AIC(mod1,mod2), DeltaAIC = AIC(mod1) - AIC(mod2))
       })

# Refitting without total PLFA as a predictor
library(pscl)

MetabResults <- lapply(seq_along(pmb[,14:length(pmb)])
                       , FUN = function(METABOLITE){
                         # Fit a model for each metabolite
                         mod <- glm(unlist(pmb[,14:length(pmb)][,METABOLITE]) ~ 
                                      FB_Ratio
                                    , data = pmb
                                    , family = Gamma)
                         
                         Result <- summary(mod)
                         
                         # Extract metabolite name, p value, and r value
                         data.frame(Metabolite = names(pmb[,14:length(pmb)][METABOLITE])
                                    , Predictor = dimnames(drop1(mod))[[1]][2]
                                    , p = drop1(mod, test = "F")[5][2,]
                                    , r = Result$coefficients[2,1]
                                    , R2 = pR2(mod)[4]
                                    , FStat = drop1(mod, test = "F")[4][2,]
                                    , df = paste(drop1(mod)[1][2,], Result$df.residual, sep = ","))
                       })
MetabResults <- rbindlist(MetabResults)

# Prettier results :)
MetabResults 

#### Plot relationships between these variables ####

#### Data tidying for plotting ####
pm <- data.frame(pm)
names(pm)
pma <- pm[,names(pm) %in% c("Sample"#, "Total_PLFA", "Percent_Actinomycetes", "Percent_AM_Fungi"
                            #, "Percent_All_Fungi", "Percent_All_Bacteria" 
                            , "FB_Ratio"
                            , "xylose", "vanillic.acid", "glucose"
                            , "fucose", "benzoic.acid", "X4.hydroxybenzoic.acid"
                            , "X3.6.anhydro.D.galactose"
)]

# Make the data long for ggplot
pma <- data.table(pma)
pma <- melt(pma, id = c("Sample", "FB_Ratio")
            , variable.name = "Metabolite", value.name = "Metabolite_Abundance")

# Create a treament variable
pma$Treatment <- "a"
pma[Sample %in% c("c1", "c2", "c3"),]$Treatment <- "10-Year Grassland"
pma[Sample %in% c("c4", "c5", "c6", "c7"),]$Treatment <- "1-Year Grassland"
pma[Sample %in% c("b1", "b2", "b3"),]$Treatment <- "10-Year Bare"
pma[Sample %in% c("b4", "b5", "b6", "b7"),]$Treatment <- "1-Year Bare"

# Tidy the names of the metabolites
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
pma$Metabolite <- as.character(pma$Metabolite)
pma[Metabolite == "X4.hydroxybenzoic.acid",]$Metabolite <- "4-Hydroxybenzoic Acid"
pma[Metabolite == "X3.6.anhydro.D.galactose",]$Metabolite <- "3,6-Anhydro-D-Galactose"
pma[Metabolite %in% c("xylose", "glucose", "fucose"),]$Metabolite <- firstup(pma[Metabolite %in% c("xylose", "glucose", "fucose"),]$Metabolite)
pma[Metabolite == "benzoic.acid",]$Metabolite <- "Benzoic Acid"
pma[Metabolite == "vanillic.acid",]$Metabolite <- "Vanillic Acid"

MetabResults$Metabolite <- as.character(MetabResults$Metabolite)
MetabResults[Metabolite == "X4.hydroxybenzoic.acid",]$Metabolite <- "4-Hydroxybenzoic Acid"
MetabResults[Metabolite == "X3.6.anhydro.D.galactose",]$Metabolite <- "3,6-Anhydro-D-Galactose"
MetabResults[Metabolite %in% c("xylose", "glucose", "fucose"),]$Metabolite <- firstup(MetabResults[Metabolite %in% c("xylose", "glucose", "fucose"),]$Metabolite)
MetabResults[Metabolite == "benzoic.acid",]$Metabolite <- "Benzoic Acid"
MetabResults[Metabolite == "vanillic.acid",]$Metabolite <- "Vanillic Acid"

# Change again for tidier plots
pma[Metabolite == "3,6-Anhydro-D-Galactose",]$Metabolite <- "Galactose"
MetabResults[Metabolite == "3,6-Anhydro-D-Galactose",]$Metabolite <- "Galactose"
pma[Metabolite == "4-Hydroxybenzoic Acid",]$Metabolite <- "HBA"
MetabResults[Metabolite == "4-Hydroxybenzoic Acid",]$Metabolite <- "HBA"

# Custom order of levels
pma$Metabolite <- factor(pma$Metabolite, levels = c("Vanillic Acid", "HBA", "Benzoic Acid"
                                                    , "Fucose", "Xylose", "Galactose"
                                                    , "Glucose"))

# Create a polymer variable
pma$Polymer <- ""
pma[Metabolite %in% c("HBA", "Benzoic Acid", "Vanillic Acid")]$Polymer <- "Lignin"
pma[Metabolite %in% c("Galactose", "Xylose", "Fucose")]$Polymer <- "Hemicellulose"
pma[Metabolite %in% c("Glucose")]$Polymer <- "Cellulose"

# Plot R2 values in nice places
head(pma)
PlotResults <- MetabResults[,c(1,3, 6:7)]

Metab_Max <- data.frame(tapply(pma$Metabolite_Abundance, pma$Metabolite, function(x){
  min(x) + ((max(x) - min(x))*1.1)}
  ))

Metab_Max$Metabolite <- rownames(Metab_Max)
Metab_Max$FB_Ratio <- tapply(pma$FB_Ratio, pma$Metabolite, function(x){
  min(x) + ((max(x) - min(x))*0.5)
  })

names(Metab_Max)[1] <- "Metabolite_Abundance"
PlotResults <- merge(PlotResults, Metab_Max)

levels(pma$Polymer) <- c("Lignin", "Hemicellulose", "Cellulose")
PolymerPalette <- c("#6C1706" # Lignin
                    ,"#066C0D" # Hemicellulose
                    ,"#0D066C" # Cellulose
                    ) 

# Tidy up statistics for plotting
PlotResults$pval <- PlotResults$p  
PlotResults$p <- with(PlotResults, ifelse(p < 0.001
                                          , yes = "< 0.001"
                                          , no = paste0(" = ", round(p, 3))
                                          ))
# Get parseable stats outputs
Labels <- c(by(PlotResults, PlotResults$Metabolite, function(x){
  bquote("F"["1,12"]~"="~.(round(x$FStat,2))~", p"~.(x$p)) # Produces something useful
}))

PlotResults$Label <- unlist(lapply(Labels, deparse))
PlotResults$Metabolite <- factor(PlotResults$Metabolite, levels = c("Vanillic Acid", "HBA", "Benzoic Acid"
                                                                    , "Fucose", "Xylose", "Galactose"
                                                                    , "Glucose"))
#### Plot the data ####

p <- ggplot(pma, aes(x = FB_Ratio, y = Metabolite_Abundance)) + 
  geom_point(aes(colour = Polymer, pch = Treatment)) + # Add points coloured by polymer
                                                       #, shapes according to sample
  
  facet_wrap(~ Metabolite, scales = 'free') + # One plot for each metabolite
  
  stat_smooth(data = pma[Metabolite %in% PlotResults[pval < 0.05,]$Metabolite]
              , method = "glm"
              , method.args = list(family = Gamma()) # Wonder-function
              ) + # Add linear models to the plots where the regression was significant
 
  # now weird pseudo-R^2 values
  geom_text(data = PlotResults, aes(x = FB_Ratio, y = Metabolite_Abundance
                                    , label = Label)
            , parse = TRUE
            , size = 2.5
            ) +
  # Tidying
  xlab("Fungal:Bacterial PLFA") +
  ylab("Metabolite Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  scale_colour_manual(breaks = levels(pma$Polymer)
                      , values = PolymerPalette) +
  geom_blank(data = pma, aes(x = FB_Ratio, y = Metabolite_Abundance*1.12), )
            
p

ggsave(filename = "Figures/Metabolomics/2020-08-29_PLFA_Composition_Metabolite_Abundance.pdf"
       , p
       , width = 7
       , height = 5)

#### Were the total PLFA values different between treatments? ####
head(pma)
head(pm)
setDT(pm)
pm$Treatment <- "a"
pm[Sample %in% c("c1", "c2", "c3"),]$Treatment <- "10-Year Grassland"
pm[Sample %in% c("c4", "c5", "c6", "c7"),]$Treatment <- "1-Year Grassland"
pm[Sample %in% c("b1", "b2", "b3"),]$Treatment <- "10-Year Bare"
pm[Sample %in% c("b4", "b5", "b6", "b7"),]$Treatment <- "1-Year Bare"


plot(Total_PLFA ~ factor(Treatment), data = pm)
