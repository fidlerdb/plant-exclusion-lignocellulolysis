rm(list = ls())

library(data.table)
library(ggplot2)
library(cowplot)
library(ggnewscale)

#### Data import ####

# Colour palette for the treatment variable
BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

# Function for 95% confidence intervals
mean_ci <- function(x){
  mean_se(x, mult = 1.96)
}

# Fibre analysis data
ff <- readRDS("Fibre_Analysis/Output_Data/Fibre_Data_For_Figure_2021-05-05.rds")
ff

# Metabolite data
mf <- readRDS("Fibre_Analysis/Output_Data/Metabolite_Data_For_Figure.rds")
mf

#### Data cleaning ####

## Fibre analysis data
ff

ff_d <- ff[[1]]
colsToKeep <- c("Sample", "Treatment","Proportion_Cellulose","Proportion_Hemicellulose", "Proportion_Lignin")
ff_d <- ff_d[,..colsToKeep]
setnames(ff_d, c("Sample", "Treatment","Cellulose","Hemicellulose", "Lignin"))
ff_d[, Measurement := "Fibre Analysis"]

ff_d <- melt(ff_d, id = c("Sample", "Treatment", "Measurement"), variable.name = "Polymer", value.name = "Abundance")
ff_d[, Abundance := Abundance * 100]

ff_d[,BinBase.name := NA] # For compatibility with the metabolite data

# Simple plot of the data in this format
ggplot(ff_d, aes(x = Treatment, y = Abundance)) + geom_point() + facet_grid(Polymer~Measurement
                                                                            #, nrow = 3
                                                                            , scales = "free_y"
                                                                            
)
levels(ff_d$Treatment) <- sub(" ", "\n" , levels(ff_d$Treatment)) # Two lines for plotting

## Metabolite data
mf

mf_d <- rbindlist(mf[1:3], idcol = "Polymer")

# Sort out a polymer variable
mf_d$Polymer <- as.character(mf_d$Polymer)
mf_d[Polymer == '1', Polymer := "Cellulose"]
mf_d[Polymer == '2', Polymer := "Hemicellulose"]
mf_d[Polymer == '3', Polymer := "Lignin"]

# Sort out the treatment variable
mf_d$Treatment
levels(mf_d$Treatment) <- c("1-Year\nBare", "10-Year\nBare", "1-Year\nGrassland", "10-Year\nGrassland")
mf_d$Treatment <- relevel(mf_d$Treatment, ref = "10-Year\nBare")

setnames(mf_d, c("Polymer", "Chemical", "Sample", "RawAbundance", "Treatment",  "Abundance"))
# Add a Measurement column
mf_d[, Measurement := "Metabolome"]

# Change the chemicals for tidier plots
mf_d[Chemical == "glucose", Chemical := "Glucose"]
mf_d[Chemical == "xylose", Chemical := "Xylose"]
mf_d[Chemical == "fucose", Chemical := "Fucose"]
mf_d[Chemical == "3,6-anhydro-D-galactose", Chemical := "Galactose"]
mf_d[Chemical == "vanillic acid", Chemical := "Vanillic Acid"]
mf_d[Chemical == "benzoic acid", Chemical := "Benzoic Acid"]
mf_d[Chemical == "4-hydroxybenzoic acid", Chemical := "HBA"]

#### Fibre analysis plot to add letters to  ####

pf <- ggplot(ff_d, aes(x = Treatment, y = Abundance
                       , fill = Treatment)) + 
  stat_summary(fun = "median", geom = "bar") +
  geom_point(position = position_jitter(0.15)) +
  ylab("Percentage of biomass") + 
  scale_fill_manual(values = BlackoutPalette
                    , breaks = levels(ff_d$Treatment)) +
  theme_bw() +
  theme(legend.position = "none"
        , panel.grid = element_blank()
        , axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(Polymer~Measurement, scales = "free_y")

#### Metabolite plot to add letters to ####

pm <- ggplot(mf_d, aesthetics = aes(x = Treatment, y = Abundance)) + 
  # Data
  geom_point(data = mf_d[Polymer == "Cellulose",]
             , aes(x = Treatment, y = Abundance, shape = Chemical)
             , position = position_jitter(0.05)
             , alpha = 0.5
  ) + 
  new_scale("shape") +
  geom_point(data = mf_d[Polymer == "Hemicellulose",]
             , aes(x = Treatment, y = Abundance, shape = Chemical)
             , position = position_jitter(0.05)
             , alpha = 0.5
  ) +
  new_scale("shape") +
  geom_point(data = mf_d[Polymer == "Lignin",]
             , aes(x = Treatment, y = Abundance, shape = Chemical)
             , position = position_jitter(0.05)
             , alpha = 0.5
  ) +
  
  # Mean values
  stat_summary(fun.data = 'mean_ci'
               , aes(x = Treatment, y = Abundance, colour = Treatment), size = 1.1) +
  scale_color_manual(values = BlackoutPalette
                     , breaks = levels(mf_d$Treatment)) +
  
  # Aesthetics
  facet_grid(Polymer ~ Measurement) +
  theme_bw() +
  theme(panel.grid = element_blank()
        , axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Relative abundance")  +
  guides(colour = "none")

pm

pl <- get_legend(pm)
pm <- pm + theme(legend.position = "none")

#### Plot legend ####

# Put the grobs in the correct order, also spaces them nicely
pl2 <- grid.arrange(grobs = list(pl$grobs[[2]], pl$grobs[[3]], pl$grobs[[1]])
                    , ncol = 1
)

library(gridExtra)
grid.arrange(pl)

#### Add significance letters ####

format_Letters <- function(x, polymer){
  setnames(x, c("Treatment", "Letters", "Abundance"))
  x[,Polymer := polymer]
  x[, Treatment := sub(" ", "\n", Treatment)]
  return(x)
}

fib_cel_letters <- format_Letters(ff[[2]], "Cellulose")
fib_hem_letters <- format_Letters(ff[[3]], "Hemicellulose")
fib_lig_letters <- format_Letters(ff[[4]], "Lignin")
fib_lig_letters$Abundance <- 12

pf2 <- pf + geom_text(data = fib_cel_letters, aes(x = Treatment, y = Abundance, label = Letters), colour = 'black') +
  geom_text(data = fib_hem_letters, aes(x = Treatment, y = Abundance, label = Letters), colour = 'black') +
  geom_text(data = fib_lig_letters, aes(x = Treatment, y = Abundance, label = Letters), colour = 'black') #+
# ylim(-1,100)

pm

mf[[4]] # Cel
mf[[5]] # Hem
mf[[6]] # Lig

format_letters_met <- function(x, polymer){
  setDT(x)
  x[,Polymer := polymer]
  return(x)
}

met_cel_letters <- format_letters_met(mf[[4]], "Cellulose")
met_hem_letters <- format_letters_met(mf[[5]], "Hemicellulose")
met_lig_letters <- format_letters_met(mf[[6]], "Lignin")

met_cel_letters$Letters <- c("b", "a", "ab", "ab") # Re-jig these so that they work the normal way round
met_lig_letters$Letters <- c("b", "a", "c", "ac") # Re-jig these so that they work the normal way round

pm2 <- pm + geom_text(data = met_cel_letters, aes(x = Treatment, y = NormAbund, label = Letters)) +
  geom_text(data = met_hem_letters, aes(x = Treatment, y = NormAbund, label = Letters)) + 
  geom_text(data = met_lig_letters, aes(x = Treatment, y = NormAbund, label = Letters))

#### Create the finalised plot ####
Figure <- plot_grid(pf2, pm2, pl2, ncol = 3, rel_widths = c(1,1,0.6))

Figure

#### PLFA - metabolite abundance ####

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
         mod2 <- lm(unlist(pmb[,14:length(pmb)][,METABOLITE]) ~ FB_Ratio + Total_PLFA
                    , data = pmb)
         data.frame(AIC(mod1,mod2), DeltaAIC = AIC(mod1) - AIC(mod2))
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
  min(x) + ((max(x) - min(x))*1.4)}
))

Metab_Max$Metabolite <- rownames(Metab_Max)
Metab_Max$FB_Ratio <- tapply(pma$FB_Ratio, pma$Metabolite, function(x){
  min(x) + ((max(x) - min(x)) * 0.5
            )
})

names(Metab_Max)[1] <- "Metabolite_Abundance"
PlotResults <- merge(PlotResults, Metab_Max)

pma$Polymer <- factor(pma$Polymer, levels = c("Cellulose","Hemicellulose","Lignin"))
#levels(pma$Polymer) <- c("Lignin", "Hemicellulose", "Cellulose")
PolymerPalette <- c("#0D066C" # Cellulose
                    ,"#066C0D" # Hemicellulose
                    ,"#6C1706" # Lignin
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
            , size = 1.75
  ) +
  geom_point(data = PlotResults, aes(x = FB_Ratio, y = Metabolite_Abundance*1.15
                                    )
            , size = 0.1
            , colour = "white"
  ) +
  # Tidying
  #xlab("Fungal:Bacterial PLFA") +
  scale_x_continuous(name = "Fungal:Bacterial PLFA"
                     , breaks = c(0.05, 0.07, 0.09)
                     #, n.breaks = 3
                     , limits = c(0.045, 0.09)
                     ) +
  
  ylab("Metabolite Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  scale_colour_manual(breaks = levels(pma$Polymer)
                      , values = PolymerPalette) +
  geom_blank(data = pma, aes(x = FB_Ratio, y = Metabolite_Abundance*1.12))

p


#### Put them together ####

p_final <- plot_grid(Figure, p
                     , labels = c("a", "b"))

save_plot(filename = "Final_figure_making/Multipanel_figures_Oct_2021/Figure_1_Fibre_metabolome_PLFA.pdf"
          , plot = p_final
          , base_height = 4.15
          , base_width = 11.7
          )
