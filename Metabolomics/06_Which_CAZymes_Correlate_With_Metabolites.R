rm(list = ls()) # Clear the workspace

# The point of this script is to find which CAZy families are correlated with 
#     the abundance of certain metabolites which are the breakdown products of 
#     elements of lignocellulose

# Load useful packages
library(data.table)
library(ggplot2)
library(ggforce)
library(ggpubr)
source("Functions/range01.R")

# Load the data: CAZyme reads and metabolome data
df <- readRDS("Metabolomics/Cleaned Data/Henfaes_CAZy_Metabolome_Reads.rds")

# Check the data
head(df)

#### Look at the relationships #####

# Create a plot to visualize the correlations between 
# enzymes and metabolites
MetabolitePlot <- function(metabolite, pageNumber){
  
  if(missing(pageNumber)){pageNumber <- 1} # Deal with lazy users
  
  Enzymes <- df[, 4:304] # Keep the enzymes
  Metabolite <- df[,names(df) %in% metabolite] # Keep the metabolite of interest
  
  working_df <- cbind(Metabolite, Enzymes) # Put the two together

  working_df <- melt(working_df, id = "Metabolite") # Make ggplot2 able to handle the data
  
  # Create the plots
  p <- ggplot(data = working_df, aes(x = value, y = Metabolite)
              ) + 
    geom_point() + 
    facet_wrap_paginate(~ variable, scales = 'free_x'
                        , ncol = 6
                        , nrow = 6
                        , page = pageNumber) + 
    theme_bw() +
    ylab(metabolite) #+
    #geom_smooth(method = 'lm')
  
  # And give the user the output
  print(p)
}

#MetabolitePlot(metabolite = "vanillic acid")

# Plots for vanillic acid
#for(i in 1:9){MetabolitePlot(metabolite = "vanillic acid", pageNumber = i)}

# There are so many plots. I cannae make sense of all that data visually.
# Let's get some regression coefficients.


#### Set up a data frame of the metabolites and the CAZymeswith scaled ####
#### predictors                                                        ####

# These are the metabolites we want to explore
metabolites <- c('xylose'
,'vanillic acid'
,'glucose'
,'fucose'
,'benzoic acid'
,'4-hydroxybenzoic acid'
,'3,6-anhydro-D-galactose')

# These are all of the CAZymes
Enzymes <- df[, 4:304] 
Metabolite <- df[,names(df) %in% metabolites] 


# Now scale the predictor and response data between 0 and 1 to make   
# the outputs of the models (with a single response variable) comparable

# Run the models on all predictors individually to assess significance and 
# correlation coefficient.

# Create the dataframe with a nested structure
regression_df <- data.frame(Metabolites = I(Metabolite)
           , Enzymes = I(Enzymes))

# Ensure that the structure of the data frame is correct
str(regression_df$Enzymes)
str(regression_df$Metabolites)
regression_df$Metabolites <- data.frame(lapply(regression_df$Metabolites
                                               , FUN = function(x){
                                                 as.numeric(as.character(x))
                                               }
))
# Scale predictor variables between 0 and 1
regression_df$Enzymes <- apply(regression_df$Enzymes+1, 2, log2)
regression_df$Enzymes <- data.frame(apply(regression_df$Enzymes, 2, range01))
regression_df$Metabolites <- data.frame(apply(regression_df$Metabolites, 2, range01))

names(regression_df$Metabolites) <- names(Metabolite)
regression_df$Metabolites

regression_df$Enzymes[1:9]

#### Run the regressions as a first screening process ####
# Now loop over each metabolite, creating a model for the relationship 
# between each predictor and the metabolite. Store the p and r values.
# Chose quasibinomial because the standard binomial distribution + logit 
# gave dispersion values of 0.1 - 0.2

Gene_Function_Results <- apply(regression_df$Metabolites, 2
                               , FUN = function(metabolite){
                                 
                                 # For each metabolite...
                                 
                                 apply(regression_df$Enzymes, 2
                                       , FUN = function(CAZyme){
                                         
                                         # Fit n models with each CAZyme as a  
                                         # seperate predictor
                                         
                                         model <- glm(metabolite ~ CAZyme
                                                      , data = df
                                                      #, family = quasibinomial(link = "logit")
                                                      )

                                         # And keep the results
                                         
                                         Result <- drop1(model, test = "F"#'Chisq'
                                                         )
    

                                         data.frame(p = #Result$`Pr(>Chi)`[2]
                                                    Result$`Pr(>F)`[2]
                                           , r = model$coefficients[2]
                                           )
                                         
                                         })
  })

# Create a vector of CAZymes used as predictors
CAZymes <- unlist(lapply(Gene_Function_Results, names))

# Make this easier to work with by sticking all results  
# within each metabolite end-on-end
Gene_Function_Results <- lapply(Gene_Function_Results
                                , FUN = function(x){
                                  do.call(rbind, x)
                                  })

# Now adjust the p-values using FPR adjustment for each
p_adj <- lapply(Gene_Function_Results
                , FUN = function(x){
                  x$p_adj <- p.adjust(x$p, method = 'fdr')
                  })


# Add adjusted p-values to the list
Gene_Function_Results <- Map(cbind # function
                             , Gene_Function_Results # data to add to
                             , p_adj = p_adj # data to add
                             ) 

str(Gene_Function_Results)
MetaboliteList <- rbindlist(Gene_Function_Results, idcol = "Metabolite")
MetaboliteList <- cbind(CAZymes, MetaboliteList)
nrow(MetaboliteList)
table(MetaboliteList[MetaboliteList$p < 0.05]$Metabolite)

#head(CAZyList)

AdjustedCorrelations <- MetaboliteList[MetaboliteList$p_adj < 0.05,]
RawCorrelations <- MetaboliteList[MetaboliteList$p < 0.05,]

# p-adjusted hits
table(AdjustedCorrelations[grep("AA", AdjustedCorrelations$CAZymes),]$Metabolite)
table(AdjustedCorrelations[grep("GH", AdjustedCorrelations$CAZymes),]$Metabolite)
table(AdjustedCorrelations[grep("CBM", AdjustedCorrelations$CAZymes),]$Metabolite)

# Positive associations
table(AdjustedCorrelations[grepl("AA", AdjustedCorrelations$CAZymes) &
                        AdjustedCorrelations$r > 0 ,]$Metabolite)
table(AdjustedCorrelations[grepl("GH", AdjustedCorrelations$CAZymes)&
                        AdjustedCorrelations$r > 0 ,]$Metabolite)
table(AdjustedCorrelations[grepl("CBM", AdjustedCorrelations$CAZymes)&
                        AdjustedCorrelations$r > 0 ,]$Metabolite)

# unadjusted hits
table(RawCorrelations[grep("AA", RawCorrelations$CAZymes),]$Metabolite)
table(RawCorrelations[grep("GH", RawCorrelations$CAZymes),]$Metabolite)
table(RawCorrelations[grep("CBM", RawCorrelations$CAZymes),]$Metabolite)

# Positive associations
table(RawCorrelations[grepl("AA", RawCorrelations$CAZymes) &
                        RawCorrelations$r > 0 ,]$Metabolite)
table(RawCorrelations[grepl("GH", RawCorrelations$CAZymes)&
                        RawCorrelations$r > 0 ,]$Metabolite)
table(RawCorrelations[grepl("CBM", RawCorrelations$CAZymes)&
                        RawCorrelations$r > 0 ,]$Metabolite)

MetaboliteList[MetaboliteList$p < 0.05,]
length(regression_df$Enzymes)
length(RawCorrelations$CAZymes)
length(unique(RawCorrelations$CAZymes))

length(grep("GH", names(regression_df$Enzymes)))
length(grep("AA", names(regression_df$Enzymes)))
length(grep("CBM", names(regression_df$Enzymes)))

length(unique(RawCorrelations$CAZymes[grep("GH", RawCorrelations$CAZymes)]))
length(unique(RawCorrelations$CAZymes[grep("AA", RawCorrelations$CAZymes)]))
length(unique(RawCorrelations$CAZymes[grep("CBM", RawCorrelations$CAZymes)]))

unique(RawCorrelations$CAZymes)
head(RawCorrelations)
hist(RawCorrelations$r)

hist(RawCorrelations$p)
nrow(RawCorrelations)
nrow(RawCorrelations[r > 0,])
nrow(RawCorrelations[r < 0,])

hist(RawCorrelations$p_adj)
RawCorrelations[p_adj < 0.05,]

write.csv(RawCorrelations
          , file = "CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Positively_Correlated_CAZymes-Metabolites.csv"
          , row.names = FALSE)


#### Old ####

# Which genes were significantly correlated with increased  
# or decreased abundances of any metabolite?
# # (1) Keep significant results
# CorrelatedCAZy <- do.call(rbind
#                           , lapply(Gene_Function_Results
#                                    , FUN = function(x){x[x$p_adj < 0.05,]}
#                                    )
#                           )
# # (2) Make the result interpretable by  information about  which CAZyme  
# #     was the predictor
# # CorrelatedCAZy <- cbind(data.frame(
# #   Metabolite = sub("\\..*", "", row.names(CorrelatedCAZy))
# #   , CAZyme = gsub("^.*?\\.", "", row.names(CorrelatedCAZy))
# #   )
# #   , CorrelatedCAZy
# # )
#       
# # plot(Enzymes$GH43_11 ~ Metabolites$fucose, data = regression_df)
# 
# # Look at the results
# CorrelatedCAZy
# 
# unique(CorrelatedCAZy$Metabolite)
# hist(CorrelatedCAZy$p_adj)
# # Ensure all metabolites are found
# CorrelatedCAZy$Metabolite <- as.character(CorrelatedCAZy$Metabolite)
# CorrelatedCAZy$Metabolite[CorrelatedCAZy$Metabolite 
#                           == 'vanillic'] <- 'vanillic acid'
# CorrelatedCAZy$Metabolite[CorrelatedCAZy$Metabolite 
#                           == 'benzoic'] <- 'benzoic acid'
# CorrelatedCAZy$Metabolite[CorrelatedCAZy$Metabolite 
#                           == 'X4'] <- '4-hydroxybenzoic acid'
# 
# #### Set up the data for plotting ####
# 
# # Keep only the significantly correlates columns (i.e. drop both metabolites 
# # and CAZymes which were not in a significant relationship)
# plot_df <- df[,names(df) %in% 
#                 unique(c(as.character(CorrelatedCAZy$Metabolite)
#                          , as.character(CorrelatedCAZy$CAZyme)))
#               ]
# 
# # Make it ggplot2able... (very long format and named properly)
# # Aiming for Metabolite, Abundance, CAZyme, and Reads columns 
# 
# names_metab <- names(plot_df)[15:16]
# plot_df  <- melt(plot_df, id = names_metab)
# names(plot_df)[3:4] <- c("CAZyme", "Reads")
# plot_df <- melt(plot_df, id = c("CAZyme", "Reads"))
# names(plot_df)[3:4] <- c("Metabolite", "Abundance")
# plot_df$Abundance <- as.numeric(plot_df$Abundance)
# 
# # Check it
# head(plot_df)
# head(CorrelatedCAZy)
# 
# # Make an index column to filter down to which combinations are wanted
# cc_MetCaz <- paste0(CorrelatedCAZy$Metabolite, CorrelatedCAZy$CAZyme)
# 
# # Make an index column to filter down to which combinations are wanted
# plot_df$Met_Caz <- paste0(plot_df$Metabolite, plot_df$CAZyme)
# 
# # Keep only these combinations
# plot_df <- plot_df[plot_df$Met_Caz %in% cc_MetCaz,]
# 
# # Get the rounded (2 decimal places) model R^2 values
# #CorrelatedCAZy$adj_r_squared <- round(CorrelatedCAZy$adj_r_squared, digits = 2)
# str(plot_df)
# unique(plot_df$Metabolite)
# # For plotting with Hadley's lameness
# CorrelatedCAZy$Reads <- 0.5
# CorrelatedCAZy$Abundance <- 35000
# # plot_df$NormReads <- tapply(plot_df$Reads, plot_df$Met_Caz
# #                             , range01)
# plot_df$NormReads <- with(plot_df, ave(Reads, Met_Caz, FUN = range01))
# 
# #### Plot the relationships ####
# 
#  p <- ggplot(plot_df, aes(x = NormReads, y = Abundance
#                           , colour = Metabolite)) + 
#   geom_point() + 
#    
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   geom_smooth(method = 'lm')
# 
#   p + facet_wrap_paginate( ~ CAZyme + Metabolite, scales = 'free'
#                        , ncol = 4, nrow = 3, page = 1)
#   p + facet_wrap_paginate( ~ CAZyme + Metabolite, scales = 'free'
#                            , ncol = 4, nrow = 3, page = 2) 
#   p + facet_wrap_paginate( ~ CAZyme + Metabolite, scales = 'free'
#                            , ncol = 4, nrow = 3, page = 3)
#   p + facet_wrap_paginate( ~ CAZyme + Metabolite, scales = 'free'
#                            , ncol = 4, nrow = 3, page = 4)
#   p + facet_wrap_paginate( ~ CAZyme + Metabolite, scales = 'free'
#                            , ncol = 4, nrow = 3, page = 5)
#   p + facet_wrap_paginate( ~ CAZyme + Metabolite, scales = 'free'
#                            , ncol = 4, nrow = 3, page = 6)
#   p + facet_wrap_paginate( ~ CAZyme + Metabolite, scales = 'free'
#                            , ncol = 4, nrow = 3, page = 7)
# #### Now run better analyses considering which CAZymes may ####
# #### actually have an effect                               ####
# 
# # Which enzymes are statistically correlated with xylose abundance?
# # Use a gaussian GLM as all predictors and response variables are 
# # scaled to allow comparison of relative strengths of the different 
# # models' and predictors' relative effects
# 
# #### Analyse the corrlations with xylose ####
# 
# # Fit a multiple regression 
# hist(regression_df$Metabolites$xylose)
# mod_xyl <- glm(regression_df$Metabolites$xylose
#   #cbind(regression_df$Metabolites$xylose
#    #                  , 1- regression_df$Metabolites$xylose)
#                ~ 
#                 regression_df$Enzymes$CBM50 
#               + regression_df$Enzymes$GH43_11 
#               + regression_df$Enzymes$GH73
#               , family = binomial(link = 'logit')
#               )
# 
# # car::vif test for colinearity, cutoff of 5 commonly used
# car::vif(mod_xyl) # Keep all terms
# 
# # Check assumptions
# par(mfrow=c(2,2));plot(mod_xyl);par(mfrow=c(1,1))
# car::qqPlot(residuals(mod_xyl, type = 'deviance'))
# deviance(mod_xyl)/df.residual(mod_xyl) # Underdispersion = conservative model
# 
# # Basic model results
# summary(mod_xyl)
# mod_xyl
# 
# # Good p-values
# drop1(mod_xyl, test = 'Chisq')
# 
# # Refit model by stepwise deletion until all terms are siggy.
# mod_xyl <- update(mod_xyl, formula = "regression_df$Metabolites$xylose ~ 
#                 regression_df$Enzymes$CBM50 + regression_df$Enzymes$GH43_11")
# 
# # Re-check assumptions
# car::qqPlot(residuals(mod_xyl, type = 'pearson'))
# par(mfrow=c(2,2));plot(mod_xyl);par(mfrow=c(1,1))
# 
# # model Results
# summary(mod_xyl)
# drop1(mod_xyl, test = 'Chisq')
# 
# # CBM50 binds to chitin or peptidoglycan. LysM domain.
# # GH43_11 activities: EC3.2.1.37  removes D-xylose from nonreducing end, 
# #                     EC3.2.1.55  removes α-L-arabinofuranoside from nonreducing end, 
# #                     EC3.2.1.8   Endohydrolysis of (1→4)-β-D-xylosidic linkages in 
# #                                 xylans
# #                     EC3.2.1.99  Endohydrolysis of (1→5)-α-arabinofuranosidic 
# #                                 linkages in (1→5)-arabinans
# #                     EC3.2.1.145 Hydrolysis of terminal, non-reducing β-D-galactose 
# #                                 residues in (1→3)-β-D-galactopyranans
# #                     EC3.2.1.1   Endohydrolysis of (1→4)-α-D-glucosidic linkages in 
# #                                 polysaccharides containing three or more 
# #                                 (1→4)-α-linked D-glucose units
# 
# ##### Fucose ####
# 
# # Which enzymes are statistically correlated with xylose abundance?
# # Use a gaussian GLM as all predictors and response variables are 
# # scaled to allow comparison of relative strengths of the different 
# # models' and predictors' relative effects
# 
# p # Look at the results visually
# 
# # Fit a multiple regression considering all main effects
# mod_fuc <- lm(regression_df$Metabolites$fucose ~ 
#                 #regression_df$Enzymes$AA10.CBM73 
#               regression_df$Enzymes$CBM20 
#               + regression_df$Enzymes$CBM50
# #              + regression_df$Enzymes$GH153 
#               + regression_df$Enzymes$GH37
#               + regression_df$Enzymes$GH43_11 
#               + regression_df$Enzymes$GH35
#               + regression_df$Enzymes$GH73 
#               )
# 
# car::vif(mod_fuc) # Doing this repeatedly gets rid of AA10, GH153 (vif > 10)
# 
# # Check assumptions
# par(mfrow=c(2,2));plot(mod_fuc);par(mfrow=c(1,1))
# car::qqPlot(residuals(mod_fuc, type = 'pearson'))
# 
# # modelcoefficients
# summary(mod_fuc)
# 
# # Good p-values
# drop1(mod_fuc, test = 'Chisq')
# 
# # Refit model by stepwise deletion.
# # First deletion - remove AA10.CBM73
# mod_fuc <- update(mod_fuc, formula = "regression_df$Metabolites$fucose ~ 
#                regression_df$Enzymes$CBM20  + regression_df$Enzymes$CBM50 + regression_df$Enzymes$GH153  + regression_df$Enzymes$GH37 + regression_df$Enzymes$GH43_11  + regression_df$Enzymes$GH35 + regression_df$Enzymes$GH73 " )
# # Check results
# drop1(mod_fuc, test = 'Chisq')
# 
# # Second deletion - remove GH35
# mod_fuc <- update(mod_fuc, formula = "regression_df$Metabolites$fucose ~ 
#                regression_df$Enzymes$CBM20  + regression_df$Enzymes$CBM50 + regression_df$Enzymes$GH153  + regression_df$Enzymes$GH37 + regression_df$Enzymes$GH43_11 + regression_df$Enzymes$GH73 " )
# # Check results
# drop1(mod_fuc, test = 'Chisq')
# 
# # Third deletion - remove GH37
# mod_fuc <- update(mod_fuc, formula = "regression_df$Metabolites$fucose ~ 
#                regression_df$Enzymes$CBM20  + regression_df$Enzymes$CBM50 + regression_df$Enzymes$GH153 + regression_df$Enzymes$GH43_11 + regression_df$Enzymes$GH73 " )
# # Check results
# drop1(mod_fuc, test = 'Chisq')
# 
# # Fourth deletion - remove GH43_11
# mod_fuc <- update(mod_fuc, formula = "regression_df$Metabolites$fucose ~ 
#                regression_df$Enzymes$CBM20  + regression_df$Enzymes$CBM50 + regression_df$Enzymes$GH153 + regression_df$Enzymes$GH73 " )
# # Check results
# drop1(mod_fuc, test = 'Chisq')
# 
# # Fifth deletion - remove CBM20
# mod_fuc <- update(mod_fuc, formula = "regression_df$Metabolites$fucose ~ 
#                regression_df$Enzymes$CBM50 + regression_df$Enzymes$GH153 + regression_df$Enzymes$GH73 " )
# # Check results
# drop1(mod_fuc, test = 'Chisq')
# 
# # Fifth deletion - remove GH73
# mod_fuc <- update(mod_fuc, formula = "regression_df$Metabolites$fucose ~ 
#                regression_df$Enzymes$CBM50 + regression_df$Enzymes$GH153")
# # Check results
# drop1(mod_fuc, test = 'Chisq') # Final model. GH153 and CBM50 well correlated with fucose abundance in our soil
# summary(mod_fuc)
# # Re-check assumptions
# par(mfrow=c(2,2));plot(mod_fuc);par(mfrow=c(1,1))
# car::qqPlot(residuals(mod_fuc, type = 'pearson'))
# 
# # CBM50 binds to chitin or peptidoglycan. LysM domain.
# # GH153 poly-β-1,6-D-glucosamine hydrolase. EC3.2.1.3 works faster on 
# #       polysaccharides than oligosaccharides and cleaves 1,6-beta-glycosidic bonds  
# #       when the next bond is a 1,4-beta-glycosidic bond.
# 
# 
# 
# 
# summary(mod_fuc)
# summary(mod_xyl)
# 
# 
# 
