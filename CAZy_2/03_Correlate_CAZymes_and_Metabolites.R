rm(list = ls()) # Clear the workspace

# The point of this script is to find which CAZy families are correlated with 
#     the abundance of certain metabolites which are the breakdown products of 
#     elements of lignocellulose

# Load useful packages
library(data.table)
source("Functions/range01.R")

# Load the data: CAZyme reads and metabolome data
df <- readRDS("Metabolomics/Cleaned Data/Henfaes_CAZy_Metabolome_Reads.rds")

# Check the data
head(df)

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

# Perform p-adjustment
AdjustedCorrelations <- MetaboliteList[MetaboliteList$p_adj < 0.05,]
RawCorrelations <- MetaboliteList[MetaboliteList$p < 0.05,]

#### Get the results ####

# How many of each CAZy type
length(unique(RawCorrelations$CAZymes[grep("GH", RawCorrelations$CAZymes)]))
length(unique(RawCorrelations$CAZymes[grep("AA", RawCorrelations$CAZymes)]))
length(unique(RawCorrelations$CAZymes[grep("CBM", RawCorrelations$CAZymes)]))

# How many significant correlations overall
unique(RawCorrelations$CAZymes)
hist(RawCorrelations$r)

# Look at the p-values
hist(RawCorrelations$p)
nrow(RawCorrelations)
nrow(RawCorrelations[r > 0,])
nrow(RawCorrelations[r < 0,])

# Look at the adjusted p-values
hist(RawCorrelations$p_adj)
RawCorrelations[p_adj < 0.05,]

RawCorrelations[r > 0,]

PositiveCorrelations <- RawCorrelations[r > 0,]

# check encoding to ensure special characters are preserved
Encoding(PositiveCorrelations$CAZymes)

#### Write the correlated CAZymes to a file ####

fwrite(PositiveCorrelations, file = "CAZy_2/outputData/Positively_Correlated_CAZymes-Metabolites.csv")
fwrite(RawCorrelations, file = "CAZy_2/outputData/All_Correlated_CAZymes-Metabolites.csv")

