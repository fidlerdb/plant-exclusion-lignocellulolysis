#### Discover abolute abundances of the metabolome data for supplementary ####

# Clean the workspace
rm(list = ls())

# Load packages
library(data.table)
library(magrittr)

#### Data Import and Cleaning #### 

# Load in the data
mf <- data.frame(fread("Metabolomics/Data/Metabolome_Blackout_Raw.csv"))
mf <- mf[, c(1,7:20)] # Keep only the most useful columns

setDT(mf) # Make mf a data.table

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

# Subset to keep only the metabolites that we are interested in
mfl <- mfl[BinBase.name %in% c('glucose'
                        , 'xylose', 'fucose', '3,6-anhydro-D-galactose'
                        , 'vanillic acid', '4-hydroxybenzoic acid'
                        , 'benzoic acid'),]

# Check it
mfl

# Find the minimum and maximum abundance of each metabolite that we 
# are interested in
absAbundance <- lapply(c('glucose', 'xylose', 'fucose'
                       , '3,6-anhydro-D-galactose'
                       , 'vanillic acid', '4-hydroxybenzoic acid'
                       , 'benzoic acid')

                       , FUN = function(x){
                         out <- data.table(Metabolite = x
                         , Min = min(mfl[BinBase.name == x,]$Abundance)
                         , Max = max(mfl[BinBase.name == x,]$Abundance)
                         )
                         return(out)
                         }) %>% rbindlist(.)

absAbundance

# Write the output to file
fwrite(absAbundance
       , file = "Metabolomics/Results/Metabolite_Absolute_Abundance.csv")
fwrite(absAbundance
       , file = "Tables/Metabolite_Absolute_Abundance.csv")
