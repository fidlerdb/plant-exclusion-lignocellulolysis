rm(list = ls())

# This script investigates the relative abundances of 
# breakdown products of cellulose hemicellulose and liugnin

# Load packages
library(data.table)
library(ggplot2)
library(cowplot)
library(randomForest)

#### Data Import and Cleaning #### 

# Load in the data
mf <- data.frame(fread("Metabolomics/Data/Metabolome_Blackout_Raw.csv"))
# of <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_Orders_CAZymes.rds")
of <- data.frame(fread("CAZy/Cleaned_Data/Henfaes_CAZyFamilyReadsPerKbp.csv"))

# Scale CAZyme reads by numer of reads in the library (billions)
of[,4:length(of)] <- of[,4:length(of)] / (of$Reads/1e9)

head(mf)
head(of)
head(sf[1:5])

mf <- mf[, c(1,7:20)] # Keep only the most useful columns

# Lignocellulose breakdown products
BreakdownProducts <- c('glucose'
                       , 'xylose', 'fucose', '3,6-anhydro-D-galactose'
                       , 'vanillic acid', '4-hydroxybenzoic acid', 'benzoic acid')

# Keep only those products
mf <- mf[mf$BinBase.name %in% BreakdownProducts,]

head(mf)

head(mf)
AbundData <- data.frame(t(mf[2:length(mf)]))
names(AbundData) <- mf[,1]
AbundData$Sample <- row.names(AbundData)


ofm <- merge(of, AbundData, by = "Sample")
head(ofm[1:5])
head(ofm[,300:311])
#names(ofm)[44] <- "Bacteroidetes_Order_II"

# Save this as a file
saveRDS(ofm, file = "Metabolomics/Cleaned Data/Henfaes_CAZy_Metabolome_Reads.rds")

