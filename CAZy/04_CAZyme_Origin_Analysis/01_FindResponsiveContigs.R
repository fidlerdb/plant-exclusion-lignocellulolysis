
# Read in the data with contig labels and CAZy IDs
k <- read.csv("../dbCAN_Results/Full_Assembly/Mapping Output/Blackout_ReadCounts_CAZyAnnotations.csv")

# Check this file
head(k)

# Responsive CAZymes (> 1 log2 fold change, p_adj < 0.05, greater read abundance than abundance of 10% 
# quantile of all other CAZymes)
RC <- c("GH30_5", "GH43_1", "GH81", "GH5_4", "GH13_21", "GH114", "GH43_9", "GH5_53", "GH135", "GH5_35"
        , "CBM8", "GH43_37", "GH55|3.2.1.-", "CBM48|GH13_13", "CBM61", "GH136", "CBM4|GH9", "CBM6|GH16"
        , "CBM20|CBM34|GH13_39", "CBM22|CBM6|GH10")

# Check these are all valid CAZyme names
library(roperators)
any(RC %ni% k$CAZyme)

# How may CAZymes are we dealing with?
nrow(k[k$CAZyme %in% RC,])
# 143 over 20 CAZy families is beleivable--check this in Exploratory_Analysis/GH_intra_family_exploration.R

# Keep only contigs which were responsive to treatment
ContigsToBLAST <- k[k$CAZyme %in% RC,]

write.csv(unique(ContigsToBLAST$Chr), file = "CAZy/Cleaned_Data/ResponsiveCAZyContigsToBLAST.txt"
          , quote = FALSE
          , row.names = FALSE
          , col.names = FALSE)

unique(ContigsToBLAST$Chr)

