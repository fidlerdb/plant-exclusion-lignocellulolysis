########################################################################################
# The purpose of this script is to find out which phyla group together based on which  # 
# genes for CAZymes(which were correlated with the abundance of a lignocellulose       #
# breakdown product) they have                                                         #
########################################################################################

library(data.table)
library(analogue)
library(clustsig)

#### (1) Import, check, and set the data up ####
rm(list = ls()) # Clear the workspace

PhylumClusters <- fread("CAZy/04_CAZyme_Origin_Analysis/Results/Clusters_of_Phyla.csv")
names(PhylumClusters)[5:length(PhylumClusters)] <- c("galactose", "vanillic_acid"
                                                     , "benzoic_acid", "hydroxybenzoic_acid")

source("Functions/Import_and_Clean_Metabolite-CAZy_Taxonomy.R")

#### (2) Find statistically distinct clusters for all metabolites ####

source("Functions/CAZymeClustering.R")

# Cellulose
GlucoseResults <- CAZymeClustering("glucose")

# Hemicellulose
FucoseResults <- CAZymeClustering("fucose")
XyloseResults <- CAZymeClustering("xylose")
GalactoseResults <- CAZymeClustering("3,6-anhydro-D-galactose")

# Lignin
VanillaResults <- CAZymeClustering("vanillic acid")
BenzResults <- CAZymeClustering("benzoic acid")
HydroxResults <- CAZymeClustering("4-hydroxybenzoic acid")

#### (3) Look at the results ####

source("Functions/ClusterSummarizer.R")
source("Functions/PlotClusterResults.R")

# trbl
par(mar = c(12, 0,0,0))

## Cellulose
simprof.plot(GlucoseResults$simprofResults) # View clusters
GlucoseResults$simprofResults$significantclusters # CAZymes in the clusters
GlucoseResults$Total_Genes # Total number of genes in each cluster
sum(GlucoseResults$Total_Genes)
GlucoseResults$Gene_Richness # Number of genes by phylum
GlucoseCluster <- ClusterSummarizer(PhylumClusters$glucose, GlucoseResults)
GlucoseCluster$MeanGenes
rownames(GlucoseCluster$PercentageGenes)
round(GlucoseCluster$PercentageGenes,2)
PlotClusterResults(GlucoseCluster$PercentageGenes, nrow = 2, ncol = 3) # What CAZyme clusters were most abundant for the phylum clusters?
# How many phyla in the cluster?
sapply(strsplit(rownames(GlucoseCluster$PercentageGenes), " ", ), length
       , USE.NAMES = TRUE)

## Hemicellulose
# Fucose
simprof.plot(FucoseResults$simprofResults) # View clusters
FucoseResults$simprofResults$significantclusters  # CAZymes in the clusters
FucoseResults$Total_Genes # Total number of genes in each cluster
FucoseResults$Gene_Richness # Number of genes by phylum
FucoseCluster <- ClusterSummarizer(PhylumClusters$fucose, FucoseResults)
FucoseCluster$MeanGenes
row.names(FucoseCluster$PercentageGenes)
round(FucoseCluster$PercentageGenes,2)
PlotClusterResults(FucoseCluster$PercentageGenes, nrow = 2, ncol = 3) # What CAZyme clusters were most abundant for the phylum clusters?
sapply(strsplit(rownames(FucoseCluster$PercentageGenes), " ", ), length
       , USE.NAMES = TRUE)
# Xylose
simprof.plot(XyloseResults$simprofResults)
XyloseResults$simprofResults$significantclusters
XyloseResults$Total_Genes
XyloseResults$Gene_Richness

XyloseCluster <- ClusterSummarizer(PhylumClusters$xylose, XyloseResults)
round(XyloseCluster$MeanGenes,2) # Mean number of genes the phyla have

row.names(XyloseCluster$PercentageGenes)
XyloseCluster$PercentageGenes
PlotClusterResults(XyloseCluster$PercentageGenes, nrow = 3, ncol = 3) # What CAZyme clusters were most abundant for the phylum clusters?

#### Problem
# 3,6-Anhydro-D-Galactose
simprof.plot(GalactoseResults$simprofResults)
GalactoseResults$simprofResults$significantclusters
GalactoseResults$Total_Genes
GalactoseResults$Gene_Richness
GalactoseCluster <- ClusterSummarizer(PhylumClusters$galactose, GalactoseResults)
round(GalactoseCluster$MeanGenes,2) # Mean number of genes the phyla have

row.names(GalactoseCluster$PercentageGenes)
GalactoseCluster$PercentageGenes
PlotClusterResults(GalactoseCluster$PercentageGenes, nrow = 2, ncol = 3) # What CAZyme clusters were most abundant for the phylum clusters?


## Lignin
# Vanillic Acid
simprof.plot(VanillaResults$simprofResults)
VanillaResults$simprofResults$significantclusters
VanillaResults$Total_Genes
VanillaResults$Gene_Richness
VanillaCluster <- ClusterSummarizer(PhylumClusters$vanillic_acid , VanillaResults)

round(VanillaCluster$MeanGenes,2) # Mean number of genes the phyla have

row.names(VanillaCluster$PercentageGenes)
VanillaCluster$PercentageGenes
PlotClusterResults(VanillaCluster$PercentageGenes, nrow = 3) # What CAZyme clusters were most abundant for the phylum clusters?

# Benzoic Acid
simprof.plot(BenzResults$simprofResults)
BenzResults$simprofResults$significantclusters
BenzResults$Total_Genes
BenzResults$Gene_Richness
BenzCluster <- ClusterSummarizer(PhylumClusters$benzoic_acid, BenzResults)
round(BenzCluster$MeanGenes,2) # Mean number of genes the phyla have

row.names(BenzCluster$PercentageGenes)
BenzCluster$PercentageGenes

PlotClusterResults(BenzCluster$PercentageGenes, nrow = 2, ncol = 3) # What CAZyme clusters were most abundant for the phylum clusters?

# 4-Hydroxybenzoic Acid
simprof.plot(HydroxResults$simprofResults)
HydroxResults$simprofResults$significantclusters
HydroxResults$Total_Genes
HydroxResults$Gene_Richness
HydroxCluster <- ClusterSummarizer(PhylumClusters$hydroxybenzoic_acid, HydroxResults)
round(HydroxCluster$MeanGenes,2) # Mean number of genes the phyla have

row.names(HydroxCluster$PercentageGenes)
round(HydroxCluster$PercentageGenes, 2)
PlotClusterResults(HydroxCluster$PercentageGenes, nrow = 3) # What CAZyme clusters were most abundant for the phylum clusters?

#### Create a list of which CAZymes are in which cluster ####

# Name the clusters 
GluClusters <- GlucoseResults$Clusters
FucClusters <- FucoseResults$Clusters 
XylClusters <- XyloseResults$Clusters 
GalClusters <- GalactoseResults$Clusters
VanClusters <- VanillaResults$Clusters 
BenClusters <- BenzResults$Clusters    
HydClusters <- HydroxResults$Clusters 

# Make the cluster names unambiguous
GluClusters$Cluster <- paste0("glucose-CazC", GlucoseResults$Clusters$Cluster)
FucClusters$Cluster  <- paste0("fucose-CazC", FucoseResults$Clusters$Cluster)
XylClusters$Cluster  <- paste0("xylose-CazC", XyloseResults$Clusters$Cluster)
GalClusters$Cluster  <- paste0("galactose-CazC", GalactoseResults$Clusters$Cluster)
VanClusters$Cluster  <- paste0("vanillic acid-CazC",VanillaResults$Clusters$Cluster)
BenClusters$Cluster    <- paste0("benzoic acid-CazC", BenzResults$Clusters$Cluster)
HydClusters$Cluster  <- paste0("hydroxybenzoic acid-CazC", HydroxResults$Clusters$Cluster)

# Stick em all together
AllClusters <- rbindlist(list(GluClusters, FucClusters, XylClusters
               , GalClusters, VanClusters, BenClusters
               , HydClusters))

# Save this as a file for functional analysis later 
fwrite(AllClusters, "CAZy/05_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")

#### Do certain CAZymes occur in clusters for multiple metabolites?
dev.off()
AllClustCAZymes <- rbindlist(list(
  Glucose = GlucoseResults$Clusters
, Fucose = FucoseResults$Clusters
, Xylose = XyloseResults$Clusters
, Galactose = GalactoseResults$Clusters
, Vanilla = VanillaResults$Clusters
, Hydroxybenz = HydroxResults$Clusters
, Benz = BenzResults$Clusters)
, idcol = TRUE
)

CAZymeCount <- table(AllClustCAZymes$CAZyme)
hist(CAZymeCount, col = 'grey50', breaks = 20)
RecurrentCAZymes <- CAZymeCount[CAZymeCount > 3]
RecurrentCAZymes

AllClustCAZymes[AllClustCAZymes$CAZyme %in% names(RecurrentCAZymes),]

names(RecurrentCAZymes)[-c(1,6,7)]

AllClustCAZymes[AllClustCAZymes$CAZyme %in% names(RecurrentCAZymes)[c(1,6,7)],]
AllClustCAZymes[AllClustCAZymes$CAZyme %in% names(RecurrentCAZymes)[-c(1,6,7)],]


