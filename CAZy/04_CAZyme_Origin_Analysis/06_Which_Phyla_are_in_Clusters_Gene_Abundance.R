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

source("Functions/Import_and_Clean_Metabolite-CAZy_Taxonomy.R")

#### (2) Find statistically distinct clusters for all metabolites ####

source("Functions/PhylumClustering.R")

# Cellulose
GlucoseResults <- PhylumClustering("glucose")

# Hemicellulose
FucoseResults <- PhylumClustering("fucose")
XyloseResults <- PhylumClustering("xylose")
GalactoseResults <- PhylumClustering("3,6-anhydro-D-galactose")

# Lignin
VanillaResults <- PhylumClustering("vanillic acid")
BenzResults <- PhylumClustering("benzoic acid")
HydroxResults <- PhylumClustering("4-hydroxybenzoic acid")

#### (3) Look at the results ####

# trbl
par(mar = c(12, 0,0,0))


### Think I should have an average number of CAZy genes per phylum cutoff for 
### describing clusters....

hist(log10(c(
  unlist(GlucoseResults$Cluster_Mean_Genes)
  ,unlist(FucoseResults$Cluster_Mean_Genes)
  ,unlist(XyloseResults$Cluster_Mean_Genes)
  ,unlist(GalactoseResults$Cluster_Mean_Genes)
  ,unlist(VanillaResults$Cluster_Mean_Genes)
  ,unlist(BenzResults$Cluster_Mean_Genes)
  ,unlist(HydroxResults$Cluster_Mean_Genes)
  ))
  , breaks = 100)
log10(50)

## Cellulose
simprof.plot(GlucoseResults$simprofResults)
GlucoseResults$simprofResults$significantclusters
GlucoseResults$Cluster_Mean_Genes
GlucoseResults$Total_Genes
GlucoseResults$Gene_Richness

## Hemicellulose
# Fucose
simprof.plot(FucoseResults$simprofResults)
FucoseResults$simprofResults$significantclusters
FucoseResults$Cluster_Mean_Genes
FucoseResults$Total_Genes
FucoseResults$Gene_Richness

# Xylose
simprof.plot(XyloseResults$simprofResults)
XyloseResults$simprofResults$significantclusters
XyloseResults$Cluster_Mean_Genes
XyloseResults$Total_Genes
XyloseResults$Gene_Richness

# 3,6-Anhydro-D-Galactose
simprof.plot(GalactoseResults$simprofResults)
GalactoseResults$simprofResults$significantclusters
GalactoseResults$Cluster_Mean_Genes
GalactoseResults$Total_Genes
GalactoseResults$Gene_Richness

## Lignin
# Vanillic Acid
simprof.plot(VanillaResults$simprofResults)
VanillaResults$simprofResults$significantclusters
VanillaResults$Cluster_Mean_Genes
VanillaResults$Total_Genes
VanillaResults$Gene_Richness

# Benzoic Acid
simprof.plot(BenzResults$simprofResults)
BenzResults$simprofResults$significantclusters
BenzResults$Cluster_Mean_Genes
BenzResults$Total_Genes
BenzResults$Gene_Richness

# 4-Hydroxybenzoic Acid
simprof.plot(HydroxResults$simprofResults)
HydroxResults$simprofResults$significantclusters
HydroxResults$Cluster_Mean_Genes
HydroxResults$Total_Genes
HydroxResults$Gene_Richness

#### Save the phylum clusters for all of the metabolites ####

PhylumClusters <- list(glucose = GlucoseResults$Clusters,
                  fucose = FucoseResults$Clusters,
                  xylose = XyloseResults$Clusters,
                  `3,6-anhydro-D-galactose` = GalactoseResults$Clusters,
                  `vanillic acid` = VanillaResults$Clusters,
                  `benzoic acid` = BenzResults$Clusters,
                  `hydroxybenzoic acid` = HydroxResults$Clusters)

# Merge all of the cluster results for the phyla 
PhylumClusters <- Reduce((function() {counter = 0
function(x, y) {
  counter <<- counter + 1
  d = merge(x, y
            , all = TRUE
            , by = 'Phylum')
  setnames(d, c(head(names(d), -1), paste0('y.', counter)))
  }})(), PhylumClusters)

names(PhylumClusters) <- c("Phylum"
,"glucose"
,"fucose "
,"xylose "
,"`3,6-anhydro-D-galactose`"
,"`vanillic acid` "
,"`benzoic acid`"
,"`4-hydroxybenzoic acid`")

# Save this for deeper analysis soon
write.csv(PhylumClusters
          , file = "CAZy/04_CAZyme_Origin_Analysis/Results/Clusters_of_Phyla.csv"
          , row.names = FALSE)

#### Compare how many genes all phyla had ####

library(tibble)
library(purrr)
library(magrittr)
library(dplyr)

TotalGeneList <- list(
  rownames_to_column(data.frame(GlucoseResults$Total_Genes))
     , rownames_to_column(data.frame(FucoseResults$Total_Genes))
     , rownames_to_column(data.frame(XyloseResults$Total_Genes))
     , rownames_to_column(data.frame(GalactoseResults$Total_Genes))
     , rownames_to_column(data.frame(VanillaResults$Total_Genes))
     , rownames_to_column(data.frame(BenzResults$Total_Genes))
     , rownames_to_column(data.frame(HydroxResults$Total_Genes))
  )


TotalGeneList <- TotalGeneList %>% purrr::reduce(left_join, by = "rowname")
TotalGeneList[is.na(TotalGeneList)] <- 0
names(TotalGeneList) <- c("rowname"
                          ,"Glucose"
                          , "Fucose"
                          , "Xylose"
                          , "3,6-Anhydro-D-Galactose"
                          , "Vanillic Acid"
                          , "Benzoic Acid"
                          , "4-Hydroxybenzoic Acid"
                          )

# Sorted Order for y axis
PhylumOrder <- TotalGeneList[order(rowSums(TotalGeneList[,-1])),]$rowname

# Ggplot2 long lameness
TotalGeneData <- melt(TotalGeneList)
TotalGeneData$rowname <- factor(TotalGeneData$rowname, levels = PhylumOrder)
colour_labels <- c(0,10,100,1000)
#colour_breaks <- c(1,11,101,1001)

# Make the plot
p <- ggplot(TotalGeneData, aes(x = variable
                          , y = rowname
                          , fill = value+1)) + 
  geom_tile() +
  scale_fill_gradient(expression("Number of Genes")
                      , high = "darkblue", low = "white"
                      , trans = 'log10'
                      , breaks = colour_labels + 1
                      , labels = colour_labels
                      ) +
  theme(axis.text.x = element_text(angle = 90
                                   , vjust = 0.5
                                   , hjust = 1)
  ) +
  xlab("Metabolite") +
  ylab("Phylum")

p

# Save it as a pdf
ggsave("Figures/CAZy_Taxonomy/Heatmaps/All_Metabolites.pdf"
       , plot = p, width = 4.53, height = 4.94, units = "in")

# What is the total number of genes which were correlated

GeneCount <- data.frame(Phylum = TotalGeneList[,1]
                        ,  Genes = rowSums(TotalGeneList[,-1]))

GeneCount[order(GeneCount$Genes),]
