########################################################################################
# The purpose of this script is to find out which orders had the genes for the CAZymes #
# which were correlated with the abundance of a lignocellulose breakdown product       #
########################################################################################
#.rs.restartR()

library(data.table)
library(doBy)

rm(list = ls()) # Clear the workspace

source("CAZy_2/Functions/prepare_Cluster_Heatmap.R")
source("Functions/plotClusterHeatmap.R")

#### (1) Read in and clean the taxonomy and CAZyme data--normalized ####

df2 <- fread("CAZy_2/outputData/Taxonomy_Data_CAZymes.csv")

head(df2[,1:24]) # This table needs summarising by phylum
# Keep only characterised contigs
df2 <- df2[df2$SuperKingdom != "Viruses"|df2$SuperKingdom != "",]
df2 <- df2[df2$Phylum != "",]

# Subset this data so we can obtain a CAZy family gene richness matrix
CAZymes <- names(df2)[24:length(df2)]
cf <- names(df2)[c(4,24:length(df2))]
df3 <- df2[, ..cf]
names(df3)[2:length(df3)] <- paste0("a_", 1:(length(df3)-1))

# Count the number of genes
CAZymeRichnessData <- summaryBy(. ~ Phylum, data = df3, FUN = sum, keep.names = TRUE)
names(CAZymeRichnessData)[2:length(CAZymeRichnessData)] <- CAZymes

# Look at the full phylum-level CAZy richness data. pretty cool.
CAZymeRichnessData

#### (2) Read in and clean the metabolite-correlated CAZyme data ####

mcc <- fread("CAZy_2/outputData/Positively_Correlated_CAZymes-Metabolites.csv")

CAZymes[CAZymes %in% unique(mcc$CAZymes)] # It's missing the weird ones for sure, and anythign with a |

# Replace . with | between CAZymes
# () around [A-Z] means that gsub stores the [A-Z] character as \\1
mcc$CAZymes <- gsub("\\.([A-Z])", "\\|\\1", mcc$CAZymes)

# Remove any remaining "." characters between CAZymes
mcc$CAZymes <- gsub("([A-Z].+?)\\.", "\\1\\|", mcc$CAZymes)

unique(mcc$CAZymes)

# Do they match up now?
CAZymes[CAZymes %in% unique(mcc$CAZymes)] # multi-domain CAZy mccassified enzyes sorted. are the EC numbers matching?

# Testing... Testing...
CAZymes[grep("\\|[0-9]",CAZymes)]
unique(mcc$CAZymes)[grep("\\|[0-9]",unique(mcc$CAZymes))]
"GH3|3.2.1.21" == "GH3|3.2.1.21"

# 3.2.1.21 β-glucosidase
# 3.2.1.8 endo-1,4-β-xylanase

# Yes.
# Good.

# Now subset the matrix to only include CAZymes which were metabolite-correlated
cf <- c("Phylum", CAZymes[CAZymes %in% unique(mcc$CAZymes)])

Corr_CAZyme_Richness <- CAZymeRichnessData[, ..cf]

#### (3) Create all of the subset variables ####

mcc_glu <- c("Phylum", mcc[Metabolite == "glucose",]$CAZymes)
mcc_fuc <- c("Phylum", mcc[Metabolite == "fucose",]$CAZymes)
mcc_xyl <- c("Phylum", mcc[Metabolite == "xylose",]$CAZymes)
mcc_gal <- c("Phylum", mcc[Metabolite == "3,6-anhydro-D-galactose",]$CAZymes)
mcc_hyd <- c("Phylum", mcc[Metabolite == "4-hydroxybenzoic acid",]$CAZymes)
mcc_van <- c("Phylum", mcc[Metabolite == "vanillic acid" ,]$CAZymes)
mcc_ben <- c("Phylum", mcc[Metabolite == "benzoic acid",]$CAZymes)

mcc_glu <- mcc_glu[mcc_glu %in% names(Corr_CAZyme_Richness)]
mcc_fuc <- mcc_fuc[mcc_fuc %in% names(Corr_CAZyme_Richness)]
mcc_xyl <- mcc_xyl[mcc_xyl %in% names(Corr_CAZyme_Richness)]
mcc_gal <- mcc_gal[mcc_gal %in% names(Corr_CAZyme_Richness)]
mcc_hyd <- mcc_hyd[mcc_hyd %in% names(Corr_CAZyme_Richness)]
mcc_van <- mcc_van[mcc_van %in% names(Corr_CAZyme_Richness)]
mcc_ben <- mcc_ben[mcc_ben %in% names(Corr_CAZyme_Richness)]

#### (4) Create heatmap data (clusters are the most important--adding functional info after) ####

# Glucose
if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_glucose.rds")){
  df_glucose <- prepare_Cluster_Heatmap(inputData = Corr_CAZyme_Richness[, ..mcc_glu]
                                        , y_var = "Phylum"
                                        , dist_method = Chi_dist
                                        , metabolite = "glucose"
                                        , simprof_expected = 10000
                                        , simprof_simulated = 10000)
  saveRDS(df_glucose, "CAZy_2/outputData/heatmap_data/heatmap_data_glucose.rds")
} else {df_glucose <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_glucose.rds")}


# Fucose
if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_fucose.rds")){
  df_fucose <- prepare_Cluster_Heatmap(inputData = Corr_CAZyme_Richness[, ..mcc_fuc]
                                       , y_var = "Phylum"
                                       , dist_method = Chi_dist
                                       , metabolite = "fucose"
                                       , simprof_expected = 10000
                                       , simprof_simulated = 10000)
  saveRDS(df_fucose, "CAZy_2/outputData/heatmap_data/heatmap_data_fucose.rds")
} else {df_fucose <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_fucose.rds")}


# Xylose
if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_xylose.rds")){
  df_xylose <- prepare_Cluster_Heatmap(inputData = Corr_CAZyme_Richness[, ..mcc_xyl]
                                       , y_var = "Phylum"
                                       , dist_method = Chi_dist
                                       , metabolite = "xylose"
                                       , simprof_expected = 10000
                                       , simprof_simulated = 10000)
  saveRDS(df_xylose, "CAZy_2/outputData/heatmap_data/heatmap_data_xylose.rds")
} else {df_xylose <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_xylose.rds")}


# Galactose
if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_galactose.rds")){
  df_galactose <- prepare_Cluster_Heatmap(inputData = Corr_CAZyme_Richness[, ..mcc_gal]
                                       , y_var = "Phylum"
                                       , dist_method = Chi_dist
                                       , metabolite = "3,6-anhydro-D-galactose"
                                       , simprof_expected = 10000
                                       , simprof_simulated = 10000)
  saveRDS(df_galactose, "CAZy_2/outputData/heatmap_data/heatmap_data_galactose.rds")
} else {df_galactose <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_galactose.rds")}


# Hydroxybenzoic acid
if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_hydroxybenzoic_acid.rds")){
  df_hyd <- prepare_Cluster_Heatmap(inputData = Corr_CAZyme_Richness[, ..mcc_hyd]
                                   , y_var = "Phylum"
                                   , dist_method = Chi_dist
                                   , metabolite = "4-hydroxybenzoic acid"
                                   , simprof_expected = 10000
                                   , simprof_simulated = 10000)
  saveRDS(df_hyd, "CAZy_2/outputData/heatmap_data/heatmap_data_hydroxybenzoic_acid.rds")
} else {df_hyd <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_hydroxybenzoic_acid.rds")}


# Vanillic acid
if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_vanillic_acid.rds")){
  df_va <- prepare_Cluster_Heatmap(inputData = Corr_CAZyme_Richness[, ..mcc_van]
                                   , y_var = "Phylum"
                                   , dist_method = Chi_dist
                                   , metabolite = "vanillic acid"
                                   , simprof_expected = 10000
                                   , simprof_simulated = 10000)
  saveRDS(df_va, "CAZy_2/outputData/heatmap_data/heatmap_data_vanillic_acid.rds")
} else {df_va <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_vanillic_acid.rds")}

# Benzoic acid
if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_benzoic_acid.rds")){
  df_ben <- prepare_Cluster_Heatmap(inputData = Corr_CAZyme_Richness[, ..mcc_ben]
                                    , y_var = "Phylum"
                                    , dist_method = Chi_dist
                                    , metabolite = "benzoic acid"
                                    , simprof_expected = 10000
                                    , simprof_simulated = 10000)
  saveRDS(df_ben, "CAZy_2/outputData/heatmap_data/heatmap_data_benzoic_acid.rds")
} else {df_ben <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_benzoic_acid.rds")}

#### (5) Plot these to check they have all worked ####

plotClusterHeatmap(df_glucose, heatmap_colour = "#0D066C"
                   , dendrogram_left_margin = -0.5)

plotClusterHeatmap(df_fucose, heatmap_colour = "#066C0D"
                   , dendrogram_left_margin = -0.5)

plotClusterHeatmap(df_xylose, heatmap_colour = "#066C0D"
                   , dendrogram_left_margin = -0.5)

plotClusterHeatmap(df_galactose, heatmap_colour = "#066C0D"
                   , dendrogram_left_margin = -0.5)

plotClusterHeatmap(df_va, heatmap_colour = "#6C1706"
                   , dendrogram_left_margin = -0.5)

plotClusterHeatmap(df_hyd, heatmap_colour = "#6C1706"
                   , dendrogram_left_margin = -0.5)

plotClusterHeatmap(df_ben, heatmap_colour = "#6C1706"
                   , dendrogram_left_margin = -0.5)


#### (6) Explore the importance of each CAZyme cluster ####

# Create a function to summarise the CAZyme clusters' gene abundance

get_n_genes <- function(data_frame){

  # (1) Find how many clusters there are and create an empty list for assignment later
  nclusters <- length(na.omit(unique(data_frame$heatmap$heatmap_data$Cluster.y)))

  clusterResults <- list(1:nclusters)

  for(i in 1:nclusters){

    # (2) Find which CAZymes are in each cluster
    CazCluster <- with(data_frame$heatmap$heatmap_data,
                       na.omit(unique(CAZyme[Cluster.y == i])))

    # (3) Summarise the total number of genes from each CAZyme
    #     within the cluster
    ClusterTotals <- unlist(with(data_frame$heatmap$heatmap_data,
                                 lapply(seq_along(CazCluster), function(i){
                                   sum(Abundance[CAZyme %in% CazCluster[i]]-1)
                                 })))

    # Output the results together
    clusterResults[[i]] <- (data.frame(Cluster = i

                                       # Total genes
                                       , nGenes = with(data_frame$heatmap$heatmap_data,
                                                       sum(Abundance[CAZyme %in% CazCluster]-1))
                                       # Mean genes per CAZy family
                                       , nGenesMean = mean(ClusterTotals)
                                       # SD of genes per cazy family
                                       , nGenesSD = sd(ClusterTotals)
    )
    )
  } # End of loop

  # Output the results
  return(rbindlist(clusterResults))
}
 
# Find out which were the most important clusters

GlucoseDensity <- setorder(get_n_genes(df_glucose), nGenesMean)
FucoseDensity <- setorder(get_n_genes(df_fucose), nGenesMean)
XyloseDensity <- setorder(get_n_genes(df_xylose), nGenesMean)
GalactoseDensity <- setorder(get_n_genes(df_galactose), nGenesMean)
VADensity <- setorder(get_n_genes(df_va), nGenesMean)
HBADensity <- setorder(get_n_genes(df_hyd), nGenesMean)
BADensity <- setorder(get_n_genes(df_ben), nGenesMean)

GlucoseDensity 
FucoseDensity 
XyloseDensity 
GalactoseDensity
VADensity
HBADensity
BADensity

#### (7) Which CAZymes were in each cluster? ####

orderbyCluster <- function(inputData){data.table(inputData$heatmap$CAZyme_pos_table)[order(Cluster), c(1,4)]}

# What CAZymes were associated with the metabolites?

AllCAZymes_Clusters <- rbindlist(list(
  data.table(Metabolite = "glucose", orderbyCluster(df_glucose))
  , data.table(Metabolite = "fucose", orderbyCluster(df_fucose))
  , data.table(Metabolite = "xylose", orderbyCluster(df_xylose))
  , data.table(Metabolite = "3,6-anhydro-D-galactose", orderbyCluster(df_galactose))
  , data.table(Metabolite = "vanillic acid", orderbyCluster(df_va))
  , data.table(Metabolite = "4-hydroxybenzoic acid", orderbyCluster(df_hyd))
  , data.table(Metabolite = "benzoic acid", orderbyCluster(df_ben))
))

AllCAZymes_Clusters

# Add a unique cluster name
AllCAZymes_Clusters$ClusterName <- paste0(AllCAZymes_Clusters$Metabolite
                                          , "-CazC"
                                          , AllCAZymes_Clusters$Cluster)

# Write this to a file
fwrite(file = "CAZy_2/05-10_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv"
       , x = AllCAZymes_Clusters)



