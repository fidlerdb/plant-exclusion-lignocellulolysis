########################################################################################
# The purpose of this script is to find out which metabolite-correlated CAZymes phyla  #
# have genomic resources in                                                           #
########################################################################################
#.rs.restartR()
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

library(data.table)


#### (1) Import, check, and set the data up ####
rm(list = ls()) # Clear the workspace

source("CAZy_2/Functions/prepare_Cluster_Heatmap.R")
source("CAZy_2/Functions/plot_Cluster_Heatmap_MetaboliteBlock.R")

# Import metabolite activity data
# ct <- fread("CAZy/05_CAZyme_Functional_information/OutputData/CAZy_Cluster_Activity_Summary.csv")
ct <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData/CAZy_Cluster_Activity_Summary.csv")

# Make the data ggplottable
ctl <- melt(ct, id = c("Cluster", "nFamilies", "Substrate"))

# For any metabolite, count the number of contigs with one of the CAZymes of interest
# Which phyla had genes related to the abundance of each metabolite?
# Which CAZymes occurred in similar phyla?

# Set up/load the data to visually answer these questions.
if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_xylose.rds")){
  df_xylose <- prepareClusterHeatmap("xylose")
  saveRDS(df_xylose, "CAZy_2/outputData/heatmap_data/heatmap_data_xylose.rds")
} else {df_xylose <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_xylose.rds")}

if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_vanillic_acid.rds")){
  df_va <- prepareClusterHeatmap("vanillic acid")
  saveRDS(df_va, "CAZy_2/outputData/heatmap_data/heatmap_data_vanillic_acid.rds")
} else {df_va <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_vanillic_acid.rds")}

if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_glucose.rds")){
  df_glucose <- prepareClusterHeatmap("glucose")
  saveRDS(df_glucose, "CAZy_2/outputData/heatmap_data/heatmap_data_glucose.rds")
} else {df_glucose <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_glucose.rds")}

if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_fucose.rds")){
  df_fucose <- prepareClusterHeatmap("fucose")
  saveRDS(df_fucose, "CAZy_2/outputData/heatmap_data/heatmap_data_fucose.rds")
} else {df_fucose <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_fucose.rds")}

if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_benzoic_acid.rds")){
  df_ba <- prepareClusterHeatmap("benzoic acid")
  saveRDS(df_ba, "CAZy_2/outputData/heatmap_data/heatmap_data_benzoic_acid.rds")
} else {df_ba <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_benzoic_acid.rds")}

if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_hydroxybenzoic_acid.rds")){
  df_hba <- prepareClusterHeatmap("4-hydroxybenzoic acid")
  saveRDS(df_hba, "CAZy_2/outputData/heatmap_data/heatmap_data_hydroxybenzoic_acid.rds")
} else {df_hba <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_hydroxybenzoic_acid.rds")}

if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_galactose.rds")){
  df_galactose <- prepareClusterHeatmap("3,6-anhydro-D-galactose")
  saveRDS(df_galactose, "CAZy_2/outputData/heatmap_data/heatmap_data_galactose.rds")
} else {df_galactose <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_galactose.rds")}

# cc = CAZyme Clusters
#cc <- fread("CAZy/05_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")
cc <- fread("CAZy_2/05-10_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")

cc

str(cc[grep("galactose", Cluster),]$Cluster)
cc[grep("galactose", Cluster),]$Cluster <- paste0("galactose-",cc[grep("galactose", Cluster),]$CAZyme)

cc$ClusterName
#### (2) Create basic function ####

####

# Make life more simple
# df <- data.table(df_xylose$heatmap$heatmap_data)
# phylum = "Actinobacteria"
# metabolite = "xylose"
# # Subset the data
# Phylum_data <- df[Phylum == phylum]
# Total_Abundance <- sum(df[Phylum == phylum]$Abundance) # Get total CAZyme abundance for the phylum
#cc_metab <- cc[grep(metabolite, ClusterName),] # Subset the clusters to only include the ones with the metabolite of choice

# # Find percentage abundance CAZy genes of that phylum belonging to each CAZyme cluster
# ClusterAbundances <- by(cc_metab, cc_metab$Cluster, function(x){
#   (sum(Phylum_data[CAZyme %in% x$CAZyme,]$Abundance) / Total_Abundance)*100
# })
# 
# # Format this nicely
# ClusterAbundances <- data.table(Cluster = names(ClusterAbundances)
#                                 , Percentage = c(ClusterAbundances)
# )
# names(ClusterAbundances)[2] <- phylum

####



# Create a function which finds the abundance of CAZymes belonging to an individual phylum, in each CAZyme cluster
GetClusterAbundances <- function(phylum, metabolite, inputdata){
  # Make life more simple
  df <- data.table(inputdata$heatmap$heatmap_data)
  
  # Subset the data
  Phylum_data <- df[Phylum == phylum]
  Total_Abundance <- sum(df[Phylum == phylum]$Abundance) # Get total CAZyme abundance for the phylum
  cc_metab <- cc[grep(metabolite, ClusterName),] # Subset the clusters to only include the ones with the metabolite of choice
  
  # Find percentage abundance CAZy genes of that phylum belonging to each CAZyme cluster
  ClusterAbundances <- by(cc_metab, cc_metab$ClusterName, function(x){
    (sum(Phylum_data[CAZyme %in% x$CAZyme,]$Abundance) / Total_Abundance)*100
  })
  
  # Format this nicely
  ClusterAbundances <- data.table(Cluster = names(ClusterAbundances)
                                  , Percentage = c(ClusterAbundances)
  )
  names(ClusterAbundances)[2] <- phylum
  return(ClusterAbundances)
}

# Test it out
#GetClusterAbundances(phylum = "Actinobacteria", "glucose", df_glucose)

#### (3) create more complex function #### 

# Now create a function to create a table showing the % of metabolite-correlated CAZymes from a phylum
# which belonged to a particular CAZyme cluster

# List of phyla
Phyla <- c("Proteobacteria",	"Actinobacteria",	"Acidobacteria"
           , "Bacteroidetes",	"Planctomycetes",	"Firmicutes"	
           , "Cyanobacteria",	"Gemmatimonadetes",	"Verrucomicrobia"	
           , "Euryarchaeota",	"Ascomycota"
)

GetMetaboliteResults <- function(inputdata, metabolite){
  
  # Get the cluster abundances for each phylum for the metabolite in question
  metabResults <- lapply(Phyla, function(x){
    GetClusterAbundances(phylum = x, metabolite, inputdata)}
  )
  
  NROW <- nrow(metabResults[[1]])
  
  NAMES <- unlist(lapply(metabResults, function(x){names(x)[2]}))
  metabResults <- cbind(metabResults[[1]][,1], matrix(unlist(lapply(metabResults, `[[`, 2)), nrow = NROW))
  names(metabResults)[2:length(metabResults)] <- NAMES
  
  return(metabResults)
}

#### (4) Get results ####
CAZyProfiles <- rbindlist(list(GetMetaboliteResults(df_glucose, "glucose")
               , GetMetaboliteResults(df_fucose, "fucose")
               , GetMetaboliteResults(df_xylose, "xylose")
               , GetMetaboliteResults(df_va, "vanillic acid")
               , GetMetaboliteResults(df_hba, "hydroxybenzoic acid")
               , GetMetaboliteResults(df_ba, "benzoic acid")
               , GetMetaboliteResults(df_galactose, "galactose")
))

CAZyProfiles
heatmap(as.matrix(CAZyProfiles[,2:length(CAZyProfiles)]))

fwrite(CAZyProfiles, "CAZy_2/05-10_CAZyme_Functional_information/OutputData/Metabolite-Correlate-CAZyme-Profiles_Phyla.csv")

#RowsToKeep <- apply(CAZyProfiles, 1, FUN = function(x){any(x[2:length(x)] > 20)})
#CAZyProfiles[RowsToKeep,]

#Abundant Clusters
#fwrite(CAZyProfiles[RowsToKeep,], "CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Metabolite-Correlate-CAZyme-Profiles_Phyla_20Percent.csv")

#### Which families are in which clusters? ####

cc <- data.table(cc)
cc[Cluster == "glucose-CazC2"]
cc[Cluster == "fucose-CazC1"]
cc[Cluster == "vanillic acid-CazC4"]
cc[Cluster == "vanillic acid-CazC2"]
cc[Cluster == "vanillic acid-CazC1"]
cc[Cluster == "hydroxybenzoic acid-CazC1"]
cc[Cluster == "hydroxybenzoic acid-CazC5"]
cc[Cluster == "hydroxybenzoic acid-CazC7"]
cc[Cluster == "benzoic acid-CazC2"]


##### (5) Richness of genes per phylum ####

GetClusterRichness <- function(phylum, metabolite, inputdata){
  # Make life more simple
  df <- data.table(inputdata$heatmap$heatmap_data)
  
  # Subset the data
  Phylum_data <- df[Phylum == phylum]
  Total_Abundance <- sum(df[Phylum == phylum]$Abundance) # Get total CAZyme abundance for the phylum
  
  # Format this nicely
  ClusterAbundances <- data.table(Phylum = phylum, Richness = Total_Abundance, Metabolite = metabolite)
  #names(ClusterAbundances)[2] <- phylum
  return(ClusterAbundances)
}

GetClusterRichness("Actinobacteria", "xylose", df_xylose)





GetMetaboliteRichness <- function(inputdata, metabolite){
  
  # Get the cluster abundances for each phylum for the metabolite in question
  metabResults <- lapply(Phyla, function(x){
    GetClusterRichness(phylum = x, metabolite, inputdata)}
  )
  
  metabResults <- rbindlist(metabResults)
  #NROW <- nrow(metabResults[[1]])
  
  #NAMES <- unlist(lapply(metabResults, function(x){names(x)[2]}))
  #metabResults <- cbind(metabResults[[1]][,1], matrix(unlist(lapply(metabResults, `[[`, 2)), nrow = NROW))
  #names(metabResults)[2:length(metabResults)] <- NAMES
  
  return(metabResults)
}

CAZyRichness <- rbindlist(list(GetMetaboliteRichness(df_glucose, "glucose")
                               , GetMetaboliteRichness(df_fucose, "fucose")
                               , GetMetaboliteRichness(df_xylose, "xylose")
                               , GetMetaboliteRichness(df_va, "vanillic acid")
                               , GetMetaboliteRichness(df_hba, "hydroxybenzoic acid")
                               , GetMetaboliteRichness(df_ba, "benzoic acid")
                               , GetMetaboliteRichness(df_galactose, "galactose")
))

CAZyRichness[order(CAZyRichness$Richness),]




