########################################################################################
# The purpose of this script is to find out which metabolite-correlated CAZymes phyla  #
# have genomic resources in                                                           #
########################################################################################

library(data.table)
library(magrittr)
library(ggplot2)
library(ggdendro)
library(egg)
library(ggrepel)

#### (1) Import, check, and set the data up ####
rm(list = ls()) # Clear the workspace
# source("Functions/customHeatmapPlot.R")
#source("Functions/customHeatmapPlot-2.R")
source("Functions/prepareClusterHeatmap_2020-07-20.R")
source("Functions/plotClusterHeatmap_MetaboliteBlock.R")

# Read in the taxonomy data--normalized
df <- fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_dbCAN_Taxonomy_ReadAbundance.csv")

# Remove viral contigs
df <- df[df$SuperKingdom != "Viruses",]
# Remove contigs which have no phylum-level assignment
df <- df[!is.na(df$Phylum),]

# Read in the CAzyme data
cl <- fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Positively_Correlated_CAZymes-Metabolites.csv")

unique(cl$CAZymes)
# Sort out the weird ones
cl[cl$CAZymes == "GH1.3.2.1..",]$CAZymes <- "GH1"
cl[cl$CAZymes == "GH10.3.2.1.8",]$CAZymes <- "GH10_endo-1,4-β-xylanase"
cl[cl$CAZymes == "GH3.3.2.1.21",]$CAZymes <- "GH3_β-Glucosidase"

# Create useful vectors for building a matrix of CAZymes vs Orders/Phyla
AllCazymes <- unique(cl$CAZymes)
AllPhyla <- unique(df$Phylum)
AllMetabolites <- unique(cl$Metabolite)
unique(df$CAZyme)

# All metabolite-cazyme combinations
Met_Caz <- paste(rep(AllMetabolites, each = length(AllCazymes)), 
                 rep(AllCazymes, 7)
                 , sep = "-"
)

# Create an empty matrix the correct size
CRM <- matrix(NA
              , nrow = length(AllPhyla)
              , ncol = length(Met_Caz))
CRM <- data.frame(CRM)
names(CRM) <- Met_Caz
rownames(CRM) <- AllPhyla

head(CRM[1:10])

# Import metabolite activity data
ct <- fread("CAZy/05_CAZyme_Functional_information/OutputData/CAZy_Cluster_Activity_Summary.csv")

# Make the data ggplottable
ctl <- melt(ct, id = c("Cluster", "nFamilies"))

# For any metabolite, count the number of contigs with one of the CAZymes of interest
# Which phyla had genes related to the abundance of each metabolite?
# Which CAZymes occurred in similar phyla?

# Set up/load the data to visually answer these questions.
if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_xylose_2020-07-20.rds")){
  df_xylose <- prepareClusterHeatmap("xylose")
  saveRDS(df_xylose, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_xylose_2020-07-20.rds")
} else {df_xylose <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_xylose_2020-07-20.rds")}

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_vanillicacid_2020-07-20.rds")){
  df_va <- prepareClusterHeatmap("vanillic acid")
  saveRDS(df_va, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_vanillicacid_2020-07-20.rds")
} else {df_va <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_vanillicacid_2020-07-20.rds")}

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_glucose.rds")){
  df_glucose <- prepareClusterHeatmap("glucose")
  saveRDS(df_glucose, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_glucose_2020-07-20.rds")
} else {df_glucose <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_glucose_2020-07-20.rds")}

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_fucose_2020-07-20.rds")){
  df_fucose <- prepareClusterHeatmap("fucose")
  saveRDS(df_fucose, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_fucose_2020-07-20.rds")
} else {df_fucose <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_fucose_2020-07-20.rds")}

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_benzoicacid_2020-07-20.rds")){
  df_ba <- prepareClusterHeatmap("benzoic acid")
  saveRDS(df_ba, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_benzoicacid_2020-07-20.rds")
} else {df_ba <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_benzoicacid_2020-07-20.rds")}

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_hydroxybenzoicacid_2020-07-20.rds")){
  df_hba <- prepareClusterHeatmap("4-hydroxybenzoic acid")
  saveRDS(df_hba, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_hydroxybenzoicacid_2020-07-20.rds")
} else {df_hba <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_hydroxybenzoicacid_2020-07-20.rds")}

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_galactose_2020-07-20.rds")){
  df_galactose <- prepareClusterHeatmap("3,6-anhydro-D-galactose")
  saveRDS(df_galactose, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_galactose_2020-07-20.rds")
} else {df_galactose <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_galactose_2020-07-20.rds")}

rm(df) # Stop it messing with the function

# cc = CAZyme Clusters
cc <- fread("CAZy/05_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")
cc

str(cc[grep("galactose", Cluster),]$Cluster)
cc[grep("galactose", Cluster),]$Cluster <- paste0("galactose-",cc[grep("galactose", Cluster),]$CAZyme)

#### (2) Create basic function ####

# Create a function which finds the abundance of CAZymes belonging to an individual phylum, in each CAZyme cluster
GetClusterAbundances <- function(phylum, metabolite, inputdata){
  # Make life more simple
  df <- data.table(inputdata$heatmap$heatmap_data)
  
  # Subset the data
  Phylum_data <- df[Phylum == phylum]
  Total_Abundance <- sum(df[Phylum == phylum]$Abundance) # Get total CAZyme abundance for the phylum
  cc_metab <- cc[grep(metabolite, Cluster),] # Subset the clusters to only include the ones with the metabolite of choice
  
  # Find percentage abundance CAZy genes of that phylum belonging to each CAZyme cluster
  ClusterAbundances <- by(cc_metab, cc_metab$Cluster, function(x){
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

fwrite(CAZyProfiles, "CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Metabolite-Correlate-CAZyme-Profiles_Phyla.csv")

RowsToKeep <- apply(CAZyProfiles, 1, FUN = function(x){any(x[2:length(x)] > 20)})
CAZyProfiles[RowsToKeep,]

#Abundant Clusters
fwrite(CAZyProfiles[RowsToKeep,], "CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Metabolite-Correlate-CAZyme-Profiles_Phyla_20Percent.csv")

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