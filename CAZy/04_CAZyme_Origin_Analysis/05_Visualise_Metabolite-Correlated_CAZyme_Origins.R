########################################################################################
# The purpose of this script is to find out which orders had the genes for the CAZymes #
# which were correlated with the abundance of a lignocellulose breakdown product       #
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
source("Functions/prepareClusterHeatmap.R")
source("Functions/plotClusterHeatmap.R")

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

# For any metabolite, count the number of contigs with one of the CAZymes of interest
# Which phyla had genes related to the abundance of each metabolite?
# Which CAZymes occurred in similar phyla?

# Set up/load the data to visually answer these questions.
if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_xylose.rds")){
   df_xylose <- prepareClusterHeatmap("xylose")
   saveRDS(df_xylose, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_xylose.rds")
} else {df_xylose <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_xylose.rds")}


if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_vanillicacid.rds")){
   df_va <- prepareClusterHeatmap("vanillic acid")
   saveRDS(df_va, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_vanillicacid.rds")
} else {df_va <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_vanillicacid.rds")}

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_glucose.rds")){
df_glucose <- prepareClusterHeatmap("glucose")
saveRDS(df_glucose, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_glucose.rds")
} else {df_glucose <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_glucose.rds")}

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_fucose.rds")){
   df_fucose <- prepareClusterHeatmap("fucose")
   saveRDS(df_fucose, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_fucose.rds")
} else {df_fucose <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_fucose.rds")}

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_benzoicacid.rds")){
   df_ba <- prepareClusterHeatmap("benzoic acid")
   saveRDS(df_ba, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_benzoicacid.rds")
} else {df_ba <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_benzoicacid.rds")}

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_hydroxybenzoicacid.rds")){
   df_hba <- prepareClusterHeatmap("4-hydroxybenzoic acid")
   saveRDS(df_hba, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_hydroxybenzoicacid.rds")
} else {df_hba <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_hydroxybenzoicacid.rds")}

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_galactose.rds")){
   df_galactose <- prepareClusterHeatmap("3,6-anhydro-D-galactose")
   saveRDS(df_galactose, "CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_galactose.rds")
} else {df_galactose <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Results/heatmap_data_galactose.rds")}



###### Create the plots ######
xyloseHeatmap <- plotClusterHeatmap(df_xylose
                                    , heatmap_colour = "#066C0D"
                                    , dendrogram_left_margin = -0.5)
#df_va$metabolite <- "vanillic acid"
VAHeatmap <- plotClusterHeatmap(df_va
                                , x_axis_size = 5
                                , dendrogram_left_margin = -0.5
                                , heatmap_colour = "#6C1706")

df_glucose$metabolite <- "glucose"
glucoseHeatmap <- plotClusterHeatmap(df_glucose
                                    , heatmap_colour = "#0D066C"
                                    , dendrogram_left_margin = -0.5
                                    )
df_fucose$metabolite <- "fucose"
fucoseHeatmap <- plotClusterHeatmap(df_fucose
                                   , dendrogram_left_margin = -0.4
                                   , colour_labels = c(0, 10, 100, 750)
                                   , heatmap_colour = "#066C0D")

df_ba$metabolite <- "benzoic acid"
BAHeatmap <- plotClusterHeatmap(df_ba, x_axis_size = 8
                               , dendrogram_left_margin = -0.5
                               , colour_labels = c(0, 5, 10, 20)
                               , heatmap_colour = "#6C1706")

df_hba$metabolite <- "4-hydroxybenzoic acid"
HBAHeatmap <- plotClusterHeatmap(df_hba, x_axis_size = 5
                                , dendrogram_left_margin = -0.5
                                , colour_labels = c(0, 10, 100, 300)
                                , heatmap_colour = "#6C1706")

df_galactose$metabolite <- "3,6-anhydro-D-galactose"
ADGHeatmap <- plotClusterHeatmap(df_galactose, x_axis_size = 8
                                , dendrogram_left_margin = -0.5
                                , colour_labels = c(0, 10, 25)
                                , heatmap_colour = "#066C0D")
# Look at the plots
# Cellulose
#glucoseHeatmap

# Hemicellulose
#xyloseHeatmap
#fucoseHeatmap
#ADGHeatmap 

# Lignin
#VAHeatmap
#BAHeatmap 
#HBAHeatmap 

# Save the plots so they are all the same size
ggsave("Figures/CAZy_Taxonomy/Heatmaps/Heatmap_Clusters_Glucose.pdf", glucoseHeatmap
       , width = 4.91,  height =  5.35, units = 'in')
ggsave("Figures/CAZy_Taxonomy/Heatmaps/Heatmap_Clusters_Xylose.pdf", xyloseHeatmap
       , width = 4.91,  height =  5.35, units = 'in')
ggsave("Figures/CAZy_Taxonomy/Heatmaps/Heatmap_Clusters_Fucose.pdf", fucoseHeatmap
       , width = 4.91,  height =  5.35, units = 'in')
ggsave("Figures/CAZy_Taxonomy/Heatmaps/Heatmap_Clusters_Galactose.pdf", ADGHeatmap
       , width = 4.91,  height =  5.35, units = 'in')
ggsave("Figures/CAZy_Taxonomy/Heatmaps/Heatmap_Clusters_Vanillic_Acid.pdf", VAHeatmap
       , width = 4.91,  height =  5.35, units = 'in')
ggsave("Figures/CAZy_Taxonomy/Heatmaps/Heatmap_Clusters_Benzoic_Acid.pdf", BAHeatmap
       , width = 4.91,  height =  5.35, units = 'in')
ggsave("Figures/CAZy_Taxonomy/Heatmaps/Heatmap_Clusters_Hydroxybenzoic_Acid.pdf", HBAHeatmap
       , width = 4.91,  height =  5.35, units = 'in')

######################
# Explore the importance of each CAZyme cluster.

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
HBADensity <- setorder(get_n_genes(df_hba), nGenesMean)
BADensity <- setorder(get_n_genes(df_ba), nGenesMean)

# Cellulose results
GlucoseDensity
hist(GlucoseDensity$nGenesMean, breaks = 100)
abline(v = max(GlucoseDensity$nGenesMean)/10, lty = 2, col = 'red')

# Glucose CAZyme Clusters to investigate further
GlucoseDensity[GlucoseDensity$nGenesMean > 
                  max(GlucoseDensity$nGenesMean)/10]

# Hemicellulose results
FucoseDensity
XyloseDensity
GalactoseDensity

hist(FucoseDensity$nGenesMean, breaks = 100)
abline(v = max(FucoseDensity$nGenesMean)/10, lty = 2, col = 'red')

hist(XyloseDensity$nGenesMean, breaks = 100)
abline(v = max(XyloseDensity$nGenesMean)/10, lty = 2, col = 'red')

hist(GalactoseDensity$nGenesMean, breaks = 100)

# Hemicellulose CAZy clusters to investigate further
FucoseDensity[FucoseDensity$nGenesMean > 
                 max(FucoseDensity$nGenesMean)/10]
XyloseDensity[XyloseDensity$nGenesMean > 
                 max(XyloseDensity$nGenesMean)/10]
# Lignin results
VADensity
HBADensity
BADensity

hist(VADensity$nGenesMean, breaks = 100)
abline(v = max(VADensity$nGenesMean)/10, lty = 2, col = 'red')
hist(HBADensity$nGenesMean, breaks = 100)
abline(v = max(HBADensity$nGenesMean)/10, lty = 2, col = 'red')
hist(BADensity$nGenesMean, breaks = 100)
abline(v = max(BADensity$nGenesMean)/10, lty = 2, col = 'red')

# Lignin CAZy clusters to investigate further
VADensity[VADensity$nGenesMean > 
                 max(VADensity$nGenesMean)/10]
HBADensity[HBADensity$nGenesMean > 
                 max(HBADensity$nGenesMean)/10]
df_hba$heatmap$CAZyme_pos_table[df_hba$heatmap$CAZyme_pos_table$Cluster == 6,]
BADensity[BADensity$nGenesMean > 
                 max(BADensity$nGenesMean)/10]

