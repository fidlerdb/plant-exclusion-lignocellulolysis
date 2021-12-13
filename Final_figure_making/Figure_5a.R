########################################################################################
# The purpose of this script is to find out which metabolite-correlated CAZymes phyla  #
# have genomic resources in                                                           #
########################################################################################
#.rs.restartR()
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

library(data.table)
library(ggplot2)

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

alldata <- list(df_ba, df_fucose, df_galactose, df_glucose, df_hba, df_va, df_xylose)

#### Get metabolite-correlated CAZy richness for each phylum and metabolite ####

GetMetabRichness <- function(x){
  Richness = c(unlist(by(x$heatmap$heatmap_data, x$heatmap$heatmap_data$Phylum, FUN = function(x){
    sum(x$Abundance)
  })))
  Metabolite = x$metabolite
  
  Out = data.table(Phylum = names(Richness)
                   , Richness = Richness
                   , Metabolite = Metabolite)
  return(Out)
}

AllMetaboliteRichness <- lapply(alldata, GetMetabRichness)
AllMetaboliteRichness <- rbindlist(AllMetaboliteRichness)

metaborder <- names(sort(tapply(AllMetaboliteRichness$Richness
                                , AllMetaboliteRichness$Metabolite
                                , FUN = sum)
                         , decreasing = TRUE))
AllMetaboliteRichness$Metabolite <- factor(AllMetaboliteRichness$Metabolite, levels = metaborder)

#### Create the plot ####

rch_plot <- ggplot(AllMetaboliteRichness, aes(x = Metabolite, y = reorder(Phylum, Richness), fill = Richness)) +
  geom_tile() +
  ylab("Phylum") +
  scale_fill_gradient("Number of genes in\ncorrelated CAZy families"
                      , high = "black", low = "white"
                      , trans = 'log10'
                      , breaks = c(1,10,100,1000)
                      #, labels =  c(0,10,100,1000)
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)
        , panel.border = element_rect(colour = "black", fill=NA, size = 0.1))

rch_plot  

ggsave("Final_figure_making/To_Put_in_Inkscape/Figure_5a_Metab.pdf"
       , plot = rch_plot
       , width = 6
       , height = 7
       , units = "in")
