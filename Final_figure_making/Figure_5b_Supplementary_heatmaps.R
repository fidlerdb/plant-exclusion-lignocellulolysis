########################################################################################
# The purpose of this script is to find out which orders had the genes for the CAZymes #
# which were correlated with the abundance of a lignocellulose breakdown product       #
########################################################################################
# .rs.restartR()

library(data.table)


#### (1) Import, check, and set the data up ####
rm(list = ls()) # Clear the workspace

source("CAZy_2/Functions/prepare_Cluster_Heatmap.R")
source("CAZy_2/Functions/plot_Cluster_Heatmap_MetaboliteBlock.R")

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

#rm(df) # Stop it messing with the function

###### Create the plots ######
library(patchwork)
layout <- "
AAAA##
AAAA##
BBBB#D
CCCCEE
CCCCEE
CCCCEE
CCCCEE
"

#####
xyloseHeatmap <- plotClusterHeatmap(df_xylose
                                    , heatmap_colour = "#066C0D"
                                    , dendrogram_left_margin = -0.5
                                    , x_axis_size = 7
                                    , y_axis_size = 7)

xyloseHeatmap <- with(xyloseHeatmap, plt_dendr_top + plt_density + plt_hmap 
                      + lgnd_density + plt_dendr_side +
                        plot_layout(design = layout))

####

VAHeatmap <- plotClusterHeatmap(df_va
                                #, x_axis_size = 5
                                , dendrogram_left_margin = -0.5
                                , heatmap_colour = "#6C1706"
                                , x_axis_size = 7
                                , y_axis_size = 7)

VAHeatmap <- with(VAHeatmap, plt_dendr_top + plt_density + plt_hmap 
                  + lgnd_density + plt_dendr_side +
                    plot_layout(design = layout))

####


glucoseHeatmap <- plotClusterHeatmap(df_glucose
                                     , heatmap_colour = "#0D066C"
                                     , dendrogram_left_margin = -0.5
                                     , x_axis_size = 7
                                     , y_axis_size = 7)

glucoseHeatmap <- with(glucoseHeatmap, plt_dendr_top + plt_density + plt_hmap 
                       + lgnd_density + plt_dendr_side +
                         plot_layout(design = layout))

####

fucoseHeatmap <- plotClusterHeatmap(df_fucose
                                    , dendrogram_left_margin = -0.4
                                    , colour_labels = c(0, 10, 100, 750)
                                    , heatmap_colour = "#066C0D"
                                    , x_axis_size = 7
                                    , y_axis_size = 7)

fucoseHeatmap <- with(fucoseHeatmap, plt_dendr_top + plt_density + plt_hmap 
                      + lgnd_density + plt_dendr_side +
                        plot_layout(design = layout))

####


BAHeatmap <- plotClusterHeatmap(df_ba#, x_axis_size = 8
                                , dendrogram_left_margin = -0.5
                                , colour_labels = c(0, 5, 10, 20)
                                , heatmap_colour = "#6C1706"
                                , x_axis_size = 7
                                , y_axis_size = 7)

BAHeatmap <- with(BAHeatmap, plt_dendr_top + plt_density + plt_hmap 
                  + lgnd_density + plt_dendr_side +
                    plot_layout(design = layout))

####

HBAHeatmap <- plotClusterHeatmap(df_hba#, x_axis_size = 5
                                 , dendrogram_left_margin = -0.5
                                 , colour_labels = c(0, 10, 100, 300)
                                 , heatmap_colour = "#6C1706"
                                 , x_axis_size = 7
                                 , y_axis_size = 7)

HBAHeatmap <- with(HBAHeatmap, plt_dendr_top + plt_density + plt_hmap 
                   + lgnd_density + plt_dendr_side +
                     plot_layout(design = layout))

####

ADGHeatmap <- plotClusterHeatmap(df_galactose#, x_axis_size = 8
                                 , dendrogram_left_margin = -0.5
                                 , colour_labels = c(0, 10, 25)
                                 , heatmap_colour = "#066C0D"
                                 , x_axis_size = 7
                                 , y_axis_size = 7)

ADGHeatmap <- with(ADGHeatmap, plt_dendr_top + plt_density + plt_hmap 
                   + lgnd_density + plt_dendr_side +
                     plot_layout(design = layout))

#### Save the plots so that they are all legible ####
get_Plot_width <- function(x){
  min_cell_width <- 0.25
  if(length(unique(x$heatmap$heatmap_data$CAZyme)) * min_cell_width < 6){
    total_width <- 6
    return(total_width)
  } else {
    total_width <- length(unique(x$heatmap$heatmap_data$CAZyme)) * min_cell_width
    return(total_width)
  }
}

ggsave("Figures/CAZy_2/Heatmap_Clusters_Glucose.pdf", glucoseHeatmap
       , width = get_Plot_width(df_glucose)
       ,  height =  8, units = 'in')
ggsave("Figures/CAZy_2/Heatmap_Clusters_Xylose.pdf", xyloseHeatmap
       , width = get_Plot_width(df_xylose)
       ,  height =  8, units = 'in')
ggsave("Final_figure_making/To_Put_in_Inkscape/Figure_5b_Fucose_clusters.pdf", fucoseHeatmap
       , width = get_Plot_width(df_fucose)
       ,  height =  8, units = 'in')
ggsave("Figures/CAZy_2/Heatmap_Clusters_Galactose.pdf", ADGHeatmap
       , width = get_Plot_width(df_galactose)
       ,  height =  8, units = 'in')
ggsave("Figures/CAZy_2/Heatmap_Clusters_Vanillic_Acid.pdf", VAHeatmap
       , width = get_Plot_width(df_va)
       ,  height =  8, units = 'in')
ggsave("Figures/CAZy_2/Heatmap_Clusters_Benzoic_Acid.pdf", BAHeatmap
       , width = get_Plot_width(df_ba)
       ,  height =  8, units = 'in')
ggsave("Figures/CAZy_2/Heatmap_Clusters_Hydroxybenzoic_Acid.pdf", HBAHeatmap
       , width = get_Plot_width(df_hba)
       ,  height =  8, units = 'in')