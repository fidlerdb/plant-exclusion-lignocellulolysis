rm(list = ls())

library(data.table)
library(clustsig)
library(ade4)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.14")
#BiocManager::install("Heatplus")
source("Functions/ClusterHeatmap_fromMatrix.R") # Yes, this is the version I want to use
library(cowplot)
library(vegan)
library(compiler)

#### Read in and clean data ####
df <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData/Metabolite-Correlate-CAZyme-Profiles_Phyla.csv")

head(df)
unique(df$Cluster) # Check all of the clusters are there
# remove fully zero clusters....
df <- df[rowSums(df[,2:length(df)]) != 0,]
par(mar = c(10.1, 4.1, 4.1, 2.1))

#### Prelim clustering analysis ####
clust_mat <- data.frame(t(df[,2:length(df)]))
names(clust_mat) <- df$Cluster
clust_dist <- dist.prop(clust_mat, method = 1)
clust_dist_clust <- hclust(clust_dist, method = "ward.D2")

plot(clust_dist_clust)
heatmap(as.matrix(clust_mat))
heatmap(log10(as.matrix(clust_mat+1)))

# Significance testing
Manly <- function(x){dist.prop(data.frame(x), method = 1)}
Manly <- cmpfun(Manly) # Make this a little faster

# Takes a while to run and we're doing this analysis with a heatmap in a moment
#clust_sig <- simprof(clust_mat, method.cluster = "ward.D2", method.distance = Manly)
#simprof.plot(clust_sig)

#### NMDS analysis ####

# Plot this as an NMDS plot
NMDS1 <- metaMDS(clust_mat, distfun = Manly)
ordiplot(NMDS1, type = "text")

coldf <- data.table(Phylum = rownames(NMDS1$points)
                    , Group = c(1,1,2,2,3,4,3,5,6,6,7))

plot(NMDS1, type = "n", xlim = c(-0.4,0.4))
text(NMDS1$species
     , col = scales::alpha("grey50", 0.5)
     , labels = rownames(NMDS1$species)
)
text(NMDS1$points
     , col = coldf$Group
     , labels = rownames(NMDS1$points)
)

#### Heatmap and clustering analysis ####
# Create a heatmap with statistically determined clusters for the x axis -- takes a while to run. Go have a cuppa.
if(!file.exists("CAZy_2/outputData/heatmap_data/heatmap_data_all_Clusters.rds")){
  
  clust_plot_object <- prepareClusterHeatmap(clust_mat
                                             , simprof_expected = 10000
                                             , simprof_simulated = 10000
                                             , distance_metric = Manly
  )
  saveRDS(clust_plot_object
          , file = "CAZy_2/outputData/heatmap_data/heatmap_data_all_Clusters.rds")
} else { 
  clust_plot_object <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_all_Clusters.rds")
}

# Check all clusters are in the data
unique(clust_plot_object$heatmap$heatmap_data$CAZyme)
unique(clust_plot_object$heatmap$CAZyme_pos_table$CAZyme)

# Simplify the cluster name
clust_plot_object$heatmap$heatmap_data$CAZyme <- factor(sub("4-", ""
           , as.character(clust_plot_object$heatmap$heatmap_data$CAZyme)))
clust_plot_object$heatmap$CAZyme_pos_table$CAZyme <- factor(sub("4-", ""
                                                            , as.character(clust_plot_object$heatmap$CAZyme_pos_table$CAZyme)))

# Remove galactose clusters
clust_plot_object$heatmap$heatmap_data <- data.frame(data.table(clust_plot_object$heatmap$heatmap_data)[-grep("galactose", CAZyme),])
clust_plot_object$heatmap$CAZyme_pos_table <- data.frame(data.table(clust_plot_object$heatmap$CAZyme_pos_table)[-grep("galactose", CAZyme),])


# data preprocessing to get a nice order for the heatmap

j <- clust_plot_object
j$heatmap$heatmap_data
unique(j$heatmap$heatmap_data$CAZyme)


# Order by number of phyla with >=20$ of genes correlated with a metabolite in that cluster
nHighProportion <- c(unlist(by(j$heatmap$heatmap_data, j$heatmap$heatmap_data$CAZyme
                               , FUN = function(x){
                                 length(x$Abundance[x$Abundance> 20]) 
                               })))

metabolite <- sub("-Caz.*", "", names(nHighProportion))


of <- data.frame(metabolite = metabolite
                 , CAZyme = names(nHighProportion)
                 , nHighProportion = nHighProportion)
of <- of[metabolite != "3,6-anhydro-D-galactose",] # remove galactose "clusters"
of$metabolite <- sub("4-", "", of$metabolite)
of$metabolite <- factor(of$metabolite, levels = c("glucose", "fucose", "xylose", "vanillic acid", "hydroxybenzoic acid", "benzoic acid"))

of <- rbindlist(by(of, of$metabolite, function(x){x[order(x$nHighProportion, decreasing = TRUE),]}))
of$column_order <- 1:length(of$metabolite)
column_order <- as.character(of$CAZyme) # input to function

# Make the plot
cp <- plotClusterHeatmap(clust_plot_object
                         , side_dend_width = 0.5
                         , dendrogram_left_margin = -0.7
                         , side_dendrogram_margin = c(0, -0.7, 0.2, 0)
                         , heatmap_margin = c(0.2, -0.7, 0.2, 0)
                         , x_axis_hjust = 1
                         , column_order = column_order
                         , heatmap_ylab = "NULL"
                         , colour_scale_title = "% Genes in Cluster\n(per Metabolite)")
cp$metabolitePlot

#### Edit the heatmap ####

cp$heatmap <- cp$heatmap + geom_vline(xintercept = c(5.5, 9.5, 14.5, 18.5, 25.5))

#### Investigate which clusters were abundant across all phyla of interst ####
unique(of$CAZyme)

of[order(of$nHighProportion),]
unique(of$CAZyme)

# 11 Phyla
#benzoic acid-CazC2
#glucose-CazC5

# 10 phyla
# 4-hydroxybenzoic acid-CazC7
# vanillic acid-CazC3

# 9 phyla
# xylose-CazC5
# fucose-CazC3
# fucose-CazC2

#### Evenness analysis ####
head(clust_mat)
head(df)

# Calculate evenness of CAZymes across CAZy clusters for each phylum
ef <- data.frame(Phylum = rownames(clust_mat)
                 , J = diversity(clust_mat)/log(specnumber(clust_mat)) # Pielou's evenness
)
# Plot this -- fits with clusters
plot(J ~ Phylum, data = ef, las = 2)

# Join this to the dataframe which is used to create the heatmap
ef <- dplyr::left_join(clust_plot_object$heatmap$Phylum_pos_table, ef, by = "Phylum")

ef

# Create a plot of CAZyme evenness by phylum
ep <- ggplot(ef, aes(x = J, y = y_center)) + 
  geom_bar(stat = 'identity', orientation = "y") +
  scale_y_continuous(#breaks = ef$y_center, 
    #labels = element_blank(),#ef$Phylum, 
    expand = c(0, 0)
  ) +
  geom_point(aes(x = 1, y = 1), colour = "white") +
  scale_x_continuous(expand = c(0, 0), name = "Pielou's J'"
                     , breaks = c(0,0.5,1)) +
  theme_classic() +
  theme(axis.text.y = element_blank()
        , axis.title.y = element_blank()
        , axis.text.x = element_text(size = 8)
        , axis.title.x = element_text(size = 10)
        , plot.margin = unit(c(0.2, 0.5, 0.2, 0.7), "cm")
        , axis.ticks.y = element_blank()
  )
ep

ef$J_rnd <- round(ef$J,2) 

ef <- ef[order(ef$J),]

# How evenly were the CAZymes distributed across CAZyme clusters for each phylum?
ef[ef$J_rnd <= 0.8,]
ef[ef$J_rnd <= 0.85 & ef$J_rnd > 0.8,]
ef[ef$J_rnd <= 0.89 & ef$J_rnd > 0.85,]
ef[ef$J_rnd >= 0.9,]

ef[order(ef$J),]
# 1     Acidobacteria       10      1       7 0.8879918  0.89
# 7        Firmicutes        9      1       6 0.8897453  0.89
# 11  Verrucomicrobia        3      1       3 0.8915203  0.89
# 4     Bacteroidetes       11      1       7 0.8968887  0.90
# 6     Euryarchaeota        4      1       3 0.8976399  0.90
# 3        Ascomycota        2      1       2 0.9158676  0.92

ggplot(ef, aes(x = J, y = y_center)) + 
  geom_bar(stat = 'identity', orientation = "y"
           , fill = alpha("grey50", 0.5)) +
  scale_y_continuous(breaks = ef$y_center, 
                     labels = ef$Phylum, 
                     expand = c(0, 0)
  ) +
  geom_point(aes(x = 1, y = 1), colour = "white") +
  scale_x_continuous(expand = c(0, 0), name = "J'"
                     #, breaks = c(0,1)
  ) +
  theme_bw() + 
  theme(axis.text.y = element_blank()
        , axis.title.y = element_blank()
        , axis.text.x = element_text(size = 8)
        , axis.title.x = element_text(size = 10)
        , plot.margin = unit(c(0.2, 0.5, 0.2, 0.7), "cm")
  ) + coord_flip()

#### Combine all of the plots together into a figure ####
cp_final <- plot_grid(cp$dendrogram + scale_x_reverse(), cp$heatmap, ep
                      , align = "h"
                      , axis = "tb"
                      , rel_widths = c(0.4, 1, 0.5)
                      , ncol = 3)
cp_final

# # Save the plot
# 
# save_plot("Final_figure_making/Final Figures/Figure_6_CAZy_cluster_richness_evenness.pdf"
#           , cp_final
#           , base_height = 4, base_width = 9
# )


# Read in plotting data
#clust_plot_object <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Final_Heatmap_Data_heatmap.rds")
clust_plot_object <- readRDS("CAZy_2/outputData/Final_Heatmap_Data_heatmap.rds")


column_order <- as.character(clust_plot_object$column_order) 
df <- clust_plot_object


x_axis_size = 8
y_axis_size = 8
x_axis_vjust = 1
x_axis_hjust = 0
heatmap_margin = c(0, 0.2, 0.2, 0.2) # trbl


dendrogram_left_margin = -1

colour_labels =  c(0,0.5,1)
heatmap_colour = "#152736"
colour_scale_labels = c(0, 0.5, 1)
colour_scale_title = "Number of Genes"
side_dend_width = 0.5
top_dend_height = 0.5
side_dendrogram_margin = c(0.2, 0.2, 0.2, dendrogram_left_margin) # trbl
heatmap_ylab = "Phylum"

library(ggrepel)

#### Heatmap plot ####

# Sort column order of the heatmap
if(!is.null(column_order)){
  column_order <- data.frame(CAZyme = column_order, column_order = 1:length(column_order))
  # A quick check
  
  # Join the two datasets
  df$heatmap$heatmap_data <- dplyr::left_join(df$heatmap$heatmap_data, column_order, by = "CAZyme")
  df$heatmap$heatmap_data$x_center <- df$heatmap$heatmap_data$column_order
  
  # Now work on CAZyme_pos_table as well
  df$heatmap$CAZyme_pos_table <- dplyr::left_join(df$heatmap$CAZyme_pos_table, column_order, by = "CAZyme")
  df$heatmap$CAZyme_pos_table$x_center <- df$heatmap$CAZyme_pos_table$column_order
}

# plt_hmap <- with(df$heatmap, (ggplot(heatmap_data, 
#                                      aes(x = factor(x_center), y = y_center
#                                          , fill = Abundance, 
#                                          height = height, width = width)) + 
#                                 geom_tile() +
#                                 scale_fill_gradient("Number of Genes"
#                                                     , high = "#152736", low = "white"
#                                                     , breaks = c(0,0.5,1)
#                                                     , labels = c(0,0.5,1)) +
#                                 scale_x_discrete(breaks = factor(CAZyme_pos_table$x_center), 
#                                                  labels = CAZyme_pos_table$CAZyme, 
#                                                  expand = c(0, 0)) + 
#                                 # For the y axis, alternatively set the labels as: Phylum_position_table$Phylum
#                                 scale_y_continuous(breaks = Phylum_pos_table[, "y_center"], 
#                                                    labels = Phylum_pos_table$Phylum,
#                                                    limits = Phylum_axis_limits, 
#                                                    expand = c(0, 0)) + 
#                                 labs(x = "CAZyme Cluster", y = "Phylum") +
#                                 theme_bw() +
#                                 theme(axis.text.x = element_text(size = x_axis_size
#                                                                  #, hjust = x_axis_hjust
#                                                                  , angle = 90
#                                                                  , vjust = x_axis_vjust
#                                                                  , hjust = x_axis_hjust),
#                                       axis.text.y = element_text(size = y_axis_size
#                                                                  , colour = Phylum_pos_table$Colour),
#                                       # margin: top, right, bottom, and left
#                                       plot.margin = unit(heatmap_margin, "cm"), 
#                                       panel.grid = element_blank(),
#                                       legend.position = 'bottom',
#                                       panel.spacing = unit(0, "lines")
#                                 ))
# )

#####

plt_hmap2 <- with(df$heatmap, (ggplot(heatmap_data, 
                                      aes(x = y_center, y = factor(x_center), 
                                          fill = Abundance, 
                                          height = height, width = width)) + 
                                 geom_tile() +
                                 scale_fill_gradient("Proportion of genes\nper metabolite"
                                                     , high = "#152736", low = "white"
                                                     #, trans = 'log10'
                                                     , breaks = c(0,0.5,1)
                                                     , labels =  c(0,0.5,1)) +
                                 scale_y_discrete(breaks = factor(CAZyme_pos_table$x_center) 
                                                  , labels = CAZyme_pos_table$CAZyme 
                                                  , expand = c(0, 0)) + 
                                 # For the y axis, alternatively set the labels as: Phylum_position_table$Phylum
                                 scale_x_continuous(breaks = Phylum_pos_table[, "y_center"] 
                                                    , labels = Phylum_pos_table$Phylum
                                                    , limits = Phylum_axis_limits 
                                                    , expand = c(0, 0)
                                                    , position = "top"
                                                    , sec.axis = dup_axis()
                                 ) + 
                                 labs(x = "Phylum"  #"CAZyme Cluster"
                                      , y = "CAZyme cluster"
                                 ) +
                                 theme_bw() +
                                 theme(axis.text.y = element_text(size = y_axis_size)
                                       , axis.text.x = element_text(size = x_axis_size
                                                                    , colour = Phylum_pos_table$Colour
                                                                    , angle = 90
                                                                    , vjust = x_axis_vjust
                                                                    , hjust = x_axis_hjust)
                                       , axis.title.x.top = element_blank()
                                       , axis.text.x.bottom = element_blank()
                                       # margin: top, right, bottom, and left
                                       , plot.margin = unit(c(-0.7,0.2, 0.5, 0.2), "cm") 
                                       , panel.grid = element_blank()
                                       #, legend.position = 'bottom'
                                       , panel.spacing = unit(0, "lines")
                                 ))
)

plt_hmap2

plt_hmap2 <- plt_hmap2 + geom_hline(yintercept = c(5.5, 9.5, 14.5
                                                   , 18.5, 25.5))


plt_dendr2 <- with(df$phylum_dend, ggplot(segment_data) + 
                     geom_segment(aes(x = y, y = x
                                      , xend = yend, yend = xend
                                      , col = col)) + 
                     #scale_x_reverse(expand = c(0, 0.5)) + 
                     scale_y_continuous(expand = c(0, 0.5)) +
                     scale_x_continuous(breaks = Phylum_pos_table$y_center, 
                                        limits = Phylum_axis_limits, 
                                        expand = c(0, 0)) + 
                     labs(x = "", y = "", colour = "", size = "") +
                     theme_bw() + 
                     theme(panel.grid = element_blank()
                           , axis.text = element_blank()
                           , rect = element_blank()
                           , axis.ticks = element_blank()
                           # trbl
                           , plot.margin = unit(c(0.2,0.2,-0.7,0.2), "cm")
                           , legend.position = "none"
                     ) 
)

plt_dendr2

plot_grid(plt_dendr2, plt_hmap2
          , ncol = 1
          , align = "v"
          , axis = "lr")

#####

#### Cluster gene richness evenness analysis -- all phyla ####
pc_evenness_all <- fread("CAZy_2/outputData/Metabolite-Correlated-CAZyme_Evenness_AllPhyla.csv")

pc_evenness_all.m <- melt(pc_evenness_all, id = "Sample")
pc_evenness_all.m 
pc_evenness_all.m <- pc_evenness_all.m[complete.cases(pc_evenness_all.m),]


sortedPhyla <- names(sort(tapply(pc_evenness_all.m$value
                                 , pc_evenness_all.m$variable, mean)
                          , decreasing = TRUE))
pc_evenness_all.m$variable <- factor(pc_evenness_all.m$variable
                                     , levels = sortedPhyla
)
names(pc_evenness_all.m) <- c("Sample", "Phylum", "Evenness")

####

df$heatmap
pc <- c("Colour", "Phylum")
phylumColours <- unique(data.table(df$heatmap$Phylum_pos_table)[,..pc])

phylumColours2 <- data.table(Colour = "grey30"
                             , Phylum = unique(pc_evenness_all.m$Phylum)[unique(pc_evenness_all.m$Phylum) %ni% phylumColours$Phylum])

phylumColours <- rbind(phylumColours, phylumColours2)
####

#pc_evenness_all.m$Colour <- "grey30"
pc_evenness_all.m <- merge(pc_evenness_all.m, phylumColours, by = "Phylum", all.x = TRUE)
#pc_evenness_all.m[is.na(Colour),]$Colour <- "grey30"

ea <- ggplot(pc_evenness_all.m, aes(x = Phylum, y = Evenness
                                    , fill =  Phylum
)) + 
  stat_summary(aes(y = Evenness
                   #, fill = Colour
  )
  , fun = mean, geom = "bar"
  ) +
  scale_fill_manual(#breaks = unique(pc_evenness_all.m$Colour), values = unique(pc_evenness_all.m$Colour)
    breaks = as.character(phylumColours$Phylum), values = phylumColours$Colour) +
  
  stat_summary(fun.data = mean_cl_boot
               , geom = "errorbar"
               , width = 0.3
               , size = 1
               #, fun.args = list(mult = 1)
  ) +
  geom_point(position = position_jitter(0.3), colour = alpha("black", 0.3), size = 1.5) + 
  ylim(0,1) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12)
        , axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.25)
        , axis.title.x = element_text(size = 14)
        , axis.title.y = element_text(size = 14)
        , plot.margin = unit(c(0.5, 0.5, 0.2, 0.2), "cm")
        , axis.ticks.y = element_blank()
        , legend.position = "none"
  ) +
  geom_point(data = data.frame(y_center = 16
                               , Evenness = 1.02
                               , Phylum = "Actinobacteria")
             , aes(x = y_center, y = Evenness)
             , colour = "white") +
  scale_y_continuous(expand = c(0, 0), name = "Pielou's J'"
                     , labels = scales::number_format(accuracy = 0.1)
  ) 

ea

#### Put them all together ####
final_Plot <- plot_grid(plt_dendr2, plt_hmap2, ea
                        , ncol = 1
                        , align = "v"
                        , axis = "lr"
                        , labels = c("(a)", "", "(b)")
                        , rel_heights = c(0.8, 1.5, 1.5))

final_Plot

# A4 plot
# save_plot(path  = "Final_figure_making/Final Figures/"
#           , filename = "Figure_7_CAZy_cluster_read_evenness.pdf"
#           , plot = final_Plot
#           , base_height = 11.75
#           , base_width = 8.25
# )


pl <- plot_grid(cp_final, ea, ncol = 1, rel_heights = c(1,1)
                , labels = c("a", "c"), label_fontface = "bold")

pr1 <- plot_grid(plt_dendr2, plt_hmap2, ncol = 1, rel_heights = c(0.8, 1.5)
                , align = "v", axis = "lr"
                #, labels = "b", label_fontface = "bold"
                )

BLANK <- ggplot() + cowplot::theme_nothing() +
  #trbl
  theme(plot.margin = unit(c(0.2, 0.2, -0.7, -0.7)
                           , "cm"))

pr <- plot_grid(pr1, BLANK 
                , labels = c("b", "")
                , label_fontface = "bold"
                , ncol = 1
                , rel_heights = c(1, 0.25)
                )

final_plot <- plot_grid(pl, pr)

save_plot(path  = "Final_figure_making/Final Figures/"
          , filename = "Figure_6_CAZy_clusters_evenness.pdf"
          , plot = final_plot
          , base_height = 8.25
          , base_width = 11.75
          
)
