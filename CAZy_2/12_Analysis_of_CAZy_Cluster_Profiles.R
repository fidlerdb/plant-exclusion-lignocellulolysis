
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
of$metabolite <- factor(of$metabolite, levels = c("glucose", "fucose", "xylose", "vanillic acid", "4-hydroxybenzoic acid", "benzoic acid"))

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

# Save the plot

save_plot("Figures/CAZy_2/Cluster_Overview_Heatmap_Richness.pdf"
          , cp_final
          , base_height = 4, base_width = 7
          )

#### Evenness Analysis V2 (from script 12 in this directory) ####

# Read in the sample-wise CAZyme cluster evenness data (phylum-level)
pc_evenness <- fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Metabolite-Correlated-CAZyme_Evenness_Phylum.csv")
pc_evenness.m <- melt(pc_evenness, id = "Sample")
names(pc_evenness.m) <- c("Sample", "Phylum", "Evenness")

pce <- dplyr::left_join(clust_plot_object$heatmap$Phylum_pos_table, pc_evenness.m, by = "Phylum")

ef

# Create the plot using means and 95% CIs
em <- ggplot(pce, aes(x = y_center, y = Evenness)) + 
  stat_summary(aes(y = Evenness)
               , fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_boot
               , geom = "errorbar"
               #, fun.args = list(mult = 1)
               , width = 0.3
               , size = 0.5) +
  geom_point(position = position_jitter(0.3), colour = alpha("black", 0.3), size = 0.8) + 
  ylim(0,1) +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  #) +
  theme(axis.text.y = element_blank()
        , axis.title.y = element_blank()
        , axis.text.x = element_text(size = 8)
        , axis.title.x = element_text(size = 10)
        , plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm")
        , axis.ticks.y = element_blank()
  ) +
  geom_point(data = data.frame(y_center = 1
                               , Evenness = 1)
             , aes(x = y_center, y = Evenness)
             , colour = "white") +
  scale_y_continuous(expand = c(0, 0), name = "J'"
                     , labels = scales::number_format(accuracy = 0.1)
                     #, breaks = c(0,1)
  ) +
  scale_x_continuous(breaks = pce$y_center, 
                     labels = pce$Phylum, 
                     expand = c(0, 0)
  ) #+
 # coord_flip()

em # Read abundance evenness plot

# Put the plots together
cp_final <- plot_grid(cp$dendrogram + scale_x_reverse(), cp$heatmap, em
                      , align = "h"
                      , axis = "tb"
                      , rel_widths = c(0.3, 1, 0.5)
                      , ncol = 3)

cp_final

# Save the figure
# save_plot("Figures/CAZy_Taxonomy/Heatmaps/Cluster_Overview_Heatmap_2020-09-28.pdf"
#           , cp_final
#           , base_height = 4, base_width = 7
# )

