rm(list = ls())

source("Functions/ClusterHeatmap_fromMatrix.R")


rf <- fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Phylum_Level_CAZyme_NormalizedReadAbundance.csv")
pf <- by(rf, rf$Phylum, FUN = function(x){
  dcast(x, Sample  ~ CAZyme, value.var = "Abundance")
})

pcpd <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Phylum_Level_CAZyme_Cluster_Percentage_Abundance_perMetabolite.rds")
heatmap(as.matrix(pcpd[[1]]))

acd <- rbindlist(pcpd, idcol = "Phylum")
acd$Sample <- rep(pf[[1]]$Sample, times = length(pcpd))

head(acd)
acd_mt <- as.matrix(acd[,2:(length(acd)-1)])
rownames(acd_mt) <- paste(acd$Phylum, acd$Sample, sep = " ")
heatmap(acd_mt)


if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/heatmap_data_all_Clusters_Reads_Replicated.rds")){
  
  clust_plot_object <- prepareClusterHeatmap(data.frame(acd[,2:(length(acd)-1)])
                                             , simprof_expected = 10000
                                             , simprof_simulated = 10000
                                             , distance_metric = Manly
  )
  saveRDS(clust_plot_object
          , file = "CAZy/04_CAZyme_Origin_Analysis/heatmap_data_all_Clusters_Reads_Replicated.rds")
} else { 
  clust_plot_object <- readRDS("CAZy/04_CAZyme_Origin_Analysis/heatmap_data_all_Clusters_Reads_Replicated.rds")
}

names(acd[,1:(length(acd)-1)])
names(acd[,2:(length(acd)-1)])

# No phylum names, however, the table is just a melted version of acd[,2:length-1]
length(melt(acd[,1:(length(acd)-1)], id.vars = "Phylum")$Phylum) == nrow(clust_plot_object$heatmap$heatmap_data)

clust_plot_object$heatmap$heatmap_data$Phylum <- melt(acd[,1:(length(acd)-1)]
                                                      , id.vars = "Phylum")$Phylum
clust_plot_object$heatmap$heatmap_data$Sample <- melt(acd[,1:(length(acd))]
                                                      , id.vars = c("Phylum","Sample"))$Sample

clust_plot_object$heatmap$Phylum_pos_table

clust_plot_object$heatmap$Phylum_pos_table <- merge(with(clust_plot_object$heatmap$heatmap_data
                                                         , data.frame(Sample = as.character(Sample)
                                                                      , Phylum = Phylum
                                                                      , y_center = y_center
                                                                      , stringsAsFactors = FALSE))
                                                    , clust_plot_object$heatmap$Phylum_pos_table
                                                    , by = "y_center"
                                                    )
names(clust_plot_object$heatmap$Phylum_pos_table)[3] <- "Phylum"

PhylumPalette <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c"
                   ,"#fb9a99","#e31a1c","#fdbf6f","#ff7f00"
                   ,"#cab2d6","#6a3d9a","#ffff99")
clust_plot_object$heatmap$Phylum_pos_table$Colour <- PhylumPalette[factor(clust_plot_object$heatmap$Phylum_pos_table$Phylum)]

# Remove the random galactose cluster left in
clust_plot_object$heatmap$heatmap_data <- data.table(clust_plot_object$heatmap$heatmap_data)[!grep("galactose",CAZyme),]

clust_plot_object$heatmap
clust_plot_object$phylum_dend
acd$Phylum


j <- clust_plot_object
j$heatmap$heatmap_data
unique(j$heatmap$heatmap_data$CAZyme)

# Order by number of phyla with >=20$ of genes correlated with a metabolite in that cluster
nHighProportion <- c(unlist(by(j$heatmap$heatmap_data, j$heatmap$heatmap_data$CAZyme
                               , FUN = function(x){
                                 length(x$Abundance[x$Abundance> 0.2]) 
                               })))
metabolite <- sub("\\.Caz.*", "", names(nHighProportion))
metabolite <- sub("\\.", " ", metabolite)
of <- data.frame(metabolite = metabolite
                 , CAZyme = names(nHighProportion)
                 , nHighProportion = nHighProportion)
of$metabolite <- factor(of$metabolite, levels = c("glucose", "fucose", "xylose", "vanillic acid", "hydroxybenzoic acid", "benzoic acid"))

of <- rbindlist(by(of, of$metabolite, function(x){x[order(x$nHighProportion, decreasing = TRUE),]}))
of$column_order <- 1:length(of$metabolite)
#of$CAZyme <- sub("\\.Caz", "-Caz", as.character(of$CAZyme))
#of$CAZyme <- sub("\\.", " ", as.character(of$CAZyme))

column_order <- as.character(of$CAZyme) # input to function

# Save the function input so we can make the final figure in a seperate script
clust_plot_object$column_order <- column_order

saveRDS(clust_plot_object, file = "CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Final_Heatmap_Data_heatmap.rds")

# Make the plot
cp <- plotClusterHeatmap(clust_plot_object
                         , side_dend_width = 0.5
                         , dendrogram_left_margin = -0.7
                         , side_dendrogram_margin = c(0, -0.7, 0.2, 0)
                         , heatmap_margin = c(0.2, -0.7, 0.2, 0)
                         , x_axis_hjust = 1
                         , column_order = column_order
                         , heatmap_ylab = "NULL"
                         , colour_scale_title = "Proportion of normalized reads\n in cluster (per Metabolite)"
                         , colour_labels = c(0, 0.5, 1))
cp$metabolitePlot

unique(clust_plot_object$heatmap$Phylum_pos_table$Phylum)

cp$heatmap + coord_flip()



#### Evenness plot ####
pc_evenness <- fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Metabolite-Correlated-CAZyme_Evenness_Phylum.csv")
pc_evenness.m <- melt(pc_evenness, id = "Sample")
names(pc_evenness.m) <- c("Sample", "Phylum", "Evenness")
pce <- pc_evenness.m
sortedPhyla <- names(sort(tapply(pce$Evenness, pce$Phylum, mean)
                          , decreasing = TRUE))
pce$Phylum <- factor(pce$Phylum, levels = sortedPhyla)

#pce <- dplyr::left_join(clust_plot_object$heatmap$Phylum_pos_table, pc_evenness.m, by = "Phylum")

#ef

# Create the plot using means and 95% CIs
em <- ggplot(pce, aes(x = Phylum, y = Evenness)) + 
  stat_summary(aes(y = Evenness)
               , fun = mean, geom = "bar", colour = "grey30") +
  stat_summary(fun.data = mean_cl_boot
               , geom = "errorbar"
               #, fun.args = list(mult = 1)
               , width = 0.3
               , size = 1) +
  geom_point(position = position_jitter(0.3), colour = alpha("black", 0.3), size = 1.5) + 
  ylim(0,1) +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  #) +
  theme(axis.text.y = element_text(size = 12)
        #, axis.title.y = element_text(label = )
        , axis.text.x = element_text(size = 12, angle = 90)
        , axis.title.x = element_text(size = 14)
        , axis.title.y = element_text(size = 14)
        , plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm")
        , axis.ticks.y = element_blank()
  ) +
  geom_point(data = data.frame(y_center = 1
                               , Evenness = 1)
             , aes(x = y_center, y = Evenness)
             , colour = "white") +
  scale_y_continuous(expand = c(0, 0), name = "Pielou's J'"
                     , labels = scales::number_format(accuracy = 0.1)
                     #, breaks = c(0,1)
  ) 

#### Cluster gene richness evenness analysis -- all phyla ####
pc_evenness_all <- fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Metabolite-Correlated-CAZyme_Evenness_AllPhyla.csv")

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

ea <- ggplot(pc_evenness_all.m, aes(x = Phylum, y = Evenness)) + 
  stat_summary(aes(y = Evenness)
               , fun = mean, geom = "bar", colour = "grey30") +
  stat_summary(fun.data = mean_cl_boot
               , geom = "errorbar"
               , width = 0.3
               , size = 1) +
  geom_point(position = position_jitter(0.3), colour = alpha("black", 0.3), size = 1.5) + 
  ylim(0,1) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12)
        , axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.25)
        , axis.title.x = element_text(size = 14)
        , axis.title.y = element_text(size = 14)
        , plot.margin = unit(c(0.5, 0.5, 0.2, 0.2), "cm")
        , axis.ticks.y = element_blank()
  ) +
  geom_point(data = data.frame(y_center = 16
                               , Evenness = 1.02)
             , aes(x = y_center, y = Evenness)
             , colour = "white") +
  scale_y_continuous(expand = c(0, 0), name = "Pielou's J'"
                     , labels = scales::number_format(accuracy = 0.1)
  ) 

ea

#### Final results ####

# (1) Read abundance evenness plot
em 

# Save this figure
save_plot(path = "Figures/CAZy_Taxonomy/Heatmaps/"
          , filename = "CAZyCluster_ReadAbundance_Evenness.pdf"
          , plot = em
          , base_height = 4
          , base_width = 4)

# (2) Big heatmap and dendrogram of read abundances across metabolite-correlated  
#     CAZyme clusters

dend_plot <- cp$dendrogram + scale_x_reverse() + 
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))

heatmap_plot <- cp$heatmap + theme(plot.margin = unit(c(0.2, 0.8, 0.2, 0), "cm")
                                   , axis.text.y = element_text(size = 5))
readHeatmap <- plot_grid(dend_plot
                         , heatmap_plot

                         , align = "h"
                         , axis = "tb"
                         , rel_widths = c(0.8, 1)
                         , ncol = 2
                         )

readHeatmap

# Save this figure as a4 so all labels can be read
save_plot(path = "Figures/CAZy_Taxonomy/Heatmaps/"
          , filename = "CazyCluster_ReadAbundance_Heatmap.pdf"
          , plot = readHeatmap
          , base_height = 11.75
          , base_width = 8.25)

#### (3) Richness evenness plot ####
# Save this figure
ea

save_plot(path = "Figures/CAZy_Taxonomy/Heatmaps/"
          , filename = "CAZyCluster_Richness_Evenness.pdf"
          , plot = ea
          , base_height = 4
          , base_width = 6)
