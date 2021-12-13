
rm(list = ls())

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
save_plot(path  = "Figures/CAZy_2/"
          , filename = "Multiple_CazyCluster_ReadAbundance_Heatmap_barColours.pdf"
          , plot = final_Plot
          , base_height = 11.75
          , base_width = 8.25
          )
  