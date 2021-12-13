
# Reference data for the plot
at <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData/CAZy_Family_Activity_Summary.csv")

plot_Metabolite_Activities <- function(metabolite, heatmapdata){
  # ensure the subsetting works properly
  if(metabolite == "benzoic acid"){metabolite <- paste0("^", metabolite)}

  # Subsetting of the data
  af <- merge(at[grep(metabolite, at$Cluster),], heatmapdata$heatmap$CAZyme_pos_table[,1:3]
              , by = "CAZyme")

  # Put it in ggplot2able format
  afl <- melt(af, id = c("CAZyme", "Cluster", "x_center", "width"))
  afl
  afl$value <- as.character(afl$value)
  afl[value == 1,]$value <- "Known activity"
  afl[value == 0,]$value <- "No known activity"
  afl$variable <- factor(afl$variable, levels = c("Oligosaccharides"
                                                  , "Cellulose"
                                                  , "Hemicellulose"
                                                  , "Lignin"))
  
  # Create the plot
  p <- ggplot(data = afl
              , aes(x = x_center, y = variable, fill = value)
              , height = 1, width = 1) + 
    geom_tile(colour = "black") +
    scale_fill_manual(values = c("black","white")) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.text.x = element_blank()
          , axis.title = element_blank()
          , panel.grid = element_blank()
          # trbl
          , plot.margin = unit(c(-1.1, -0.7, -0.5, 0.2), "cm")
          , panel.spacing = unit(0, "lines")
          #, legend.position = "none"
          , panel.border = element_rect(fill = NA)
    ) +
    
    guides(fill = guide_legend(""))
  
  return(p) 
}


#plot_Metabolite_Activities("vanillic acid", df_va)
