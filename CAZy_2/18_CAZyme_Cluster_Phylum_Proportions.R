rm(list = ls())

# Read in plotting data
#clust_plot_object <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Final_Heatmap_Data_heatmap.rds")
df <- readRDS("CAZy_2/outputData/Final_Heatmap_Data_heatmap.rds")

df <- df$heatmap$heatmap_data 

head(df)
df <- data.table(df[,c(1:3,10)])
df

abcl <- c("glucose-CazC5",
          "fucose-CazC2",
          "fucose-CazC3",
          "xylose-CazC5",
          "vanillic acid-CazC3",
          "hydroxybenzoic acid-CazC7",
          "benzoic acid-CazC2"
          )
unique(df$CAZyme)[unique(df$CAZyme) %in% abcl]

df <- df[CAZyme %in% abcl,]

# Calculate mean and sd
df_extra <- data.table(CAZyme = names(tapply(df$Abundance,  df$CAZyme, mean))
           , Mean_Abundance = c(tapply(df$Abundance,  df$CAZyme, mean))
           , SD_Abundance = c(tapply(df$Abundance,  df$CAZyme, sd))
           , Percentile_25 = c(tapply(df$Abundance,  df$CAZyme, FUN = function(x){quantile(x, 0.25)}))
           , Percentile_75 = c(tapply(df$Abundance,  df$CAZyme, FUN = function(x){quantile(x, 0.75)}))
           , Percentile_16 = c(tapply(df$Abundance,  df$CAZyme, FUN = function(x){quantile(x, 0.159)}))
           , Percentile_84 = c(tapply(df$Abundance,  df$CAZyme, FUN = function(x){quantile(x, 0.841)}))
           )

ggplot(df, aes(y = Abundance, x = Phylum)) + 
  #geom_point() + 
  facet_wrap(~CAZyme, ncol = 4) +
  stat_summary(fun.data = mean_cl_boot) + coord_flip() + 
  theme_classic() +
  theme(panel.border = element_rect(fill = NA)) +
  geom_hline(data = df_extra, aes(yintercept = Mean_Abundance), colour = "red", size = 1.2) +
  geom_hline(data = df_extra, aes(yintercept = Mean_Abundance + SD_Abundance), linetype = "dashed", colour = "red", size = 1) +
  geom_hline(data = df_extra, aes(yintercept = Mean_Abundance - SD_Abundance), linetype = "dashed", colour = "red", size = 1) +
  geom_hline(data = df_extra, aes(yintercept = Percentile_16), linetype = "dashed", colour = "darkgreen", size = 1) +
  geom_hline(data = df_extra, aes(yintercept = Percentile_84), linetype = "dashed", colour = "darkgreen", size = 1)



lapply(seq_along(unique(df$CAZyme)), FUN = function(j){
  
  n <- unique(df$CAZyme)[[j]]
  dw <- df[df$CAZyme == n,]
  
  LessThanLowerSD = c(tapply(dw$Abundance,  dw$Phylum, function(x){
    mean(x) < (mean(dw$Abundance) - sd(dw$Abundance))
  }))
  
  Phyla = names(LessThanLowerSD)
  
  GreaterThanUpperSD = c(tapply(dw$Abundance,  dw$Phylum, function(x){
    mean(x) > (mean(dw$Abundance) + sd(dw$Abundance))
  }))
  
  return(data.table(Phylum = Phyla
                    , LessThanLowerSD
                    , GreaterThanUpperSD
                    , Cluster = n
  ))
  
})



#### Now repeat this for the richness data ####

rm(list = ls())

# Read in plotting data
#clust_plot_object <- readRDS("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Final_Heatmap_Data_heatmap.rds")
df <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_all_Clusters.rds")

df <- df$heatmap$heatmap_data 

head(df)
df <- data.table(df[,1:3])
df

abcl <- c("glucose-CazC5",
          "fucose-CazC2",
          "fucose-CazC3",
          "xylose-CazC5",
          "vanillic acid-CazC3",
          "4-hydroxybenzoic acid-CazC7",
          "benzoic acid-CazC2"
)
unique(df$CAZyme)[unique(df$CAZyme) %in% abcl]

df <- df[CAZyme %in% abcl,]

# Calculate mean and sd
df_extra <- data.table(CAZyme = names(tapply(df$Abundance,  df$CAZyme, mean))
                       , Mean_Abundance = c(tapply(df$Abundance,  df$CAZyme, mean))
                       , SD_Abundance = c(tapply(df$Abundance,  df$CAZyme, sd))
                       , Percentile_25 = c(tapply(df$Abundance,  df$CAZyme, FUN = function(x){quantile(x, 0.25)}))
                       , Percentile_75 = c(tapply(df$Abundance,  df$CAZyme, FUN = function(x){quantile(x, 0.75)}))
                       , Percentile_16 = c(tapply(df$Abundance,  df$CAZyme, FUN = function(x){quantile(x, 0.159)}))
                       , Percentile_84 = c(tapply(df$Abundance,  df$CAZyme, FUN = function(x){quantile(x, 0.841)}))
)

ggplot(df, aes(y = Abundance, x = Phylum)) + 
  geom_point() + 
  facet_wrap(~CAZyme, ncol = 4) + 
  coord_flip() + 
  theme_classic() +
  ylim (0,100) +
  theme(panel.border = element_rect(fill = NA)) +
  geom_hline(data = df_extra, aes(yintercept = Mean_Abundance), colour = "red", size = 1.2) +
  geom_hline(data = df_extra, aes(yintercept = Mean_Abundance + SD_Abundance), linetype = "dashed", colour = "red", size = 1) +
  geom_hline(data = df_extra, aes(yintercept = Mean_Abundance - SD_Abundance), linetype = "dashed", colour = "red", size = 1) +
geom_hline(data = df_extra, aes(yintercept = Percentile_16), linetype = "dashed", colour = "darkgreen", size = 1) +
geom_hline(data = df_extra, aes(yintercept = Percentile_84), linetype = "dashed", colour = "darkgreen",size = 1)

unique(df$CAZyme)[[6]]

# Which phyla had the lowest and highest percentage of genes in each cazyme cluster?

lapply(seq_along(unique(df$CAZyme)), FUN = function(j){
  
  n <- unique(df$CAZyme)[[j]]
  dw <- df[df$CAZyme == n,]
  
  LessThanLowerSD = c(tapply(dw$Abundance,  dw$Phylum, function(x){
    mean(x) < (mean(dw$Abundance) - sd(dw$Abundance))
  }))
  
  Phyla = names(LessThanLowerSD)
  
  GreaterThanUpperSD = c(tapply(dw$Abundance,  dw$Phylum, function(x){
    mean(x) > (mean(dw$Abundance) + sd(dw$Abundance))
  }))
  
  return(data.table(Phylum = Phyla
                    , LessThanLowerSD
                    , GreaterThanUpperSD
                    , Cluster = n
                    ))
  
})

