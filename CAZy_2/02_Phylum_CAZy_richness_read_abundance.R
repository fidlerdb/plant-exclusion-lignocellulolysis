rm(list = ls())

library(data.table)
library(ggplot2)
library(scales)
library(doBy)
library(cowplot)

#### CAZy family richness plot ####

# Read in data
df <- fread("CAZy_2/outputData/Taxonomy_Data_CAZymes.csv")
df <- df[Phylum != "",] # Keep only data with phylum-level information

#  Isolate the CAZyme matrix for the phylum in question
CAZymes <- names(df)[24:length(df)]

# Calculate phylum-level CAZy family richness values
rf <- c(by(df, df$Phylum, FUN = function(x){
  # Calculate richness
  sum(ifelse(colSums(x[, ..CAZymes]) > 0, yes = 1, no = 0))
  }, simplify = TRUE
  )
)
rf <- data.table(Phylum = names(rf), Richness = rf) # Make the data mergeable
tf <- unique(df[,c("SuperKingdom", "Phylum")]) # Data to merge to
rf <- merge(tf, rf, by = "Phylum") # Merge the data

rf
rf[order(rf$Richness,)]

# Palette
DomainPalette <- c('#e41a1c'# Archaea
                   , '#377eb8'# Bacteria
                   , '#4daf4a'# Eukaryotes
)

# Create a plot
r <- ggplot(rf, aes(y = reorder(Phylum, Richness), x = Richness
                    , fill = SuperKingdom)) + 
  geom_vline(xintercept = c(1, 10, 100)
             , colour = "grey90") +
  geom_bar(stat = 'identity')+ 
  scale_x_continuous(trans = 'log10') +
  annotation_logticks(sides = "tb") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()
        , panel.grid.minor.x = element_blank()) +
  xlab("CAZy family richness (log scale)") +
  ylab("Phylum") +
  scale_fill_manual(values = DomainPalette) + 
  guides(fill = 'none'#guide_legend(title = "Domain")
  )

r

#### Read abundance plot ####

# Read in data
df <- fread("CAZy_2/outputData/Taxonomy_Data_CAZymes.csv")
df <- df[Phylum != "",] # Keep only data with phylum-level information

head(df[,1:30])

# Summarise the normalised read abundances at the phylum level
ad <- summaryBy(b1+b2+b3+b4+b5+b6+b7+
                  c1+c2+c3+c4+c5+c6+c7 ~  Phylum + SuperKingdom + CAZyme
                , data = df, FUN = sum, keep.names = TRUE)

# Make ggplot understand the data
ad <- melt(ad, id = c("Phylum", "SuperKingdom")
           , variable.name = "Sample"
           , value.name = "Reads")

head(ad)

# Palette
DomainPalette <- c('#e41a1c'# Archaea
                   , '#377eb8'# Bacteria
                   , '#4daf4a'# Eukaryotes
)

# Make the plot
p <- ggplot(ad, aes(y = reorder(Phylum, Reads), x = Reads+1
                    , fill = SuperKingdom)) + 
  geom_vline(xintercept = c(0.1, 1, 10, 100, 1000, 10000, 100000, 1e6)
             , colour = "grey90") +
  stat_summary(fun = mean, geom = "bar", ) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar"
               , width = 0.3, fun.args = list(mult = 1)) +
  #geom_bar(stat = 'identity') + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank()
        , panel.grid.minor.x = element_blank()
        , axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Normalized read abundance + 1 (log scale)") +
  ylab("Phylum") +
  scale_fill_manual(values = DomainPalette) + 
  guides(fill = guide_legend(title = "Domain")) +
  scale_x_continuous(
     breaks = c(0.01 ,0.1, 1, 10, 100, 1000, 10000, 100000, 1e6)
    , trans = "log10"
    #, labels = number_format(accuracy = 0.1)
    , labels = formatC(c(0.01 ,0.1, 1, 10, 100, 1000, 10000, 100000, 1e6))
  ) 

p

# Create a legend for both plots
legend <- get_legend(p)

# Strip the legend from a plot
p <- p + guides(fill = "none")
p

#### Combine the plots ####

# Put the two plots together
pp1 <- plot_grid(r, p
                 , labels = c("(a)", "(b)")
                 , nrow = 1, rel_widths = c(1,1)
                 , align = "h"
                 , axis = "tb"
                 )

# Add the legend
pp <- plot_grid(pp1, legend
                , nrow = 1, rel_widths = c(1,0.15))

pp

# Save the figure
ggsave("Figures/CAZy_2/Phylum_CAZy_Richness+Abundance_2020-10-05.pdf"
       , pp
       , device = "pdf"
       , width = 8.5
       , height = 5
)

#### Values for writing about ####

# Richness of CAZy families belonging to each phylum
r
rf[order(rf$Richness,)]

# Percentage of normalized reads belonging to each phylum
p
ad_means <- summaryBy(Reads ~ Phylum + SuperKingdom, data = ad, FUN = c(mean, sd))
ad_means <- ad_means[order(ad_means$Reads.mean),]
ad_means$Percentage <- round(100*(ad_means$Reads.mean / sum(ad_means$Reads.mean)),2)
ad_means$sd.Percentage <- round(100*(ad_means$Reads.sd / sum(ad_means$Reads.mean)),2)

####

SampleTotals <- tapply(X = ad$Reads
                       , INDEX = ad$Sample
                       , sum)
SampleTotals <- data.table(Sample = names(SampleTotals), SampleTotal = c(SampleTotals))
ad <- merge(ad, SampleTotals, by = "Sample")

# Calculate the percentage each species made up of each sample
ad$PercentageAbundance <- 100 * (ad$Reads/ ad$SampleTotal)

# Calculate mean percentage (and SD) each species made up of each sample
MeanPercentageSD  <- summaryBy(PercentageAbundance ~ Phylum, data = ad, FUN = c(mean, sd))

MeanPercentageSD <- MeanPercentageSD[order(MeanPercentageSD$PercentageAbundance.mean),]
MeanPercentageSD$MeanPercentage <- round(MeanPercentageSD$PercentageAbundance.mean, 2)
MeanPercentageSD$SDPercentage <- round(MeanPercentageSD$PercentageAbundance.sd, 2)

MeanPercentageSD

####

##### Origins of CAZymes ####

sp <- c("Species", "Phylum")

# 10-year fallow treatment
# Actinobacteria
df[AA10 == 1, ..sp]
df[CBM3 == 1 , ..sp]
df[GH101 == 1 , ..sp]
df[GH101 == 1 , ..sp]
df[GH43_4 == 1 , ..sp]
df[`CBM22|CBM6|GH10` == 1 , ..sp]

# The others
df[`CBM8|GH44` == 1, ..sp]$S
df[GH43_35 == 1, ..sp]
df[GH45 == 1, ..sp]
df[GH5_38 == 1, ..sp]

# Other treatments
df[GH43_32 == 1, ..sp]
#                        Species              Phylum
# 1: Deinococcus peraridilitoris Deinococcus-Thermus
# 2:       Hymenobacter swuensis       Bacteroidetes
# 3:          Nostoc punctiforme       Cyanobacteria
df[GH135  == 1, ..sp]
#                             Species         Phylum
# 1: Candidatus Nitrosotenuis cloacae Thaumarchaeota
# 2:           Unknown Actinobacteria Actinobacteria
# 3:             Salinimonas sp. N102 Proteobacteria
# 4:            Nitrosospira briensis Proteobacteria
df[`CBM48|CBM68|GH13_14` == 1, ..sp]
#     Species     Phylum
# 1: Bacillus asahii Firmicutes
df[CBM68 == 1, ..sp]
#           Species            Phylum
# 1: Geobacillus sp. WCH70 Firmicutes








