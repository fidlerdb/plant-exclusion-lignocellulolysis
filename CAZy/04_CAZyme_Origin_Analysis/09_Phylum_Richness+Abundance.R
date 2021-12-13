# Summarize and explore the link between taxonomy and lignocellulase genes

# Clear the workspace
.rs.restartR()
rm(list = ls())

# Load useful packages
library(data.table)
library(doBy) #
library(ggplot2) #
library(cowplot)

# Steal old code to make prettier figures.


########### Richness  ############

# Read in the taxonomy data--normalized
df <- data.frame(fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_dbCAN_Taxonomy_ReadAbundance.csv"))

# Keep only contigs which have both CAZy genes and a taxonomic assignment
df <- df[!is.na(df$CAZyme),]
# Remove viral contigs
df <- df[df$SuperKingdom != "Viruses",]

df[df$n1 != 0,]$Contig # data has been cleaned properly

##### Check for reasons behind disparity in results ####
unique(gsub("\\|", ".", unique(df[df$Phylum == "Proteobacteria",]$CAZyme)))

cf <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZymes_NotGrouped_Boolean.csv")
cf[, 4:length(cf)] <- cf[, 4:length(cf)] / (cf$Reads/1e9) # Normalise per Gbp


head(cf[1:5])

AllCAZymes <- unique(gsub("\\.[0-9]+", "", names(cf)[4:length(cf)]))
length(AllCAZymes)
head(AllCAZymes)
length(grep("GH", AllCAZymes))
length(grep("CBM", AllCAZymes))
length(grep("AA", AllCAZymes))

AllCAZymes_richness <- unique(gsub("\\|", ".", unique(df$CAZyme)))

library(roperators)

AllCAZymes_richness[AllCAZymes_richness %ni% AllCAZymes]


#####

# Which phyla had CAZymes, and how many?
sf <- summaryBy(CAZyme ~ Phylum + SuperKingdom, data = df, FUN = function(x)length(unique(x)))
sf <- sf[order(sf$`CAZyme.function(x) length(unique(x))`),]
names(sf)[3] <- "CAZymes"
sf$Phylum[is.na(sf$Phylum)] <- paste("Unknown", sf$SuperKingdom[is.na(sf$Phylum)], sep = " ")
sf$Phylum[sf$Phylum == "Unknown NA"] <- "Unknown"
sf <- sf[-grep("Unknown", sf$Phylum),] # Remove "unknown whatever"

sf

# What CAZymes did each phylum have?
by(data = df, INDICES = df$Phylum
   , FUN = function(x){
     g <- x$CAZyme
     return(unique(g))
   }
   , simplify = TRUE)

unique(df[df$Phylum == "Proteobacteria",]$CAZyme)

# Palette
DomainPalette <- c('#e41a1c'# Archaea
                   , '#377eb8'# Bacteria
                   , '#4daf4a'# Eukaryotes
)

# Create a plot
r <- ggplot(sf, aes(y = reorder(Phylum, CAZymes), x = CAZymes
                     , fill = SuperKingdom))+ 
  geom_bar(stat = 'identity')+ 
  scale_x_continuous(trans = 'log10') +
  annotation_logticks(sides = "tb") +
  theme_bw() +
  xlab("CAZy family richness (log scale)") +
  ylab("Phylum") +
  scale_fill_manual(values = DomainPalette) + 
  guides(fill = 'none'#guide_legend(title = "Domain")
         )

############ Read Abundance ##############

# Read in the taxonomy data--normalized
rf <- data.frame(fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_dbCAN_Taxonomy_ReadAbundance.csv"))

# Keep only contigs which have both CAZy genes and a taxonomic assignment
rf <- rf[!is.na(rf$CAZyme),]
#rf <- rf[!is.na(rf$Phylum),]

# Summarise the data at the phylum level -- CAZyme reads in each sample annotated by phylum
rf <- summaryBy(b1+b2+b3+b4+b5+b6+b7+ 
                  c1+c2+c3+c4+c5+c6+c7 ~  Phylum + SuperKingdom + CAZyme
                , data = rf, FUN = sum, keep.names = TRUE)
rf <- data.table(rf)

# Make it long
rf <- melt(rf, id.vars = c("Phylum", "SuperKingdom", "CAZyme")
           , variable.name = "Sample", value.name = "Abundance")

# Add a treatment column
rf$Treatment <- "a"
rf[Sample %in% c("b1","b2","b3"),]$Treatment <- "10-Year Fallow"
rf[Sample %in% c("b4","b5","b6", "b7"),]$Treatment <- "1-Year Fallow"
rf[Sample %in% c("c1","c2","c3"),]$Treatment <- "10-Year Grassland"
rf[Sample %in% c("c4","c5","c6", "c7"),]$Treatment <- "1-Year Grassland"

# Sum values from all samples together to figure out % of CAZyme reads coming from each phylum
#rf <- summaryBy(Abundance ~  Phylum + SuperKingdom + CAZyme
#                , data = rf, FUN = sum, keep.names = TRUE)
#rf <- rbindlist(by(rf, INDICES = rf$CAZyme
#                  , FUN = function(x){
#                     TotAbund <- sum(x$Abundance)
#                     x$Abundance <- (x$Abundance / TotAbund) * 100
#                     if(TotAbund == 0){x$Abundance <- 0}
#                     return(x)
#                   }))

# Where did the majority of CAZyme reads come from?
tcf <- summaryBy(Abundance ~ Phylum + SuperKingdom, data = rf, FUN = sum)
tcf$Phylum <- factor(tcf$Phylum, levels = tcf$Phylum[order(tcf$Abundance.sum)])
tcf$PercentAbundance <- (tcf$Abundance.sum / sum(tcf$Abundance.sum))*100

# Palette
DomainPalette <- c('#e41a1c'# Archaea
                   , '#377eb8'# Bacteria
                   , '#4daf4a'# Eukaryotes
)

# Make the plot
p <- ggplot(tcf, aes(y = Phylum, x = PercentAbundance
                     , fill = SuperKingdom)) + 
  geom_bar(stat = 'identity') + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank()
        , panel.grid.minor.x = element_blank()) +
  xlab("Percentage normalized read abundance (log scale)") +
  scale_fill_manual(values = DomainPalette) + 
  guides(fill = guide_legend(title = "Domain")) +
  scale_x_continuous(breaks = c(0.00001, 0.0001, 0.001
                                , 0.01, 0.1, 1
                                , 10, 100)
  , trans = "log10"
  , labels = number_format(accuracy = 0.001)
  ) 

p

legend <- get_legend(p)

p <- p + guides(fill = "none")
p

########## Combine the two plots ############

pp <- plot_grid(r, p, legend
          , labels = c("a)", "b)", "")
          , nrow = 1, rel_widths = c(1,1,0.3)
          )

ggsave("Figures/CAZy_Taxonomy/Phylum_CAZy_Richness+Abundance_2020-10-03.pdf"
       , pp
       , device = "pdf"
       , width = 8.5
       , height = 5
       )
