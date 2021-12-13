# What were the origins of the most abundant CAZy families?
.rs.restartR()
rm(list = ls())

library(data.table)
library(doBy)
library(ggplot2)

#### Import, check, and set the data up ####

# Read in the taxonomy data--normalized
df <- data.frame(fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_dbCAN_Taxonomy_ReadAbundance.csv"))

# See how this looks
head(df)

# Keep only contigs which have both CAZy genes and a taxonomic assignment

df <- df[!is.na(df$CAZyme),]
df <- df[!is.na(df$Phylum),]

head(df)

# Summarise the data at the phylum level -- CAZyme reads in each sample annotated by phylum
rf <- summaryBy(b1+b2+b3+b4+b5+b6+b7+ 
            c1+c2+c3+c4+c5+c6+c7 ~  Phylum + SuperKingdom + CAZyme
          , data = df, FUN = sum, keep.names = TRUE)
rf <- data.table(rf)
head(rf)

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

rf <- summaryBy(Abundance ~  Phylum + SuperKingdom + CAZyme
                , data = rf, FUN = sum, keep.names = TRUE)

rf <- rbindlist(by(rf, INDICES = rf$CAZyme
                   , FUN = function(x){
                     TotAbund <- sum(x$Abundance)
                     x$Abundance <- (x$Abundance / TotAbund) * 100
                     if(TotAbund == 0){x$Abundance <- 0}
                     return(x)
                   }))

# Keep only the most abundant cazymes or the ones which were 
# shifted to in the different treatments

AbundantCAZymes <- c("GH15",'GH13_11','GH23','GH3','GH103')

IndicativeCAZymes_b10 <- c('CBM68', 'GH43_32', 'CBM48|CBM68|GH13_14'
                           , 'GH135')
IndicativeCAZymes_gb1 <- c('GH43_4', 'AA10', 'CBM3', 'GH43_35', 'GH101'
                           , 'GH45', 'CBM8|GH44', 'CBM22|CBM6|GH10', 'GH5_38')

# Subset the data to make plotting easier
#rf_ib <- rf[CAZyme %in% IndicativeCAZymes_b10]
#rf_ig <- rf[CAZyme %in% IndicativeCAZymes_gb1]
#rf_a <- rf[CAZyme %in% AbundantCAZymes]

# Get taxonomic origins of CAZy families from these plots.
ggplot(rf[CAZyme %in% IndicativeCAZymes_b10,]
       , aes(x = Phylum, y = Abundance)) + 
  geom_bar(stat = 'identity') +
  facet_wrap(~ CAZyme)

ggplot(rf[CAZyme %in% IndicativeCAZymes_gb1,]
       , aes(x = Phylum, y = Abundance)) + 
  geom_bar(stat = 'identity') +
  facet_wrap(~ CAZyme)


head(rf)

# Where did the majority of CAZyme reads come from?
tcf <- summaryBy(Abundance ~ Phylum + SuperKingdom, data = rf, FUN = sum)
tcf$Phylum <- factor(tcf$Phylum, levels = tcf$Phylum[order(tcf$Abundance.sum)])
tcf$PercentAbundance <- (tcf$Abundance.sum / sum(tcf$Abundance.sum))*100
#tcf$X = 1
#tcf <- tcf[order(-tcf$Abundance.sum),]

DomainPalette <- c('#e41a1c'# Archaea
                   , '#377eb8'# Bacteria
                   , '#4daf4a'# Eukaryotes
)

p <- ggplot(tcf, aes(y = Phylum, x = Abundance.sum
                , fill = SuperKingdom
                #, label = Phylum
                )) + 
  geom_bar(stat = 'identity') + 
  scale_x_continuous(trans = 'log10') +
  annotation_logticks(sides = "tb") +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5
  #                                 , hjust  = 1)) +
  xlab("Normalized Read Abundance") +
  scale_fill_manual(values = DomainPalette) + 
  guides(fill = guide_legend(title = "Domain"))

ggsave("Figures/CAZy_Taxonomy/Phylum_CAZy_Read_Abundance.pdf"
       , p
       , width = 6
       , height = 5
       , units = "in"
       )

tcf[Phylum %in% c("Proteobacteria", "Actinobacteria", "Acidobacteria"),]
tcf[Phylum %in% c("Caldiserica"),]
tcf[Abundance.sum %between% c(100,1000),]
tcfo <- tcf[order(PercentAbundance),]

tcfo[PercentAbundance > 10,]
tcfo[PercentAbundance %between% c(2,10),]
tcfo[PercentAbundance %between% c(1,2),]
tcfo[PercentAbundance >= 0.1,]

#, position = position_stack()) +
#  geom_text(position = position_stack()) 

sum(c(1.255094
,1.140103
,1.092831
))

rf$CAZyme[grep('CBM48.CBM68', rf$CAZyme)]
any(rf$CAZyme == 'CBM48.CBM68.GH13_14')
##############################

# Analyse as percentage abundance

rf <- rbindlist(by(rf, INDICES = list(rf$CAZyme, rf$Sample)
                   , FUN = function(x){
  TotAbund <- sum(x$Abundance)
  x$Abundance <- (x$Abundance / TotAbund) * 100
  if(TotAbund == 0){x$Abundance <- 0}
  return(x)
}))

# c("1-Year Fallow"
#   , "10-Year Fallow"
#   , "1-Year Grassland"
#   , "10-Year Grassland")

# Keep only the most abundant cazymes or the ones which were 
# shifted to in the different treatments

AbundantCAZymes <- c("GH15",'GH13_11','GH23','GH3','GH103')

IndicativeCAZymes_b10 <- c('CBM68', 'GH43_32', 'CBM48|CBM68|GH13_14'
                           , 'GH135')
IndicativeCAZymes_gb1 <- c('GH43_4', 'AA10', 'CBM3', 'GH43_35', 'GH101'
, 'GH45', 'CBM8|GH44', 'CBM22|CBM6|GH10', 'GH5_38')

# Subset the data to make plotting easier
rf_ib <- rf[CAZyme %in% IndicativeCAZymes_b10]
rf_ig <- rf[CAZyme %in% IndicativeCAZymes_gb1]
rf_a <- rf[CAZyme %in% AbundantCAZymes]

#### Plot the CAZyme origin data ####

# Plot of origins of CAZymes which were indicative of 10-year fallow treatment 
pb <- ggplot(rf_ib, aes(x = Phylum, y = Abundance, group = Treatment
                        , colour = Treatment)) + 
  stat_summary(fun.data = mean_se, position = position_dodge(0.7)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90
                                   , vjust = 0.5
                                   , hjust = 1)
        , panel.grid = element_blank()
  ) + 
  facet_wrap(~ CAZyme, ncol = 2) + 
  geom_vline(xintercept = seq(1.5, length(unique(rf_ib$Phylum))-0.5, by = 1)
  )#+

pb

# Plot of origins of CAZymes which were indicative of all other treatments
pg <- ggplot(rf_ig, aes(x = Phylum, y = Abundance, group = Treatment
                       , colour = Treatment)) + 
  stat_summary(fun.data = mean_se, position = position_dodge(0.7)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90
                                   , vjust = 0.5
                                   , hjust = 1)
        , panel.grid = element_blank()
        ) + 
  facet_wrap(~ CAZyme, drop = TRUE) + 
  geom_vline(xintercept = seq(1.5, length(unique(rf_ig$Phylum))-0.5, by = 1)
             )#+

pg

# Find the mean value for each treatment
mf <- summaryBy(Abundance ~ Phylum + CAZyme + Treatment, data = rf, FUN = mean
                , keep.names = TRUE)
# Find the range of means for each treatment
mf <- summaryBy(Abundance ~ Phylum + CAZyme, data = mf, FUN = range)
mf <- split(mf, by = "CAZyme")
lapply(mf, FUN = function(x){
  FullTable <- x[order(x$Abundance.FUN1),]
  FullTable[Abundance.FUN1 > 1,]
  
  })




hist(rf[CAZyme == "GH15"]$Abundance)
hist(rf[CAZyme == "GH13_11"]$Abundance)
hist(rf[CAZyme == "GH23"]$Abundance)
hist(rf[CAZyme == "GH3"]$Abundance)
hist(rf[CAZyme == "GH103"]$Abundance)

