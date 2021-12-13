#### Import, check, and set the data up ####
# Clear the workspace
rm(list = ls())

library(data.table)
library(future.apply)
library(roperators)
library(doBy)
library(ggplot2)
library(MASS)
library(lme4)
library(vegan)
library(factoextra)

# Read in the taxonomy data--normalized
df <- data.frame(fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_dbCAN_Taxonomy_ReadAbundance.csv"))

#### Data cleaning and manipulation ####

# Remove viral contigs
df <- df[df$SuperKingdom != "Viruses",]

head(df)

# Create a variable which allows for classification of all unique taxonomyfindings
df$Lineage <- paste(df[,3],df[,4],df[,5],df[,6],df[,7],df[,8], sep = "_")

# Now summarise how many unique taxonomies and CAZymes were found in each phylum
sj <- summaryBy(Lineage + CAZyme ~ Order + Phylum + SuperKingdom, data = df
                , FUN = function(x){
                  #x <- na.omit(x)
                  length(unique(x))
                })
# Make this data nice to look at
# names(sj)[3:4] <- c("Species_Richness", "CAZyme_Richness")
names(sj)[4:5] <- c("Species_Richness", "CAZyme_Richness")

# Deal with NAs
sj$Phylum[is.na(sj$Phylum) & sj$SuperKingdom == 'Archaea'] <- "Unclassified Archaea"
sj$Phylum[is.na(sj$Phylum) & sj$SuperKingdom == 'Bacteria'] <- "Unclassified Bacteria"
sj$Phylum[is.na(sj$Phylum) & sj$SuperKingdom == 'Eukaryota'] <- "Unclassified Eukaryote"
sj$Phylum[is.na(sj$Phylum) & is.na(sj$SuperKingdom)] <- "Unclassified"

head(df)

# CAZyme rich orders
CRO <- c("Bryobacterales"       , "Acidothermales"       , "Acidithiobacillales"                   
         ,"Calditrichales"      , "Jiangellales"         , "Caldilineales"                         
         ,"Fimbriimonadales"    , "Solirubrobacterales"  , "Gloeobacterales"                       
         , "Catenulisporales"   , "Sedimentisphaerales"  , "Ardenticatenales"                      
         , "Glycomycetales"     , "Immundisolibacterales", "Kineosporiales"                        
         , "Limnochordales"     , "Phycisphaerales"      , "Nakamurellales"                        
         , "Rubrobacterales"    , "Sphaerobacterales"    , "Unclassified"                          
         , "Chloroflexales"     , "Deinococcales"        , "Bacteroidetes Order II. Incertae sedis"
         , "Unclassified"       , "Gemmatimonadales"     , "Unclassified"                          
         , "Unclassified"       , "Streptomycetales"     , "Frankiales"                            
         , "Nitrospirales"      , "Opitutales"           , "Geodermatophilales"                    
         , "Caulobacterales"    , "Acidobacteriales"     , "Unclassified"                          
         , "Bifidobacteriales"  , "Micromonosporales"    , "Sphingobacteriales"                    
         , "Chitinophagales"    , "Pseudomonadales"      , "Streptosporangiales"                   
         , "Planctomycetales"   , "Pseudonocardiales"    , "Xanthomonadales"                       
         , "Corynebacteriales"  , "Myxococcales"         , "Sphingomonadales"                      
         , "Propionibacteriales", "Unclassified"         , "Burkholderiales"                       
         , "Rhizobiales")

Many_CAZymes <- df[df$Order %in% CRO,]
head(Many_CAZymes)


mc <- melt(Many_CAZymes, id = c("Contig", 'Species', 'SuperKingdom','Phylum','Class'
                                ,'Order','Family','Genus','CAZyme','Lineage'))

rm(Many_CAZymes)

head(mc)
# Count the number of reads matching to each species
mcs <- summaryBy(value ~ variable + Species + Genus + Family + Order + Class + Phylum + Superkingdom + CAZyme
                 , data = mc
                 , FUN = sum)

# check no CAZymes were found in the negative
sum(mcs$value[mcs$variable == 'n1'])

# check what we have
head(mcs)

# Remove the negative
mcs <- mcs[mcs$variable != 'n1',]
mcs$variable <- factor(mcs$variable)

# Now sort out duplicated species because of CAZymes
mcs <- by(mcs, INDICES = mcs$Species, FUN = function(x){
  
  # Get all CAZymes for each species
  CAZyme <- paste(unique(x$CAZyme), collapse = "_")
  
  # Get sum of all reads
  NormReads <- tapply(x$value.sum, INDEX = x$variable, FUN = function(y){
    sum(y)
  })
  
  # Put out a dataframe of all uniqueinformation for each species in the CAZyme rich orders
  return(data.frame(Sample = unique(x$variable)
                    , Species =  unique(x$Species)
                    , Genus = unique(x$Genus)
                    , Family = unique(x$Family)
                    , Order = unique(x$Order)
                    , Class = unique(x$Class)
                    , Phylum  = unique(x$Phylum)
                    , CAZyme = CAZyme
                    , NormReads = NormReads
  )
  )
  
  
})

mcs <- do.call(rbind, mcs)
rm(mc) # Save memory

# Create a treatment variable for analysis
mcs$Treatment <- NA
mcs$Treatment[grep("b1|b2|b3", mcs$Sample)] <- "b-Old"
mcs$Treatment[grep("b4|b5|b6|b7", mcs$Sample)] <- "b-New"
mcs$Treatment[grep("c1|c2|c3", mcs$Sample)] <- "C-Old"
mcs$Treatment[grep("c4|c5|c6|c7", mcs$Sample)] <- "c-New"

mcs$Treatment <- factor(mcs$Treatment)
mcs$Treatment <- relevel(mcs$Treatment, ref = 'b-Old')

# Now summarise by order--ignore CAZymes, these can be figured out after
## --using CAZyList
mco <- summaryBy(NormReads ~ Sample + Treatment 
                 + Order + Class + Phylum + Superkingdom
                 , data = mcs
                 , FUN = sum
                 , keep.names = TRUE)
mco$Treatment <- factor(mco$Treatment)
mco$Treatment <- relevel(mco$Treatment, ref = 'b-Old')

# Summarise at the family level. How many reads were there per family
mcf <- summaryBy(NormReads ~ Sample + Treatment + Family + 
                   Order + Class + Phylum + Superkingdom
                 , data = mcs
                 , FUN = sum
                 , keep.names = TRUE)
mcf$Treatment <- factor(mcf$Treatment)
mcf$Treatment <- relevel(mcf$Treatment, ref = 'b-Old')



# Which CAZymes were present in each order?
CAZyList <- tapply(mcs$CAZyme, INDEX = mcs$Order, FUN = function(x){
  paste(unique(x), collapse = '-')
})

head(CAZyList)
head(mco)

#### Abundances of each CAZyme rich order in the different treatments ####

ggplot(data = mco
       , aes(x = Order, y = NormReads + 1
             , colour = Treatment)) +
  stat_summary(fun.data = 'mean_se'
               , position = position_dodge(0.9)) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))


#### Which of the CAZyme rich orders had different abundances under different treatments? ####

# Fit models for each of the orders--then compute the adjusted p-value 
# for the explained differences
OrderResults <- c(by(mco, INDICES = mco$Order
                     , FUN = function(x){
                       m <- glm(NormReads ~ Treatment
                                , data = x
                                , family = gaussian()
                       )
                       Result <- drop1(m, test = 'F')
                       return(p.adjust(Result[2,5]
                                       , n = length(unique(mco$Order))
                                       , method = 'BH'))
                     }))

# Fit models for each of the families--then compute the adjusted p-value 
# for the explained differences
FamilyResults <- c(by(mcf, INDICES = mcf$Family
                      , FUN = function(x){
                        m <- glm(NormReads ~ Treatment
                                 , data = x
                                 , family = gaussian())
                        Result <- drop1(m, test = 'F')
                        return(p.adjust(Result[2,5]
                                        , n = length(unique(mcf$Family))
                                        , method = 'BH'))
                      }))

# Fit models for each of the species--then compute the adjusted p-value 
# for the explained differences
SpeciesResults <- c(by(mcs, INDICES = mcs$Species
                       , FUN = function(x){
                         m <- glm(NormReads ~ Treatment
                                  , data = x
                                  , family = gaussian())
                         Result <- drop1(m, test = 'F')
                         return(
                           p.adjust(Result[2,5]
                                    , n = length(unique(mcs$Species))
                                    , method = 'BH'))
                       }))

# View the distribution of the p-values for each of the taxonomic ranks analysed
hist(OrderResults, breaks = 200);abline(v = 0.05, lty = 2)
hist(FamilyResults, breaks = 200);abline(v = 0.05, lty = 2)
hist(SpeciesResults, breaks = 200);abline(v = 0.05, lty = 2)

# Keep only significant results
SignificantSpecies <- SpeciesResults[SpeciesResults < 0.05]
SignificantFamilies <- FamilyResults[FamilyResults < 0.05]
SignificantOrders <- OrderResults[OrderResults < 0.05]

# Which orders, families, and species had significantly 
# different "numbers" of reads mapping to their contigs?
# (And what are the adjusted p-values 
SignificantSpecies 
SignificantFamilies
SignificantOrders

# Look at the differences in mean numbers of reads mapping back 
# to each order/family/species with high degradative potential
# which responded significantly

## order
# Colour palette
BlackoutPalette <- c('#c2a5cf'
                     ,'#7b3294'
                     ,'#a6dba0'
                     ,'#008837')
OrderPlotData <- mco[mco$Order %in% names(SignificantOrders),]
OrderPlotData$Phylum <- factor(OrderPlotData$Phylum)
OrderPlotData$Order <- factor(OrderPlotData$Order)

ggplot(data = OrderPlotData
       , aes(x = Order, y = NormReads
             , colour = Treatment)
) +
  geom_point(alpha = 0.3, position = position_dodge(0.5), shape = 16) +
  
  stat_summary(fun.data = 'mean_se'
               , position = position_dodge(0.5), shape = 17, size = 1) +
  
  theme_classic() +
  theme(text = element_text(size = 16)) +
  ylab(expression("Normalized Read Abundance (Mean  " %+-% " SE)")) +
  scale_color_manual(values = BlackoutPalette
                     , labels = c("Fallow Old", "Fallow New"
                                  , "Grassland Old", "Grassland New")) +
  
  facet_wrap(facets = vars(Phylum), scales = 'free_x' 
  )


OrderPlotData
head(df)

# Maybe relevel this analysis to check if fallow new is also different?
# OrderResults2 <- c(by(mco, INDICES = mco$Order
#                      , FUN = function(x){
#                        m <- glm(NormReads ~ Treatment
#                                 , data = x
#                                 , family = gaussian()
#                        )
#                        Result <- drop1(m, test = 'F')
#                        return(p.adjust(Result[2,5]
#                                        , n = length(unique(mco$Order))
#                                        , method = 'BH'))
#                      }))

#### Which CAZymes did the differentially abundant CAZy rich orders have? ####

# How many CAZymes did each order have?
# Actinobacteria
df[df$Order == "Propionibacteriales",]

length(unique(df[df$Order == "Propionibacteriales",]$CAZyme)) -1 # Remove the NA
# How many AAs
length(grep("AA", unique(df[df$Order == "Propionibacteriales",]$CAZyme)))
length(grep("GH", unique(df[df$Order == "Propionibacteriales",]$CAZyme)))
length(grep("CBM", unique(df[df$Order == "Propionibacteriales",]$CAZyme)))

length(unique(df[df$Order == "Nakamurellales",]$CAZyme)) -1 # Remove the NA
# How many AAs
length(grep("AA", unique(df[df$Order == "Nakamurellales",]$CAZyme)))
length(grep("GH", unique(df[df$Order == "Nakamurellales",]$CAZyme)))
length(grep("CBM", unique(df[df$Order == "Nakamurellales",]$CAZyme)))

length(unique(df[df$Order == "Acidobacteriales",]$CAZyme)) -1 # Remove the NA
# How many AAs
length(grep("AA", unique(df[df$Order == "Acidobacteriales",]$CAZyme)))
length(grep("GH", unique(df[df$Order == "Acidobacteriales",]$CAZyme)))
length(grep("CBM", unique(df[df$Order == "Acidobacteriales",]$CAZyme)))

length(unique(df[df$Order == "Bryobacterales",]$CAZyme)) -1 # Remove the NA
# How many AAs
length(grep("AA", unique(df[df$Order == "Bryobacterales",]$CAZyme)))
length(grep("GH", unique(df[df$Order == "Bryobacterales",]$CAZyme)))
length(grep("CBM", unique(df[df$Order == "Bryobacterales",]$CAZyme)))

#### Bit more data cleaning -- how did CAZyme composition change between treatments ####

# Group by species, keep all individual CAZymes, see changes in total 
# read numbers mapping to contigs from species with AAs, GHs and CBMs

# Keep the Orders we care about
rf <- df[df$Order %in% c("Propionibacteriales", "Bryobacterales", "Acidobacteriales", "Nakamurellales"),]

# Check the data
head(rf)
max(rf$n1)

# Get all the reads from a sample which map to a particular species
TaxAbundance <- summaryBy(b1 + b2 + b3 + b4 + b5 + b6 + b7 +
                            c1 + c2 + c3 + c4 + c5 + c6 + c7 ~
                            Species + Genus + Family + Order + Class + Phylum + Superkingdom
                          , data = rf
                          , FUN = sum
                          , keep.names = TRUE)

head(TaxAbundance)

# Put it in a useable format
sf <- melt(TaxAbundance[,c(1, 4, 6, 7:length(TaxAbundance))], id = c('Species', "Phylum", "Order"))
names(sf)[4:5] <- c("Sample", "Abundance")
sf$Treatment <- NA
sf[sf$Sample %in% c("b1", "b2", "b3"),]$Treatment <- "Fallow-Old"
sf[sf$Sample %in% c("b4", "b5", "b6", "b7"),]$Treatment <- "Fallow-New"
sf[sf$Sample %in% c("c1", "c2", "c3"),]$Treatment <- "Grassland-Old"
sf[sf$Sample %in% c("c4", "c5", "c6", "c7"),]$Treatment <- "Grassland-New"
sf$Treatment <- factor(sf$Treatment)

head(sf)

# Now deal with the CAZymes
# 1) Which CAZymes were in which genomes
CAZymeList <- by(rf, INDICES = rf$Species, FUN = function(x){unique(na.omit(x$CAZyme))})

# 2) Get all of the CAZymes in the genomes of these organisms
AllCAZymes <- unique(unlist(CAZymeList))

# If the CAZyme is present in the genome give 1, if no give 0. Whack it all end-on-end
CZf <- data.frame(do.call(rbind, lapply(CAZymeList, FUN = function(x){
  ifelse(AllCAZymes %in% x, yes = 1, no = 0)
})))
# And name the columns sensibly
names(CZf) <- AllCAZymes

# Czech it aaaht
head(CZf)
CZf$Species <- row.names(CZf)

hist(sf$Abundance, breaks = 1000000, xlim = c(0,10000))

row.names(CZf) == sf$Species


nrow(sf)/ nrow(CZf)

sf <- merge(sf, CZf, by = "Species", all = TRUE)

head(sf)

#### how did CAZyme composition change between treatments ####

# How many species were in each of the differentially abundant orders?
for(i in unique(sf$Order)){print(paste0(i, " - ", length(unique(sf[sf$Order == i,]$Species))))}

# What is the one Acidobacteria in Bryobacteriales with 63? CAZymes
unique(sf[sf$Order == 'Bryobacterales',]$Species) 

# Make the CAZymes as abundant as the species
sf[,7:length(sf)] <- sf[,7:length(sf)] * sf$Abundance

# Now for each order get the abundance of certain CAZymes in each sample
cf <- by(sf[,7:length(sf)], INDICES = list(sf$Order, sf$Treatment, sf$Sample), FUN = colSums)

any(is.na(unlist(sf[7:length(sf)])))
any(rowSums(sf[7:length(sf)]) == 0)
which(rowSums(sf[7:length(sf)]) == 0)


#nmds1 <- metaMDS(sf[7:length(sf)])
pca1 <- prcomp(log2(sf[7:length(sf)]+1)#, scale. = TRUE, center = TRUE
)

plot(pca1)

# fviz_pca(pca1
#          , habillage = paste(sf$Order,sf$Treatment)
#          , addEllipses = TRUE
#          , ellipse.type = "convex"
#          , mean.point = FALSE
#          , label = 'var'
#          , select.var = list(contrib = 25)
#          #, xlim = c(-10,10)
# )

# Visualise how different factors affect CAZyme composition.

# RColorBrewer::brewer.pal(8, "Paired") # Chose the distinct hard colours from the paired palette

# Which factors influenced CAZyme community composition?

# Statistical test
PERMANOVA_Results <- adonis(log2(sf[7:length(sf)]+1) ~ sf$Treatment * sf$Order
                            , method = 'euclidean'
                            , permutations = 10000)
PERMANOVA_Results

# Make the plots


p1 <- fviz_pca(pca1
         , habillage = paste(sf$Order)
         , addEllipses = TRUE
         , ellipse.type = "convex"
         , mean.point = FALSE
         , label = 'var'
         , select.var = list(contrib = 30)
         , palette = OrderPalette
         , ggtheme = theme_classic()
         , title = "a)"
         ) +
  geom_text(mapping = aes(x = -55, y = 90
                          , label = "Treatment x Order: p = 1")
            , size = 6) +
  geom_text(mapping = aes(x = -55, y = 80
                          , label = "Treatment: p = 1")
            , size = 6) +
  geom_text(mapping = aes(x = -55, y = 70
                          , label = "Order: p = 0.001")
            , size = 6) +
  theme(legend.position = "bottom"
        ,legend.text = element_text(size = 15)
        , axis.text = element_text(size = 15)
        , axis.title = element_text(size = 15)
        , title = element_text(size = 15)) +
  labs(color = "Order" 
      , shape = "Order"
      , fill = "Order") +
  guides(color = guide_legend(title.position = "top"
                              , nrow = 2, byrow = TRUE)
         , shape = guide_legend(title.position = "top"
                                , nrow = 2, byrow = TRUE)
         , fill = guide_legend(title.position = "top"
                               , nrow = 2, byrow = TRUE))
p1

p2 <- fviz_pca(pca1
         , habillage = paste(sf$Treatment)
         , addEllipses = TRUE
         , ellipse.type = "convex"
         , mean.point = FALSE
         , label = 'var'
         , select.var = list(contrib = 30)
         , palette = BlackoutPalette
         , ggtheme = theme_classic()
         , title = "b)"
         ) +
  theme(legend.position = "bottom"
        , legend.text = element_text(size = 15)
        , axis.title = element_text(size = 15)
        , title = element_text(size = 15)) +
  labs(color = "Treatment" 
       , shape = "Treatment"
       , fill = "Treatment") +
  guides(color = guide_legend(title.position = "top"
                              , nrow = 2, byrow=TRUE)
         , shape = guide_legend(title.position = "top"
                                , nrow = 2, byrow=TRUE)
         , fill = guide_legend(title.position = "top"
                               , nrow = 2, byrow=TRUE))

p2

# CAZyme_Community_Composition
cowplot::plot_grid(p1, p2)


