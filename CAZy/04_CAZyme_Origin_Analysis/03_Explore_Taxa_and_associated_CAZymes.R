# Summarize and explore the link between taxonomy and lignocellulase genes

# Clear the workspace
rm(list = ls())

library(data.table)
library(future.apply)
library(roperators)
library(doBy)
library(ggplot2)
library(MASS)
library(lme4)

#### Import, check, and set the data up ####

# Read in the taxonomy data--normalized
df <- data.frame(fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_dbCAN_Taxonomy_ReadAbundance.csv"))

# See how this looks
head(df)

# Keep only contigs which have both CAZy genes and a taxonomic assignment

df <- df[!is.na(df$CAZyme),]
head(df)

# Remove viral contigs
df <- df[df$SuperKingdom != "Viruses",]

# Which phyla had CAZymes, and how many?
sf <- summaryBy(CAZyme ~ Phylum + SuperKingdom, data = df, FUN = function(x)length(unique(x)))
sf <- sf[order(sf$`CAZyme.function(x) length(unique(x))`),]
names(sf)[3] <- "CAZymes"
sf$Phylum[is.na(sf$Phylum)] <- paste("Unknown", sf$SuperKingdom[is.na(sf$Phylum)], sep = " ")
sf$Phylum[sf$Phylum == "Unknown NA"] <- "Unknown"
head(sf)
sf <- sf[-grep("Unknown", sf$Phylum),] # Remove "unknown whatever"

# What CAZymes did each phylum have?
by(data = df, INDICES = df$Phylum
   , FUN = function(x){
     g <- x$CAZyme
     return(unique(g))
     }
   , simplify = TRUE)

#### What was the total CAZy Family Richness for each Phylum? ####

# Make plot of total CAZy family richness for each phylum

# Palette
DomainPalette <- c('#e41a1c'# Archaea
                   , '#377eb8'# Bacteria
                   , '#4daf4a'# Eukaryotes
                   )
levels(factor(sf$SuperKingdom))

library(ggplot2)
ggplot(sf, aes(x = reorder(Phylum, CAZymes), y = log10(CAZymes)
               , fill = SuperKingdom)) + 
  geom_bar(stat = 'identity', position = position_dodge(0.9)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  scale_fill_manual(values = DomainPalette) +
  ylab(expression("Log"[10]*"(CAZy Family Richness)")) +
  xlab("Phylum") + coord_flip()

sf[sf$CAZymes >= 50,]
sf[sf$Phylum == "Nitrospirae",]

##### How many reads mapped to each of the phyla in each treatment? ####

# Wrangle the data for a plot exploring this:
# Summarise total number of reads by phylum, then take mean +- 95% CI values

jf <- melt(summaryBy(b1+b2+b3+b4+b5+b6+b7+ 
          c1+c2+c3+c4+c5+c6+c7 ~  Phylum
          , data = df, FUN = sum, keep.names = TRUE), id = c("Phylum"))
names(jf)[2:3] <- c("Sample", "CAZymes")
jf$Treatment <- paste(ifelse(grepl("c", jf$Sample), yes = "Grassland", no = "Fallow")
      , ifelse(grepl("1|2|3", jf$Sample), yes = "Old", no = "New"), sep = "-")
jf$Phylum <- as.character(jf$Phylum)
jf$Phylum[is.na(jf$Phylum)] <- "Unknown"
source("Functions/boot_ci.R")
jf <- na.omit(jf)
jf2 <- summaryBy(CAZymes ~ Treatment + Phylum, data = jf, FUN = c(mean,boot_ci))
jf2
# summarise the data further and make weird bits for the plot:
# Abundance categories
jf3 <- tapply(jf2$CAZymes.ans, jf2$Phylum, FUN = mean)
Rare <- names(jf3[jf3 < 500])
Abundant <- names(jf3[jf3 > 500 & jf3 < 5000])
Dominant <- names(jf3[jf3 > 5000])

jf2$Abundance <- NA
jf2$Abundance[jf2$Phylum %in% Rare] <- "Rare"
jf2$Abundance[jf2$Phylum %in% Abundant] <- "Abundant"
jf2$Abundance[jf2$Phylum %in% Dominant] <- "Dominant"

jf2$Abundance <- factor(jf2$Abundance
                        , levels = c("Dominant", "Abundant", "Rare"))

# Nicer name
names(jf2)[3] <- "CAZymes.mean"

# Colour palette
BlackoutPalette <- c('#c2a5cf'
                     ,'#7b3294'
                     ,'#a6dba0'
                     ,'#008837')

# Axis start values for each facet
jf2$y_start <- NA
jf2$y_start[jf2$Abundance == "Rare"] <- 1
jf2$y_start[jf2$Abundance == "Abundant"] <- 1
jf2$y_start[jf2$Abundance == "Dominant"] <- 1e3

# Nicer name for a plot element
jf2$Phylum[jf2$Phylum == "Deinococcus-Thermus"] <- "D-T"

# Create a dataframe of x values at which to add vertical lines
lineframe <- data.frame(
  Abundance = c(rep("Dominant", 13), rep("Abundant", 10)
                , rep("Rare", 10))
  , linevals = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5
                 , 1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5
                 , 1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5))

# Extra point for each facet so that you can actually see the data...
ExtraPoints <- data.frame(
  Abundance = c("Dominant", "Abundant","Rare")
  , Phylum = c("Acidobacteria", "Chlorobi", "Fibrobacteres")
  , CAZymes.mean = c(3e6, 20000, 1e4)
  , y_start = c(1e5, 1, 1)
  , Treatment = factor("Fallow-old", levels = levels(jf2$Treatment))
  )


# Finally... 

# Create the basic plot structure
ggplot(jf2, aes(x = reorder(Phylum, CAZymes.mean)
                , ymax = CAZymes.mean + 1
                , ymin = y_start
                , fill = Treatment
                , colour = Treatment)
       ) + 
  
  # We want to add lines as segments (because Hadley is an idiot)
  geom_linerange(position = position_dodge(0.9), size = 4) +
  
  # Make the y axis on a log10 scale, and put the bars all 
  # the way to the bottom
  scale_y_continuous(trans = 'log10'
                     , expand = c(0,0)) +
  
  # Put abundant and rarer community members etc. in different boxes
  facet_wrap(. ~ Abundance, scales = "free"
             , nrow = 3) +
  
  # Add errorbars--these are bootstrapped CIs. Gaussian CIs suggest 
  # negative numbers of  reads mapping to some phyla
  geom_errorbar(data = jf2
                , mapping = aes(ymax = CAZymes.u_ci + 1
                                , ymin = CAZymes.l_ci + 1 )
                , width = 0.1
                , position = position_dodge(0.9),
                colour = 'black') +
  
  # Do some pretties
  theme_classic() +
  theme(panel.border = element_rect(fill = NA)) +
  scale_fill_manual(values = BlackoutPalette) + 
  scale_colour_manual(values = BlackoutPalette) + 
  ylab("Normalized Read Abundance Mapping to Contigs which Contain CAZymes + 1") +
  xlab("Phylum") +
  
  # Add lines which seperate each phylum from each other--makes it easier to look at
  geom_vline(data = lineframe, aes(xintercept = linevals))  +
  
  # Now add dummy points so that you can see the values of each CAZyme
  geom_point(data = ExtraPoints, aes(x = Phylum, y = CAZymes.mean))

#### Correlate Richness with Abundance ####

#### Import, check, and set the data up ####
rm(list = ls()) # Clear the workspace
# Read in the taxonomy data--normalized
df <- data.frame(fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_dbCAN_Taxonomy_ReadAbundance.csv"))

# Remove viral contigs
df <- df[df$SuperKingdom != "Viruses",]

head(df)
unique(df$CAZyme)

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

# Quick and dirty plot--NB some low diversity taxa have > 100 CAZymes i.e. log 2
plot(CAZyme_Richness ~ 
     Species_Richness
     , data = sj
     #, pch = 16
     , type = 'n'
     )
text(CAZyme_Richness ~ Species_Richness
     , data = sj
     , labels = sj$Order
     , cex = 0.8)


# Fit a model and then plot it properly highlighting the
# orders with the most degradative potential. First tried lm, awful
# Then tried poisson, awful, then tried negative binomial, less bad.
# Tried negative binomial lmer with RE of SuperKingdom/Phylum. Good.
# Getting rid of Superkingdom did not change the result so went with 
# the more parsimonious option. Implement this now

# Scaled x variables allow model convergence
sj$Log_SR <- log10(sj$Species_Richness)

# Fit the model
lme_nb <- glmer.nb(CAZyme_Richness ~ Log_SR + (1|Phylum) #* SuperKingdom
      , data = sj
      )

plot(lme_nb) # The model fits well
car::qqPlot(residuals(lme_nb, type = 'deviance')) # The model fits well
drop1(lme_nb, test = 'Chisq') # Terms are siggy
summary(lme_nb)
#### Predictions ####

# Try for a prediction interval
newdat <- data.frame(Log_SR=log10(1:74))
mm <- model.matrix(~ Log_SR, newdat)
y <- mm %*% fixef(lme_nb)
pvar1 <- diag(mm %*% tcrossprod(vcov(lme_nb),mm))
tvar1 <- pvar1+VarCorr(lme_nb)$Phylum[1]#+VarCorr(lme_nbv)$f2[1]  ## must be adapted for more complex models
newdat <- data.frame(
  Log_SR=newdat$Log_SR
  , Species_Richness = 10^newdat$Log_SR # Actual species richness values
  , y=exp(y) # Predicted values
  , plo = exp(y-1.96*sqrt(pvar1)) # Lower 95% CI
  , phi = exp(y+1.96*sqrt(pvar1)) # Upper 95% CI
  , tlo_95 = exp(y-1.96*sqrt(tvar1)) # Lower 95% prediction interval
  , thi_95 = exp(y+1.96*sqrt(tvar1)) # Upper 95% prediction interval
  , tlo_68 = exp(y-sqrt(tvar1)) # Lower 1 SD prediction interval
  , thi_68 = exp(y+sqrt(tvar1)) # Upper 1 SD prediction interval
)



#### Make a pretty picture of this! :D ####
# Name unidentified Orders
sj[is.na(sj$Order),]$Order <- "Unclassified"

# Orders above the prediction interval
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

# Stick prediction interval and observed data together
I_N <- merge(sj, newdat, id = "Species_Richness")
H <- I_N[I_N$CAZyme_Richness > I_N$thi_68,]

# Make the basic plot
plot(CAZyme_Richness ~ Species_Richness
     , data = sj
     , pch = 16
     , type = 'n'
     , ylim = c(0, 220)
     , xlim = c(-5, 80)
     , bty = 'l'
     , xlab = "Species Richness"
     , ylab = "CAZy Family Richness"
)

# Add model predictions
lines(y ~ Species_Richness, data = newdat, col = 'black', lwd = 2.5)

# Add a prediction interval (1 SD)
with(newdat,
     polygon(x = c(Species_Richness, rev(Species_Richness))
             , y = c(tlo_68, rev(thi_68))
             , border = NA
             , col = scales::alpha('black', 0.2))
)

# Make subsetting look nicer
BL <- sj[sj$Order %ni% CRO,]
BL$SuperKingdom <- factor(BL$SuperKingdom)
points(CAZyme_Richness ~ Species_Richness
     , data = BL
     #, labels = paste(Phylum, Order, sep = "-")
     , cex = 0.8
     , pch = 16
     # , col = scales::alpha(
     #   c('#DE5E4E' # Archaea
     #     , '#0C7BDC' # Bacteria
     #     , '#37802F')
     #   , 0.8)[factor(BL$SuperKingdom)]
     , col = scales::alpha(1,0.5)
)

AC <- sj[sj$Order %in% CRO[CRO != 'Unclassified'],]

text(CAZyme_Richness ~ Species_Richness
     , data = AC#sj[sj$Order %in% CRO[CRO != 'Unclassified'],]
     , labels = Order# paste(Phylum, Order, sep = "-")
     , cex = 0.6
     # , col = c('#DE5E4E' # Archaea
     #           , '#0C7BDC' # Bacteria
     #           , '#37802F' # Eukarya
     #           )[factor(AC$SuperKingdom, levels = levels(BL$SuperKingdom))]
     )  

AU <- H[H$Order == 'Unclassified',]
text(CAZyme_Richness ~ Species_Richness
     , data = AU[AU$Phylum != 'Unclassified Bacteria',]
     , labels = paste(Order, Phylum, sep = " ")
     , cex = 0.6
     # , col = c('#DE5E4E' # Archaea
     #           , '#0C7BDC' # Bacteria
     #           , '#37802F' # Eukarya
     #)[factor(AU$SuperKingdom, levels = levels(BL$SuperKingdom))]
)

text(CAZyme_Richness ~ Species_Richness
     , data = AU[AU$Phylum == 'Unclassified Bacteria',]
     , labels = Phylum
     , cex = 0.6
)

#### Woo! Now find out which Orders had what CAZymes ####

nrow(I_N[I_N$CAZyme_Richness > I_N$thi_68,])

write.csv(I_N[I_N$CAZyme_Richness > I_N$thi_68,]
          , "CAZy/04_CAZyme_Origin_Analysis/Results/CAZyme_Rich_Orders.csv"
          , row.names = FALSE)

# Orders which contained a single species
I_N[I_N$CAZyme_Richness > I_N$thi_68 & I_N$Species_Richness == 1,]
nrow(I_N[I_N$CAZyme_Richness > I_N$thi_68 & I_N$Species_Richness == 1,])
nrow(I_N[I_N$CAZyme_Richness > I_N$thi_68 & I_N$Species_Richness == 2,])
nrow(I_N[I_N$CAZyme_Richness > I_N$thi_68 & I_N$Species_Richness == 3,])
nrow(I_N[I_N$CAZyme_Richness > I_N$thi_68 & I_N$Species_Richness == 4,])
nrow(I_N[I_N$CAZyme_Richness > I_N$thi_68 & I_N$Species_Richness == 5,])

# Orders which contained multiple species
I_N[I_N$CAZyme_Richness > I_N$thi_68 & I_N$Species_Richness > 4,]

sort(-table(I_N[I_N$CAZyme_Richness > I_N$thi_68,]$Phylum))


H[H$Phylum == "Proteobacteria",]

paste(H$Order)

H[H$Order == 'Unclassified' ,]

H2 <- summaryBy(Species_Richness ~ Phylum, data = H, FUN = sum)

sort(-table(H$Phylum)) # How many CAZyme rich orders were in each phylum
sort(-table(H[H$Species_Richness == 1,]$Phylum)) # How many single-species orders were in each phylum

H2[order(-H2$Species_Richness.sum),]

unique(sj$Species_Richness)

save.image(file = "CAZy/04_CAZyme_Origin_Analysis/RData_Files/03_Explore_Taxa_and_CAZymes_Enrivonment.RData")

