rm(list = ls())

library(ggplot2)
library(data.table)
library(multcomp)

BlackoutPalette <- c('#c2a5cf'
                     , '#7b3294'
                     , '#a6dba0'
                     , '#008837')

#df <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZymes_NotGrouped_ReadsPerKbp.csv")

#####

# CAZy family read data
df <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZyFamilyReadsPerKbp.csv")
df[, 4:length(df)] <- df[, 4:length(df)] / (df$Reads/1e9) 

# List of CAZy families which were responsive to treatment
ResponsiveCAZymes <- read.csv("CAZy/Cleaned_Data/ResponsiveCAZymes_FurtherExploration.csv")

#####

# Keep only the CAZymes which had 
# large significant changes in richness

rf <- cbind(df[,1:4]
            , df[,names(df) %in% ResponsiveCAZymes$x]
)

head(rf)

# Set up the data for a ggplot2 bar chart
rfl <- melt(rf, id = c("Treatment", "Sample", "Reads")
            , variable.name = "CAZyme"
            , value.name = "Abundance"
     )


head(rfl)

#### Did treatment cause a difference in read abundance ####
####  for the responsive CAZy families   
####
CAZymeResults <- c(by(rfl, INDICES = rfl$CAZyme
                       , FUN = function(x){
                         m <- glm(Abundance ~ Treatment
                                  , data = x
                                  , family = gaussian())
                         Result <- drop1(m, test = 'F')
                         return(
                           p.adjust(Result[2,5]
                                    , n = length(unique(rfl$CAZyme))
                                    , method = 'BH'))
                       }))

CAZymeResults[CAZymeResults < 0.05]
SignificantCAZymes <- names(CAZymeResults[CAZymeResults <= 0.05])

rfl2 <- rfl[rfl$CAZyme %in% SignificantCAZymes,]
rfl2$CAZyme <- factor(rfl2$CAZyme)
levels(rfl2$Treatment) <- list(`1-Year Fallow` = "b-New"
                               , `10-Year Fallow` = "b-Old"
                               , `1-Year Grassland` = "c-New"
                               , `10-Year Grassland` = "c-Old")

TukeyLabels <- by(rfl2, INDICES = rfl2$CAZyme, FUN = function(x){
  m <- glm(Abundance ~ Treatment
           , data = x
           , family = gaussian())
  
  TukeyLabels <- cld(glht(m, linfct = mcp(Treatment = "Tukey"))
                     , level = 0.05)
  TukeyLabels <- data.frame(Treatment = names(TukeyLabels$mcletters$Letters)
                            , Letters = TukeyLabels$mcletters$Letters
                            , Abundance = max(x$Abundance) + 0.1*max(x$Abundance)
                        )
  })

TukeyLabels <- rbindlist(TukeyLabels, idcol = "CAZyme")

#### Misc things for plotting ####
mean_ci <- function(x){
  mean_se(x, mult = 1.96)
}


# Adjust colour values of the palette
#devtools::install_github("briandconnelly/colormod")
library(raster)
library(colormod)

BlackoutPalette2 <- adjust_hsv(col = BlackoutPalette
                               , v = -0.25)

#### Plot the data ####

p <- ggplot(rfl2, aes(x = Treatment, y = Abundance
                , fill = Treatment
                , colour = Treatment)) +
  
  # Add mean and errorbars
  stat_summary(fun.y = "mean"
               #, position = position_dodge(1.2)
               , geom = "bar") +
  stat_summary(fun.data = "mean_cl_boot"
               #, position = position_dodge(0.9)
               , geom = "errorbar"
               , colour = 'black'
               , width = 0.3) +
  # Show raw data
  geom_point(position = position_jitterdodge(0.9)
             #, alpha = 0.3
             , size = 1.2) +

  # Prettiness in ggplot2-shitty format
  theme_classic() + 
  theme(legend.position = "bottom"
        , axis.text.x = element_blank()
        , strip.text = element_text(size = 10)
        ) + 
  scale_colour_manual(breaks = levels(rfl2$Treatment)
                      , values = BlackoutPalette2) +
  scale_fill_manual(breaks = levels(rfl2$Treatment)
                      , values = BlackoutPalette) +
  ylab("Normalized read abundance") +

  # Give every CAZy family its own plot
  facet_wrap(~ CAZyme, scales = "free_y") +
  geom_text(data = TukeyLabels
            , aes(label = Letters
                  , x = Treatment
                  , y = Abundance), colour = 'black')
#

p

ggsave("Figures/CAZy/Differentially_Rich_CAZyme_Abundance.pdf", plot = p
       , device = "pdf", width = 6.95, height = 5.35, units = "in")
