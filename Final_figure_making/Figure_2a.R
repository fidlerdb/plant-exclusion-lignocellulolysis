#### (1) Taxonomy NMDS plot ####
rm(list = ls())
library(vegan)
library(scales)


# All Contigs
dfa <- data.frame(fread("Taxonomy/Multi-Method/Cleaned_Data/Henfaes_Phylum_Multi-Method_NormalizedReadAbundance.csv"))
dfa <- dfa[, !names(dfa) %in% "UNASSIGNED"] # Remove unassigned contigs as they will belong to many phyla
dfa <- dfa[, -grep("Unknown", names(dfa))] # Remove Domain-level classifications

# Percentage/proportion abundance
dfa.rel <-  decostand(dfa[,4:length(dfa)], method = "total")
nmds1 <- vegan::metaMDS(dfa.rel, distance = 'bray')
# Create useful vectors
BlackoutPalette <- c('#7b3294','#c2a5cf','#a6dba0','#008837')
Treatment <- factor(dfa$Treatment)
levels(Treatment) <- c("1-Year Bare", "10-Year Bare", "1-Year Grassland", "10-Year Grassland")
Treatment <- relevel(Treatment, ref = "10-Year Bare")
colours <- BlackoutPalette[Treatment]

source("Functions/range01.R")
magnitude <- function(x){x <- sqrt(x^2); return(x)}

speciesAlpha <- range01(log(magnitude(nmds1$species[,1]) + magnitude(nmds1$species[,2])))
#hist(speciesAlpha)
#speciesAlpha <- range01(magnitude(nmds1$species[,1]) + magnitude(nmds1$species[,2])) 

dev.off()
pdf(file = "Final_figure_making/To_Put_in_Inkscape/Figure_2a_Taxonomic_composition.pdf"
    , height =  4, width = 5)

# Create the plot
plot(nmds1$points, type = 'n'
     , xlim = c(-0.17, 0.07)
     #, ylim = c(-0.04, 0.06)
     , bty = 'L'
)

orditorp(nmds1, display = "species"
         , col = alpha("black", speciesAlpha)
         , air = 0.1)

points(nmds1$points
       , pch = 16
       , col = colours)
# Add convex hulls
for(i in levels(Treatment)){
  ordihull(nmds1$points[grep(i, Treatment),]
           , groups = Treatment[Treatment == i]
           , draw = "polygon"
           , col = colours[grep(i, Treatment)]
           , border = FALSE
           , alpha = 0.7
  )
}

text(x = 0.04, y = 0.05
     , labels = paste0("Stress = ",round(nmds1$stress, 3))
     , cex = 0.8)

dev.off()
# legend("topleft", legend = levels(Treatment)
#        , col = BlackoutPalette
#        , pch = 16
#        , cex = 0.8
#        , bty = "n") # Legend in second plot should work fine

