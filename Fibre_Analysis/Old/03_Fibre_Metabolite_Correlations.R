
rm(list = ls())

ff <- fread("Fibre_Analysis/Output_Data/Fibre_Analysis_Data_Processed.csv")
ff <- ff[,c("Sample", "Percentage_NDF",  "Percentage_ADF", "Treatment")  ]

mf <- fread("Metabolomics/Data/Metabolome_Blackout_Raw.csv")

head(mf)
mf <- mf[, c(1,7:20)] # Keep only the most useful columns
mf <- mf[`BinBase name` %in% c("glucose"
                         , "xylose", "fucose", '3,6-anhydro-D-galactose'
                         , 'vanillic acid', '4-hydroxybenzoic acid', 'benzoic acid'),]
head(mf)           

# Transpose this so each metabolite has its own column
mft <- data.frame(t(mf)[-1,]) 
names(mft) <- mf$`BinBase name`
mft <- data.table(mft, keep.rownames = TRUE)
names(mft)[1] <- "Sample"




mft

af <- merge(ff, mft, by = "Sample")
af$Treatment <- factor(af$Treatment)
af$xylose           <- as.numeric(af$xylose          )
af$`vanillic acid`  <- as.numeric(af$`vanillic acid` )
af$glucose          <- as.numeric(af$glucose         )
af$fucose           <- as.numeric(af$fucose          )
af$`benzoic acid`   <- as.numeric(af$`benzoic acid`  )
af$`4-hydroxybenzoic acid` <- as.numeric(af$`4-hydroxybenzoic acid`)
af$`3,6-anhydro-D-galactose` <- as.numeric(af$`3,6-anhydro-D-galactose`)


# Do metabolites and fibre content correlate?

## Cellulose-- need to change predictor variable
cor_glucose <- cor.test(x = af$Percentage_ADF, y = af$glucose)
plot(glucose ~ Percentage_ADF, data = af, pch = 16, col = af$Treatment
     , main = 'glucose') # Not for this one--mostly expected
legend("topright", legend = levels(af$Treatment), pch = 16, col = 1:4)
abline(lm(glucose ~ Percentage_ADF, data = af))
text(x = 50, y = 14000, labels = paste0("r = ", round(cor_glucose[[4]],2)
                                        , ", t = ", round(cor_glucose[[1]],2)
                                        , "\ndf = ", cor_glucose[[2]]
                                        , ", p = ", round(cor_glucose[[3]],2)
))

## Hemicellulose -- need to change predictor variable
cor_xylose <- cor.test(x = af$Percentage_NDF, y = af$xylose)
plot(xylose ~ Percentage_NDF, data = af, pch = 16, col = af$Treatment
     , main = "xylose")
legend("topright", legend = levels(af$Treatment), pch = 16, col = 1:4)
abline(lm(xylose ~ Percentage_NDF, data = af))
text(x = 65, y = 20000, labels = paste0("r = ", round(cor_xylose[[4]],2)
                                        , ", t = ", round(cor_xylose[[1]],2)
                                        , "\ndf = ", cor_xylose[[2]]
                                        , ", p = ", round(cor_xylose[[3]],2)
                                        ))

cor_fucose <- cor.test(x = af$Percentage_NDF, y = af$fucose)
plot(fucose ~ Percentage_NDF, data = af, pch = 16, col = af$Treatment
     , main = "fucose")
legend("topright", legend = levels(af$Treatment), pch = 16, col = 1:4)
abline(lm(fucose ~ Percentage_NDF, data = af))
text(x = 68, y = 10000, labels = paste0("r = ", round(cor_fucose[[4]],2)
                                        , ", t = ", round(cor_fucose[[1]],2)
                                        , "\ndf = ", cor_fucose[[2]]
                                        , ", p = ", round(cor_fucose[[3]],2)
))

cor_galactose <- cor.test(x = af$Percentage_NDF, y = af$`3,6-anhydro-D-galactose`)
plot(`3,6-anhydro-D-galactose` ~ Percentage_NDF, data = af, pch = 16, col = af$Treatment
     , main = "3,6-anhydro-D-galactose") # Not really for this one
legend("topright", legend = levels(af$Treatment), pch = 16, col = 1:4)
abline(lm(`3,6-anhydro-D-galactose` ~ Percentage_NDF, data = af))
text(x = 68, y = 750, labels = paste0("r = ", round(cor_galactose[[4]],2)
                                        , ", t = ", round(cor_galactose[[1]],2)
                                        , "\ndf = ", cor_galactose[[2]]
                                        , ", p = ", round(cor_galactose[[3]],2)
))



## lignin -- need to change predictor variable
plot(Percentage_ADF ~ `vanillic acid`, data = af, pch = 16, col = af$Treatment)


plot(Percentage_ADF ~ `benzoic acid`, data = af, pch = 16, col = af$Treatment)


plot(Percentage_ADF ~ `4-hydroxybenzoic acid`, data = af, pch = 16, col = af$Treatment)



af$Proportion_NDF <- af$Percentage_NDF/100
betamod <- betareg(Proportion_NDF ~ fucose + xylose , data = af, link = "loglog")




