
rm(list = ls())

ff <- fread("Fibre_Analysis/Output_Data/Fibre_Analysis_Data_Processed_2021-04-29.csv")

head(ff)

ff[,Cellulose := Cellulose / (Sample_Mass - Ash)]
ff[,Hemicellulose := Hemicellulose / (Sample_Mass - Ash)]
ff[,Lignin := Lignin / (Sample_Mass - Ash)]

ff <- ff[,c("Sample", "Hemicellulose",  "Cellulose", "Lignin", "Treatment")  ]
#ff <- ff[,c("Sample", "NDF",  "ADF", "ADL", "Treatment")  ]



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
cor_glucose <- cor.test(x = af$Cellulose, y = af$glucose)
plot(glucose ~ Cellulose, data = af, pch = 16, col = af$Treatment
     , main = 'glucose') # Not for this one--mostly expected
legend("topright", legend = levels(af$Treatment), pch = 16, col = 1:4)
abline(lm(glucose ~ Cellulose, data = af))
text(x = 0.22, y = 14000, labels = paste0("r = ", round(cor_glucose[[4]],2)
                                          , ", t = ", round(cor_glucose[[1]],2)
                                          , "\ndf = ", cor_glucose[[2]]
                                          , ", p = ", round(cor_glucose[[3]],2)
))

## Hemicellulose -- need to change predictor variable
cor_xylose <- cor.test(x = af$Hemicellulose, y = af$xylose)
plot(xylose ~ Hemicellulose, data = af, pch = 16, col = af$Treatment
     , main = "xylose")
legend("topright", legend = levels(af$Treatment), pch = 16, col = 1:4)
abline(lm(xylose ~ Hemicellulose, data = af))
text(x = 0.2, y = 20000, labels = paste0("r = ", round(cor_xylose[[4]],2)
                                        , ", t = ", round(cor_xylose[[1]],2)
                                        , "\ndf = ", cor_xylose[[2]]
                                        , ", p = ", round(cor_xylose[[3]],2)
))

cor_fucose <- cor.test(x = af$Hemicellulose, y = af$fucose)
plot(fucose ~ Hemicellulose, data = af, pch = 16, col = af$Treatment
     , main = "fucose")
legend("topright", legend = levels(af$Treatment), pch = 16, col = 1:4)
abline(lm(fucose ~ Hemicellulose, data = af))
text(x = 0.2, y = 10000, labels = paste0("r = ", round(cor_fucose[[4]],2)
                                        , ", t = ", round(cor_fucose[[1]],2)
                                        , "\ndf = ", cor_fucose[[2]]
                                        , ", p = ", round(cor_fucose[[3]],2)
))

cor_galactose <- cor.test(x = af$Hemicellulose, y = af$`3,6-anhydro-D-galactose`)
plot(`3,6-anhydro-D-galactose` ~ Hemicellulose, data = af, pch = 16, col = af$Treatment
     , main = "3,6-anhydro-D-galactose") # Not really for this one
legend("topright", legend = levels(af$Treatment), pch = 16, col = 1:4)
abline(lm(`3,6-anhydro-D-galactose` ~ Hemicellulose, data = af))
text(x = 0.2, y = 850, labels = paste0("r = ", round(cor_galactose[[4]],2)
                                      , ", t = ", round(cor_galactose[[1]],2)
                                      , "\ndf = ", cor_galactose[[2]]
                                      , ", p = ", round(cor_galactose[[3]],2)
))



## lignin -- Results will change when ash is removed
cor_va <- cor.test(x = af$Lignin, y = af$`vanillic acid`)
plot(`vanillic acid` ~ Lignin, data = af, pch = 16, col = af$Treatment
     , main = "vanillic acid")
abline(lm(`vanillic acid` ~ Lignin, data = af))
text(x = 0.042, y = 2250, labels = paste0("r = ", round(cor_va[[4]],2)
                                       , ", t = ", round(cor_va[[1]],2)
                                       , "\ndf = ", cor_va[[2]]
                                       , ", p = ", round(cor_va[[3]],2)
))

cor_ba <- cor.test(x = af$Lignin, y = af$`benzoic acid`)
plot(`benzoic acid` ~ Lignin, data = af, pch = 16, col = af$Treatment
     , main = "benzoic acid")
abline(lm(`benzoic acid` ~ Lignin, data = af))
text(x = 0.042, y = 50000, labels = paste0("r = ", round(cor_ba[[4]],2)
                                         , ", t = ", round(cor_ba[[1]],2)
                                         , "\ndf = ", cor_ba[[2]]
                                         , ", p = ", round(cor_ba[[3]],2)
))

cor_hba <- cor.test(x = af$Lignin, y = af$`4-hydroxybenzoic acid`)
plot(`4-hydroxybenzoic acid` ~ Lignin, data = af, pch = 16, col = af$Treatment
     , main = "4-hydroxybenzoic acid")
abline(lm(`4-hydroxybenzoic acid` ~ Lignin, data = af))
text(x = 0.042, y = 25000, labels = paste0("r = ", round(cor_hba[[4]],2)
                                          , ", t = ", round(cor_hba[[1]],2)
                                          , "\ndf = ", cor_hba[[2]]
                                          , ", p = ", round(cor_hba[[3]],2)))

# Zero correlations of interest