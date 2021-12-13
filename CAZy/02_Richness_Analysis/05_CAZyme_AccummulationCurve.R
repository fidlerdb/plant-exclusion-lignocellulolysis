
rm(list = ls())

#### Data import  and checking ####
# df <- read.csv("../Blackout_R/CAZy/Cleaned_Data/Henfaes_CAZyFamilyReads_KbpGene_NoNegative.csv")
df <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZymes_NotGrouped_Boolean.csv")

head(df[1:10])

ncol(df) -3
names(df)[3:length(df)]

#############

df[, 4:length(df)] <- df[, 4:length(df)] / (df$Reads/1e9) # Normalise per Gbp

df[,1:10]



hist(unlist(df[, 4:length(df)]), breaks = 100) # Range of read counts

#head(rft)

library(vegan)

# Species accummula
df2 <- df[4:length(df)]# tion curve
sp1 <- specaccum(df2)
sp2 <- specaccum(df2, "random")
sp2
head(df2)[1:10]
summary(sp2)
plot(sp1
     , ci.type="poly"
     , col="blue"
     , lwd=2
     , ci.lty=0
     , ci.col="lightblue"
     , ylab = "Cumulative number of CAZymes"
     , xlab = "Number of samples"
     #, xaxs = "i"
     , yaxs = "i"
     , bty = 'l')
boxplot(sp2#, col="yellow"
        , add=TRUE
        , pch="+")


#### How many CAZy families were there really? ####

head(df)
unique(names(df))

AllCAZymes <- unique(gsub("\\.([0-9])*", "", names(df)[4:length(df)]))
length(AllCAZymes)
length(grep("GH", AllCAZymes))
length(grep("CBM", AllCAZymes))
length(grep("AA", AllCAZymes))

# -- Reported numbers seem correct according to this methodology
#       assume this is to do with removal of whole families found 
#       in the negative rather than just removing the contigs.

# ############ Rarefaction curve
# 
# df3 <- read.csv("../Blackout_R/CAZy/Cleaned_Data/Henfaes_CAZymes_NotGrouped_Boolean.csv")
# 
# df3[,1:10]
# 
# # Identify fully 0 columns
# ZeroCols <- names(which(colSums(df3[4:length(df3)]) == 0))
# 
# # identify things present in negative
# NegCols <- names(df3[, which(df3[4:length(df3), 15] != 0)+3])
# 
# # Remove genes which are found in the negative
# df3 <- df3[names(df3) %ni% c(ZeroCols, NegCols)]
# 
# # Remove the negative
# df3 <- df3[-15,]
# # df3 <- apply(df3,2, as.numeric)
# 
# # Take the matrix
# df4 <- df3[4:length(df3)]
# df4[1:10]
# 
# # Do the rarefaction
# S <- specnumber(df4) # observed number of species
# (raremax <- min(rowSums(df4)))
# Srare <- rarefy(df4, raremax)
# plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# abline(0, 1)
# rarecurve(df4, step = 20, sample = raremax, col = "blue", cex = 0.6)
# 
# 
