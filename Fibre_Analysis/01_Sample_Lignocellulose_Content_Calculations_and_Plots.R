# Clean workspace 
rm( list = ls())
gc(full = TRUE)

library(data.table)
library(ggplot2)

#### Import data ####
df <- fread("../Data/Fibre Analysis/Ankom_Raw_Data.csv")

df

#### Calculate probable actual mass of c7 sample ####
# added into bag, based on difference of others in the control group

# View the distribution of differences for all samples and the control group
par(mfrow = c(2,1))
with(df, hist(na.omit(Mass_added_minus_actual)))
with(df, abline(v = mean(na.omit(Mass_added_minus_actual)), col = "red"))
with(df, abline(v = median(na.omit(Mass_added_minus_actual)), col = "blue"))

with(df[grep("c", Sample)], hist(na.omit(Mass_added_minus_actual)))
with(df[grep("c", Sample)], abline(v = mean(na.omit(Mass_added_minus_actual)), col = "red"))
with(df[grep("c", Sample)], abline(v = median(na.omit(Mass_added_minus_actual)), col = "blue"))
par(mfrow = c(1,1))

# Ratio of sample mass to sample mass added
hist(df[grep("c", Sample), na.omit(Sample_Mass / Sample_Mass_Added)])
with(df[grep("c", Sample)], abline(v = mean(na.omit(Sample_Mass / Sample_Mass_Added)), col = "red"))
with(df[grep("c", Sample)], abline(v = median(na.omit(Sample_Mass / Sample_Mass_Added)), col = "blue"))

# Calculate the mean proportion of sample mass lost between sample weighing and sample+bag weighing
df[grep("c", Sample), mean(na.omit(Sample_Mass / Sample_Mass_Added))]

Mass_added_correction <- df[grep("c", Sample), median(na.omit(Sample_Mass / Sample_Mass_Added))]

df[Sample == "c7", Sample_Mass := Sample_Mass_Added * Mass_added_correction]

#### Calculate blank bag correction factor and ash content ####

bf <- df[grep("blank", Sample),]

bf

# What happened to bag mass between procedures?
plot(value ~ variable, data = melt(bf[,c("Sample", "Bag_Mass", "NDF", "ADF", "ADL")], id = "Sample"))

bf$NDF / bf$Bag_Mass

NDF_blank_Correction <- mean(bf$NDF/ bf$Bag_Mass)
NDF_blank_Correction_min <- min(bf$NDF / bf$Bag_Mass)
ADF_blank_Correction <- mean(bf$ADF / bf$Bag_Mass) # Use NDF as t0 as the analysis is sequential
ADF_blank_Correction_min <- min(bf$ADF / bf$Bag_Mass) # Use NDF as t0 as the analysis is sequential
ADL_blank_Correction <- mean(bf$ADL / bf$Bag_Mass) # Use ADF as t0 as the analysis is sequential

bf[, Ash := Crucible_and_Ash - Crucible]

hist(bf$Ash)

Ash_Correction <- mean(bf$Ash)

#bf$Bag_Mass / (bf$Crucible_and_Ash - bf$Crucible)

# Keep the samples only
dfs <- df[grep("^[a-z]{1}[0-9]", Sample),]

#### c1 is a bit funny ####
dfs[Sample == "c1",]
# NDF                      bag mass
dfs[Sample == "c1",]$NDF - (dfs[Sample == "c1",]$Bag_Mass * NDF_blank_Correction) 
dfs[Sample == "c1",]$NDF - (dfs[Sample == "c1",]$Bag_Mass * NDF_blank_Correction_min) 
dfs[Sample == "c1",]$NDF - dfs[Sample == "c1",]$Bag_Mass

bf$ADF / bf$NDF
dfs[Sample == "c1",]$ADF - (dfs[Sample == "c1",]$Bag_Mass * ADF_blank_Correction) 


#### Calculate Ash mass ####

dfs[, Ash := Crucible_and_Ash - Crucible - Ash_Correction]

dfs$Ash

#### Calculation of polymer abundance ####

# Create a new data table which keeps only the samples where the useful values can be calculated
df2 <- dfs[, c("Sample", "Sample_Mass", "Bag_Mass")]

# Bagless sample masses
df2[, Ash := dfs$Ash]
df2[, ADL := dfs$ADL - (ADL_blank_Correction * dfs$Bag_Mass)] # Sample mass after ADL (bag mass removed)
df2[, ADF := dfs$ADF - (ADF_blank_Correction_min * dfs$Bag_Mass)] # Sample mass after ADF (bag mass removed)
df2[, NDF := dfs$NDF - (NDF_blank_Correction_min * dfs$Bag_Mass)] # Sample mass after NDF (bag mass removed)

df2
with(df2, ADL < ADF)
with(df2, ADF < NDF)

# Distribution of ash-containing digestion masses
hist(df2$NDF)
hist(df2$ADF)
hist(df2$ADL)

# Mass of lignocellulosic polymers
df2[, Lignin := ADL - Ash]
# df2[, Lignin := ADL] # Testing
df2[, Cellulose := ADF - ADL]
df2[, Hemicellulose := NDF - ADF]
df2[, Labile := Sample_Mass - NDF]

# Are the values possible -- twiddle with correction favtors if not

hist(df2$Lignin) # Possible
hist(df2$Cellulose) # Two impossible values using mean -- c1 and c4 - possible using ADF_blank_Correction_min
hist(df2$Hemicellulose) # Two imossible values using mean -- possible using NDF_blank_Correction_min
hist(df2$Labile) 

df2

# Proportion mass of lignocellulosic polymers
# df2[, Ash_prop := Ash / Sample_Mass]
# df2[, Lignin_prop := Lignin / Sample_Mass]
# df2[, Cellulose_prop := Cellulose / Sample_Mass]
# df2[, Hemicellulose_prop := Hemicellulose / Sample_Mass]


#### Plotting what we have ####

colsToKeep <- c("Sample", "Sample_Mass", "Hemicellulose", "Cellulose", "Lignin", "Ash", "Labile")
dfl <- melt(df2[,..colsToKeep], id = c("Sample", "Sample_Mass"))

#dfl$variable <- factor(dfl$variable, levels = c("Ash", "ADL", "ADF", "NDF"), ordered = TRUE)

ggplot(dfl
       , aes(x = Sample, y = value, fill = variable)) + geom_bar(stat = "identity")

ggplot(dfl
       , aes(x = Sample, y = value/Sample_Mass, fill = variable)) + geom_bar(stat = "identity")

# All adds up to 1. Nice.

ggplot(dfl
       , aes(x = Sample, y = value/Sample_Mass, fill = variable)) + geom_bar(stat = "identity", position = "dodge")

#### Save file to .csv ####
df2$Treatment <- as.factor(paste(sub("[0-9]", "", df2$Sample)
                                , ifelse(sub("[a-z]", "", df2$Sample) > 3
                                         , yes = "New", no = "Old"
                                )
                                , sep= "-")
)
df2$Treatment <- relevel(df2$Treatment, ref = "b-Old")

df2

# Final check that values are correct
data.frame(Total = df2$Labile + df2$Hemicellulose + df2$Cellulose + df2$Lignin + df2$Ash
           , Sample_Mass = df2$Sample_Mass)

df2$Labile + df2$Hemicellulose + df2$Cellulose + df2$Lignin + df2$Ash == df2$Sample_Mass


# Write this to an output file
fwrite(df2, "Fibre_Analysis/Output_Data/Fibre_Analysis_Data_Processed_2021-04-29.csv")

