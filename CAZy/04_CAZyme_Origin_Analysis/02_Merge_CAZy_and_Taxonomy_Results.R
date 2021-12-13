rm(list = ls())

library(data.table)
library(future.apply)
library(roperators)

#### Import, check, and set the data up -- Scale data ####
####     by contig length ####

# Read in the taxonomy data
df <- data.frame(fread("Taxonomy/Multi-Method/Cleaned_Data/Henfaes_Taxonomy_Multi-Method_ReadCounts.csv"))

# See how this looks
head(df)

# Divide by Kbp contig length, remove hits in the negative. 
# Control for library size (billion reads) later

df[, 14:length(df)] <- df[, 14:length(df)]/(df$ContigLength/1000)

# Work with the species matrix - call this sf
sf <- (df[, 14:length(df)])
head(sf[,1:5]) # Check the data
sf
#### Change NA species to UNASSIGNED ####

class(df$Species) # Should be character

# Reassign the unidentified contigs so they can be analysed as a unit
df$Species[is.na(df$Species)] <- "UNASSIGNED"

#### Remove unassigned contigs mapping to the negative ####

# Select rows where the contig could not be assigned taxonomy, and where 
# there are reads in the negative
NegFound <- sf[df$Multi.Assignment_TaxID == 0 &  sf$n1 > 0,] # THE ZERO HERE WANTS TO BE CHANGED TO account for read hopping?

RowsToRemove <- as.numeric(
  row.names(NegFound)
) # Create a vector of unassigned rows found in the negative

# Create a new full data frame without unassigned species found in the negative

if(length(RowsToRemove) == 0){sn <- df} else {sn <- df[-RowsToRemove,]}
head(sn)

# Now make a set of taxonomies
Lineage <- paste(sn$Species, sn$Genus, sn$Family
                 , sn$Order, sn$Class, sn$Phylum
                 , sn$SuperKingdom, sep = ';')
head(Lineage)

# Take the species matrix of sn 
# (full data minus unassigned contigs which occur in the negative)
sn_s <- sn[, 14:length(sn)]

#### Group contigs by sequence (create grouping factor)              ####
####     of reads mapping to each species--allows analysis of number ####

# sn[7:13] is the species and taxonomy. Create a vector of useable names
# e.g. Unknown Proteobacteria
plan(multiprocess)
Labels <- future_apply(sn[7:13], 1, FUN = function(x){
  
  if(x[1] # Species 
     == 'UNASSIGNED' & !is.na(x[2] # SuperKingdom
     )){
    
    First_Non_NA <- which(!is.na(rev(x)))[1] # Take fisrt non NA classification
    
    paste("Unknown", rev(x)[First_Non_NA], sep = " ")
    
  } else {x[1] # Species
  }
}) 
plan(sequential) # Go back to sequential task processing

head(Labels)
head(sn)
sn_tree <- sn
sn_tree$Species <- Labels

df <- sn_tree
head(df)

# Subset to the columns we want (contig name, taxonomy, and read abundance)
df <- df[, c(1,7:length(df))]
head(df)

df <- df[df$n1 == 0,] # Remove all contigs found in the negative
head(df)

max(df$n1)

#### Now add CAZy data ####

# Import thre CAZyme data
cf <- data.frame(fread("../dbCAN_Results/Full_Assembly/Mapping Output/Blackout_ReadCounts_CAZyAnnotations.csv"))
head(cf)

# Keep only the columns of CAZyme and contig
cf <- cf[, names(cf) %in% c("CAZyme", "Chr")]

# Keep only CBM, AA, and GH CAZymes as these are the ones involved in
# processing of lignocellulose
cf <- cf[grep("GH|AA|CBM", cf$CAZyme),]

# Put all CAZymes which occur on the same contig together
CAZyList <- tapply(cf$CAZyme, INDEX = cf$Chr, FUN = function(x)paste(x, collapse = "."))
CAZyList <- data.frame(Contig = names(CAZyList), CAZyme = CAZyList)

# Pufec
head(CAZyList)
head(df)

#### Bind together the CAZyme, taxonomy, and read abundance data
df2 <- merge(df, CAZyList, all.x = TRUE, by = "Contig")

head(df2)
max(df$n1)

# Control for library size
dr <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZyFamilyReadsPerKbp.csv")

df2$b1 <- df2$b1 / (dr$Reads[dr$Sample == 'b1']/1e9)
df2$b2 <- df2$b2 / (dr$Reads[dr$Sample == 'b2']/1e9)
df2$b3 <- df2$b3 / (dr$Reads[dr$Sample == 'b3']/1e9)
df2$b4 <- df2$b4 / (dr$Reads[dr$Sample == 'b4']/1e9)
df2$b5 <- df2$b5 / (dr$Reads[dr$Sample == 'b5']/1e9)
df2$b6 <- df2$b6 / (dr$Reads[dr$Sample == 'b6']/1e9)
df2$b7 <- df2$b7 / (dr$Reads[dr$Sample == 'b7']/1e9)

df2$c1 <- df2$c1 / (dr$Reads[dr$Sample == 'c1']/1e9)
df2$c2 <- df2$c2 / (dr$Reads[dr$Sample == 'c2']/1e9)
df2$c3 <- df2$c3 / (dr$Reads[dr$Sample == 'c3']/1e9)
df2$c4 <- df2$c4 / (dr$Reads[dr$Sample == 'c4']/1e9)
df2$c5 <- df2$c5 / (dr$Reads[dr$Sample == 'c5']/1e9)
df2$c6 <- df2$c6 / (dr$Reads[dr$Sample == 'c6']/1e9)
df2$c7 <- df2$c7 / (dr$Reads[dr$Sample == 'c7']/1e9)

head(df2) # Check the format
max(df2$n1) # Check we removed the contaminants

write.csv(df2
          , file = "CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_dbCAN_Taxonomy_ReadAbundance.csv"
          , row.names = FALSE)





