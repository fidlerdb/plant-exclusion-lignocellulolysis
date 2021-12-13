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
cf <- fread("../dbCAN_Results/Full_Assembly/Mapping Output/Blackout_ReadCounts_CAZyAnnotations.csv")
head(cf)

cf[n1 != 0,]
cf <- cf[n1 == 0,] # remove contigs mapped to in the negative

unique(cf$CAZyme) # Same as before.

# Keep only CBM, AA, and GH CAZymes as these are the ones involved in
# processing of lignocellulose
cf <- cf[grep("GH|AA|CBM", cf$CAZyme),]

# Get all unique CAZymes. Want to create a presence-absence matrix for all families.
AllCazymes <- unique(cf$CAZyme)

# Create a contig-wise CAZyme matrix
cm <- by(cf, cf$Chr, FUN = function(x){data.table(CAZymes = AllCazymes %in% x$CAZyme)})

# Format it nicely
cmm <- do.call("rbind", cm)
cmm <- data.table(cmm)
cmm <- cbind(names(cm), cmm)
names(cmm) <- c("Contig", AllCazymes)

# Converting to numerical values
(to.replace <- names(which(sapply(cmm, is.logical))))
for (var in to.replace) cmm[, (var):= as.numeric(get(var))]

head(cmm)

# largely clear the workspace
rm(list=setdiff(ls(), c("cmm", "df")))


df <- data.table(df)
head(df)

#### Merge the two datasets -- can take a while ~ 3 Gb of data###

# Perform the merge according to contig.

tcf <- merge(df, cmm, by = "Contig")
#tcf <- tcf[!is.na(GH116),]

#### Control for library size ####

dr <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZyFamilyReadsPerKbp.csv")

tcf$b1 <- tcf$b1 / (dr$Reads[dr$Sample == 'b1']/1e9)
tcf$b2 <- tcf$b2 / (dr$Reads[dr$Sample == 'b2']/1e9)
tcf$b3 <- tcf$b3 / (dr$Reads[dr$Sample == 'b3']/1e9)
tcf$b4 <- tcf$b4 / (dr$Reads[dr$Sample == 'b4']/1e9)
tcf$b5 <- tcf$b5 / (dr$Reads[dr$Sample == 'b5']/1e9)
tcf$b6 <- tcf$b6 / (dr$Reads[dr$Sample == 'b6']/1e9)
tcf$b7 <- tcf$b7 / (dr$Reads[dr$Sample == 'b7']/1e9)

tcf$c1 <- tcf$c1 / (dr$Reads[dr$Sample == 'c1']/1e9)
tcf$c2 <- tcf$c2 / (dr$Reads[dr$Sample == 'c2']/1e9)
tcf$c3 <- tcf$c3 / (dr$Reads[dr$Sample == 'c3']/1e9)
tcf$c4 <- tcf$c4 / (dr$Reads[dr$Sample == 'c4']/1e9)
tcf$c5 <- tcf$c5 / (dr$Reads[dr$Sample == 'c5']/1e9)
tcf$c6 <- tcf$c6 / (dr$Reads[dr$Sample == 'c6']/1e9)
tcf$c7 <- tcf$c7 / (dr$Reads[dr$Sample == 'c7']/1e9)

# Write this as a file -- all taxonomic data with CAZyme data appended.
#fwrite(tcf, file = "CAZy_2/outputData/Taxonomy_Data_CAZymes_Full.csv")

# Write this as a file. All taxonomic data with associated CAZyme data.
fwrite(tcf, file = "CAZy_2/outputData/Taxonomy_Data_CAZymes.csv")





