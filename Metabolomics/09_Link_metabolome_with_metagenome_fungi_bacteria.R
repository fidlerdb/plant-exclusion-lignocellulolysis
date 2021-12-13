
# Clear the workspace
rm(list = ls())

# Data cleaning--produce a long and wide file of read abundance for each species

#### Load the data ####
# Load the taxonomy data
library(data.table)
library(magrittr)
tf <- data.frame(fread("Taxonomy/Cleaned_Data/Henfaes_Taxonomy.csv"))

head(tf)
tf$Species[!is.na(tf$Species)]

#### Scale the data by contig size and read length ####
# Scale reads by contig size
tf[, 5:length(tf)] <- tf[, 5:length(tf)]/(tf$ContigLength/1000) 

# Now scale by library size
dr <- read.csv("CAZy/Cleaned_Data/Henfaes_CAZyFamilyReadsPerKbp.csv")
tf$b1 <- tf$b1 / (dr$Reads[dr$Sample == 'b1']/1e9)
tf$b2 <- tf$b2 / (dr$Reads[dr$Sample == 'b2']/1e9)
tf$b3 <- tf$b3 / (dr$Reads[dr$Sample == 'b3']/1e9)
tf$b4 <- tf$b4 / (dr$Reads[dr$Sample == 'b4']/1e9)
tf$b5 <- tf$b5 / (dr$Reads[dr$Sample == 'b5']/1e9)
tf$b6 <- tf$b6 / (dr$Reads[dr$Sample == 'b6']/1e9)
tf$b7 <- tf$b7 / (dr$Reads[dr$Sample == 'b7']/1e9)

tf$c1 <- tf$c1 / (dr$Reads[dr$Sample == 'c1']/1e9)
tf$c2 <- tf$c2 / (dr$Reads[dr$Sample == 'c2']/1e9)
tf$c3 <- tf$c3 / (dr$Reads[dr$Sample == 'c3']/1e9)
tf$c4 <- tf$c4 / (dr$Reads[dr$Sample == 'c4']/1e9)
tf$c5 <- tf$c5 / (dr$Reads[dr$Sample == 'c5']/1e9)
tf$c6 <- tf$c6 / (dr$Reads[dr$Sample == 'c6']/1e9)
tf$c7 <- tf$c7 / (dr$Reads[dr$Sample == 'c7']/1e9)

# See how this looks
head(tf)

# Add other taxonomic information
class(tf)
tf$Domain <- sub(";.*", "",tf$Taxonomy)
tf$Phylum <- sub("^[A-Za-z]+;", "",tf$Taxonomy, perl = TRUE) %>% sub(";.*", "", .)

head(tf)
tf$Phylum

sort(unique(tf$Phylum))
#fungi <- tf[tf$Phylum %in% c("Basidiomycota", "Ascomycota"),]
#bacteria <- tf[tf$Domain == "Bacteria",]

# Create a new data frame for further analysis
df <- tf
#### Get rid of unassigned contigs, and remove species  ####
####     which were found in the negative ##################

# How many responsive CAZymes 
# are there, and how many of these match to species?
nrow(df)                      # 6.7 million contigs
nrow(df[complete.cases(df),]) # 1.2 million contigs  
#     which could be matched by descriminative k-mers at  
#     the species level to genomes on RefSeq by CLARK

# 17% of all contigs

# Get rid of unmatched species
#df <- na.omit(df)
head(df)

# Find out what species were mapped to in the negative
sum(df$n1)
NegSpecies <- unique(df$Species[df$n1 != 0])

# Remove these species from the analysis
df <- df[!df$Species %in% NegSpecies,]

# How many species are left
length(unique(df$Species))

# Remove the negative
df <- df[!names(df) %in% 'n1'] 

head(df)

setDT(df)

# Keep only bacterial or fungal contigs

fungi <- df[Phylum %in% c("Basidiomycota", "Ascomycota"),]
bacteria <- df[df$Domain == "Bacteria",]

df <- rbind(fungi, bacteria)

Samples <- c("b1", "b2", "b3", "b4", "b5", "b6", "b7"
             , "c1", "c2", "c3", "c4", "c5", "c6", "c7")
df <- df[, lapply(.SD, sum), by = Domain, .SDcols = Samples]

df2 <- df[1, ..Samples]/df[2, ..Samples]

df3 <- melt(df2, variable.name = "Sample", value.name = "FB")

df3

df <- df3

rm(list = ls()[!ls() %in% "df"])

#### Work on the metabolites ####
mf <- data.frame(fread("Metabolomics/Data/Metabolome_Blackout_Raw.csv"))

mf <- mf[, c(1,7:20)] # Keep only the most useful columns

# Lignocellulose breakdown products
BreakdownProducts <- c('glucose'
                       , 'xylose', 'fucose', '3,6-anhydro-D-galactose'
                       , 'vanillic acid', '4-hydroxybenzoic acid', 'benzoic acid')

# Keep only those products
mf <- mf[mf$BinBase.name %in% BreakdownProducts,]

head(mf)
AbundData <- data.frame(t(mf[2:length(mf)]))
names(AbundData) <- mf[,1]
AbundData$Sample <- row.names(AbundData)

setDT(mf)
dt <- merge(df, AbundData, by = "Sample")
dt <- melt(dt, id = c("Sample", "FB"), variable.name = "Chemical"
     , value.name = "Abundance")
dt

ggplot(dt, aes(x = FB, y = Abundance)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~Chemical)

lapply(unique(dt$Chemical), FUN = function(x){
  mod <- lm(Abundance ~ FB
            , data = dt[Chemical == x,])
  Test <- drop1(mod, test = "F")
  P <- Test$`Pr(>F)`[[2]]
  return(data.table(Chemical = x, p = P))
}) %>% rbindlist(.)





