#### A more robust evenness analysis

rm(list = ls())

#library(doBy)
library(data.table)
library(magrittr)

#### Read in and clean the taxonomy+CAZy data--normalized ####

# Read in data
df <- fread("CAZy_2/outputData/Taxonomy_Data_CAZymes.csv")
df <- df[Phylum != "",] # Keep only data with phylum-level information

df

ccf <- fread("CAZy_2/outputData/All_Correlated_CAZymes-Metabolites.csv")
# Make sure that these will match with the CAZymes in the dataset
ccf$CAZymes <- gsub("\\.([A-Z])", "\\|\\1", ccf$CAZymes)
# Remove any remaining "." characters between CAZymes
ccf$CAZymes <- gsub("([A-Z].+?)\\.", "\\1\\|", ccf$CAZymes)
unique(ccf$CAZymes)
AllCorrelatedCAZymes <- unique(ccf[p < 0.05,]$CAZymes)


# Keep only correlated CAZy families
ColumnsToKeep <- c(names(df)[1:23], names(df)[names(df) %in% unique(ccf$CAZymes)])
df <- df[,..ColumnsToKeep]

# Delete fully zero columns
colsToDelete <- names(which(colSums(df[,24:length(df)]) == 0))
data.table::set(df, , colsToDelete, NULL)

# Get a list of the CAZymes
CAZymes <- names(df)[24:length(df)]

# Sum abundance from each smple for contigs where the CAZyme was present
FindOneCAZymeAbundance <- function(i){
  summaryCols <- c(c("Phylum", "b1", "b2", "b3", "b4"
                     , "b5" ,"b6" ,"b7" ,"c1" ,"c2" ,"c3" ,"c4" ,"c5", "c6", "c7"), CAZymes[[i]])

  # Individual CAZyme abundance for each phylum
  CAZy_col <- CAZymes[[i]] # Keep only contigs where that CAZyme was present
  rowsToKeep <- which(df[, ..CAZy_col] == 1)
  df_s <- df[rowsToKeep,]
  
  OneCAZymeAbundance <- by(df_s[,..summaryCols], df_s$Phylum, FUN = function(x){
    Total_abundance <- colSums(x[,2:(ncol(x)-1)])

    data.table(Abundance = Total_abundance
               , Sample = names(Total_abundance)
               , CAZyme = CAZy_col)
  })

  # All data for one CAZyme
  return(OneCAZymeAbundance)
}

rf <- lapply(seq_along(CAZymes), FindOneCAZymeAbundance)

# Stick them all together
rf <- lapply(rf, FUN = function(x){return(rbindlist(x, idcol = "Phylum"))})
rf <- rbindlist(rf)

unique(rf$CAZyme)
# Add a treatment column
rf$Treatment <- "a"
rf[Sample %in% c("b1","b2","b3"),]$Treatment <- "10-Year Fallow"
rf[Sample %in% c("b4","b5","b6", "b7"),]$Treatment <- "1-Year Fallow"
rf[Sample %in% c("c1","c2","c3"),]$Treatment <- "10-Year Grassland"
rf[Sample %in% c("c4","c5","c6", "c7"),]$Treatment <- "1-Year Grassland"

head(rf)

fwrite(rf, file = "CAZy_2/outputData/Phylum_Level_CAZyme_NormalizedReadAbundance.csv")

#### Read in the CAZyme cluster data ####

# Read in the CAZymes in the clusters
czf <- fread("CAZy_2/05-10_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")
czf

AllClusters <- unique(czf$ClusterName)

#### For each Phylum, create a matrix of read abundances across all CAZymes (columns) for each sample (rows) ####

# Create matrices for all phyla
pf <- by(rf, rf$Phylum, FUN = function(x){
  dcast(x, Sample  ~ CAZyme, value.var = "Abundance")
  })

# The ones we're really interested in
TaxaOfInterest <- c("Proteobacteria"
                    ,"Actinobacteria"
                    ,"Acidobacteria"
                    ,"Bacteroidetes"
                    ,"Firmicutes"
                    ,"Planctomycetes"
                    ,"Cyanobacteria"
                    ,"Verrucomicrobia"
                    ,"Euryarchaeota"
                    ,"Gemmatimonadetes"
                    ,"Ascomycota")
pf[TaxaOfInterest]

#### For each phylum, find the proportion of reads which belongs to each CAZyme cluster ####

# Create a function to sum reads from the same CAZyme cluster
get_cluster_Abundance <- function(x, cluster){
  # Get a list of CAZymes in the clusters
  CC_families <- czf[ClusterName == cluster,]$CAZyme
  
  # Make sure data.table works
  CC_families <- names(x) %in% CC_families
  
  # Find columns of CAZy families in the cluster for this phylum, sum the abundances
  rowSums(x[,..CC_families])
  
  Result <- data.table(cluster = rowSums(x[,..CC_families]))
  if(length(Result$cluster) != 14){
    Result <- data.table(cluster = rep(0,14))
  }
  colnames(Result) <- cluster
  return(Result)
}

# Prove it works
get_cluster_Abundance(pf[[1]], "glucose-CazC1")
get_cluster_Abundance(pf[[2]], "glucose-CazC1")
get_cluster_Abundance(pf[[3]], "4-hydroxybenzoic acid-CazC1")

# Create a function to replace NA or NaN with 0. Many 0/0 in this dataset 
DT_NA_0 = function(DT) {
  # either of the following for loops
  
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    data.table::set(DT,which(is.na(DT[[j]])),j,0)
}

# For every phylum, sum abundance of reads from each sample
#  loop along every cluster
pf_ca <- lapply(pf, FUN = function(x){
  Cluster_Abundances <- lapply(seq_along(AllClusters), FUN = function(i){get_cluster_Abundance(pf[[1]], AllClusters[i])})
  Cluster_Abundances <- do.call(cbind, Cluster_Abundances)
  return(Cluster_Abundances)
})

pf_ca
heatmap(as.matrix(pf[[1]][,-1]))
heatmap(as.matrix(pf[[2]][,-1]))
heatmap(as.matrix(pf[[3]][,-1]))
heatmap(as.matrix(pf[[4]][,-1]))

#####

# create function to find proportional abundance of reads mapping to genes which 
#   correlated with an individual metabolite, for all metabolites

get_phylum_cluster_proportions <- function(data){
  # Loop over all metabolites and calculate the 
  #   total number of reads mapping to contigs containing 
  #   these CAZymes
  Phylum_Cluster_Abundance <- lapply(AllClusters, FUN = function(cluster){
    get_cluster_Abundance(data, cluster)
  })
  
  Phylum_Cluster_Abundance <- do.call(cbind, Phylum_Cluster_Abundance)
  
  # Turn this into % of reads from CAZy genes which were correlated with a particular metabolite
  Phylum_Cluster_Abundance <- lapply(unique(sub("-.*", "", names(Phylum_Cluster_Abundance))),
                                     FUN = function(metabolite){
                                       
                                       # Keep daa.table happy
                                       NamesToKeep <- names(Phylum_Cluster_Abundance)[sub("-.*", "", names(Phylum_Cluster_Abundance))
                                                                                      == metabolite]
                                       # Select all columns which deal with the same metabolite
                                       mf <- Phylum_Cluster_Abundance[,..NamesToKeep]
                                       return(mf/rowSums(mf))
                                     }
  )
  Phylum_Cluster_Abundance <- do.call(cbind, Phylum_Cluster_Abundance)
  DT_NA_0(Phylum_Cluster_Abundance) #Remove NaNs
  return(Phylum_Cluster_Abundance)
}

# Find cluster proportions for all phyla
Phylum_Cluster_Proportions <- lapply(pf, get_phylum_cluster_proportions)
Phylum_Cluster_Proportions[[2]]


#### Begin the evenness analysis ####

# See how genes are distributed in a few phyla
heatmap(as.matrix(Phylum_Cluster_Proportions[[1]]))
heatmap(as.matrix(Phylum_Cluster_Proportions[[2]]))
heatmap(as.matrix(Phylum_Cluster_Proportions[[3]]))
heatmap(as.matrix(Phylum_Cluster_Proportions[[4]]))

names(Phylum_Cluster_Proportions)

# Abundant degradative taxa
TaxaOfInterest <- c("Proteobacteria"
,"Actinobacteria"
,"Acidobacteria"
,"Bacteroidetes"
,"Firmicutes"
,"Planctomycetes"
,"Cyanobacteria"
,"Verrucomicrobia"
,"Euryarchaeota"
,"Gemmatimonadetes"
,"Ascomycota")

pcpd <- Phylum_Cluster_Proportions[names(Phylum_Cluster_Proportions) %in% TaxaOfInterest]

saveRDS(pcpd
        , "CAZy_2/outputData/Phylum_Level_CAZyme_Cluster_Percentage_Abundance_perMetabolite.rds")

heatmap(as.matrix(pcpd[[8]]))

library(vegan)
library(ggplot2)

# Find CAZy cluster evenness
n = names(pcpd)
pc_evenness <- mapply(function(list.elem, names) {diversity(list.elem)/log(specnumber(list.elem))}
       , list.elem = pcpd, names = n) %>%
  data.table()
pc_evenness$Sample <- pf[[1]]$Sample
pc_evenness.m <- melt(pc_evenness, id = "Sample")
sortedPhyla <- names(sort(tapply(pc_evenness.m$value
                                 , pc_evenness.m$variable, mean)
                          , decreasing = TRUE))
pc_evenness.m$variable <- factor(pc_evenness.m$variable
                                     , levels = sortedPhyla
)

# quickly check the species numbers
n = names(pcpd)
pc_spec_number <- mapply(function(list.elem, names) {specnumber(list.elem)}
                      , list.elem = pcpd, names = n) %>%
  data.table()
pc_spec_number # All the same :( evenness is the same as diversity

unique(pc_evenness.m$variable)

# Commit this to a file-- now use it in script 11 in this directory
fwrite(pc_evenness
       , "CAZy_2/outputData/Metabolite-Correlated-CAZyme_Evenness_Phylum.csv"
       )

n = names(Phylum_Cluster_Proportions)
pc_evenness_all <- mapply(function(list.elem, names) {diversity(list.elem)/log(specnumber(list.elem))}
                      , list.elem = Phylum_Cluster_Proportions, names = n) %>%
  data.table()
pc_evenness_all$Sample <- pf[[1]]$Sample
fwrite(pc_evenness_all
       , "CAZy_2/outputData/Metabolite-Correlated-CAZyme_Evenness_AllPhyla.csv"
       )
pc_evenness_all.m <- melt(pc_evenness_all, id = "Sample")
pc_evenness_all.m 
pc_evenness_all.m <- pc_evenness_all.m[complete.cases(pc_evenness_all.m),]


sortedPhyla <- names(sort(tapply(pc_evenness_all.m$value
                                 , pc_evenness_all.m$variable, mean)
                          , decreasing = TRUE))
pc_evenness_all.m$variable <- factor(pc_evenness_all.m$variable
                                     , levels = sortedPhyla
)

#### Plot the results ####

# Look at all of the evenness data
ggplot(pc_evenness_all.m , aes(x = variable, y = value)) + 
  geom_point(position = position_jitter(0.1)) + 
  stat_summary(fun.data = mean_sdl, colour = 'red', size = 1, fun.args = list(mult = 1)) +
  ylim(0,1) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Ignavibacteridae
# Spirochaetes
# Calditrichaeota
# Nitrospirae
# Thermotogae

# Look only at the phyla in the plot:
ggplot(pc_evenness.m, aes(x = variable, y = value)) + 
  stat_summary(aes(y = value)
               , fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_sdl
               , geom = "errorbar"
               , fun.args = list(mult = 1)
               , width = 0.3
               , size = 1.5) +
  geom_point(position = position_jitter(0.3), colour = alpha("black", 0.5)) + 
  ylim(0,1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        ) +
  geom_point(data = data.frame(variable = "Actinobacteria"
                               , value = 1)
                               , aes(x = variable, y = value)
             , colour = "white") +
  scale_y_continuous(expand = c(0, 0), name = "J'"
                               #, breaks = c(0,1)
  ) 

#head(pc_evenness.m)

# Get the values for writing in the results/discussion
pce_msd <- data.table(summaryBy(value ~ variable
                     , data = pc_evenness.m
                     , FUN = c(mean,sd)))

pce_msd$value.mean <- round(pce_msd$value.mean, 2)
pce_msd$value.sd <- round(pce_msd$value.sd, 2)

pce_msd[order(value.mean),]

# Gemmatimonadetes  0.76     0.02
# Actinobacteria    0.79     0.01
# Planctomycetes    0.79     0.01

# Firmicutes        0.80     0.03
# Proteobacteria    0.80     0.01
# Verrucomicrobia   0.81     0.03
# Cyanobacteria     0.83     0.02
# Acidobacteria     0.84     0.02
# Euryarchaeota     0.84     0.04
# Bacteroidetes     0.85     0.02

# Ascomycota        0.92     0.03

#### One final bit of exploration, mean number of reads mapping to each cluster  ####

# pc_evenness <- mapply(function(list.elem, names) {diversity(list.elem)/log(specnumber(list.elem))}
#                       , list.elem = pcpd, names = n) %>%
#   data.table()
# pc_evenness$Sample <- pf[[1]]$Sample
# pc_evenness.m <- melt(pc_evenness, id = "Sample")

acd <- rbindlist(pcpd, idcol = "Phylum")
acd$Sample <- rep(pf[[1]]$Sample, times = length(pcpd))
names(acd)

# `glucose-CazC1`                 + `glucose-CazC2`                + 
# `glucose-CazC3`                 + `glucose-CazC4`                + `glucose-CazC5`                
# `fucose-CazC1`                  + `fucose-CazC2`                 + `fucose-CazC3`                 
# `fucose-CazC4`                  + `xylose-CazC1`                 + `xylose-CazC2`                 
# `xylose-CazC3`                  + `xylose-CazC4`                 + `xylose-CazC5`                 
# `3,6-anhydro-D-galactose-CazC0` + `3,6-anhydro-D-galactose-CazC1`+ `vanillic acid-CazC1`          
# `vanillic acid-CazC2`           + `vanillic acid-CazC31`         + `vanillic acid-CazC4`          
# `4-hydroxybenzoic acid-CazC1`   + `4-hydroxybenzoic acid-CazC2`  + `4-hydroxybenzoic acid-CazC3`  
# `4-hydroxybenzoic acid-CazC4`   + `4-hydroxybenzoic acid-CazC5`  + `4-hydroxybenzoic acid-CazC6`  
# `4-hydroxybenzoic acid-CazC7`   + `benzoic acid-CazC1`           + `benzoic acid-CazC2` 
``
acd_m <- summaryBy(`glucose-CazC1` + `glucose-CazC2` 
                   + `glucose-CazC3`                 + `glucose-CazC4`                + `glucose-CazC5`                
                   + `fucose-CazC1`                  + `fucose-CazC2`                 + `fucose-CazC3`                 
                   + `fucose-CazC4`                  + `xylose-CazC1`                 + `xylose-CazC2`                 
                   + `xylose-CazC3`                  + `xylose-CazC4`                 + `xylose-CazC5`                 
                   + `3,6-anhydro-D-galactose-CazC0` + `3,6-anhydro-D-galactose-CazC1`+ `vanillic acid-CazC1`          
                   + `vanillic acid-CazC2`           + `vanillic acid-CazC3`         + `vanillic acid-CazC4`          
                   + `4-hydroxybenzoic acid-CazC1`   + `4-hydroxybenzoic acid-CazC2`  + `4-hydroxybenzoic acid-CazC3`  
                   + `4-hydroxybenzoic acid-CazC4`   + `4-hydroxybenzoic acid-CazC5`  + `4-hydroxybenzoic acid-CazC6`  
                   + `4-hydroxybenzoic acid-CazC7`   + `benzoic acid-CazC1`           + `benzoic acid-CazC2` 
          ~ Phylum
          , data = acd
          , FUN = mean
          , keep.names = TRUE)

# Basic read abundance data
heatmap(as.matrix(acd_m[,2:length(acd_m)]))

sub("-.*", "", names(acd_m))

rowSums(acd_m[, .SD, .SDcols = names(acd_m) %like% "fucose"])
acd_m[, .SD, .SDcols = names(acd_m) %like% "fucose"]
acd_m[, .SD, .SDcols = names(acd_m) %like% "^benzoic"]

# Manly distance matrix (for percentage abundances)
Manly <- function(x){dist.prop(data.frame(x), method = 1)}

source("Functions/ClusterHeatmap_fromMatrix.R")

if(!file.exists("CAZy/04_CAZyme_Origin_Analysis/heatmap_data_all_Clusters_Reads.rds")){
  
  clust_plot_object <- prepareClusterHeatmap(data.frame(acd_m[,2:length(acd_m)])
                                             , simprof_expected = 100
                                             , simprof_simulated = 100
                                             , distance_metric = Manly
  )
  saveRDS(clust_plot_object
          , file = "CAZy/04_CAZyme_Origin_Analysis/heatmap_data_all_Clusters_Reads.rds")
} else { 
  clust_plot_object <- readRDS("CAZy/04_CAZyme_Origin_Analysis/heatmap_data_all_Clusters_Reads.rds")
}

# Check all clusters are in the data: they are not. Why is this...?
unique(clust_plot_object$heatmap$heatmap_data$CAZyme)

# data preprocessing to get a nice order for the heatmap

j <- clust_plot_object
j$heatmap$heatmap_data
unique(j$heatmap$heatmap_data$CAZyme)

# Order by number of phyla with >=20$ of genes correlated with a metabolite in that cluster
nHighProportion <- c(unlist(by(j$heatmap$heatmap_data, j$heatmap$heatmap_data$CAZyme
                               , FUN = function(x){
                                 length(x$Abundance[x$Abundance> 0.2]) 
                               })))
metabolite <- sub("\\.Caz.*", "", names(nHighProportion))
metabolite <- sub("\\.", " ", metabolite)
of <- data.frame(metabolite = metabolite
                 , CAZyme = names(nHighProportion)
                 , nHighProportion = nHighProportion)
of$metabolite <- factor(of$metabolite, levels = c("glucose", "fucose", "xylose", "vanillic acid", "hydroxybenzoic acid", "benzoic acid"))

of <- rbindlist(by(of, of$metabolite, function(x){x[order(x$nHighProportion, decreasing = TRUE),]}))
of$column_order <- 1:length(of$metabolite)
#of$CAZyme <- sub("\\.Caz", "-Caz", as.character(of$CAZyme))
#of$CAZyme <- sub("\\.", " ", as.character(of$CAZyme))

column_order <- as.character(of$CAZyme) # input to function


# Make the plot
cp <- plotClusterHeatmap(clust_plot_object
                         , side_dend_width = 0.5
                         , dendrogram_left_margin = -0.7
                         , side_dendrogram_margin = c(0, -0.7, 0.2, 0)
                         , heatmap_margin = c(0.2, -0.7, 0.2, 0)
                         , x_axis_hjust = 1
                         , column_order = column_order
                         , heatmap_ylab = "NULL"
                         , colour_scale_title = "% Genes in Cluster\n(per Metabolite)")
cp$metabolitePlot



# # Abundant degradative taxa
# TaxaOfInterest <- c("Proteobacteria"
#                     ,"Actinobacteria"
#                     ,"Acidobacteria"
#                     ,"Bacteroidetes"
#                     ,"Firmicutes"
#                     ,"Planctomycetes"
#                     ,"Cyanobacteria	"
#                     ,"Verrucomicrobia"
#                     ,"Euryarchaeota"
#                     ,"Gemmatimonadetes"
#                     ,"Ascomycota"
#                     
#                     , "Thaumarchaeota"
#                     , "Nitrospirae"
#                     , "Crenarchaeota"
#                     , "Calditrichaeota"
#                     , "Spirochaetes"
#                     , "Ignavibacteridae"
#                     )
# # Ignavibacteridae
# # Spirochaetes
# # Calditrichaeota
# # Nitrospirae
# # Thaumarchaeota
# 
# ggplot(pc_evenness_all.m[variable %in% TaxaOfInterest] 
#        , aes(x = variable, y = value)) + 
#   stat_summary(aes(y = value)
#                , fun = mean, geom = "bar") +
#    stat_summary(fun.data = mean_cl_boot
#                , geom = "errorbar"
#                #, fun.args = list(mult = 1)
#                , width = 0.3
#                , size = 1.5) +
#   geom_point(position = position_jitter(0.1)) + 
#   ylim(0,1) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

