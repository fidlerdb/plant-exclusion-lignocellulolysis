

# Clusters which are abundant in abundant degradative phyla:
# 11 Phyla
# benzoic acid-CazC2
# glucose-CazC5

# 10 phyla
# 4-hydroxybenzoic acid-CazC7
# vanillic acid-CazC3

# 9 phyla
# xylose-CazC5
# fucose-CazC3
# fucose-CazC2

## (1) Which phyla are they in?
## (2) What CAZymes are in them? 
## (3) known activities of those CAZymes

rm(list = ls())

library(data.table)

cc <- fread("CAZy_2/05-10_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")
cc

UbiqClusters <- c("benzoic acid-CazC2","glucose-CazC5","4-hydroxybenzoic acid-CazC7"
                  ,"vanillic acid-CazC3","xylose-CazC5","fucose-CazC3","fucose-CazC2")

ADP <- c("Proteobacteria"
, "Actinobacteria"
, "Acidobacteria"
, "Planctomycetes"
, "Cyanobacteria"
, "Gemmatimonadetes"
, "Firmicutes"
, "Verrucomicrobia"
, "Bacteroidetes"
, "Euryarchaeota"
, "Ascomycota")

cc <- cc[ClusterName %in% UbiqClusters,]

#### (1) Which phyla have high abundances of the near-ubiquitous clusters in? ####

clust_plot_object <- readRDS("CAZy_2/outputData/heatmap_data/heatmap_data_all_Clusters.rds")

# Remove galactose clusters
clust_plot_object$heatmap$heatmap_data <- data.frame(data.table(clust_plot_object$heatmap$heatmap_data)[-grep("galactose", CAZyme),])
clust_plot_object$heatmap$CAZyme_pos_table <- data.frame(data.table(clust_plot_object$heatmap$CAZyme_pos_table)[-grep("galactose", CAZyme),])

j <- data.table(clust_plot_object$heatmap$heatmap_data)
j <- j[CAZyme %in% UbiqClusters & Phylum %in% ADP,]

head(j)

lapply(UbiqClusters, FUN = function(x){
  obj <- j[CAZyme == x, c(1:3)]
  obj$Cluster <- x
  obj$Abundance <- round(obj$Abundance, 2)
  return(obj[order(Abundance),])
})

#### (2) What CAZymes are in the shared clusters? ####

Percent_matches <- function(cluster1, cluster2){

#cluster1 = "xylose-CazC5"
#cluster2 = "fucose-CazC3"
 
  # Boolean data frame of CAZymes in abundant clusters
  comp_df <- cc[ClusterName %in% c(cluster1, cluster2),]
  comp_mat <- comp_df[,c(2,4)]
  
  comp_mat <- dcast(comp_mat, ClusterName ~ CAZyme, fill = 0)
  Clusters <- comp_mat$ClusterName
  comp_mat <- as.matrix(comp_mat[,2:length(comp_mat)])
  comp_mat[comp_mat != 0] <- 1
  comp_mat <- matrix(as.numeric(comp_mat), nrow = 2)

  # Remove non-occuring CAZy families
  #comp_mat <- comp_mat[,colSums(comp_mat) > 0]
  
  # Find shared CAZy families (colSums == 2)
  if(is.vector(comp_mat[,colSums(comp_mat) == 2])){nmatches <- 1} else {
    nmatches <- length(which(colSums(comp_mat) == 2))
  }
  
  # Nice interpretable output
  return(data.table(Cluster = Clusters
                    , Compared_With = rev(Clusters)
                    , Percent_Matches = 100*(nmatches / rowSums(comp_mat))
                    , nMatches = paste0(nmatches, " / ", rowSums(comp_mat))
  ))
  
}

## Clusters shared by all abundant degradative phyla ##

cc[ClusterName == "benzoic acid-CazC2",]$CAZyme
cc[ClusterName == "glucose-CazC5",]$CAZyme

## Clusters shared by all but one abundant degradative phyla ##

# Very compositionally similar
cc[ClusterName == "4-hydroxybenzoic acid-CazC7",]$CAZyme # >= 20% for all
cc[ClusterName == "vanillic acid-CazC3",]$CAZyme # Low abundance in Gemmatimonadetes - basically the same as 4-hydroxybenzoic acid-CazC7

## Clusters shared by all but two abundant degradative phyla ##

cc[ClusterName == "xylose-CazC5",]$CAZyme # Ascomycota 19.61%
cc[ClusterName == "fucose-CazC3",]$CAZyme # Gemmatimonadetes and Ascomycota ~18% - similar to 4-hydroxybenzoic acid-CazC7
cc[ClusterName == "fucose-CazC2",]$CAZyme # Lower abundance in Proteobacteria (15%), slightly lower abundance in Planctomycetes (19.88%)

# Create a comparison function

comparisons <- rbindlist(lapply(UbiqClusters, FUN = function(x){
  out <- lapply(UbiqClusters, function(y){
    Percent_matches(x, y)
    
  })
  return(rbindlist(out))
}))

comparisons

comparisons[Cluster == "benzoic acid-CazC2",] # no similarities :)
comparisons[Cluster == "glucose-CazC5",] # no similarities :)

comparisons[Cluster == "4-hydroxybenzoic acid-CazC7",] # Similar to vanillic acid-CazC3, xylose-CazC5, fucose-CazC3 
comparisons[Cluster == "vanillic acid-CazC3",] # Similar to 4-hydroxybenzoic acid-CazC7, xylose-CazC5, fucose-CazC3
comparisons[Cluster == "xylose-CazC5",] # Similar to 4-hydroxybenzoic acid-CazC7, xylose-CazC5, fucose-CazC3
comparisons[Cluster == "fucose-CazC3",] # Similar to xylose-CazC5, vanillic acid-CazC3, 4-hydroxybenzoic acid-CazC7

comparisons[Cluster == "fucose-CazC2",] # no similarities :)

# So for the nearly ubiquitous CAZyme clusters, "fucose-CazC3", "vanillic acid-CazC3", and "4-hydroxybenzoic acid-CazC7"
# Were made of the same things.

#### (4) Evenness of reads in ubiquitous clusters ####

library(doBy)

pc_evenness_all <- fread("CAZy_2/outputData/Metabolite-Correlated-CAZyme_Evenness_AllPhyla.csv")
pc_evenness_all.m <- melt(pc_evenness_all, id = "Sample")
pc_evenness_all.m 
pc_evenness_all.m <- pc_evenness_all.m[complete.cases(pc_evenness_all.m),]
sortedPhyla <- names(sort(tapply(pc_evenness_all.m$value
                                 , pc_evenness_all.m$variable, mean)
                          , decreasing = TRUE))
pc_evenness_all.m$variable <- factor(pc_evenness_all.m$variable
                                     , levels = sortedPhyla)
names(pc_evenness_all.m) <- c("Sample", "Phylum", "Evenness")

cre <- summaryBy(Evenness ~ Phylum, pc_evenness_all.m, keep.names = TRUE)
cre <- cre[Phylum %in% ADP,]
cre[order(Evenness),]

#              Phylum  Evenness
# 1:   Actinobacteria 0.8012103
# 2:   Planctomycetes 0.8060255
# 3: Gemmatimonadetes 0.8109070
# 4:   Proteobacteria 0.8131088
# 5:       Firmicutes 0.8176928
# 6:  Verrucomicrobia 0.8177666
# 7:    Cyanobacteria 0.8512564
# 8:    Bacteroidetes 0.8637262
# 9:    Euryarchaeota 0.8742126
# 10:   Acidobacteria 0.8832011
# 11:      Ascomycota 0.9270765

#### (5) Know activities from CAZyme clusters ####

ca <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData/All_Correlated_CAZy_Families_Activities.csv")
head(ca)

cc[ClusterName == "benzoic acid-CazC2",] 
# All but one have documented activities relating to hemicellulose degradation
# CBM11 -- coniferyl binding? Lignin degradation
ca[CAZyme %in% cc[ClusterName == "benzoic acid-CazC2",]$CAZyme,][1,]
ca[CAZyme %in% cc[ClusterName == "benzoic acid-CazC2",]$CAZyme,][2,]
ca[CAZyme %in% cc[ClusterName == "benzoic acid-CazC2",]$CAZyme,][3,]
ca[CAZyme %in% cc[ClusterName == "benzoic acid-CazC2",]$CAZyme,][4,]
ca[CAZyme %in% cc[ClusterName == "benzoic acid-CazC2",]$CAZyme,][5,]
ca[CAZyme %in% cc[ClusterName == "benzoic acid-CazC2",]$CAZyme,][6,]
ca[CAZyme %in% cc[ClusterName == "benzoic acid-CazC2",]$CAZyme,][7,]
ca[CAZyme %in% cc[ClusterName == "benzoic acid-CazC2",]$CAZyme,][8,]


cc[ClusterName == "glucose-CazC5",]
# GH6 - hemicellulose and cellulose
# GH44 - hemicellulose and cellulose
# GH8 - hemicellulose and cellulose
# GH43_12 - hemicellulose
ca[CAZyme %in% cc[ClusterName == "glucose-CazC5",]$CAZyme,][1,]
ca[CAZyme %in% cc[ClusterName == "glucose-CazC5",]$CAZyme,][2,]
ca[CAZyme %in% cc[ClusterName == "glucose-CazC5",]$CAZyme,][3,]
ca[CAZyme %in% cc[ClusterName == "glucose-CazC5",]$CAZyme,][4,]
ca[CAZyme %in% cc[ClusterName == "glucose-CazC5",]$CAZyme,][5,]
ca[CAZyme %in% cc[ClusterName == "glucose-CazC5",]$CAZyme,][6,]
ca[CAZyme %in% cc[ClusterName == "glucose-CazC5",]$CAZyme,][7,]
ca[CAZyme %in% cc[ClusterName == "glucose-CazC5",]$CAZyme,][8,]


cc[ClusterName == "fucose-CazC2",]
# AA6 - lignin
# GH130 - hemicellulose
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][1,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][2,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][3,] #  -- clearly a lignin degradation capability
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][4,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][5,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][6,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][7,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][8,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][9,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][10,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][11,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][12,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][13,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][14,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][15,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][16,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][17,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][18,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][19,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][20,]
ca[CAZyme %in% cc[ClusterName == "fucose-CazC2",]$CAZyme,][21,]

# 5 cellulolytic
# 11 hemicellulolytic
cc[ClusterName == "fucose-CazC2",]

# Similar clusters
cc[ClusterName == "vanillic acid-CazC3",] 
cc[ClusterName == "xylose-CazC5",]
cc[ClusterName == "fucose-CazC3",]
cc[ClusterName == "4-hydroxybenzoic acid-CazC7",] 

similarClusterCAZymes <- unique(cc[ClusterName == "vanillic acid-CazC3" |
     ClusterName == "xylose-CazC5" |
     ClusterName == "fucose-CazC3" |
     ClusterName == "4-hydroxybenzoic acid-CazC7",]$CAZyme)

similarClusterCAZymes

ca[ca$CAZyme %in% similarClusterCAZymes,]
ca[CAZyme %in% similarClusterCAZymes,][1,]
ca[CAZyme %in% similarClusterCAZymes,][2,]
ca[CAZyme %in% similarClusterCAZymes,][3,] #  -- clearly a lignin degradation capability
ca[CAZyme %in% similarClusterCAZymes,][4,]
ca[CAZyme %in% similarClusterCAZymes,][5,]
ca[CAZyme %in% similarClusterCAZymes,][6,]
ca[CAZyme %in% similarClusterCAZymes,][7,]
ca[CAZyme %in% similarClusterCAZymes,][8,]
ca[CAZyme %in% similarClusterCAZymes,][9,]
ca[CAZyme %in% similarClusterCAZymes,][10,]
