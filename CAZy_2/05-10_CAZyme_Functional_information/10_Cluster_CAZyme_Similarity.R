rm(list = ls()) 

# Load required packages
library(data.table)

#cc <- fread("CAZy/05_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")
cc <- fread("CAZy_2/05-10_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")

unique(cc$ClusterName)

cc$ClusterName <- sub("4-", "", cc$ClusterName)
cc$ClusterName <- sub("3,6-anhydro-D-", "", cc$ClusterName)

ac <- c("glucose-CazC1",
"glucose-CazC2",
"glucose-CazC3",
"fucose-CazC1",
"fucose-CazC6",
"fucose-CazC7",
"xylose-CazC10",
"xylose-CazC9",
"vanillic acid-CazC1",
"vanillic acid-CazC2",
"vanillic acid-CazC3",
"vanillic acid-CazC4",
"hydroxybenzoic acid-CazC1",
"hydroxybenzoic acid-CazC5",
"hydroxybenzoic acid-CazC6",
"hydroxybenzoic acid-CazC7",
"benzoic acid-CazC1",
"benzoic acid-CazC2")
ac

# Make a matrix to compare the CAZyme clusters
cca <- cc[ClusterName %in% ac,]
cca <- dcast(cca, ClusterName ~ CAZyme)
cca[is.na(cca)] <- 0
cm <- as.matrix(cca[,2:length(cca)])
cm[cm!=0] <- 1
class(cm) <- "numeric"
rownames(cm) <- cca$Cluster
cca <- data.table(Cluster = cca$Cluster, cm)

heatmap(cm)
ccl <- melt(cca, id = "Cluster")

# Binary similarity
jacc <- dist(cm, method = "binary", diag = TRUE, upper = TRUE
             )
jacc <- as.data.table(as.matrix(jacc))
jaccsim <- 1 - jacc

#jaccsim$names <- rownames(cm)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

lwr.tri <- get_lower_tri(jaccsim)
rownames(lwr.tri) <- rownames(cm)
lwr.tri$names <- rownames(cm)
lwr.tri.m <- melt(lwr.tri, id = "names", na.rm = TRUE)

rownames(jaccsim) <- rownames(cm)
jaccsim$names <- rownames(cm)
jaccsim.m <- melt(jaccsim, id = "names", na.rm = TRUE)




#jacc.m <- melt(jaccsim, id = "names")
#jacc.m$names <- factor(jacc.m$names, rownames(cm))
#jacc.m$variable <- factor(jacc.m$variable, rev(rownames(cm)))
#sort the data frame
#jacc.m <- plyr::arrange(jacc.m, variable, plyr::desc(names))

# ggplot(jacc.m, aes(x = names, y = variable, fill = value)) + 
#   geom_tile() +
#   geom_text(aes(label = round(value*100, 0))) +
#   scale_fill_gradient(low = "white", high = "black") +
#   theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1))

ggplot(lwr.tri.m, aes(x = names, y = variable, fill = value)) + 
  geom_tile() +
  geom_text(aes(label = round(value*100, 0))) +
  scale_fill_gradient(low = "white", high = "black") +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1))


# jacc.m[value > 0 & value != 1,]
# jacc.m %>% 
#   as.matrix %>%
#   cor %>%
#   as.data.frame %>%
#   rownames_to_column(var = 'var1') %>%
#   gather(var2, value, -var1) %>%
#   filter(value > 0.8 | value < -0.8) %>%
#   filter(value != 1) %>%
#   filter(var1 < var2)
# 
# lwr.tri.m[value > 0]

#######

# Percentage of matching CAZy families

cca[Cluster %in% ac[1:2],]

Percent_matches <- function(cluster1, cluster2){
  
  # Boolean data frame of CAZymes in abundant clusters
  comp_df <- cca[Cluster %in% c(cluster1, cluster2),]
  comp_mat <- data.frame(comp_df[,2:length(comp_df)])
  
  # Remove non-occuring CAZy families
  comp_mat <- comp_mat[,colSums(comp_mat) > 0]
  
  # Find shared CAZy families (colSums == 2)
  if(is.vector(comp_mat[,colSums(comp_mat) == 2])){nmatches <- 1} else {
    nmatches <- length(comp_mat[,colSums(comp_mat) == 2])
  }
  
  # Nice interpretable output
  return(data.table(Cluster = comp_df$Cluster
             , Compared_With = rev(comp_df$Cluster)
             , Percent_Matches = 100*(nmatches / rowSums(comp_mat))
             , nMatches = paste0(nmatches, " / ", rowSums(comp_mat))
  ))
  
}

#Percent_matches("fucose-CazC6", "glucose-CazC1") # For demonstration (% matches is matches )
Percent_matches("fucose-CazC1", "glucose-CazC1") # For demonstration (% matches is matches )


# Take all cluster comparisons with a shared CAZyme
ltm <- data.frame(lwr.tri.m[value > 0 & value != 1, 1:2])
jtm <- data.frame(jaccsim.m[value > 0 & value != 1, 1:2])

# What were all of the similarities between abundant CAZyme clusters which shared a CAZy family?

SimilarityResults <- rbindlist(apply(jtm, 1, FUN = function(j){Percent_matches(j[1], j[2])}))

SimilarityResults
SimilarityResults <- SimilarityResults[order(Percent_Matches, decreasing = TRUE)]
SimilarityResults <- unique(SimilarityResults)
SimilarityResults

# Most widely distributed groups
library(dplyr)
Percent_matches("glucose-CazC3", "hydroxybenzoic acid-CazC1") # No matches
Percent_matches("glucose-CazC3", "benzoic acid-CazC2") # No matches
Percent_matches("benzoic acid-CazC2 ", "hydroxybenzoic acid-CazC1") # No matches


cca[Cluster == "hydroxybenzoic acid-CazC1", ] %>% select_if(~ !is.numeric(.) || sum(.) != 0)
cca[Cluster == "benzoic acid-CazC2", ] %>% select_if(~ !is.numeric(.) || sum(.) != 0)
cca[Cluster == "vanillic acid-CazC3", ] %>% select_if(~ !is.numeric(.) || sum(.) != 0)


##### START FROM HERE ####

# Import CAZy family and cluster activity information
at <- fread("CAZy/05_CAZyme_Functional_information/OutputData/CAZy_Family_Activity_Summary.csv")
atc <- fread("CAZy/05_CAZyme_Functional_information/OutputData/CAZy_Cluster_Activity_Summary.csv")


at

#  Investigate glu-CC3
cca[Cluster == "glucose-CazC3", ] %>% select_if(~ !is.numeric(.) || sum(.) != 0)
atc[Cluster == "glucose-CazC3",]
SimilarityResults[Cluster == "glucose-CazC3",]

# Investigate hba-CC1
cca[Cluster == "hydroxybenzoic acid-CazC1", ] %>% select_if(~ !is.numeric(.) || sum(.) != 0)
atc[Cluster == "hydroxybenzoic acid-CazC1",]
at[CAZyme == "GH130",][1,]
at[CAZyme == "GH13_11",][1,]
at[CAZyme == "GH77",][1,]
SimilarityResults[Cluster == "hydroxybenzoic acid-CazC1",]
cca[Cluster == "xylose-CazC10", ] %>% select_if(~ !is.numeric(.) || sum(.) != 0)

# Investigate bza-CC2
cca[Cluster == "benzoic acid-CazC2", ] %>% select_if(~ !is.numeric(.) || sum(.) != 0)
atc[Cluster == "benzoic acid-CazC2",]
at[CAZyme == "CBM38",][1,]
at[CAZyme == "GH141",][1,]
at[CAZyme == "GH147",][1,]
at[CAZyme == "GH30_2",][1,]
at[CAZyme == "GH30_4",][1,]
at[CAZyme == "GH43",][1,]
at[CAZyme == "GH5_26",][1,]
at[CAZyme == "GH5_27",][1,]
SimilarityResults[Cluster == "benzoic acid-CazC2",]
cca[Cluster == "glucose-CazC2", ] %>% select_if(~ !is.numeric(.) || sum(.) != 0)
Percent_matches("glucose-CazC2", "benzoic acid-CazC2") # The other way round

#test <- "AA9" # LPMOs
#cca[,..test] # Just understanding data.table a little better

###############################
# Create a coocurance matrix
# coocCounts <- t(cm) %*% cm
# require(igraph)
# 
# net <- graph_from_adjacency_matrix(coocCounts, weighted = TRUE, mode = "undirected")
# net <- simplify(net, remove.loops = TRUE
#                 , remove.multiple = TRUE
#                 , edge.attr.comb = c(weight ="sum") 
#                 )
# E(net)$width <- E(net)$weight^2
# E(net)$color <- "orange"
# V(net)$color <- "grey50"
# V(net)$frame.color <- "#ffffff"
# V(net)$label.color <- "black"
# V(net)$size <- rowSums(coocCounts)/2
# 
# lf <- layout.fruchterman.reingold(net)
# lg <- layout.grid(net)
# lc<- layout.circle(net)
# 
# set.seed(10091)
# plot(net)
# plot(net, layout = lf, vertex.label.cex = 0.8,
#      vertex.shape = "circle",
#      vertex.label.dist = 0.5)
# plot(net, layout = lg)
# plot(net, layout = lc)
# 
# 
# links <- E(net)
# cut.off <- mean(links$weight)
# cut.off.size <- mean(V(net)$size)
# 
# 
# net.sp <- delete.edges(net, E(net)[weight<cut.off])
# net.sp <- delete.vertices(net.sp, V(net.sp)[size<cut.off.size])
# #net.sp <- delete.vertices(net, V(net)[size<cut.off.size])
# 
# #net.sp <- delete.vertices(net, E(net)[weight<cut.off])
# 
# lf.sp <- layout.fruchterman.reingold(net.sp)
# lg.sp <- layout.grid(net.sp)
# lc.sp <- layout.circle(net.sp)
# ls.sp <- layout.star(net.sp)
# 
# 
# set.seed(100)
# plot(net.sp, layout=lf.sp)
# plot(net.sp, layout = lg.sp)
# plot(net.sp, layout = lc.sp)
# plot(net.sp, layout = ls.sp)


