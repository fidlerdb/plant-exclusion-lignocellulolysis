# Look at the activities of all of the CAZy families in the 
# most abundant CAZyme clusters

rm(list = ls())

# Import the data
#cf <- fread("Metabolomics/Data/CAZyme_Cluster_Activities.csv")

cf <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData-old/All_Correlated_CAZy_Families_Activities.csv")

str(cf)
cf$CAZyme
#####
#AllClusters <- fread("CAZy/05_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")
AllClusters <- fread("CAZy_2/05-10_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")

AllCAZymes <- data.table(Family = unique(as.character(AllClusters$CAZyme)))
AllCAZymes <- data.table(Family = unique(unlist(strsplit(AllCAZymes$Family, split = "\\|"))))

# Unincluded families
AllCAZymes$Family[!AllCAZymes$Family %in% unique(cf$CAZyme)] 

# These are CAZy families which need that activities 
# finding for, with these being added to the list of activities in the handmade file Metabolomics/Data/CAZyme_Cluster_Activities.csv

# (1) get unique CAZy families
AllCAZymes$Family <- enc2utf8(AllCAZymes$Family)

fwrite(AllCAZymes, "CAZy_2/05-10_CAZyme_Functional_information/OutputData/All_Correlated_CAZy_Families.csv"
       , bom = TRUE)

# (2) go on CAZy.org and scrape all recorded activities (added to unique families file)
  # 02_Scrape_CAZy_Functional_Information.R

# (3) Merge with "CAZy/05_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv"
# (4) Update "CAZy/05_CAZyme_Functional_information/InputData/Unique_Activities_Substrates.csv"
# (5) Look at what the 


# #####
# 
# # Make the table horribly long. The idea is to then find unique 
# # enzyme classifications, identify meaningful broad classifications
# # and to link them back to the GH families and clusters
# 
# splitoutput <- strsplit(as.character(cf[4][,3]), ";")
# unlist(splitoutput)
# 
# # Create a data frame of all CAZymes and the activities they have
# # (One activity per row) 
# cf2 <- rbindlist(apply(cf, 1, FUN = function(x){
#   # Find all individual activities
#   ActivityList <- strsplit(x[3], ";")
#   n <- length(ActivityList)
#   
#   # Create a long table of these
#   return(data.table(CAZyme_Cluster = rep(x[1], n)
#              , CAZy_Family = rep(x[2], n)
#              , Activity = unlist(ActivityList)
#   ))
# }))
# 
# # Strip leading white space
# cf2$Activity <- sub("^\\s+", "", cf2$Activity)
# 
# # Extract the EC number
# ActivityEC <- data.table(Activity = unique(cf2$Activity)
#                          , EC = paste0("(EC", sub("^.*\\(EC"
#                                                , ""
#                                                , unique(cf2$Activity)
#                          ))
#                          )
# 
# # Save all of the unique enzyme activities for 
# # manual "what does it mean for my dataset"-ing
# fwrite(ActivityEC
#        , "CAZy/05_CAZyme_Functional_information/OutputData/Unique_Activities.csv")
