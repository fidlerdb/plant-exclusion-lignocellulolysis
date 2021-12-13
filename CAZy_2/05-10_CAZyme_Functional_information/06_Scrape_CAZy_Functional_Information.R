# Look at the activities of all of the CAZy families in the 
# all CAZyme clusters

# CAZyme clusters
rm(list = ls())

# install.packages("rvest")
library(data.table)
library(rvest) # Use this to scrape CAZy.org for all of the CAZy families in my dataset 

##### (1) Import the data #####
cf <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData-old/All_Correlated_CAZy_Families_Activities.csv")

AllClusters <- fread("CAZy_2/05-10_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")
AllCAZymes <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData/All_Correlated_CAZy_Families.csv")
AllCAZymes$Family # Check there are no compound families
AllCAZymes <- AllCAZymes[-grep("^[0-9]", AllCAZymes$Family),] # Remove CAZymes only with EC. numbers


# Test
# address <- paste0("http://www.cazy.org/", "GH5_46", ".html")
# 
# simple <- read_html(address)
# 
# CAZyHTML <- simple %>%
#   html_nodes(":contains(Activities in Family)") %>%
#   html_text()
# 
# Activities <- CAZyHTML[9]
# Activities <- sub("Activities in Family", "", Activities)


# Scrape CAZy
AllActivities <- lapply(AllCAZymes$Family, FUN = function(x){
  
  address <- paste0("http://www.cazy.org/", x, ".html")
  
  simple <- read_html(address)
  
  CAZyHTML <- simple %>%
    html_nodes(":contains(Activities in Family)") %>%
    html_text()
  
  Activities <- CAZyHTML[9]
  Activities <- sub("Activities in Family", "", Activities)
  
  return(Activities)
})

FamilyActivities <- data.table(CAZyme = AllCAZymes$Family, Activity = unlist(AllActivities))

FamilyActivities$Activity <- enc2utf8(FamilyActivities$Activity) # Ensure special characters don't go weird

# Commit this to a file
fwrite(FamilyActivities, "CAZy_2/05-10_CAZyme_Functional_information/OutputData/All_Correlated_CAZy_Families_Activities.csv", bom = TRUE)

#### (3) IMPORTANT BY-HAND STEP ####

# Now search BRENDA to add columns indicating activity on cellulose, hemicellulose, lignin
# , oligosaccharides -- do this by hand.

# To make my life much easier, these are the CAZy families that weren't in the previous version of the file 
# [1] "CBM22"    "GH2"      "CBM6"     "GH10"     "3.2.1.8"  "GH3"      "3.2.1.21" "CBM48"    "GH13_14"  "GH14"    
# [11] "CE4"      "GT2"      "GH5_1"    "GH5_22"   "GH62"     "CBM35"    "CBM51"    "GH35"     "CBM11" 

##### (4) Make files which the updated file just created can be merged to #####

rm(list = ls())

FamilyActivities <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData/All_Correlated_CAZy_Families_Activities.csv", encoding = "UTF-8")
AllClusters <- fread("CAZy_2/05-10_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")
AllCAZymes <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData/All_Correlated_CAZy_Families.csv")

FamilyActivities
AllClusters
AllCAZymes

# For merging with:
# "CAZy/05_CAZyme_Functional_information/OutputData/All_Correlated_CAZy_Families_Activities.csv"

# Get all activities relating to the compound families, and put them all in one place
CompoundFamilies <- AllClusters$CAZyme[grep("\\|", AllClusters$CAZyme)]

CompoundActivities <- lapply(CompoundFamilies, FUN = function(i){
  apply(FamilyActivities, 1, FUN = function(x){ifelse(grep(x[1], i), yes = x[2], no = NULL)})
  })

CompoundActivities <- lapply(CompoundActivities, FUN = function(x){
  x <- paste(unlist(x))
  paste0(x, collapse = ";")
})


# Now merge all functional information with the cluster and CAZyme information
CAZyActivities <- rbind(FamilyActivities
      , data.table(CAZyme = CompoundFamilies, Activity = unlist(CompoundActivities))
           )
CAZyActivities <- merge(AllClusters, CAZyActivities, by = "CAZyme", all.x = TRUE)
CAZyActivities$Activity <- enc2utf8(CAZyActivities$Activity)

# Commit this to a file
fwrite(CAZyActivities, "CAZy_2/05-10_CAZyme_Functional_information/OutputData/CAZy_Clusters_Activities.csv", bom = TRUE)

