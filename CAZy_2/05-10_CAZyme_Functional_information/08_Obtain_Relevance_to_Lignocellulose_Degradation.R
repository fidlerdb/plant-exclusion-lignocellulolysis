# Clear the workspace
# .rs.restartR() # New R session in case something isn't
rm(list = ls()) 

# Load required packages
library(data.table)

##### Import data #####

# ca = CAZyme Activities -- a list of which enzyme type has activities on the substrates I care about
# Encoding not right. Had to save as .xls and use dirty tidyverse code ;)
ca <- readxl::read_excel("CAZy_2/05-10_CAZyme_Functional_information/InputData/All_Unique_Activities_Substrates_WithUTF-8.xlsx")
ca <- data.table(ca)

colSums(ca[,3:length(ca)])

# cf = CAZyme Functions. List of functions CAZy identiies for each CAZy family 
#cf <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData/All_Correlated_CAZy_Families_Activities.csv"
#            , encoding = "UTF-8")
#cf

ccf <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData/CAZy_Clusters_Activities.csv"
            , encoding = "UTF-8")
ccf


# cc = CAZyme Clusters
#cc <- fread("CAZy_2/05-10_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")

#cc <- fread("CAZy/05_CAZyme_Functional_information/InputData/CAZyme_Clusters_All_Metabolites.csv")
#cc

# Put functional and cluster and cazyme info into one data frame
#ccf <- merge(cc, cf , all.x = TRUE, by = "CAZyme")
#ccf

#### Create a table for further analysis ####
# Keep grep happy
ca$Activity <- gsub("[[:space:]]", "", ca$Activity)
ccf$Activity <- gsub("[[:space:]]", "", ccf$Activity)

#ca$Activity
#ccf$Activity
any(ca$Activity %in% ccf$Activity)

####

ccf[grep("\\|", ccf$CAZyme),]
ccf[1,]
####

# Create an activity table for each CAZyme in each cluster
at <- rbindlist(
  # For every row of ccf (CAZyme, Cluster, Function)
  apply(ccf, 1, FUN = function(x){
  # Get all of the inidvidual activities belonging to this CAZy family
  CAZy_Activities <- unlist(strsplit(x[5], ";"))
  
  # Find the rows of the reference table where the activities are matched
  # Return the sum of how many activities were noted for the CAZy family
  Activity_Totals <- colSums(ca[ca$Activity %in% CAZy_Activities, 3:6])
  
  
  # Create a table of: CAZyme cluster, CAZyme, number of activities 
  # related to one of the constituents of lignocellulose that I care about
  return(data.table(Cluster = x[4], CAZyme = x[1]
                    , Cellulose = Activity_Totals[1]
                    , Hemicellulose = Activity_Totals[2]
                    , Lignin = Activity_Totals[3]
                    , Oligosaccharides = Activity_Totals[4]
                    ))
}))

at
colSums(at[,3:6])

# Make all of the columns Boolean
at$Cellulose <- ifelse(at$Cellulose != 0, yes = 1, no = 0)
at$Hemicellulose <- ifelse(at$Hemicellulose != 0, yes = 1, no = 0)
at$Lignin <- ifelse(at$Lignin != 0, yes = 1, no = 0)
at$Oligosaccharides <- ifelse(at$Oligosaccharides != 0, yes = 1, no = 0)

at
# at[grep("galactose", Cluster),]

fwrite(at, "CAZy_2/05-10_CAZyme_Functional_information/OutputData/CAZy_Family_Activity_Summary.csv")

# Find proportion of CAZymes with each of the important activities in each cluster
ct <- rbindlist(
  by(at, INDICES = at$Cluster, FUN = function(x){
    boolSums <- data.table(Cluster = x$Cluster[1]
                           , nFamilies = nrow(x)
                           , Percentage = round(c(colSums(x[,3:6]*100)/length(x$Cluster)),2))
    }))
ct$Substrate <- rep(names(at)[3:6], length.out = nrow(ct))

ct # previously documented activities (summarised)
ct[grep("glucose", Cluster),]
ct[grep("fucose", Cluster),]
ct[grep("xylose", Cluster),]
ct[grep("galactose", Cluster),]
ct[grep("van", Cluster),]
ct[grep("hydroxybenz", Cluster),]
ct[grep("^benzoic", Cluster),]

at[grep("fucose-CazC6", Cluster),]
at[grep("xylose-CazC10", Cluster),]
at[grep("galactose", Cluster),]
at[grep("vanillic acid-CazC3", Cluster),]
at[grep("hydroxybenz", Cluster),]
at[grep("hydroxybenzoic acid-CazC5", Cluster),]
at[grep("hydroxybenz", Cluster),]
at[grep("hydroxybenz", Cluster),]


# Commit this to a file
fwrite(ct, "CAZy_2/05-10_CAZyme_Functional_information/OutputData/CAZy_Cluster_Activity_Summary.csv")

