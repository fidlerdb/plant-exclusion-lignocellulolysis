
rm(list = ls()) # Clear the workspace

library(data.table)
library(roperators)

cf <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData/All_Correlated_CAZy_Families_Activities.csv"
                          , encoding = "UTF-8")

cf_old <- fread("CAZy_2/05-10_CAZyme_Functional_information/OutputData-old/All_Correlated_CAZy_Families_Activities.csv"
            , encoding = "UTF-8")

head(cf)

#####

# Make the table horribly long. The idea is to then find unique 
# enzyme classifications, identify meaningful broad classifications
# and to link them back to the GH families and clusters

# splitoutput <- strsplit(as.character(cf$Activity), ";")
# unlist(splitoutput)

# Create a data frame of all CAZymes and the activities they have
# (One activity per row) 
cf2 <- rbindlist(apply(cf, 1, FUN = function(x){
  # Find all individual activities
  ActivityList <- strsplit(x[2], ";")
  n <- length(ActivityList)
  
  # Create a long table of these
  return(data.table(CAZy_Family = rep(x[1], n)
                    , Activity = unlist(ActivityList)
  ))
}))

# Strip leading white space
cf2$Activity <- sub("^\\s+", "", cf2$Activity)

# Extract the EC number
EC <- sub("^.*\\(EC", "", unique(cf2$Activity))
EC <- sub("\\)$", "", EC)
EC <- sub("[[:space:]]", "", EC)

ActivityEC <- data.table(Activity = unique(cf2$Activity)
                         , EC = EC
                         )

#### Check if there are any new activities with the uncombining of the CAZy domains ####
cf2_old <- rbindlist(apply(cf_old, 1, FUN = function(x){
  # Find all individual activities
  ActivityList <- strsplit(x[2], ";")
  n <- length(ActivityList)
  
  # Create a long table of these
  return(data.table(CAZy_Family = rep(x[1], n)
                    , Activity = unlist(ActivityList)
  ))
}))

# Strip leading white space
cf2_old$Activity <- sub("^\\s+", "", cf2_old$Activity)

# Extract the EC number
EC_old <- sub("^.*\\(EC", "", unique(cf2_old$Activity))
EC_old <- sub("\\)$", "", EC_old)
EC_old <- sub("[[:space:]]", "", EC_old)

ActivityEC_old <- data.table(Activity = unique(cf2_old$Activity)
                         , EC = EC_old
)

# Add these to the file All_Unique_Activities_Substrates_WithUTF-8.xlsx


OnesToAdd <- ActivityEC[ActivityEC$Activity %ni% ActivityEC_old$Activity]
#ActivityEC[ActivityEC$EC %ni% ActivityEC_old$EC]

ExtraRows <- data.table(Activity = c("3.2.1.8", "3.2.1.21")
           , EC = c("3.2.1.8", "3.2.1.21"))

RowsToAdd <- rbind(OnesToAdd, ExtraRows)
RowsToAdd

RowsToAdd$Activity <- enc2utf8(RowsToAdd$Activity)
RowsToAdd$EC <- enc2utf8(RowsToAdd$EC)

fwrite(RowsToAdd, "CAZy_2/05-10_CAZyme_Functional_information/OutputData/NewActivities_UTF-8.csv"
       , bom = TRUE)
#####
# ActivityEC$Activity <- enc2utf8(ActivityEC$Activity) # Ensure special characters don't go weird
# ActivityEC$EC <- enc2utf8(ActivityEC$EC) # Ensure special characters don't go weird
# 
# str(ActivityEC)
# Encoding(ActivityEC$Activity)
# ActivityEC[1,]
# # Commit this to a file
# # Save all of the unique enzyme activities for manual "what does it mean for my dataset"-ing
# fwrite(ActivityEC, "CAZy/05_CAZyme_Functional_information/OutputData/All_Unique_Activities.csv", bom = TRUE)

