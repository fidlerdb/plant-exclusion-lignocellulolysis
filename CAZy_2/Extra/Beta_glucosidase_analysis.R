
#### The purpose of this script is to understand if there was a treatment effect 
#### on beta glucosidases. We would expect higher beta glucosidase abundance
#### (but probably not richness) in the 10-year grasslands.

library(data.table)
library(rvest) # Use this to scrape CAZy.org for all of the CAZy families with beta glucosidase activity

#### Import the activity data and keep only beta glucosidase ECs ####
ac <- readxl::read_xlsx("CAZy_2/05-10_CAZyme_Functional_information/InputData/All_Unique_Activities_Substrates_WithUTF-8.xlsx")
head(ac)
setDT(ac)

# Keep only glucosidases
bglu <- ac[grep("-glucosidase", Activity, fixed = TRUE),]

# ECs to remove (not beta glucosidases)
aglus <- c("3.2.1.84", "3.2.1.48", "3.2.1.33", "3.2.1.20"
  , "3.2.1.107", "3.2.1.106", "3.2.1.10")

# Keep only the relevant enzyme activities
bglu <- bglu[!EC %in% aglus,]


#### Scrape CAZy.org ####

address <- paste0("http://www.cazy.org//search?page=recherche&lang=en&recherche=", bglu$EC[1], "&tag=9")

simple <- read_html(address)
html_nodes()

AllFamilies <- lapply(bglu$Activity[1], FUN = function(x){
  
  #address <- paste0("http://www.cazy.org/", x, ".html")
  address <- paste0("http://www.cazy.org//search?page=recherche&lang=en&recherche=", x, "&tag=9")
  
  simple <- read_html(address)
  
  CAZyHTML <- simple %>%
    html_nodes("body") %>%
    html_text()
  
  Families <- CAZyHTML#[9]
  #Activities <- sub("Activities in Family", "", Activities)
  
  return(Families)
})

FamilyActivities <- data.table(CAZyme = AllCAZymes$Family, Activity = unlist(AllActivities))

FamilyActivities$Activity <- enc2utf8(FamilyActivities$Activity) # Ensure special characters don't go weird
