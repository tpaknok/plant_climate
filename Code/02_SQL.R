library(inborutils)
library(dplyr)
library(DBI)
#can replace Magnoliopsidae with liliopsida
csv_to_sqlite(csv_file = "GBIF_plants.csv", "GBIF_plants.sqlite", "Plants", delim="\t",pre_process_size = 1000, chunk_size = 50000, show_progress_bar = T)

con <- DBI::dbConnect(RSQLite::SQLite(), "C:/Users/pakno/OneDrive - University of Toronto/GLONAF/GBIF_plants.sqlite")
plants <- tbl(con, "Plants")

### SQL code - need to fetch it here before running some data cleaning function
subset_plants <- dbGetQuery(con,
                             'SELECT decimalLongitude,decimalLatitude,genus,species,issue,basisOfRecord
                             FROM Plants
                             WHERE coordinateUncertaintyInMeters <= 1000
                             AND (basisOfRecord = "HUMAN_OBSERVATION" 
                             OR basisOfRecord = "PRESERVED_SPECIMEN"
                             OR basisOfRecord = "OCCURRENCE"
                             OR basisOfRecord = "MATERIAL_SAMPLE"
                             OR basisOfRecord = "MATERIAL_CITATION")')

### dplyr version (just a part of it) ###
#subset_plants1 <- plants %>% 
  #filter(coordinateUncertaintyInMeters <= 1000) %>% 
  #select(decimalLongitude,decimalLatitude,genus,species,issue,basisOfRecord)

write.csv(subset_plants,"plants.csv")

