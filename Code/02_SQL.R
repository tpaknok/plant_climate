library(inborutils)
library(dplyr)

#can replace Magnoliopsidae with liliopsida
csv_to_sqlite(csv_file = "GBIF_plants.csv", "GBIF_plants.sqlite", "Plants", delim="\t",pre_process_size = 1000, chunk_size = 50000, show_progress_bar = T)

con <- DBI::dbConnect(RSQLite::SQLite(), "C:/Users/pakno/OneDrive - University of Toronto/GLONAF/GBIF_plants.sqlite")
plants <- tbl(con, "Plants")

subset_plants <- plants %>% 
  group_by(species) %>% 
  filter(coordinateUncertaintyInMeters <= 1000) %>% 
  select(decimalLongitude,decimalLatitude,genus,species,issue,basisOfRecord)

write.csv(subset_plants,"plants.csv")

plants <- read.csv("plants.csv")
