#03 - clean data
library(CoordinateCleaner)
plants <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/plants.csv")

plants_db_cleaned <- clean_coordinates(plants,
                                       lon="decimalLongitude",
                                       lat="decimalLatitude",
                                       tests = c("capitals","centroids","institutions","sea","zero","duplicates"),
                                       capitals_rad=100,
                                       centroids_rad=100,
                                       inst_rad = 100)


head(plants_db_cleaned)
plants_db_cleaned <- plants_db_cleaned[plants_db_cleaned$.summary,]
plants_db_cleaned <- plants_db_cleaned[,2:11]
plants_db_cleaned <- plants_db_cleaned[!is.na(plants_db_cleaned$species),]
write.csv(plants_db_cleaned,"plants_cleaned.csv")
