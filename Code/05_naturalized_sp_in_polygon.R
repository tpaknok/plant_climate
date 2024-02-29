#05 - check number of naturalized species in each polygon (useful for excluding unnecessary poylgon to minimize computation time)
#40144227
library(WorldFlora)
library(tidyverse)
library(terra)
plants_db_cleaned <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/plants_cleaned.csv")
label <- paste0(plants_db_cleaned$decimalLongitude,"_",plants_db_cleaned$decimalLatitude)
coord_df <- plants_db_cleaned[!duplicated(label),c("decimalLongitude","decimalLatitude")]

region <- vect("regions2.shp")
region_df <- read.csv("regionID.csv")
GLONAF <- read.csv("GLONAF.csv")
tdwg3 <- vect("tdwg3/level3.shp")

tdwg3_list <- NULL
region_list <- NULL
seg <-10000

for (i in 1: ceiling(nrow(plants_db_cleaned)/seg)) {
  message(i,"in", ceiling(nrow(plants_db_cleaned)/seg))
  start<- Sys.time()
  lower = 1+(i-1)*seg
  upper = ifelse(nrow(plants_db_cleaned) < seg*i,nrow(plants_db_cleaned),seg*i)
  region_list[[i]]<- terra::extract(region,plants_db_cleaned[lower:upper,c("decimalLongitude","decimalLatitude")])
  tdwg3_list[[i]]<- terra::extract(tdwg3,plants_db_cleaned[lower:upper,c("decimalLongitude","decimalLatitude")])
  end <- Sys.time()
  print(end-start)
}

###################################################################################\
plants_db_cleaned <- plants_db_cleaned[!is.na(plants_db_cleaned$polygon_number) & !is.na(plants_db_cleaned$tdwg3),]
plants_db_cleaned$new_species <- plants_db_cleaned$species

WFO.remember("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/classification.csv")
GLONAF <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/GLONAF.csv")

GLONAF_name_list <- unique(word(unique(GLONAF$standardized_name),1,2," "))
GLONAF_name_list_unmatched <- GLONAF_name_list[!GLONAF_name_list %in% WFO.data$scientificName]
GLONAF_name_list_not_accpeted <- WFO.data[WFO.data$scientificName %in% GLONAF_name_list & WFO.data$taxonomicStatus != "Accepted",]
GLONAF_not_accepted <- WFO.match(GLONAF_name_list_not_accpeted$scientificName,WFO.data=WFO.data,no.dates=T,counter=1,Fuzzy=0)
GLONAF_match <- WFO.match(GLONAF_name_list_unmatched,WFO.data=WFO.data,no.dates=T,counter=1,Fuzzy=0.1)
GLONAF_not_accepted_one <- WFO.one(GLONAF_not_accepted)
GLONAF_unmatched_one <- WFO.one(GLONAF_match)

########################################################################
GBIF_name_list <- unique(word(unique(plants_db_cleaned$species),1,2," ")) #faster
GBIF_name_list_unmatched <- GBIF_name_list[!GBIF_name_list %in% WFO.data$scientificName]
GBIF_list_not_accpeted <- WFO.data[WFO.data$scientificName %in% GBIF_name_list & WFO.data$taxonomicStatus != "Accepted",]
GBIF_match_not_accepted <- WFO.match(unique(GBIF_list_not_accpeted$scientificName),WFO.data=WFO.data,no.dates=T,counter=1,Fuzzy=0)
GBIF_match_unmatched <- WFO.match(GBIF_name_list_unmatched,WFO.data=WFO.data,no.dates=T,counter=1,Fuzzy=0.1)

#GBIF_match_not_accepted <- subset(GBIF_match_not_accepted,taxonRank == "species" | taxonRank == "form" | taxonRank == "subspecies" | taxonRank == "variety")
pos <- GBIF_match_not_accepted$taxonRank == "form" | GBIF_match_not_accepted$taxonRank == "subspecies" | GBIF_match_not_accepted$taxonRank == "variety"
GBIF_match_not_accepted$scientificName[pos] <- word(GBIF_match_not_accepted$scientificName[pos],1,2," ")
GBIF_not_accepted_one <- WFO.one(GBIF_match_not_accepted)

GBIF_unmatched_one <- WFO.one(GBIF_match_unmatched)
GBIF_unmatched_one <- subset(GBIF_match_unmatched,taxonRank == "species" | taxonRank == "form")
GBIF_unmatched_one$scientificName[GBIF_unmatched_one$taxonRank == "form"] <- word(GBIF_unmatched_one$scientificName[GBIF_unmatched_one$taxonRank == "form"],1,2," ")

pos <- plants_db_cleaned$new_species %in% GBIF_unmatched_one$spec.name
plants_db_cleaned[plants_db_cleaned$new_species %in% GBIF_unmatched_one$spec.name,"new_species"] <- GBIF_unmatched_one[na.omit(match(plants_db_cleaned$new_species,GBIF_unmatched_one$spec.name)),"scientificName"]

pos <- plants_db_cleaned$new_species %in% GBIF_match_not_accepted$spec.name
plants_db_cleaned[plants_db_cleaned$new_species %in% GBIF_match_not_accepted$spec.name,"new_species"] <- GBIF_match_not_accepted[na.omit(match(plants_db_cleaned$new_species,GBIF_match_not_accepted$spec.name)),"scientificName"]
plants_db_cleaned[!plants_db_cleaned$new_species %in% WFO.data$scientificName,"new_species"] <- NA

####################################################################
GLONAF$new_species <- word(GLONAF$standardized_name,1,2," ")
pos <- GLONAF$new_species %in% GLONAF_unmatched_one$spec.name
GLONAF[GLONAF$new_species %in% GLONAF_unmatched_one$spec.name,"new_species"] <- GLONAF_unmatched_one[na.omit(match(GLONAF$new_species,GLONAF_unmatched_one$spec.name)),"scientificName"]

pos <- GLONAF$new_species %in% GLONAF_not_accepted_one$spec.name
GLONAF[GLONAF$new_species %in% GLONAF_not_accepted_one$spec.name,"new_species"] <- GLONAF_not_accepted_one[na.omit(match(GLONAF$new_species,GLONAF_not_accepted_one$spec.name)),"scientificName"]
GLONAF[!GLONAF$new_species %in% WFO.data$scientificName,"new_species"] <- NA
GLONAF$new_species <- word(GLONAF$new_species,1,2, " ")
write.csv(GLONAF,"GLONAF_new.csv")
write.csv(plants_db_cleaned,"plants_cleaned_standardized.csv")
####################################### no need run these

library(sf)
unique_exotic <- unique(plants_db_cleaned$new_species)[unique(plants_db_cleaned$new_species) %in% unique(GLONAF$new_species)]
count <- table(GLONAF[GLONAF$status == "naturalized","new_species"])
count <- names(count[count >= 20])
unique_exotic <- unique_exotic[unique_exotic %in% count]

polygon_df <- species_df <-  list()
record_threshold <- 30

region <- read_sf("regions2.shp")
region_df <- read.csv("regionID.csv")

for (i in 1:length(unique_exotic)) {
  message(i,"/",length(unique_exotic))
  sp <- unique_exotic[[i]]
  start <- Sys.time()
  GLONAF_plants <- subset(GLONAF, new_species == sp & status == "naturalized")
  subset_plants <- subset(plants_db_cleaned,new_species==sp)
  
  coords <- data.frame(decimalLongitude=subset_plants$decimalLongitude,decimalLatitude=subset_plants$decimalLatitude,species=subset_plants$new_species)
  points <- st_as_sf(coords, coords = c("decimalLongitude","decimalLatitude"))
  st_crs(points) <- st_crs(region)
  
  region_p <- st_transform(region, 2163) 
  point_p <- st_transform(points, 2163) 
  
  intersect_df <- st_intersects(point_p,region_p,sparse=F)
  intersection <- apply(intersect_df,1,which)
  polygon_seq <- unlist(lapply(intersection,function(x) ifelse(length(x) == 0,NA,x)))
  
  polygon_associated <- region_p[polygon_seq, ] 
  above_threshold <- which(table(polygon_associated$OBJIDsic) >= record_threshold)
  
  at_region_id <- region_df[region_df$OBJIDsic %in% names(above_threshold),"region_id"]
  at_naturalized_polygon <- unique(GLONAF_plants$region_id) %in%  at_region_id
  at_naturalized_polygon_id <- unique(GLONAF_plants$region_id)[at_naturalized_polygon]
  num_n_p <- length(which(at_naturalized_polygon))
  
  species_df[[i]] <- data.frame(sp,num_n_p,length(unique(polygon_seq)),nrow(coords))
  if (length(at_naturalized_polygon_id) > 0){
    polygon_df[[i]] <- data.frame(sp=sp,region_id=at_naturalized_polygon_id)
  } else {
    polygon_df[[i]] <- data.frame(sp=sp,region_id=NA)
    
  }
  
  end <- Sys.time()
  print(polygon_df[[i]])
  print(end-start)
}

sp_df <- do.call(rbind,species_df)
polygon_df <- do.call(rbind,polygon_df)
write.csv(polygon_df,"polygon_df_plants.csv")
