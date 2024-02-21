final_plants_db <- read.csv("plants_final.csv")
library(terra)
library(taxize)
library(tidyverse)
library(radiant.data)
region <- vect("regions2.shp")
region_df <- read.csv("regionID.csv")
GLONAF <- read.csv("GLONAF_new.csv")

library(WorldFlora)
treelist <- read.csv("global_tree_search_trees_1_7.csv")
#WFO.remember("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/classification.csv")
#not_found <- treelist[!treelist$TaxonName %in% WFO.data$scientificName,]
#match_final <- WFO.match(not_found$TaxonName,WFO.data=WFO.data,no.dates=T,counter=1,Fuzzy=0.1)
#match_final_one <- WFO.one(match_final)
#treelist$new_species <- treelist$TaxonName
#treelist[treelist$TaxonName %in% match_final_one$spec.name, "new_species" ] <- match_final_one$scientificName
final_plants_db$tree <- ifelse(final_plants_db$new_species %in% treelist$TaxonName,"Tree","Non-tree")
#colnames(final_plants_db)[13:14] <- c("Agr","Urb")
#write.csv(treelist,"global_tree_search_trees_1_7.csv")
############
threshold = 20 #at least 20 records for the plant species
LC_pref <- list()
i = 1

#all_effort <- final_plants_db %>% group_by(tree) %>% summarise(mean.Agr = mean(Agr,na.rm=T),
                                                               #mean.Urb = mean(Urb,na.rm=T),
                                                               #sd.Agr = sd(Agr,na.rm=T), 
                                                               #sd.Urb=sd(Urb,na.rm=T))

all_effort <- final_plants_db %>% summarise(mean.Agr = mean(Agr,na.rm=T),
                                                               mean.Urb = mean(Urb,na.rm=T),
                                                               sd.Agr = sd(Agr,na.rm=T), 
                                                               sd.Urb=sd(Urb,na.rm=T))

for (n in na.omit(unique(final_plants_db$polygon_number))) { #polygon_number is OBJIDsic
  message(n)
  subset_plants_df <- subset(final_plants_db,polygon_number == n)

  if (nrow(subset_plants_df) == 0) {
    LC_pref[[i]] <- data.frame(species=NA,
                               count=NA,
                               Agr_obs_mean = NA, Urb_obs_mean =NA ,Agr_obs_sd = NA, Urb_obs_sd =NA, Tree = NA,
                               Agr_Effort_obs_mean = NA,Urb_Effort_obs_mean = NA,Agr_Effort_obs_sd =NA,Urb_Effort_obs_sd =NA,
                               Agr_Effort_obs_mean_all = NA,Urb_Effort_obs_mean_all = NA,Agr_Effort_obs_sd_all =NA,Urb_Effort_obs_sd_all =NA,
                               OBJIDsic = n,
                               Reason = "No naturalized sp in GBIF")
    i = i+1
    
    next
  } else {
  
  record_n <- table(subset_plants_df$new_species)
  
  analyzed_sp <- names(record_n)[record_n >= threshold]
  region_id <- region_df[region_df$OBJIDsic == n,"region_id"]
  regional_naturalized_sp_list <- GLONAF[GLONAF$status=="naturalized" & GLONAF$region_id == region_id,]
  analyzed_sp <- analyzed_sp[analyzed_sp %in% regional_naturalized_sp_list$new_species] #only include species that are naturalized in the region!
  
  if (length(analyzed_sp) == 0) {
    LC_pref[[i]] <- data.frame(species=NA,
                               count=NA,
                               Agr_obs_mean = NA, Urb_obs_mean =NA ,Agr_obs_sd = NA, Urb_obs_sd =NA, Tree=NA,
                               Agr_Effort_obs_mean = NA,Urb_Effort_obs_mean = NA,Agr_Effort_obs_sd =NA,Urb_Effort_obs_sd =NA,
                               Agr_Effort_obs_mean_all = NA,Urb_Effort_obs_mean_all = NA,Agr_Effort_obs_sd_all =NA,Urb_Effort_obs_sd_all =NA,
                               OBJIDsic = n,
                               Reason = "No species / no species with sufficient records")
    i = i+1
    
    next
  } else {
    
    subset_analyzed_sp <- subset_plants_df[subset_plants_df$new_species %in% analyzed_sp,] #records for species > threshold
    
    ##################no longer relevant weighted mean part
    #Agr_fun <- approxfun(density(subset_plants_df$Agr)) #create a distirubtion of records along agr gradient
    #Urb_fun <- approxfun(density(subset_plants_df$Urb))
    
    #subset_analyzed_sp$Urb_d <- Urb_fun(subset_analyzed_sp$Urb) #calculate density of each record - 
    #subset_analyzed_sp$Agr_d <- Agr_fun(subset_analyzed_sp$Agr)
    
    #Urb_naturalized_sp <- subset_analyzed_sp %>% group_by(species) %>% summarise(Urb = weighted.mean(Urb,1/Urb_d)) #1/d = weight 
    #Agr_naturalized_sp <- subset_analyzed_sp %>% group_by(species) %>% summarise(Agr = weighted.mean(Agr,1/Agr_d))
    #Urb_naturalized_sp_sd <- subset_analyzed_sp %>% group_by(species) %>% summarise(Urb = weighted.sd(Urb,1/Urb_d))
    #Agr_naturalized_sp_sd <- subset_analyzed_sp %>% group_by(species) %>% summarise(Agr = weighted.sd(Agr,1/Agr_d))
    #####################
    
    Species_obs <- subset_analyzed_sp %>% dplyr::group_by(new_species) %>% dplyr::summarise(count=n(),mean.Agr=mean(Agr),mean.Urb=mean(Urb),sd.Agr=sd(Agr),sd.Urb=sd(Urb))
    Effort_obs<- subset_plants_df %>% dplyr::group_by(tree) %>% dplyr::summarise(mean.Agr=mean(Agr),mean.Urb=mean(Urb),sd.Agr=sd(Agr),sd.Urb=sd(Urb))
    
    LC_pref_data <- Species_obs
    
    #pos <- ifelse(LC_pref_data$new_species %in% treelist$TaxonName,2,1)
    
    LC_pref_data <- cbind(LC_pref_data,Effort_obs[pos,],all_effort[,-1])
    
    LC_pref_data <- cbind(LC_pref_data,n,"Fine")
  
    colnames(LC_pref_data) <- c("species","count","Agr_obs_mean","Urb_obs_mean","Agr_obs_sd","Urb_obs_sd",
                              "Tree","Agr_Effort_obs_mean","Urb_Effort_obs_mean","Agr_Effort_obs_sd","Urb_Effort_obs_sd",
                              "Agr_Effort_obs_mean_all","Urb_Effort_obs_mean_all","Agr_Effort_obs_sd_all","Urb_Effort_obs_sd_all",
                              "OBJIDsic","Reason")
  
  LC_pref[[i]] <- LC_pref_data
  i = i+1
    }
  }
}

LC_pref_df <- do.call(rbind,LC_pref)
LC_pref_df <- na.omit(LC_pref_df)
write.csv(LC_pref_df,"plants_LC_pref_df.csv")

####Polygon-level invaded region climate
library(terra)
LC_pref_df <- read.csv("plants_LC_pref_df.csv")
#################
Tmax <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/CHELSA/CHELSA_bio10_05.tif")
Tmin <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/CHELSA/CHELSA_bio10_06.tif")
Tsd <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/CHELSA/CHELSA_bio10_04.tif")

#Tmax <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/Terraclim/TerraClimate19812010_tmax.nc")
#Tmin <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/Terraclim/TerraClimate19812010_tmin.nc")
#P <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/Terraclim/TerraClimate19812010_ppt.nc")
#Soil <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/Terraclim/TerraClimate19812010_soil.nc")
#PET <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/Terraclim/TerraClimate19812010_pet.nc")
#Tmonth <- mean(Tmin,Tmax)
#Tsd <- stdev(Tmonth)
#Tmin <- min(Tmin)
#Tmax <- max(Tmax)

#Tmin <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/WorldClim/wc2.1_2.5m_bio_6.tif") #worldclim
#Tmax <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/WorldClim/wc2.1_2.5m_bio_5.tif") #worldclim
#Tsd <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/WorldClim/wc2.1_2.5m_bio_4.tif") #worldclim
temp <- c(Tmin,Tmax,Tsd)
temp <- aggregate(temp,5)
#swc_dir <- list.files("C:/Users/pakno/OneDrive - University of Toronto/Raster/swc_fr")
#swc_dir <- paste0("C:/Users/pakno/OneDrive - University of Toronto/Raster/swc_fr/",paste0(swc_dir),"/w001001.adf")
#Soil <- rast(swc_dir)
#Soil <- terra::aggregate(Soil,5,na.rm=T)

SM_max <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/SM_agg/max_sm_5km_2020.tif")
SM_min <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/SM_agg/min_sm_5km_2020.tif")
SM_sd <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/SM_agg/sd_sm_5km_2020.tif")

#SM_min <- min(Soil)
#SM_max <- max(Soil)
#SM_sd <- stdev(Soil)
#Pmin <- min(P)
#Pmax <- max(P)
#Psd <- stdev(P)
#Pmean <- mean(P)
#Psd <- stdev(P)
#PET_mean <- mean(PET)
#PET_sd <- stdev(PET)
#SMmean <- mean(Soil)
#SMsd <- stdev(Soil)
#CMI <- P-PET
#CMI_min <- min(CMI)
#CMI_max <- max(CMI)
#CMI_sd <- stdev(CMI)

#clim <- c(Tmin,Tmax,Tsd,SM_min,SM_max,SM_sd)
#clim
#cor(as.data.frame(clim))

water <- c(SM_min,SM_max,SM_sd)
water

temp <- project(temp,water)

pca_water <- prcomp(na.omit(as.data.frame(water)),scale=T)
summary(pca_water)
pca_water
BiodiversityR::PCAsignificance(vegan::rda(scale(na.omit(as.data.frame(water)))~1))
water_df <- na.omit(as.data.frame(water,xy=T))[,1:2]
PCA1_water <- predict(pca_water)[,1]
PCA1_water <- cbind(water_df,PCA1_water)
PCA1_water <- rast(PCA1_water)
writeRaster(PCA1_water,"PCA1_water_0810.tif",overwrite=T)

pca_temp <- prcomp(na.omit(as.data.frame(temp)),scale=T)
summary(pca_temp)
pca_temp
BiodiversityR::PCAsignificance(vegan::rda(scale(na.omit(as.data.frame(temp)))~1))
temp_df <- na.omit(as.data.frame(temp,xy=T))[,1:2]
PCA1_temp <- predict(pca_temp)[,1]
PCA1_temp <- cbind(temp_df,PCA1_temp)
PCA1_temp <- rast(PCA1_temp)
writeRaster(PCA1_temp,"PCA1_temp_0810.tif",overwrite=T)

############################################
PCA1_temp <- rast("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/PCA1_temp_0810.tif")
PCA1_water <- rast("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/PCA1_water_0810.tif")

PCA1_temp <- project(PCA1_temp,PCA1_water)
clim <- c(PCA1_temp,PCA1_water)
polygon_clim <- terra::extract(clim,region,fun=mean,na.rm=T)
polygon_clim$OBJIDsic <- region$OBJIDsic
colnames(polygon_clim) <- c("ID", "PCA1temp","PCA1water","OBJIDsic")
#colnames(polygon_clim) <- c("ID", "BIO1","BIO4","BIO12","BIO15","PET_mean","PET_sd","SMmean","SMsd","CMI_mean","CMI_sd","OBJIDsic")

Agr <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/EarthEnv/consensus_full_class_7.tif")
Urb <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/EarthEnv/consensus_full_class_9.tif")
LC <- c(Agr,Urb)

polygon_max_LC <- terra::extract(LC,region,fun=max,na.rm=T)
colnames(polygon_max_LC) <- c("ID", "Agr_max","Urb_max")
polygon_min_LC <- terra::extract(LC,region,fun=min,na.rm=T)
colnames(polygon_min_LC) <- c("ID", "Agr_min","Urb_min")

polygon_clim <- cbind(polygon_clim,polygon_max_LC[,-1],polygon_min_LC[,-1])
for_combine <- polygon_clim[match(LC_pref_df$OBJIDsic,polygon_clim$OBJIDsic),]
LC_pref_df <- cbind(LC_pref_df,for_combine)
write.csv(LC_pref_df,"plants_LC_pref_df.csv")

###don't run this
LC_pref_df <- read.csv("Magnoliopsidae_LC_pref_df.csv")

library(sf)
library(taxize)
library(tidyverse)
library(radiant.data)
library(GIFT)
load("Code/GIFT.RData")

points <- st_as_sf(data.frame(lon=final_plants_db$decimalLongitude,lat=final_plants_db$decimalLatitude),coords=c("lon","lat"))
st_crs(points) <- st_crs(gift_shape)
sf_use_s2(FALSE)
gift_shape <- st_make_valid(gift_shape)
segment <- nrow(points)/200000

for (i in 1:segment) {
  message(i/segment)
  start <- Sys.time()
  first_data <- (1+(i-1)*200000)
  last_data <- ifelse(i*200000 > nrow(points),nrow(points),i*200000)
  pts_gift_poly <- st_intersects(points[first_data:last_data,],gift_shape)
  end <- Sys.time()
  print(end-start)
  write.csv(pts_gift_poly,paste0("gift/",i,".csv"))
}

LC_pref_df_unique_sp <- unique(LC_pref_df$species)
distr_threshold = 15
clim_niche <- list()
i = 1
GIFT_list <- GIFT_species()
GIFT_parsed <- gbif_parse(GIFT_list$work_species)

for (sp in LC_pref_df_unique_sp) { #polygon_number is OBJIDsic
  message(sp)
  subset_plants_df <- subset(final_plants_db,species == sp)
  
  GIFT_sp <- GIFT_parsed[match(sp,GIFT_parsed$canonicalname),"canonicalname"]
  
  genus <- unlist(lapply(strsplit(GIFT_sp," "),function(x) x[[1]]))
  epithet <- unlist(lapply(strsplit(GIFT_sp," "),function(x) x[[2]]))
  
  lookup <- GIFT_species_lookup(genus = genus[[1]], epithet = epithet[[1]])
  
  lookup_distribution <- GIFT_species_distribution(
    genus = genus[[1]],epithet = epithet[[1]], aggregation = TRUE)
  statuses <- lookup_distribution %>% 
    mutate(native = ifelse(native == 1, "native", "non-native"),
    naturalized = ifelse(naturalized == 1, "naturalized","non-naturalized"),
    endemic_list = ifelse(endemic_list == 1, "endemic_list",
    "non-endemic_list")) %>%
    dplyr::select(entity_ID, native, naturalized, endemic_list)
  native_polygon <- subset(statuses,native == "native")
  shape <- gift_shape[which(gift_shape$entity_ID %in% 
                              unique(native_polygon$entity_ID)), ]
  
  sp_points <- st_as_sf(data.frame(lon=subset_plants_df$decimalLongitude,lat=subset_plants_df$decimalLatitude),coords=c("lon","lat"))
  st_crs(sp_points) <- st_crs(shape)

  
  if (nrow(subset_naturalized_df) == 0) {
    LC_pref[[i]] <- data.frame(species=NA,
                               count=NA,
                               Agr=NA,
                               Urb=NA,
                               Agr_sd=NA,
                               Urb_sd = NA,
                               Tree = NA,
                               OBJIDsic = n,
                               Reason = "No naturalized sp in GBIF")
    i = i+1
    
    next
  } else {
    record_n <- table(subset_naturalized_df$species)
    
    analyzed_sp <- names(record_n)[record_n >= threshold]
    if (length(analyzed_sp) == 0) {
      LC_pref[[i]] <- data.frame(species=NA,
                                 count=NA,
                                 Agr=NA,
                                 Urb=NA,
                                 Agr_sd=NA,
                                 Urb_sd = NA,
                                 Tree= NA,
                                 OBJIDsic = n,
                                 Reason = "No species with sufficient records")
      i = i+1
      
      next
    } else {
      
      Agr_fun <- approxfun(density(subset_plants_df$Agr))
      Urb_fun <- approxfun(density(subset_plants_df$Urb))
      
      subset_analyzed_sp <- subset_naturalized_df[subset_naturalized_df$species %in% analyzed_sp,]
      subset_analyzed_sp$Urb_d <- Urb_fun(subset_analyzed_sp$Urb)
      subset_analyzed_sp$Agr_d <- Agr_fun(subset_analyzed_sp$Agr)
      
      Urb_naturalized_sp <- subset_analyzed_sp %>% group_by(species) %>% summarise(Urb = weighted.mean(Urb,1/Urb_d))
      Agr_naturalized_sp <- subset_analyzed_sp %>% group_by(species) %>% summarise(Agr = weighted.mean(Agr,1/Agr_d))
      Urb_naturalized_sp_sd <- subset_analyzed_sp %>% group_by(species) %>% summarise(Urb = weighted.sd(Urb,1/Urb_d))
      Agr_naturalized_sp_sd <- subset_analyzed_sp %>% group_by(species) %>% summarise(Agr = weighted.sd(Agr,1/Agr_d))
      
      
      LC_pref_data <- data.frame(table(subset_analyzed_sp$species),
                                 Agr_naturalized_sp$Agr,
                                 Urb_naturalized_sp$Urb,
                                 Agr_naturalized_sp_sd$Agr,
                                 Urb_naturalized_sp_sd$Urb,
                                 n,
                                 "Fine")
      
      colnames(LC_pref_data) <- c("species","count","Agr","Urb","Agr_sd","Urb_sd","OBJIDsic","Reason")
      LC_pref[[i]] <- LC_pref_data
      i = i+1
    }
  }
}

LC_pref_df <- do.call(rbind,LC_pref)
LC_pref_df <- na.omit(LC_pref_df)
