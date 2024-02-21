library(kewr)
library(purrr)
library(tidyverse)
library(rnaturalearth)
library(terra)
library(radiant.data)
library(plyr)
library(WorldFlora)

#final_plants_db <- read.csv("plants_final.csv")
LC_pref_df <- read.csv("plants_LC_pref_df.csv")
tdwg3 <- vect("tdwg3/level3.shp")
#tdwg3_assoc <- read.csv("plants_pts_in_tdwg3.csv")

region <- vect("regions2.shp")
region_df <- read.csv("regionID.csv")
GLONAF <- read.csv("GLONAF_new.csv")
WFO.remember("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/classification.csv")

#tdwg3_assoc <- terra::extract(tdwg3,final_plants_db[,c("decimalLongitude","decimalLatitude")])
#write.csv(tdwg3_assoc,"tdwg3_pts.csv")

#tdwg3_assoc <- read.csv("tdwg3_pts.csv")
#final_plants_db$tdwg3 <- NULL
#tdwg3_assoc <- tdwg3_assoc[!duplicated(tdwg3_assoc$row.id),]
#final_plants_db[tdwg3_assoc$row.id,"tdwg3"] <- tdwg3[tdwg3_assoc$col.id,"LEVEL3_COD"]

#remove(tdwg3_assoc)
gc()
#region$tdwg3 <- region_df[match(region$OBJIDsic,region_df$OBJIDsic),"tdwg3"]

#################
#Tmax <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/Terraclim/TerraClimate19812010_tmax.nc")
#Tmin <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/Terraclim/TerraClimate19812010_tmin.nc")
#P <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/Terraclim/TerraClimate19812010_ppt.nc")
#Soil <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/Terraclim/TerraClimate19812010_soil.nc")
#PET <- rast("C:/Users/pakno/OneDrive - University of Toronto/Raster/Terraclim/TerraClimate19812010_pet.nc")

#Tmonth <- (Tmax+Tmin)/2 #no Tmean in Terraclimate

#Tsd <- stdev(Tmonth)
#Tmean <- mean(Tmonth)
#Pmean <- mean(P)
#Psd <- stdev(P)
#PET_mean <- mean(PET)
#PET_sd <- stdev(PET)
#SMmean <- mean(Soil)
#SMsd <- stdev(Soil)
#CMI_mean <- Pmean - PET_mean
#CMI_sd <- stdev(P-PET)

#clim <- c(Tsd,Tmean,Pmean,Psd,PET_mean,PET_sd,SMmean,SMsd,CMI_mean,CMI_sd)
####################
alien_list <- unique(LC_pref_df$species[LC_pref_df$species %in% GLONAF$new_species]) #only get climatic niche for alien species included for analyses later
#clim_var <- c("BIO4","BIO1","BIO12","BIO15","PET_mean","PET_sd","SMmean","SMsd","CMI_mean","CMI_sd")
#clim_niche <- as.data.frame(matrix(nrow=1,ncol=13))
#colnames(clim_niche) <- c("sp",clim_var,paste0(clim_var,".sd"))

clim_niche_tdwg3 <- terra::extract(clim,tdwg3,fun=mean,na.rm=T)
clim_niche_tdwg4 <- terra::extract(clim,region,fun=mean,na.rm=T)
#colnames(clim_niche_tdwg3) <- c("ID","BIO4","BIO1","BIO12","BIO15","PET_mean","PET_sd","SMmean","SMsd","CMI_mean","CMI_sd")
clim_niche_tdwg4 <- cbind(clim_niche_tdwg4,region_df[match(region$OBJIDsic,region_df$OBJIDsic),"region_id"])
colnames(clim_niche_tdwg3) <-  c("ID","PCA1temp","PCA1water")
colnames(clim_niche_tdwg4) <- c("ID","PCA1temp","PCA1water","region_id")

clim_niche_tdwg3<-cbind(as.data.frame(tdwg3),clim_niche_tdwg3)
clim_niche_tdwg3$area <- expanse(tdwg3)

clim_niche <- NULL
i=1
for (sp in alien_list) {
  message(sp,";",i,"/",length(alien_list))
  i=i+1
  start <- Sys.time()
  
  powo <- search_powo(sp)
  powo_info <- tidy(powo)
  
  if (nrow(powo_info) == 0) {
    clim_niche <- plyr::rbind.fill(clim_niche,data.frame(sp=sp))
    message("no match found")
    next()
  }
  
  accepted_match <- subset(powo_info,(name == paste0(strsplit(sp," ")[[1]][1], " × ",strsplit(sp," ")[[1]][2]) | name == sp) & accepted == T) #x for Erysimum cheiri??
  
  if (nrow(accepted_match) == 0) {
    clim_niche <- plyr::rbind.fill(clim_niche,data.frame(sp=sp))
    message("no match found")
    next()
  }
  
  synonym <- 1
  
  if (!("synonymOf" %in% colnames(powo_info))) {
    synonym <- 0
    message("no synonym found")
    powo_info <- subset(powo_info,(name == paste0(strsplit(sp," ")[[1]][1], " × ",strsplit(sp," ")[[1]][2]) | name == sp) & accepted == T)
    powo_ID <- sub(".*names:","",powo_info$fqId)
  }
  
  if ("synonymOf" %in% colnames(powo_info)) {
    if (nrow(accepted_match) != 0 & is.null(nrow(accepted_match$synonymOf))) {
      synonym <- 0
      message("found synonym, but the name is correct")
      powo_ID <- sub(".*names:","",accepted_match$fqId)
    }
  }
  
  if ("synonymOf" %in% colnames(powo_info)) {
    if (!is.null(powo_info[powo_info$name == sp,"synonymOf"]) & synonym == 1) {
      message("found synonym")
      synonym <- do.call(rbind.fill,powo_info$synonymOf)
      powo <- lapply(synonym$name,search_powo)
      powo_info <- do.call(plyr::rbind.fill,lapply(powo,tidy))
      powo_info <- powo_info[powo_info$name %in% synonym$name & powo_info$accepted == T,]
      powo_ID <- sub(".*names:","",powo_info$fqId)
    } 
  }
  
  tdwg3_code <- NULL
  for (j in 1:length(powo_ID)){
    record <- lookup_powo(powo_ID[[j]], distribution=TRUE)
    tidied <- tidy(record)
    if (!"natives" %in% colnames(as.data.frame(tidied$distribution))){
      message("no native distribution")
      next()
    }  
    
    if ("distribution" %in% colnames(tidied)) {
      tidy_distribution <- tidied %>%
        select(fqId, distribution) %>%
        unnest(cols=distribution) %>%
        unnest(cols=natives)
      
      tdwg3_code <- c(tdwg3_code,tidy_distribution$tdwgCode)
    }
  }
  
  if (is.null(tdwg3_code)) {
    clim_niche <- plyr::rbind.fill(clim_niche,data.frame(sp=sp))
    next()
  }
  
  native_polygon <- clim_niche_tdwg3[clim_niche_tdwg3$LEVEL3_COD %in% unique(tdwg3_code),]
  
  exotic_polygon_tdwg4 <- subset(GLONAF,new_species == sp & status == "naturalized")
  exotic_polygon_clim <- clim_niche_tdwg4[clim_niche_tdwg4$region_id %in% exotic_polygon_tdwg4$region_id,]
  range_polygon_exotic <- t(data.frame(range(exotic_polygon_clim$PCA1temp,na.rm=T),range(exotic_polygon_clim$PCA1water,na.rm=T)))
  
  tdwg3_exotic_polygon_code <- region_df[match(exotic_polygon_tdwg4$region_id, region_df$region_id),"tdwg3"]
  exotic_polygon <- clim_niche_tdwg3[clim_niche_tdwg3$LEVEL3_COD %in% tdwg3_exotic_polygon_code,]
  
  all_polygon <- rbind(native_polygon,exotic_polygon)
  range_polygon_native <- t(data.frame(range(native_polygon$PCA1temp,na.rm=T),range(native_polygon$PCA1water,na.rm=T)))
  range_polygon_all <- t(data.frame(range(all_polygon$PCA1temp,na.rm=T),range(all_polygon$PCA1water,na.rm=T)))
  n_naturalized <- nrow(exotic_polygon)
  w <- native_polygon$area
  clim_niche_sp_native <- apply(native_polygon[,6:7],2,function(x) mean(x,w=w,na.rm=T))
  clim_niche_sp_exotic <- apply(exotic_polygon[,6:7],2,function(x) mean(x,w=w,na.rm=T))
  
  clim_niche_sp <- data.frame(sp=sp,t(clim_niche_sp_native),t(clim_niche_sp_exotic),t(range_polygon_native[1,]),t(range_polygon_native[2,]),t(range_polygon_all[1,]),t(range_polygon_all[2,]),
                              t(range_polygon_exotic[1,]),t(range_polygon_exotic[2,]),n_naturalized)
  
  clim_niche <- rbind.fill(clim_niche,clim_niche_sp)
  end <- Sys.time()
  print(end-start)
}

colnames(clim_niche)<-c("sp","PCA1temp.mean.native","PCA1water.mean.native","PCA1temp.mean.exotic","PCA1water.mean.exotic",
                        "PCA1temp.min","PCA1temp.max","PCA1water.min","PCA1water.max",
                        "all.PCA1temp.min","all.PCA1temp.max","all.PCA1water.min","all.PCA1water.max",
                        "exotic.PCA1temp.min","exotic.PCA1temp.max","exotic.PCA1water.min","exotic.PCA1water.max","exotic_n_region")
colnames(clim_niche)
write.csv(clim_niche,"clim_niche.csv")

########################################################### no need ...point-based mean
########################################################### point-based

for (sp in alien_list) {
  message(sp)
  start <- Sys.time()
  
  powo <- search_powo(sp)
  powo_info <- tidy(powo)
  
  if (nrow(powo_info) == 0) {
    clim_niche <- plyr::rbind.fill(clim_niche,data.frame(sp=sp))
    message("no match found")
    next()
  }
  
  accepted_match <- subset(powo_info,(name == paste0(strsplit(sp," ")[[1]][1], " × ",strsplit(sp," ")[[1]][2]) | name == sp) & accepted == T) #x for Erysimum cheiri??
  
  if (nrow(accepted_match) == 0) {
    clim_niche <- plyr::rbind.fill(clim_niche,data.frame(sp=sp))
    message("no match found")
    next()
  }
  
  synonym <- 1
  
  if (!("synonymOf" %in% colnames(powo_info))) {
    synonym <- 0
    message("no synonym found")
    powo_info <- subset(powo_info,(name == paste0(strsplit(sp," ")[[1]][1], " × ",strsplit(sp," ")[[1]][2]) | name == sp) & accepted == T)
    powo_ID <- sub(".*names:","",powo_info$fqId)
  }
  
  if ("synonymOf" %in% colnames(powo_info)) {
    if (nrow(accepted_match) != 0 & is.null(nrow(accepted_match$synonymOf))) {
      synonym <- 0
      message("found synonym, but the name is correct")
      powo_ID <- sub(".*names:","",accepted_match$fqId)
    }
  }
  
  if ("synonymOf" %in% colnames(powo_info)) {
    if (!is.null(powo_info[powo_info$name == sp,"synonymOf"]) & synonym == 1) {
      message("found synonym")
      synonym <- do.call(rbind.fill,powo_info$synonymOf)
      powo <- lapply(synonym$name,search_powo)
      powo_info <- do.call(plyr::rbind.fill,lapply(powo,tidy))
      powo_info <- powo_info[powo_info$name %in% synonym$name & powo_info$accepted == T,]
      powo_ID <- sub(".*names:","",powo_info$fqId)
    } 
  }
  
  tdwg3_code <- NULL
  for (j in 1:length(powo_ID)){
    record <- lookup_powo(powo_ID[[j]], distribution=TRUE)
    tidied <- tidy(record)
    if (!"natives" %in% colnames(as.data.frame(tidied$distribution))){
      message("no native distribution")
      next()
    }  
    
    if ("distribution" %in% colnames(tidied)) {
      tidy_distribution <- tidied %>%
        select(fqId, distribution) %>%
        unnest(cols=distribution) %>%
        unnest(cols=natives)
      
      tdwg3_code <- c(tdwg3_code,tidy_distribution$tdwgCode)
    }
  }
  
  if (is.null(tdwg3_code)) {
    clim_niche <- plyr::rbind.fill(clim_niche,data.frame(sp=sp))
    next()
  }
  
  native_polygon <- tdwg3[tdwg3$LEVEL3_COD %in% unique(tdwg3_code),]
  
  # Area of distribution with floristic status
  #orig_p <- ggplot(world) +
  #geom_sf(color = "gray70") +
  #geom_sf(data = native_polygon, color = "black") +
  #scale_fill_brewer("Status", palette = "Set2") +
  #labs(title = paste("Distribution map of ",
  #sp),
  #subtitle = "Unprojected (GCS: WGS84)") +
  #theme_void()
  
  #plot(orig_p)
  #ggsave(paste0("POWO_map/",sp,".tiff"),width=8,height=8)
  
  subset_all_plants <- final_plants_db[tdwg3_assoc$id.y[tdwg3_assoc$LEVEL3_COD %in% native_polygon$LEVEL3_COD],]
  subset_sp_occ <- subset(subset_all_plants,species==sp)
  
  clim_niche_sp <- NULL
  for (i in 1:length(clim_var)) {
    fun <- approxfun(density(na.omit(subset_all_plants[,clim_var[[i]]])))
    clim_density <- fun(subset_sp_occ[,clim_var[[i]]])
    clim_niche_sp[[i]] <- weighted.mean(subset_sp_occ[,clim_var[[i]]],1/clim_density,na.rm=T)
    clim_niche_sp[[length(clim_var)+i]] <- weighted.sd(subset_sp_occ[,clim_var[[i]]],1/clim_density,na.rm=T)
  }
  
  clim_niche_sp <- data.frame(sp=sp,do.call(cbind,clim_niche_sp))
  names(clim_niche_sp) <- c("sp",clim_var,paste0(clim_var,".","sd"))
  
  clim_niche <- rbind.fill(clim_niche,clim_niche_sp)
  end <- Sys.time()
  print(end-start)
}

########################################################### point-based

for (sp in alien_list[1:1117]) {
  message(sp)
  start <- Sys.time()
  
  powo <- search_powo(sp)
  powo_info <- tidy(powo)
  
  if (nrow(powo_info) == 0) {
    clim_niche <- plyr::rbind.fill(clim_niche,data.frame(sp=sp))
    message("no match found")
    next()
  }
  
  accepted_match <- subset(powo_info,(name == paste0(strsplit(sp," ")[[1]][1], " × ",strsplit(sp," ")[[1]][2]) | name == sp) & accepted == T) #x for Erysimum cheiri??
  
  if (nrow(accepted_match) == 0) {
    clim_niche <- plyr::rbind.fill(clim_niche,data.frame(sp=sp))
    message("no match found")
    next()
  }
  
  synonym <- 1
  
  if (!("synonymOf" %in% colnames(powo_info))) {
    synonym <- 0
    message("no synonym found")
    powo_info <- subset(powo_info,(name == paste0(strsplit(sp," ")[[1]][1], " × ",strsplit(sp," ")[[1]][2]) | name == sp) & accepted == T)
    powo_ID <- sub(".*names:","",powo_info$fqId)
  }
  
  if ("synonymOf" %in% colnames(powo_info)) {
    if (nrow(accepted_match) != 0 & is.null(nrow(accepted_match$synonymOf))) {
      synonym <- 0
      message("found synonym, but the name is correct")
      powo_ID <- sub(".*names:","",accepted_match$fqId)
    }
  }
    
  if ("synonymOf" %in% colnames(powo_info)) {
    if (!is.null(powo_info[powo_info$name == sp,"synonymOf"]) & synonym == 1) {
    message("found synonym")
    synonym <- do.call(rbind.fill,powo_info$synonymOf)
    powo <- lapply(synonym$name,search_powo)
    powo_info <- do.call(plyr::rbind.fill,lapply(powo,tidy))
    powo_info <- powo_info[powo_info$name %in% synonym$name & powo_info$accepted == T,]
    powo_ID <- sub(".*names:","",powo_info$fqId)
    } 
  }
  
  tdwg3_code <- NULL
  for (j in 1:length(powo_ID)){
    record <- lookup_powo(powo_ID[[j]], distribution=TRUE)
    tidied <- tidy(record)
    if (!"natives" %in% colnames(as.data.frame(tidied$distribution))){
      message("no native distribution")
      next()
    }  
    
    if ("distribution" %in% colnames(tidied)) {
    tidy_distribution <- tidied %>%
      select(fqId, distribution) %>%
      unnest(cols=distribution) %>%
      unnest(cols=natives)
    
    tdwg3_code <- c(tdwg3_code,tidy_distribution$tdwgCode)
    }
  }
  
  if (is.null(tdwg3_code)) {
    clim_niche <- plyr::rbind.fill(clim_niche,data.frame(sp=sp))
    next()
  }

  native_polygon <- tdwg3[tdwg3$LEVEL3_COD %in% unique(tdwg3_code),]
  
  # Area of distribution with floristic status
  #orig_p <- ggplot(world) +
    #geom_sf(color = "gray70") +
    #geom_sf(data = native_polygon, color = "black") +
    #scale_fill_brewer("Status", palette = "Set2") +
    #labs(title = paste("Distribution map of ",
                       #sp),
         #subtitle = "Unprojected (GCS: WGS84)") +
    #theme_void()
  
  #plot(orig_p)
  #ggsave(paste0("POWO_map/",sp,".tiff"),width=8,height=8)
  
  subset_all_plants <- final_plants_db[tdwg3_assoc$id.y[tdwg3_assoc$LEVEL3_COD %in% native_polygon$LEVEL3_COD],]
  subset_sp_occ <- subset(subset_all_plants,species==sp)
  
  clim_niche_sp <- NULL
  for (i in 1:length(clim_var)) {
    fun <- approxfun(density(na.omit(subset_all_plants[,clim_var[[i]]])))
    clim_density <- fun(subset_sp_occ[,clim_var[[i]]])
    clim_niche_sp[[i]] <- weighted.mean(subset_sp_occ[,clim_var[[i]]],1/clim_density,na.rm=T)
    clim_niche_sp[[length(clim_var)+i]] <- weighted.sd(subset_sp_occ[,clim_var[[i]]],1/clim_density,na.rm=T)
  }

  clim_niche_sp <- data.frame(sp=sp,do.call(cbind,clim_niche_sp))
  names(clim_niche_sp) <- c("sp",clim_var,paste0(clim_var,".","sd"))
  
  clim_niche <- rbind.fill(clim_niche,clim_niche_sp)
  end <- Sys.time()
  print(end-start)
}

write.csv(clim_niche,"clim_niche.csv")
##use GIFT, currently not in use

devtools::install_github("https://github.com/BioGeoMacro/GIFT")
library("GIFT")
library(tidyverse)
library(taxize)
library(kewr)
load("Code/GIFT.RData")

region <- read_sf("regions2.shp")
region_df <- read.csv("regionID.csv")
GLONAF <- read.csv("GLONAF.csv")

LC_pref_df_unique_sp <- unique(LC_pref_df$species)
distr_threshold = 15
clim_niche <- list()
i = 1
GIFT_sp_list <- GIFT_species()
GIFT_study_list <- GIFT_lists()
GIFT_parsed <- gbif_parse(GIFT_sp_list$work_species)
GIFT_parsed <- subset(GIFT_parsed,type=="SCIENTIFIC")

final_plants_db <- read.csv("Magnoliopsidae_final.csv")
LC_pref_df <- read.csv("Magnoliopsidae_LC_pref_df.csv")
glonaf <- GIFT_overlap(resource = "glonaf")

clim_niche <- NULL
for (sp in unique(coef_result_df$species)) {
  message(sp)
  start <- Sys.time()
  potential_name <- as.data.frame(name_lookup(sp)$data)
  potential_name <- unlist(lapply(strsplit(potential_name$canonicalName," "),function(x) paste(x[1],x[2])))
  potential_name <- unique(potential_name)
  pos <- which(GIFT_parsed$canonicalname %in% potential_name)
  parsed_sp <- GIFT_parsed$canonicalname[pos]
  
  genus <- unlist(lapply(strsplit(parsed_sp," "),function(x) x[[1]]))
  epithet <- unlist(lapply(strsplit(parsed_sp," "),function(x) x[[2]]))
  
  lookup <- GIFT_species_lookup(genus = genus, epithet = epithet)
  lookup_distribution <- GIFT_species_distribution(
    genus = genus, epithet = epithet, area_th_mainland = 0,aggregation = FALSE, by_ref_ID= T)
  lookup_distribution$native[lookup_distribution$quest_native == 1 & lookup_distribution$native == 1] <- NA
  statuses <- lookup_distribution%>%
    mutate(native = ifelse(native == 1, "native", "non-native"),
           naturalized = ifelse(naturalized == 1, "naturalized",
                                "non-naturalized"),
           endemic_list = ifelse(endemic_list == 1, "endemic_list",
                                 "non-endemic_list")) %>%
    dplyr::select(entity_ID, native, naturalized, endemic_list)
  
  anemone_shape <- gift_shape[which(gift_shape$entity_ID %in% 
                                      unique(lookup_distribution$entity_ID)), ]
  anemone_map <- dplyr::left_join(anemone_shape, statuses,
                                  by = "entity_ID")
  
  # Area of distribution with floristic status
  orig_p <- ggplot(world) +
    geom_sf(color = "gray70") +
    geom_sf(data = anemone_map, color = "black", aes(fill = as.factor(native))) +
    scale_fill_brewer("Status", palette = "Set2") +
    labs(title = paste("Distribution map of ",
                                  sp),
         subtitle = "Unprojected (GCS: WGS84)") +
    theme_void()
  
  plot(orig_p)
  ggsave(paste0("map/",sp,"_orig.tiff"),width=8,height=8)
  
  native_poly <- subset(statuses,native=="native")
  gift_subset <- glonaf[glonaf$entity_ID %in% native_poly$entity_ID,]
  glonaf_polygon_native <- subset(gift_subset,overlap21 >= 0.5)
  glonaf_polygon_native$tdwg3 <- region_df[match(glonaf_polygon_native$glonaf_ID,region_df$OBJIDsic),"tdwg3"]
    
  naturalized_polygon <- subset(GLONAF,standardized_name == sp)
  OBJIDsic_naturalized <- region_df[match(naturalized_polygon$region_id,region_df$region_id),"OBJIDsic"]
  TDWG3_naturalized <- region_df[match(naturalized_polygon$region_id,region_df$region_id),"tdwg3"]

  glonaf_polygon_native <- glonaf_polygon_native[!glonaf_polygon_native$glonaf_ID %in% OBJIDsic_naturalized,]
  glonaf_polygon_native <- glonaf_polygon_native[!glonaf_polygon_native$tdwg3 %in% TDWG3_naturalized,]
  
  native_polygon_glonaf <- region[region$OBJIDsic %in% glonaf_polygon_native$glonaf_ID,1]
  #OBJIDsic_native <- region_df[match(glonaf_polygon_native$glonaf_ID,region_df$region_id),"OBJIDsic"]
  
  revised_p <- ggplot(world) +
    geom_sf(color = "gray70") +
    geom_sf(data = native_polygon_glonaf, color = "black",fill="green") +
    scale_fill_brewer("Status", palette = "Set2") +
    labs(title = paste("Distribution map of ",
                       sp),
         subtitle = "Unprojected (GCS: WGS84)") +
    theme_void()
  
  plot(revised_p)
  
  ggsave(paste0("map/",sp,"_revised.tiff"),width=8,height=8)
  
  subset_all_plants <- final_plants_db[final_plants_db$polygon_number %in% native_polygon_glonaf$OBJIDsic,]
  subset_sp_occ <- subset(subset_all_plants,species==sp)
  
  Bio1_fun <- approxfun(density(na.omit(subset_all_plants$BIO1)))
  Arid_fun <- approxfun(density(na.omit(subset_all_plants$Arid_mean)))
  
  subset_sp_occ$BIO1_d <- Bio1_fun(subset_sp_occ$BIO1)
  subset_sp_occ$Arid_mean_d <- Arid_fun(subset_sp_occ$Arid_mean)
  
  sp_Bio1_niche <- weighted.mean(subset_sp_occ$BIO1,1/subset_sp_occ$BIO1_d,na.rm=T)
  sp_Arid_niche <- weighted.mean(subset_sp_occ$Arid_mean,1/subset_sp_occ$Arid_mean_d,na.rm=T)
  
  clim_niche <- rbind(clim_niche,data.frame(species=sp,sp_Bio1_niche=sp_Bio1_niche,sp_Arid_niche=sp_Arid_niche))
  end <- Sys.time()
  print(end-start)
}
