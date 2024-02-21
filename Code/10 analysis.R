library(phyr)
library(terra)
library(ape)
LC_pref_df <- read.csv("plants_LC_pref_df.csv")
LC_pref_df$id <- paste0(LC_pref_df$species,"_",LC_pref_df$OBJIDsic)
region_df <- read.csv("regionID.csv")
region <- vect("regions2.shp")
GLONAF <- read.csv("GLONAF_new.csv")
GLONAF$OBJIDsic <- region_df[match(GLONAF$region_id,region_df$region_id),"OBJIDsic"]
GLONAF_naturalized <- subset(GLONAF,status=="naturalized")
GLONAF_naturalized$id <- paste0(GLONAF_naturalized$new_species,"_",GLONAF_naturalized$OBJIDsic)
LC_pref_df$naturalized <- GLONAF_naturalized[match(LC_pref_df$id,GLONAF_naturalized$id),"status"]

clim_niche <- read.csv("clim_niche.csv")[,-1]
clim_niche <- na.omit(clim_niche)
econ <- read.csv("EcoUsers3.csv")
econ$econ_freq <- rowSums(econ[,4:16])
econ$species.name <- gsub(" x "," ",econ$species.name)
library(WorldFlora)
WFO.remember("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/classification.csv")
match_final <- WFO.match(unique(econ$species.name)[unique(econ$species.name) %in% GLONAF$new_species],WFO.data=WFO.data,no.dates=T,counter=1,Fuzzy=0.1)
match_final_one <- WFO.one(match_final)
econ$new_species <- econ$species.name
mismatch <- match_final_one[match_final_one$spec.name != match_final_one$scientificName,]

######################
LC_pref_df$econ <- econ[match(LC_pref_df$species,gsub(" x "," ",econ$species.name)),"econ_freq"]
LC_pref_df$Environmental <- econ[match(LC_pref_df$species,gsub(" x "," ",econ$species.name)),"Environmental"]
LC_pref_df$Human_Food <- econ[match(LC_pref_df$species,gsub(" x "," ",econ$species.name)),"Human.food"]

LC_pref_df <- cbind(LC_pref_df,clim_niche[match(LC_pref_df$species,clim_niche$sp),2:18])

LC_pref_df$native.BIO4 <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"BIO4"]
LC_pref_df$native.BIO4.sd <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"BIO4.sd"]

LC_pref_df$native.BIO1 <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"BIO1"]
LC_pref_df$native.BIO1.sd <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"BIO1.sd"]

LC_pref_df$native.BIO4 <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"BIO4"]
LC_pref_df$native.BIO4.sd <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"BIO4.sd"]

LC_pref_df$native.SMmean <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"SMmean"]
LC_pref_df$native.SMsd <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"SMsd"]

LC_pref_df$native.CMImean <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"CMI_mean"]
LC_pref_df$native.CMIsd <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"CMI_sd"]

LC_pref_df$native.BIO12 <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"BIO12"]
LC_pref_df$native.BIO12.sd <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"BIO12.sd"]

LC_pref_df$native.BIO15 <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"BIO15"]
LC_pref_df$native.BIO15.sd <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"BIO15.sd"]

LC_pref_df$native.PETmean <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"PET_mean"]
LC_pref_df$native.PETsd <- clim_niche[match(LC_pref_df$species,clim_niche$sp),"PET_sd"]

#######################
options(glmmTMB.cores=6)
library(glmmTMB)
#mean_all <- mean(final_plants_db$consensus_full_class_9,na.rm=T)
#sd_all <- mean(final_plants_db$consensus_full_class_9,na.rm=T)
#LC_pref_df$SES_urb <- (LC_pref_df$Urb_obs_mean-mean_all)/(sd_all)

#LC_pref_df$SES_urb <- (LC_pref_df$Urb_obs_mean-LC_pref_df$Urb_Effort_obs_mean_all)/(LC_pref_df$Urb_Effort_obs_sd_all)
#LC_pref_df$SES_agr <- (LC_pref_df$Agr_obs_mean-LC_pref_df$Agr_Effort_obs_mean_all)/(LC_pref_df$Agr_Effort_obs_sd_all)

LC_pref_df$SES_urb <- (LC_pref_df$Urb_obs_mean-LC_pref_df$Urb_Effort_obs_mean)/(LC_pref_df$Urb_Effort_obs_sd)
LC_pref_df$SES_agr <- (LC_pref_df$Agr_obs_mean-LC_pref_df$Agr_Effort_obs_mean)/(LC_pref_df$Agr_Effort_obs_sd)

#LC_pref_df[LC_pref_df$Urb_obs_mean-LC_pref_df$Urb_Effort_obs_mean == 0,"SES_urb"] <- 0 #force (0-0)/0 = 0

alien_df <- LC_pref_df[LC_pref_df$species %in% GLONAF$new_species,]
alien_df[is.na(alien_df$econ),c("econ","Environmental","Human_Food")] <- 0
alien_df <- na.omit(alien_df)

alien_df <- cbind(alien_df,as.data.frame(region)[match(alien_df$OBJIDsic,region$OBJIDsic),c("LON","LAT")])
alien_df$OBJIDsic <- as.factor(as.character(alien_df$OBJIDsic))

#alien_df <- alien_df[alien_df$Urb >= 100,]

#native.LV2 <- alien_df %>% group_by(OBJIDsic) %>% summarise(native.PCA1.water.LV2=mean(native.PCA1.water),
                                                       #native.PCA1.temp.LV2 = mean(native.PCA1.temp))
#region.LV2 <- alien_df %>% group_by(species) %>% summarise(PCA1water.LV2=mean(PCA1water.1),
                                                            #PCA1temp.LV2 = mean(PCA1temp.1))

#alien_df <- cbind(alien_df,native.LV2[match(alien_df$OBJIDsic,native.LV2$OBJIDsic),2:3])
#alien_df <- cbind(alien_df,region.LV2[match(alien_df$species,region.LV2$species),2:3])
#alien_df <- cbind(alien_df,alien_df[,c("PCA1temp.1","PCA1water.1")]-alien_df[,c("PCA1temp.LV2","PCA1water.LV2")])
#alien_df <- cbind(alien_df,alien_df[,c("native.PCA1.temp","native.PCA1.water")]-alien_df[,c("native.PCA1.temp.LV2","native.PCA1.water.LV2")])
#colnames(alien_df)[44:47] <- c("LV1.PCAtemp.1","LV1.PCAwater.1","LV1.native.PCA1.temp","LV1.native.PCA1.water")

#m_urb <- glmmTMB(SES_urb~PCA1temp.LV2*native.PCA1.temp+PCA1water.LV2*native.PCA1.water+
                  #LV1.PCAtemp.1*native.PCA1.temp+
                   #LV1.PCAwater.1*native.PCA1.water+econ+
                   #(LV1.PCAtemp.1+LV1.PCAwater.1||species)+
                   #(native.PCA1.temp+native.PCA1.water||OBJIDsic),data=alien_df,REML=T)
#summary(m_urb)
alien_df$native.PCA1temp.range <- alien_df$PCA1temp.max-alien_df$PCA1temp.min
alien_df$native.PCA1water.range <- alien_df$PCA1water.max-alien_df$PCA1water.min
#cor(alien_df$native.PCA1water.max,alien_df$native.PCA1water.min)
m_urb <- glmmTMB(SES_urb~scale(PCA1temp)*scale(PCA1temp.min)+scale(PCA1water)*scale(PCA1water.min)+
                   econ+(scale(PCA1temp)+scale(PCA1water)||species)+
                   (scale(PCA1water.min)+scale(PCA1temp.min)||OBJIDsic),data=alien_df,REML=T)
summary(m_urb)
car::Anova(m_urb)
performance::r2(m_urb)

m_urb <- glmmTMB(SES_urb~scale(PCA1temp)*scale(all.PCA1temp.min)+scale(PCA1water)*scale(all.PCA1water.min)+
                   econ+(scale(PCA1temp)+scale(PCA1water)||species)+
                   (scale(all.PCA1water.min)+scale(all.PCA1temp.min)||OBJIDsic),data=alien_df,REML=T)
summary(m_urb)

performance::r2(m_urb)

alien_df$pred1 <- predict(m_urb)
cor(alien_df$pred1,alien_df$SES_urb)
m_urb1 <- glmmTMB(SES_urb~1+(1|OBJIDsic),data=alien_df[order(alien_df$species),],REML=T)
summary(m_urb1)
performance::r2(m_urb1)

alien_df$pred2 <- predict(m_urb1)

####################
library(DHARMa)
res <- simulateResiduals(m_urb)
res2 <- recalculateResiduals(res,alien_df[,"OBJIDsic"])
testSpatialAutocorrelation(res2,x =  aggregate(alien_df$LON, list(alien_df$OBJIDsic), mean)$x, 
                           y = aggregate(alien_df$LAT, list(alien_df$OBJIDsic), mean)$x)

library(WorldFlora)
library(tidyverse)
WFO.remember("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/classification.csv")
match_final <- WFO.match(unique(alien_df$species)[unique(alien_df$species) %in% GLONAF$new_species],WFO.data=WFO.data,no.dates=T,counter=1,Fuzzy=0)
match_final_one <- WFO.one(match_final)
phylo_df <- data.frame(species=match_final_one$spec.name,genus=word(match_final_one$spec.name,1,1," "),family=match_final_one$family)

library(rtrees)
library(ape)
tree <- get_tree(sp_list = unique(alien_df$species),
                 taxon = "plant",
                 scenario = "at_basal_node",
                 show_grafted = TRUE)

res <- simulateResiduals(m_urb)
res2 <- recalculateResiduals(res,alien_df$species)
phyD <- cophenetic.phylo(tree)
phyD <- phyD[order(rownames(phyD)),order(colnames(phyD))]
testSpatialAutocorrelation(res2,distMat=phyD)

cor(alien_df[,20:37])
############################################


predict_data <- ggemmeans(m,terms=("scaled.diff.BIO12[-6:4,by=1]]"))
p <- ggplot(data=predict_data)+
  geom_line(aes(x=x,y=predicted,color=group))
plot(p)

max(scale(alien_df$diff.BIO12))
#########test phylo AC

library(ggeffects)
library(ggplot2)
library(interactions)

johnson_neyman(m,pred=Arid_mean,modx=native.Arid_mean)
predict.data <- ggemmeans(m,terms=c("BIO1[-10:10,by=2]","native.BIO1[-10:10],by=2]"))

p <- ggplot(data=predict.data)+
  geom_line(aes(x=x,y=predicted,color=group))
plot(p)
m <- glmmTMB(Urb~scale(BIO1)*scale(native.BIO1)+
               scale(BIO4)*scale(native.BIO4)+
               scale(BIO12)*scale(native.BIO12)+
               scale(BIO15)*scale(native.BIO15)+
               (scale(BIO1)+scale(BIO4)+scale(BIO12)+scale(BIO15)||species)+
               (scale(native.BIO1)+scale(native.BIO4)+scale(native.BIO12)+scale(native.BIO15)||OBJIDsic),data=alien_df)
summary(m)
car::Anova(m)
######################################################################
vcv_sp <- vcv(phy_tree)

LC_pref_df$species <- gsub(" ","_",LC_pref_df$species)
z <- pglmm(Urb~ BIO1+BIO12+BIO15+log(count)+
           (1|OBJIDsic)+(1|species__)+(BIO1|species__)+(BIO12|species__)+(BIO15|species__), 
           data = LC_pref_df,
           cov_ranef = list(species = phy_tree), REML = TRUE, verbose = TRUE,bayes=T)
end <- Sys.time()

###
library(pgirmess)
SAC_test <- correlog(coord_polygon[,2:3],resid)
library(DHARMa)
resid = simulateResiduals(m)
resid = recalculateResiduals(resid,group=as.factor(LC_pref_df$OBJIDsic),sum)
coord_polygon <- coord_polygon[order(coord_polygon$OBJIDsic),]
testSpatialAutocorrelation(resid, x =  coord_polygon$LON, y = coord_polygon$LAT)

resid = simulateResiduals(m)
resid = recalculateResiduals(resid,group=LC_pref_df$species)
phy_d <- cophenetic.phylo(phy_tree)
phy_d <- phy_d[order(rownames(phy_d)),order(colnames(phy_d))]
testSpatialAutocorrelation(resid, dist=phy_d)

library(PCAtest)
result<-PCAtest(phy_d, 100, 100, 0.05, varcorr=FALSE, counter=T, plot=F)
#####

#####
library(vegan)
phy_d <- cophenetic.phylo(phy_tree)
phy_pcoa <- pcoa(phy_d)
Phy_Eigen <- phy_pcoa$vectors[,phy_pcoa$values$Relative_eig >= 1/2572]

Phy_Eigen_result <- NULL
for (i in 1:ncol(Phy_Eigen)) {
  message(i)
  sp_phy_eigen <- Phy_Eigen[match(LC_pref_df$species,rownames(Phy_Eigen)),i]
  m_resid <- glmmTMB(resid~scale(sp_phy_eigen)+(1|OBJIDsic)+(1|species),data=LC_pref_df)
  p_value <- summary(m_resid)$coefficients$cond[2,4]
  Phy_Eigen_result <- rbind(Phy_Eigen_result,data.frame(i=i,p_value=p_value))
}

#####
library(adespatial)
library(SoDA)
library(sf)
library(spaMM)
library(lmerTest)
region <- read_sf("regions2.shp")
#coord_polygon <- as.data.frame(region[region$OBJIDsic %in% LC_pref_df$OBJIDsic,])
coord_polygon <- as.data.frame(region[match(LC_pref_df$OBJIDsic,region$OBJIDsic),])
coord_polygon <- coord_polygon[,1:3]
coord_polygon[,2:3]<-geoXY(longitude=coord_polygon[,2],latitude=coord_polygon[,3])
clim_niche <- read.csv("clim_niche.csv")

LC_pref_df <- cbind(LC_pref_df,as.data.frame(region[match(LC_pref_df$OBJIDsic,region$OBJIDsic),])[c("LON","LAT")])

sp_list <- names(table(LC_pref_df$species))[table(LC_pref_df$species) >= 30]
LC_pref_df_subset <- LC_pref_df[LC_pref_df$species %in% sp_list,]

coef_result <- p_value_result <- list()

for (sp in unique(LC_pref_df_subset$species)) {
  message(sp)
  LC_pref_df_subset_sp <- subset(LC_pref_df_subset,species==sp)
  LC_pref_df_subset_sp$OBJIDsic <- as.factor(LC_pref_df_subset_sp$OBJIDsic)
  LC_pref_df_subset_sp <- cbind(LC_pref_df_subset_sp,clim_niche[match(sp,clim_niche$species),c("min_Bio1","min_Arid")])
  LC_pref_df_subset_sp$Bio1_diff <- LC_pref_df_subset_sp$BIO1 - LC_pref_df_subset_sp$min_Bio1 
  LC_pref_df_subset_sp$Arid_diff <- LC_pref_df_subset_sp$Arid_mean - LC_pref_df_subset_sp$min_Arid 
  
  m_Urb <- fitme(Urb~BIO1+Arid_mean+log(count)+Matern(1|LON+LAT),data=LC_pref_df_subset_sp)
  m_Urb <- fitme(Urb~Bio1_diff+Arid_diff+log(count)+Matern(1|LON+LAT),data=LC_pref_df_subset_sp)
  
  #m_Urb1 <- fitme(Urb~scale(log(count))+Matern(1|LON+LAT),data=LC_pref_df_subset_sp)
  #p_Urb <- LRT(m_Urb,m_Urb1)
  p_Urb <- anova(m_Urb)$`Pr(>F)`[1:2]
  
  m_Agr <- fitme(Agr~BIO1+Arid_mean+log(count)+Matern(1|LON+LAT),data=LC_pref_df_subset_sp)
  m_Agr <- fitme(Agr~Bio1_diff+Arid_diff+log(count)+Matern(1|LON+LAT),data=LC_pref_df_subset_sp)
  #m_Agr1 <- fitme(Agr~scale(log(count))+Matern(1|LON+LAT),data=LC_pref_df_subset_sp)
  #p_Agr <- LRT(m_Agr,m_Agr1)
  p_Agr <- anova(m_Agr)$`Pr(>F)`[1:2]
  
  m_Urb_result <- as.data.frame(summary(m_Urb)$beta_table)[2:3,]
  m_Agr_result <- as.data.frame(summary(m_Agr)$beta_table)[2:3,]
  
  coef_result[[sp]] <- data.frame(species=sp,
                                  count=sum(LC_pref_df_subset_sp$count),
                                  mean_BIO1 = mean(LC_pref_df_subset_sp$BIO1),
                                  mean_Arid = mean(LC_pref_df_subset_sp$Arid_mean),
                                  range_BIO1 = max(LC_pref_df_subset_sp$BIO1) - min(LC_pref_df_subset_sp$BIO1),
                                  range_Arid = max(LC_pref_df_subset_sp$Arid_mean) - min(LC_pref_df_subset_sp$Arid_mean),
                                  mean_BIO1_diff = mean(LC_pref_df_subset_sp$Bio1_diff),
                                  mean_Arid_diff = mean(LC_pref_df_subset_sp$Arid_diff),
                                  range_BIO1_diff = max(LC_pref_df_subset_sp$Bio1_diff) - min(LC_pref_df_subset_sp$Bio1_diff),
                                  range_Arid_diff = max(LC_pref_df_subset_sp$Arid_diff) - min(LC_pref_df_subset_sp$Arid_diff),
                                  t(m_Urb_result$Estimate),
                                  t(m_Urb_result$`t-value`),
                                  p_Urb=t(p_Urb),
                                  "Urb",
                                  t(m_Agr_result$Estimate),
                                  t(m_Agr_result$`t-value`),
                                  p_Agr=t(p_Agr),"Agr")
}

  coef_result_df <- do.call(rbind,coef_result)
  colnames(coef_result_df) <- c("species","count",
                                "mean_BIO1","mean_Arid",
                                "range_BIO1","range_Arid",
                                "mean_BIO1_diff","mean_Arid_diff",
                                "range_BIO1_diff","range_Arid_diff",
                                "BIO1_slope_urb","Arid_mean_slope_urb",
                                "BIO1_urb_t","Arid_urb_t",
                                "p_value_BIO1_urb","p_value_Arid_urb",
                                "Urb",
                                "BIO1_slope_agr","Arid_mean_slope_agr",
                                "BIO1_agr_t","Arid_agr_t",
                                "p_value_BIO1_agr","p_value_Arid_agr",
                                "agr")
  
  check_urb_temp <- coef_result_df[which(coef_result_df$p_value_BIO1_urb < 0.05),]
  check_urb_arid <- coef_result_df[which(coef_result_df$p_value_Arid_urb < 0.05),]
  length(which(check_urb_temp$BIO1_slope_urb < 0))  
  length(which(check_urb_temp$BIO1_slope_urb > 0))  
  length(which(check_urb_arid$Arid_mean_slope_urb < 0))  
  length(which(check_urb_arid$Arid_mean_slope_urb > 0))  
  
  check_agr_temp <- coef_result_df[which(coef_result_df$p_value_BIO1_agr < 0.05),]
  check_agr_arid <- coef_result_df[which(coef_result_df$p_value_Arid_agr < 0.05),]
  length(which(check_agr_temp$BIO1_slope_agr < 0))  
  length(which(check_agr_temp$BIO1_slope_agr > 0))  
  length(which(check_agr_arid$Arid_mean_slope_agr < 0))  
  length(which(check_agr_arid$Arid_mean_slope_agr > 0))  
  