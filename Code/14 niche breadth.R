clim_niche <- read.csv("clim_niche.csv")[,-1]
clim_niche <- na.omit(clim_niche)

econ <- read.csv("EcoUsers3.csv")
econ$econ_freq <- rowSums(econ[,4:16])
econ$species.name <- gsub(" x "," ",econ$species.name)
GLONAF <- read.csv("GLONAF_new.csv")

library(WorldFlora)
WFO.remember("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/classification.csv")
match_final <- WFO.match(unique(econ$species.name)[unique(econ$species.name) %in% GLONAF$new_species],WFO.data=WFO.data,no.dates=T,counter=1,Fuzzy=0.1)
match_final_one <- WFO.one(match_final)
econ$new_species <- econ$species.name
mismatch <- match_final_one[match_final_one$spec.name != match_final_one$scientificName,]


#################
library(tidyverse)
trait <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/joswig_trait.csv")
WFO.remember("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/classification.csv")

trait_avg <- trait %>% group_by(Species) %>% dplyr::summarise(LeArea_avg = mean(LeArea),
                                                              SSD_avg = mean(SSD),
                                                              SLA_avg = mean(SLA),
                                                              LeC_avg = mean(LeC),
                                                              LeN_avg= mean(LeN),
                                                              LeP_avg = mean(LeP),
                                                              H_avg = mean(PlantHeight),
                                                              SeedMass_avg = mean(SeedMass),
                                                              SeLen_avg =mean(SeLen),
                                                              LeNArea_avg = mean(LeNArea),
                                                              LeNP_avg = mean(LeNP),
                                                              Led15N_avg = mean(Led15N),
                                                              SenbU_avg=mean(SenbU),
                                                              LeFMass_avg=mean(LeFMass),
                                                              ConduitDens_avg = mean(ConduitDens),
                                                              DispULen_avg = mean(DispULen),
                                                              VesLen_avg = mean(VesLen))

match_final <- WFO.match(unique(trait_avg$Species)[unique(trait_avg$Species) %in% GLONAF$new_species],WFO.data=WFO.data,no.dates=T,counter=1,Fuzzy=0.1)
match_final_one <- WFO.one(match_final)
trait_avg$new_species <- trait_avg$Species
mismatch <- match_final_one[match_final_one$spec.name != match_final_one$scientificName,]

clim_niche$econ_freq <- econ[match(clim_niche$sp,econ$species.name),"econ_freq"]
clim_niche <- cbind(clim_niche,trait_avg[match(clim_niche$sp,trait_avg$Species),-1])

clim_niche_cleaned <- na.omit(clim_niche)

WFO.remember("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/classification.csv")
match_final <- WFO.match(clim_niche_cleaned$sp,WFO.data=WFO.data,no.dates=T,counter=1,Fuzzy=0)
match_final_one <- WFO.one(match_final)
phylo_df <- data.frame(species=match_final_one$spec.name,genus=word(match_final_one$spec.name,1,1," "),family=match_final_one$family)

library(rtrees)
test_tree = get_tree(sp_list = phylo_df,
                     taxon = "plant",
                     scenario = "at_basal_node",
                     show_grafted = TRUE)

test_tree$tip.label <- gsub("_"," ",test_tree$tip.label)
test_tree$tip.label <- gsub("\\*","",test_tree$tip.label)

library(nlme)
pca_trait <- vegan::rda(log(trait_avg[,2:18])~1,scale=T)
summary(pca_trait)
sig_PCA <- BiodiversityR::PCAsignificance(pca_trait,axes=17)
clim_niche_cleaned <- cbind(clim_niche_cleaned,vegan::scores(pca_trait,choices=1:2)$sites[match(clim_niche_cleaned$new_species,trait_avg$new_species),])

##########################
library(nlme)
library(ape)
clim_niche_cleaned$exotic_breadth_temp <- clim_niche_cleaned$exotic.PCA1temp.max-clim_niche_cleaned$exotic.PCA1temp.min
clim_niche_cleaned$exotic_breadth_water <- clim_niche_cleaned$exotic.PCA1water.max-clim_niche_cleaned$exotic.PCA1water.min
clim_niche_cleaned$breadth_temp <- clim_niche_cleaned$all.PCA1temp.max-clim_niche_cleaned$all.PCA1temp.min
clim_niche_cleaned$breadth_water <- clim_niche_cleaned$all.PCA1water.max-clim_niche_cleaned$all.PCA1water.min
clim_niche_cleaned$native_breadth_temp <- clim_niche_cleaned$PCA1temp.max-clim_niche_cleaned$PCA1temp.min
clim_niche_cleaned$native_breadth_water <- clim_niche_cleaned$PCA1water.max-clim_niche_cleaned$PCA1water.min

m1temp <- gls(exotic_breadth_temp~PC1+PC2+econ_freq+native_breadth_temp+native_breadth_water,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))
summary(m1temp)
m1water <- gls(exotic_breadth_water~PC1+PC2+econ_freq+native_breadth_temp+native_breadth_water,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))
summary(m1water)

mPC1 <- gls(PC1~native_breadth_temp+native_breadth_water,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))
mPC1_null <- gls(PC1~1,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))

mPC2 <- gls(PC2~native_breadth_temp+native_breadth_water,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))
mPC2_null <- gls(PC2~1,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))

#mPC3 <- gls(PC3~native_breadth_temp+native_breadth_water,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))
#mPC3_null <- gls(PC3~1,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))

#mPC4 <- gls(PC4~native_breadth_temp+native_breadth_water,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))
#mPC4_null <- gls(PC4~1,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))

mEcon <- gls(econ_freq~native_breadth_temp+native_breadth_water+PC1+PC2,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))
mEcon_null <- gls(econ_freq~1,data=clim_niche_cleaned,correlation = corPagel(1,phy =test_tree,form=~sp))

library(rr2)
R2(mPC1,mPC1_null)
R2(mPC2,mPC2_null)
#R2(mPC3,mPC3_null)
#R2(mPC4,mPC4_null)
R2(mEcon,mEcon_null)
library(piecewiseSEM)
SEM <- psem(m1temp,m1water,mPC1,mPC2,mEcon,
            PC1 %~~% PC2,
            #PC1 %~~% PC3,
            #PC1 %~~% PC4,
            #PC2 %~~% PC3,
            #PC2 %~~% PC4,
            #PC3 %~~% PC4,
            exotic_breadth_temp %~~% exotic_breadth_water)
summary(SEM)

library(phyr)
start<-Sys.time()
cor_phylo(variates=~exotic_breadth_temp+exotic_breadth_water,
          data=clim_niche_cleaned,
          method= "nelder-mead-r",
          phy=test_tree,
          species=clim_niche_cleaned$sp,
          verbose=T,
          constrain_d=T)
end <- Sys.time()
end-start
cor(clim_niche_cleaned$exotic_breadth_temp,clim_niche_cleaned$breadth_temp)

library(ape)
v_e_temp <- clim_niche_cleaned$exotic_breadth_temp
v_e_water <- clim_niche_cleaned$exotic_breadth_water
names(v_e_temp) <- names(v_e_water) <- clim_niche_cleaned$sp
pic(v_,test_tree)
