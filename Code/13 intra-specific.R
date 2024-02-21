library(phyr)
library(ape)
#LC_pref_df <- read.csv("plants_LC_pref_df.csv")
GLONAF <- read.csv("GLONAF_new.csv")
clim_niche <- read.csv("clim_niche.csv")[,-1]
clim_niche <- na.omit(clim_niche)
econ <- read.csv("EcoUsers3.csv")
econ$econ_freq <- rowSums(econ[,4:16])
econ$species.name <- gsub(" x "," ",econ$species.name)

#################
sp_analysis_list <- names(table(alien_df$species))[which(table(alien_df$species) >= 20)]

result_df<-list()
library(spaMM)
for (i in 1:length(sp_analysis_list)){
  message(i,"/",length(sp_analysis_list))
  subset_alien_df <- subset(alien_df,species == sp_analysis_list[[i]])
  subset_alien_df$Urb_range <- subset_alien_df$Urb_max-subset_alien_df$Urb_min
  r_cor <- cor(subset_alien_df$PCA1temp,subset_alien_df$PCA1water)
  if(length(unique(subset_alien_df$Urb_range)) == 1){
  sp_m <- fitme(SES_urb~PCA1temp+PCA1water+log(count)+Matern(1|LON+LAT),data=subset_alien_df)
  range_df <- as.data.frame(t(data.frame(c(NA,NA),range(subset_alien_df$PCA1temp),range(subset_alien_df$PCA1water))))
  } else {
    sp_m <- fitme(SES_urb~PCA1temp+PCA1water+log(count)+Urb_range+Matern(1|LON+LAT),data=subset_alien_df)
    Urb_range <- range(subset_alien_df$Urb_range)
    range_df <- as.data.frame(t(data.frame(c(NA,NA),range(subset_alien_df$PCA1temp),range(subset_alien_df$PCA1water),Urb_range)))
  }
  result <- summary(sp_m,detail=T)
  rownames(range_df) <- NULL
  result_df[[i]] <- data.frame(sp=sp_analysis_list[[i]],
                                        Tree=unique(subset_alien_df$Tree),
                                        analyzed_gradient_PCA1temp = range_df[2,2]-range_df[2,1],
                                        analyzed_mean_PCA1temp = mean(subset_alien_df$PCA1temp),
                                        true_mean_PCA1temp =unique(subset_alien_df$PCA1temp.mean.exotic),
                                        #analyzed_mean_PCA1temp = mean(range_df[2,2],range_df[2,1]),
                                        #true_mean_PCA1temp = mean(unique(subset_alien_df$exotic.PCA1temp.min),unique(subset_alien_df$exotic.PCA1temp.max)),
                                        analyzed_gradient_PCA1water = range_df[3,2]-range_df[3,1],
                                        analyzed_mean_PCA1water = mean(subset_alien_df$PCA1water),
                                        true_mean_PCA1water = unique(subset_alien_df$PCA1water.mean.exotic),
                                        #analyzed_mean_PCA1water = mean(range_df[3,2],range_df[3,1]),
                                        #true_mean_PCA1water = mean(unique(subset_alien_df$exotic.PCA1water.min),unique(subset_alien_df$exotic.PCA1water.max)),
                                        analyzed_gradient_urb = ifelse(is.na(range_df[4,2]-range_df[4,1]),0,range_df[4,2]-range_df[4,1]),
                                        var=rownames(result$beta_table[-4,]),
                                        result$beta_table[-4,],
                                        n=nrow(subset_alien_df),
                                        corr = r_cor,
                                        econ = unique(subset_alien_df$econ),
                                        Environmental = unique(subset_alien_df$Environmental),
                                        Human_Food = unique(subset_alien_df$Human_Food))
}

intra_sp_result <- NULL
intra_sp_result <- do.call(rbind,result_df)

temp_sig <- subset(intra_sp_result,var=="PCA1temp" & `p.value` < 0.05)
water_sig <- subset(intra_sp_result,var=="PCA1water" & `p.value` < 0.05)

#################
library(tidyverse)
trait <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/joswig_trait.csv")
library(WorldFlora)
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

intra_sp_result <- cbind(intra_sp_result,clim_niche[match(intra_sp_result$sp,clim_niche$sp),2:18])
#intra_sp_result <- cbind(intra_sp_result,selected_trait[match(intra_sp_result$sp,selected_trait$SpeciesName),-1])
intra_sp_result <- cbind(intra_sp_result,trait_avg[match(intra_sp_result$sp,trait_avg$Species),-1])

#pierce <- read.csv("pierce_2017.csv")
#intra_sp_result <- cbind(intra_sp_result,pierce[match(gsub("_"," ",intra_sp_result$sp),pierce$species),c("C....","S....","R....")])

library(GIFT)
trait_meta <- GIFT_traits_meta()
trait <- GIFT_traits(trait_IDs= c("1.2.1","2.1.1"))
colnames(trait)[4:5] <- c("Form","Cycle")
intra_sp_result <- cbind(intra_sp_result,trait[match(intra_sp_result$sp,trait$work_species),4:5])
intra_sp_result$group <- paste0(intra_sp_result$Form,"_",intra_sp_result$Cycle)
library(GIFT)

WFO.remember("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/classification.csv")
match_final <- WFO.match(unique(intra_sp_result$sp)[gsub("_"," ",unique(intra_sp_result$sp)) %in% GLONAF$new_species],WFO.data=WFO.data,no.dates=T,counter=1,Fuzzy=0)
match_final_one <- WFO.one(match_final)
phylo_df <- data.frame(species=match_final_one$spec.name,genus=word(match_final_one$spec.name,1,1," "),family=match_final_one$family)

library(rtrees)
test_tree = get_tree(sp_list = phylo_df,
                     taxon = "plant",
                     scenario = "random_below_basal",
                     show_grafted = TRUE)

test_tree$tip.label <- gsub("\\*","",test_tree$tip.label)

############# metafor framework
library(metafor)
library(vegan)
library(nlme)

intra_sp_result <- na.omit(intra_sp_result)
intra_sp_result$sp <- gsub(" ","_",intra_sp_result$sp)
intra_sp_result$true_gradient_PCA1water <- intra_sp_result$exotic.PCA1water.max-intra_sp_result$exotic.PCA1water.min
intra_sp_result$true_gradient_PCA1temp <- intra_sp_result$exotic.PCA1temp.max-intra_sp_result$exotic.PCA1temp.min

test <- na.omit(subset(intra_sp_result,var=="PCA1temp"))
colnames(test)[[1]] <- "Species"
test <- subset(test,econ!=0 & corr >= -0.7 & corr <= 0.7)

#pca_niche <- rda(test[,c("all.PCA1water.min","all.PCA1water.max","all.PCA1temp.max")]~1,scale=T)
#test$pca1_niche <- vegan::scores(pca_niche,choice=1)$sites
test$breadth_temp <- test$all.PCA1temp.max-test$all.PCA1temp.min
test$breadth_water <- test$all.PCA1water.max-test$all.PCA1water.min
test$diff_gradient_temp <- test$true_gradient_PCA1temp - test$analyzed_gradient_PCA1temp
test$diff_gradient_water <- test$true_gradient_PCA1water - test$analyzed_gradient_PCA1water
test$diff_mean_temp <- test$true_mean_PCA1temp-test$analyzed_mean_PCA1temp
test$diff_mean_water<- test$true_mean_PCA1water-test$analyzed_mean_PCA1water

pca_trait <- vegan::rda(log(trait_avg[,2:18])~1,scale=T)
summary(pca_trait)
sig_PCA <- BiodiversityR::PCAsignificance(pca_trait,axes=17)

#extract scores for axis > 95% of bs model , plus observed var > 5%
#test <- cbind(test,vegan::scores(pca_trait,choices=which(as.data.frame(sig_PCA)[6,] >= 0.95 &  as.data.frame(sig_PCA)[2,] >= 10))$sites)
test <- cbind(test,vegan::scores(pca_trait,choices=1:2)$sites[match(test$Species,gsub(" ","_",trait_avg$new_species)),])
test$n <- log(test$n)

mr <- rma(Estimate~diff_mean_temp+analyzed_mean_PCA1temp+scale(analyzed_gradient_PCA1temp)+diff_gradient_temp+econ+n,sei=test$`Cond..SE`,data=test)
summary(mr)
vif(mr)

mr <- rma(Estimate~diff_mean_temp+analyzed_mean_PCA1temp+diff_gradient_temp+econ+n+PC1+PC2,sei=test$`Cond..SE`,data=test)
summary(mr)
vif(mr,data.frame())

test$analyzed_gradient_PCA1temp_non_linear <- test$analyzed_gradient_PCA1temp^2
mr <- rma(Estimate~diff_mean_temp+analyzed_mean_PCA1temp+analyzed_gradient_PCA1temp+diff_gradient_temp+econ+n+PC1+PC2,sei=test$`Cond..SE`,data=test)
summary(mr)
vif(mr)

newmods <- data.frame(analyzed_gradient_PCA1temp=seq(0.5,4,by=0.25),
                      t(apply(test[,c("diff_mean_temp","analyzed_mean_PCA1temp","diff_gradient_temp","econ","n","PC1","PC2")],2,mean)))

predict_rma <- cbind(predict(mr,newmods= as.matrix(newmods)),newmods)

plot(test$analyzed_gradient_PCA1temp,test$Estimate)
plot(test$analyzed_mean_PCA1water,test$Estimate)

library(ggplot2)
p1_temp <- ggplot(data=test)+
                  geom_point(aes(x=analyzed_gradient_PCA1temp,y=Estimate))+
                  geom_errorbar(aes(x=analyzed_gradient_PCA1temp,ymin=Estimate-`Cond..SE`,ymax=Estimate+`Cond..SE`))+
                  geom_line(data=predict_rma,aes(x=analyzed_gradient_PCA1temp,y=pred))+
                  annotate("text",-Inf,Inf,hjust=0,vjust=1,label = "(a)")+
                  ylab("Slope of urban usage and temperature relationship")+
                  xlab("Length of temperature PCA1 gradient (Higher = Wider)")+
                  theme_classic()
plot(p1_temp)
###################################################################################
library(phytools)
library(geiger)
test_tree$tip.label <- gsub(" ","_",test_tree$tip.label)
trimmed_tree <- keep.tip(test_tree,test$Species)
Estimate_c <- residuals(mr)
names(Estimate_c) <- test$Species
residual_temp_lambda <- fitContinuous(trimmed_tree,dat=Estimate_c,model=c("lambda"))

Estimate_c <- test$Estimate
se <- test$`Cond..SE`
names(Estimate_c) <- names(se) <-test$Species
#phylosig(trimmed_tree,Estimate_c,method="lambda",se=se,test=T)

temp_lambda <- fitContinuous(trimmed_tree,dat=Estimate_c,SE=se,model=c("lambda"))
temp_lambda2 <- fitContinuous(trimmed_tree,dat=Estimate_c,model=c("lambda"))

temp_an <- phytools::fastAnc(trimmed_tree, Estimate_c, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(trimmed_tree,test$Species),trait = Estimate_c)
nd <- data.frame(node = names(temp_an$ace),trait = temp_an$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
pruned.tree.temp<- full_join(trimmed_tree, d, by = 'node')

p_tree_temp <- ggtree(pruned.tree.temp,layout="circular",ladderize=T,size=1.2,hang=0.1) %<+% as.data.frame(Estimate_c)+
  geom_tree(aes(color=trait),size=1)+
  geom_tiplab(align=T, hjust = -0.1,show.legend=F,size=1.9)+
  scale_color_gradient2(high="#ca0020",mid="#ffffbf",low="#0571b0",midpoint=0,name="Slope")+
  theme(legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6), 
        legend.key.size = unit(0.35, 'cm'),legend.position="bottom",
        legend.box.spacing= unit(-1, 'cm'))

plot(p_tree_temp)  

####################################################################################
library(metafor)
intra_sp_result <- na.omit(intra_sp_result)
intra_sp_result$sp <- gsub(" ","_",intra_sp_result$sp)
intra_sp_result$true_gradient_PCA1water <- intra_sp_result$exotic.PCA1water.max-intra_sp_result$exotic.PCA1water.min
intra_sp_result$true_gradient_PCA1temp <- intra_sp_result$exotic.PCA1temp.max-intra_sp_result$exotic.PCA1temp.min

test <- na.omit(subset(intra_sp_result,var=="PCA1water"))
colnames(test)[[1]] <- "Species"
test <- subset(test,econ!=0 & corr >= -0.7 & corr <= 0.7)

#pca_niche <- rda(test[,c("all.PCA1water.min","all.PCA1water.max","all.PCA1temp.max")]~1,scale=T)
#test$pca1_niche <- vegan::scores(pca_niche,choice=1)$sites
test$breadth_temp <- test$all.PCA1temp.max-test$all.PCA1temp.min
test$breadth_water <- test$all.PCA1water.max-test$all.PCA1water.min
test$diff_gradient_temp <- test$true_gradient_PCA1temp - test$analyzed_gradient_PCA1temp
test$diff_gradient_water <- test$true_gradient_PCA1water - test$analyzed_gradient_PCA1water
test$diff_mean_temp <- test$true_mean_PCA1temp-test$analyzed_mean_PCA1temp
test$diff_mean_water<- test$true_mean_PCA1water-test$analyzed_mean_PCA1water

library(vegan)
library(nlme)
pca_trait <- vegan::rda(log(trait_avg[,2:18])~1,scale=T)
summary(pca_trait)
sig_PCA <- BiodiversityR::PCAsignificance(pca_trait,axes=17)
test <- cbind(test,vegan::scores(pca_trait,choices=1:2)$sites[match(test$Species,gsub(" ","_",trait_avg$new_species)),])

mr <- rma(Estimate~diff_mean_water+analyzed_mean_PCA1water+analyzed_gradient_PCA1water+diff_gradient_water+econ+n+PC1+PC2,sei=test$`Cond..SE`,data=test)
summary(mr)
vif(mr)

newmods <- data.frame(analyzed_gradient_PCA1water=seq(1,5,by=0.25),
                      t(apply(test[,c("diff_mean_water","analyzed_mean_PCA1water","diff_gradient_water","econ","n","PC1","PC2")],2,mean)))

predict_rma <- cbind(predict(mr,newmods= as.matrix(newmods)),newmods)

library(ggplot2)
p1_water <- ggplot(data=test)+
  geom_point(aes(x=analyzed_gradient_PCA1water,y=Estimate))+
  geom_errorbar(aes(x=analyzed_gradient_PCA1water,ymin=Estimate-`Cond..SE`,ymax=Estimate+`Cond..SE`))+
  geom_line(data=predict_rma,aes(x=analyzed_gradient_PCA1water,y=pred))+
  annotate("text",-Inf,Inf,hjust=0,vjust=1,label = "(b)")+
  ylab("Slope of urban usage and soil moisture relationship")+
  xlab("Length of soil moisture PCA1 gradient (Higher = Wider)")+
  theme_classic()
plot(p1_water)

newmods <- data.frame(analyzed_mean_PCA1water=seq(min(test$analyzed_mean_PCA1water),max(test$analyzed_mean_PCA1water),length.out=100),
                      t(apply(test[,c("diff_mean_water","analyzed_gradient_PCA1water","diff_gradient_water","econ","n","PC1","PC2")],2,mean)))

predict_rma <- cbind(predict(mr,newmods= as.matrix(newmods)),newmods)

library(ggplot2)
p2_water <- ggplot(data=test)+
  geom_point(aes(x=analyzed_mean_PCA1water,y=Estimate))+
  geom_errorbar(aes(x=analyzed_mean_PCA1water,ymin=Estimate-`Cond..SE`,ymax=Estimate+`Cond..SE`))+
  geom_line(data=predict_rma,aes(x=analyzed_mean_PCA1water,y=pred))+
  annotate("text",-Inf,Inf,hjust=0,vjust=1,label = "(c)")+
  ylab("Slope of urban usage and soil moisture relationship")+
  xlab("Average soil moisture PCA1 gradient of invaded regions (Higher = Wetter)")+
  theme_classic()
plot(p2_water)

library(ggpubr)
fig1 <- ggarrange(p1_temp,p1_water,p2_water)
ggsave("fig1.tiff",width=24,height=24,dpi=600,units="cm",compression="lzw",bg="white")

library(phytools)
trimmed_tree <- keep.tip(test_tree,test$Species)
Estimate_c <- residuals(mr)
names(Estimate_c) <- test$Species
residual_water_lambda <- fitContinuous(trimmed_tree,dat=Estimate_c,model=c("lambda"))

Estimate_c <- test$Estimate
se <- test$`Cond..SE`
names(Estimate_c) <- names(se) <- test$Species

water_lambda <- fitContinuous(trimmed_tree,dat=Estimate_c,SE=se,model=c("lambda"))

#####################
library(vegan)
library(nlme)
pca_trait <- rda(log(test[,31:47])~1,scale=T)
summary(pca_trait)
sig_PCA <- BiodiversityR::PCAsignificance(pca_trait,axes=25)

#extract scores for axis > 95% of bs model , plus observed var > 5%
test <- cbind(test,vegan::scores(pca_trait,choices=which(as.data.frame(sig_PCA)[6,] >= 0.95 &  as.data.frame(sig_PCA)[2,] >= 5))$sites)

m1 <- gls(PCA1water.min~PC1+PC3+PC4,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m1)
m2 <- gls(PCA1water.max~PC1+PC3+PC4,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m2)
m3 <- gls(PCA1temp.min~PC1+PC3+PC4,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m3)
m4 <- gls(PCA1temp.max~PC1+PC3+PC4,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m4)
m5 <- gls(econ~PC1+PC3+PC4+PCA1water.min+PCA1water.max+PCA1temp.min+PCA1temp.max,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m5)
m6 <- gls(true_gradient_PCA1temp~PC1+PC3+PC4+econ+PCA1water.min+PCA1water.max+PCA1temp.min+PCA1temp.max,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m6)
m7 <- gls(true_gradient_PCA1water~PC1+PC3+PC4+econ+PCA1water.min+PCA1water.max+PCA1temp.min+PCA1temp.max,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m7)
m8 <- gls(analyzed_gradient_PCA1temp~PC1+PC3+PC4+econ+PCA1water.min+PCA1temp.min+PCA1temp.max+PCA1water.max+true_gradient_PCA1temp,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m8)
m9 <- gls(analyzed_gradient_PCA1water~PC1+PC3+PC4+econ+PCA1water.min+PCA1temp.min+PCA1temp.max+PCA1water.max+true_gradient_PCA1water,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m9)

library(piecewiseSEM)
SEM_gradient <- psem(m1,m2,m3,m4,m5,m6,m7,m8,m9)
summary(SEM_gradient)

m6 <- gls(true_gradient_PCA1temp~econ+PCA1temp.min+PCA1temp.max+log(n),data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m6)
m7 <- gls(true_gradient_PCA1water~econ+PCA1water.min+PCA1water.max+log(n),data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m7)
m8 <- gls(analyzed_gradient_PCA1temp~PCA1temp.min+PCA1temp.max+true_gradient_PCA1temp+econ+log(n),data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m8)
m9 <- gls(analyzed_gradient_PCA1water~econ+PCA1water.min+PCA1water.max+true_gradient_PCA1water+log(n),data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m9)

###########old code
library(nlme)

intra_sp_result <- na.omit(intra_sp_result)
intra_sp_result$sp <- gsub(" ","_",intra_sp_result$sp)
test <- na.omit(subset(intra_sp_result,var=="PCA1temp"))
colnames(test)[[1]] <- "Species"

library(phytools)
trimmed_tree <- keep.tip(test_tree,test$Species)
Estimate_c <- test$Estimate
names(Estimate_c) <- test$Species

library(piecewiseSEM)
m1 <- gls(gradient~PCA1water.min+PCA1temp.min+PCA1temp.max+PCA1water.max+econ+group,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m1)
m2 <- gls(Estimate~all.PCA1water.min+all.PCA1temp.min+all.PCA1temp.max+all.PCA1water.max+gradient+econ+group,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m2)

SEM1 <- psem(m2)
summary(SEM1)

intra_sp_result$sp <- gsub(" ","_",intra_sp_result$sp)
test <- na.omit(subset(intra_sp_result,var=="PCA1water"))
colnames(test)[[1]] <- "Species"

m3 <- gls(gradient~all.PCA1water.min+all.PCA1temp.min+all.PCA1temp.max+all.PCA1water.max+econ+group,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m3)
m4 <- gls(Estimate~all.PCA1water.min+all.PCA1temp.min+all.PCA1temp.max+all.PCA1water.max+gradient+econ+group+gradient,data=test,correlation = corPagel(1,phy =test_tree,form=~Species))
summary(m4)
performance::r2(m4)
car::Anova(m4)
check_collinearity(m4)
SEM2 <- psem(m4)
summary(SEM2)
trimmed_tree <- keep.tip(test_tree,test$Species)
Estimate_c <- test$Estimate
names(Estimate_c) <- test$Species
phylosig(trimmed_tree,Estimate_c,method="lambda",test=T)

plot(SEM2)
