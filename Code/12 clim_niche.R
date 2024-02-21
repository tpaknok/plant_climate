library(ape)
library(geiger)
library(nlme)
library(phytools)
library(taxize)
trait <- read.delim("23800.txt")
clim_niche <- read.csv("clim_niche.csv")
econ <- read.csv("EcoUsers3.csv")
subset_trait <- trait[trait$SpeciesName %in% coef_result_df$species,]
unique(trait$TraitName)
wood_d <- subset(subset_trait,TraitName == "Stem specific density (SSD, stem dry mass per stem fresh volume) or wood density" &
                   (ValueKindName == "Mean" | ValueKindName == "Single"))

aggregate(StdValue~SpeciesName,data=wood_d,FUN=mean)

seed_dry_mass <- subset(subset_trait,TraitName == "Seed dry mass" &
                   (ValueKindName == "Mean" | ValueKindName == "Single"))

seed_dry_mass_aggregated <- aggregate(StdValue~SpeciesName,data=seed_dry_mass,FUN=mean)
coef_result_df$seed_mass <- seed_dry_mass_aggregated[match(rownames(coef_result_df),seed_dry_mass_aggregated$SpeciesName),"StdValue"]
coef_result_df$log_seed_mass <- log(coef_result_df$seed_mass)

leaf_N <- subset(subset_trait,TraitName == "Leaf nitrogen (N) content per leaf dry mass"   &
                          (ValueKindName == "Mean" | ValueKindName == "Single"))

aggregate(StdValue~SpeciesName,data=leaf_N,FUN=mean)
leaf_N_aggregated <- aggregate(StdValue~SpeciesName,data=leaf_N,FUN=mean)
#coef_result_df$leaf_N <- leaf_N_aggregated[match(rownames(coef_result_df),leaf_N_aggregated$SpeciesName),"StdValue"]

leaf_P <- subset(subset_trait,TraitName == "Leaf phosphorus (P) content per leaf dry mass"    &
                   (ValueKindName == "Mean" | ValueKindName == "Single"))

aggregate(StdValue~SpeciesName,data=leaf_P,FUN=mean)

P50 <- subset(subset_trait,TraitName == "Xylem hydraulic vulnerability, xylem cavitation vulnerability, embolism vulnerability, (P20, P50, P80)" &
                   (ValueKindName == "Mean" | ValueKindName == "Single")) #too few data...

aggregate(StdValue~SpeciesName,data=P50,FUN=mean)

LMA <- subset(subset_trait,TraitName == "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded" &
                (ValueKindName == "Mean" | ValueKindName == "Single")) 
LMA_aggregated <- aggregate(StdValue~SpeciesName,data=LMA,FUN=mean)

coef_result_df$LMA <- 1/LMA_aggregated[match(rownames(coef_result_df),LMA_aggregated$SpeciesName),"StdValue"]

height <- subset(subset_trait,TraitName == "Plant height vegetative"  &
                (ValueKindName == "Mean" | ValueKindName == "Single")) 
height_aggregated <- aggregate(StdValue~SpeciesName,data=height,FUN=mean)

coef_result_df$height <- height_aggregated[match(rownames(coef_result_df),height_aggregated$SpeciesName),"StdValue"]
coef_result_df$log_height <- log(coef_result_df$height)

leaf_area <- subset(subset_trait,TraitName == "Leaf area (in case of compound leaves: leaf, petiole included)"   &
                   (ValueKindName == "Mean" | ValueKindName == "Single")) 
leaf_area_aggregated <- aggregate(StdValue~SpeciesName,data=leaf_area,FUN=mean)

coef_result_df$height <- height_aggregated[match(rownames(coef_result_df),height_aggregated$SpeciesName),"StdValue"]
coef_result_df$log_height <- log(coef_result_df$height)

econ_gbif <- gbif_parse(econ$species.name)
econ$gbif_species <- econ_gbif$canonicalname
econ$econ_freq <- rowSums(econ[,4:16])
coef_result_df$econ <- econ[match(coef_result_df$species,econ$gbif_species),"econ_freq"]
coef_result_df[is.na(coef_result_df$econ),"econ"] <- 0

###
phy_tree <- read.tree("pruned_plant.tre")

phy_tree$tip.label <- gsub("_"," ",phy_tree$tip.label)

pruned_phy_tree <- drop.tip(phy_tree,phy_tree$tip.label[!phy_tree$tip.label %in% coef_result_df$species])


###
coef_result_df <- cbind(coef_result_df,clim_niche)

coef_result_df$log_count <- log(coef_result_df$count)

coef_result_df_subset<-na.omit(coef_result_df)

pglsModel_1 <- gls(log_height~min_Bio1+min_Arid+econ+log_count, correlation = corPagel(1, pruned_phy_tree),
                          data = coef_result_df_subset)
summary(pglsModel_1)

pglsModel_2 <- gls(LMA~min_Bio1+min_Arid+econ+log_count, correlation = corPagel(1, pruned_phy_tree),
                   data = coef_result_df_subset)
summary(pglsModel_2)

pglsModel_3 <- gls(log_seed_mass~min_Bio1+min_Arid+econ+log_count, correlation = corPagel(1, pruned_phy_tree),
                   data = coef_result_df_subset)
summary(pglsModel_3)

pglsModel_4 <- gls(BIO1_slope_urb~min_Bio1+min_Arid+econ+log_height+LMA+log_seed_mass+econ+log_count, correlation = corPagel(1, pruned_phy_tree),
                 data = coef_result_df_subset)
summary(pglsModel_4)

pglsModel_5 <- gls(Arid_mean_slope_urb~min_Bio1+min_Arid+econ+log_height+LMA+log_seed_mass+econ+log_count, correlation = corPagel(1, pruned_phy_tree),
                   data = coef_result_df_subset)
summary(pglsModel_5)

library(piecewiseSEM)

SEM <- psem(pglsModel_1,pglsModel_2,pglsModel_3,pglsModel_4,pglsModel_5)
summary(SEM)
####
coef_result_df$BIO1_breadth <- coef_result_df$max_Bio1 - coef_result_df$min_Bio1
coef_result_df$Arid_breadth <- coef_result_df$max_Arid - coef_result_df$min_Arid

pglsModel_arid_urb <- gls(Arid_mean_slope_urb~min_Bio1+BIO1_breadth+min_Arid+Arid_breadth+econ+log_count, correlation = corPagel(1, pruned_phy_tree),
                 data = coef_result_df)
summary(pglsModel_arid_urb)

pglsModel_arid_urb_subset <- gls(Arid_mean_slope_urb~min_Bio1+BIO1_breadth+min_Arid+Arid_breadth+econ+log_height+LMA+leaf_N+log_seed_mass+econ+log_count+range_BIO1+range_Arid, correlation = corPagel(1, pruned_phy_tree),
                 data = coef_result_df_subset)
summary(pglsModel_arid_urb_subset)

