library(GIFT)

trait <- read.delim("23800.txt")
trait_meta <- GIFT_traits_meta()

subset_trait <- trait[trait$SpeciesName %in% alien_df$species,]

trait_agg <- subset_trait %>% dplyr::group_by(SpeciesName,TraitName) %>% filter(ValueKindName == "Mean" | ValueKindName == "Single") %>%
  dplyr::summarise(Value=mean(StdValue))
trait_agg <- na.omit(trait_agg)
trait_agg <- subset(trait_agg,TraitName != "")
trait_agg <- trait_agg %>% pivot_wider(names_from = TraitName,values_from=Value)
summary(trait_agg)

selected_trait <- trait_agg[,c("SpeciesName","Seed dry mass","Plant height vegetative","Leaf nitrogen (N) content per leaf dry mass","Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded")]
selected_trait <- na.omit(selected_trait)
colnames(selected_trait) <- c("SpeciesName","Seed_dry_mass","Height","Leaf_N","SLA")

alien_df <- cbind(alien_df,selected_trait[match(alien_df$species,selected_trait$SpeciesName),-1])

m_urb_seed <- glmmTMB(SES_urb~PCA1temp.1*scale(Seed_dry_mass)+PCA1water.1*scale(Seed_dry_mass)+econ+(PCA1temp.1+PCA1water.1||species)+
                   (1|OBJIDsic),data=alien_df,REML=T)
summary(m_urb_seed)
performance::r2(m_urb_seed)

m_urb_h <- glmmTMB(SES_urb~PCA1temp.1*scale(Height)+PCA1water.1*scale(Height)+econ+(PCA1temp.1+PCA1water.1||species)+
                        (scale(Height)||OBJIDsic),data=alien_df,REML=T)
summary(m_urb_h)
performance::r2(m_urb_h)

m_urb_SLA <- glmmTMB(SES_urb~PCA1temp.1*scale(SLA)+PCA1water.1*scale(SLA)+econ+(PCA1temp.1+PCA1water.1||species)+
                     (1|OBJIDsic),data=alien_df,REML=T)
summary(m_urb_SLA)
performance::r2(m_urb_SLA)

m_urb_SLA <- glmmTMB(SES_urb~PCA1temp.1*scale(Leaf_N)+PCA1water.1*scale(Leaf_N)+econ+(PCA1temp.1+PCA1water.1||species)+
                       (1|OBJIDsic),data=alien_df,REML=T)
summary(m_urb_SLA)
performance::r2(m_urb_SLA)

m_urb_overall <- glmmTMB(SES_urb~PCA1temp.1*scale(Seed_dry_mass)+PCA1water.1*scale(Seed_dry_mass)+
                        PCA1temp.1*scale(Height)+PCA1water.1*scale(Height)+
                        PCA1temp.1*scale(SLA)+PCA1water.1*scale(SLA)+
                        econ+(PCA1temp.1+PCA1water.1||species)+
                        (scale(Height)+scale(Seed_dry_mass)+scale(SLA)|OBJIDsic),data=alien_df,REML=T)
summary(m_urb_overall)
performance::r2(m_urb_overall)
