#01- download data from gbif
library(rgbif)
library(WorldFlora)
library(tidyverse)

name <- vascular.families %>% filter(Group == "angiosperms" | Group == "gymnosperms")
matching_GBIF <- name_backbone_checklist(name,verbose=F) #All exact matches

occ_download(
  pred_in("TAXON_KEY", matching_GBIF$usageKey),
  pred_gte("year",2012),
  pred_lte("year",2022),
  pred("hasCoordinate",TRUE),
  pred("hasGeospatialIssue",FALSE),
  format = "SIMPLE_CSV",
  user="paknok",pwd="Toby5997",email="tpaknok@gmail.com"
)
