#04 - associate points with polygons - this can speed up the process for future analyses
library(sf)
plants_db_cleaned <- read.csv("plants_cleaned.csv")

region <- read_sf("regions2.shp")
region_df <- read.csv("regionID.csv")
GLONAF <- read.csv("GLONAF.csv")
tdwg3 <- read_sf("tdwg3/level3.shp")

pts <- data.frame(lon=plants_db_cleaned$decimalLongitude,lat=plants_db_cleaned$decimalLatitude)
pts <- st_as_sf(pts,coords= c("lon","lat"))
st_crs(pts) <- st_crs(region)
region <- st_make_valid(region)
tdwg3 <- st_make_valid(tdwg3)

start <- Sys.time()
pts_in_region <- st_intersects(pts,region)
end <- Sys.time()
end-start

write.csv(pts_in_region,"plants_pts_in_polygon.csv")

remove(pts_in_region)
gc()

start <- Sys.time()
pts_in_tdwg3<- st_intersects(pts,tdwg3)
end <- Sys.time()
end-start

write.csv(pts_in_tdwg3,"pts_in_tdwg3.csv")

#############################################################

### Check if polygons are nested - if yes use the largest one (and delete smaller one)
library(sf)
region <- read_sf("regions2.shp")
region_df <- read.csv("regionID.csv")
region <- st_make_valid(region)
overlap_region <- st_within(region)

overlap <- as.data.frame(overlap_region)
overlap$dummy <- overlap$row.id - overlap$col.id
overlap <- overlap[overlap$dummy != 0,]

overlap$Area_row_ID <- as.data.frame(region[overlap$row.id,])$GeodAREA
overlap$Area_col_ID <- as.data.frame(region[overlap$col.id,])$GeodAREA
overlap$Area_diff <- overlap$Area_row_ID - overlap$Area_col_ID

overlap$delete <- ifelse(overlap$Area_diff < 0, overlap$row.id,overlap$col.id)

##############################
plants_db_cleaned <- read.csv("C:/Users/pakno/OneDrive - University of Toronto/GLONAF/plants_cleaned.csv")

pts_in_polygon <- read.csv("plants_pts_in_polygon.csv")

pts_in_polygon <- pts_in_polygon[!pts_in_polygon$col.id %in% overlap$delete,] #remove nested polygons

count <- table(pts_in_polygon$row.id) #a few data points are associated with multiple polygons - delete these (238050 / 37425529)
count <- count[count > 1]

pts_in_polygon$revised_col.id <- ifelse(pts_in_polygon$row.id %in% names(count), NA,pts_in_polygon$col.id) #change col.id to NA if they are associated with multiple polygons

seq <- pts_in_polygon[match(rownames(plants_db_cleaned),pts_in_polygon$row.id),"revised_col.id"]
plants_db_cleaned$polygon_number <- as.data.frame(region)[seq,"OBJIDsic"]

###############################
pts_in_tdwg3 <- read.csv("pts_in_tdwg3.csv")
count <- table(pts_in_tdwg3$row.id) #a few data points are associated with multiple polygons - delete these (74 / 38285903)
count <- count[count > 1]

pts_in_tdwg3$revised_col.id <- ifelse(pts_in_tdwg3$row.id %in% names(count), NA,pts_in_tdwg3$col.id) #change col.id to NA if they are associated with multiple polygons
seq <- pts_in_tdwg3[match(rownames(plants_db_cleaned),pts_in_tdwg3$row.id),"revised_col.id"]

tdwg3 <- read_sf("tdwg3/level3.shp")
plants_db_cleaned$tdwg3 <- tdwg3[seq,"LEVEL3_COD"]
