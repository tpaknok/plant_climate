library(sf)
##################now incoporated into script 04. no need run these again
### Check if polygons overlapped - if overlapped use the largest one (and delete smaller one)
plants_db_cleaned <- read.csv("plants_cleaned_standardized.csv")
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

###associate points with polygon 
  
pts_in_polygon <- read.csv("plants_pts_in_polygon.csv")
polygon_df <- read.csv("polygon_df_plants.csv")

polygon_df$OBJISdic <- region_df[match(polygon_df$region_id, region_df$region_id),"OBJIDsic"]
polygon_df$col.id <- match(polygon_df$OBJISdic,region$OBJIDsic)  
pts_in_polygon <- pts_in_polygon[pts_in_polygon$col.id %in% na.omit(unique(polygon_df$col.id)),]

count <- table(pts_in_polygon$row.id)
count <- count[count > 1]

pts_in_polygon$revised_col.id <- ifelse(pts_in_polygon$row.id %in% names(count), NA,pts_in_polygon$col.id)

seq <- pts_in_polygon[match(rownames(plants_db_cleaned),pts_in_polygon$row.id),"revised_col.id"]
plants_db_cleaned$polygon_number <- as.data.frame(region)[seq,"OBJIDsic"]

write.csv(plants_db_cleaned,"plants_cleaned_with_polygon.csv")
 