#!/bin/Rscript
rm(list=ls())
# We collect the arguments supplied to the Rscript call
args <- commandArgs(TRUE)

require(pacman)
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# The Cpp code for ForestGales was compiled in its own small package for simplifying its use in R - fgRcppNorway
pacman::p_load(purrr, data.table, tidyverse, terra, fgRcppNorway)
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# We set up the path for saving the files
path_0 <- args[1]

exit <- function() {invokeRestart("abort")}
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# We get the variables from the bash arguments
print("Loading rasters")
# The tile identifier
tile_id=args[2]
# the full tile shapefile - it includes the region (fylke) as it is key for using the proper stem volume equations
sr16_tiles=vect(args[3])
# the species raster
species=rast(args[4])
# the mean tree height raster
mean_ht=rast(args[5])
# the mean DBH raster
mean_dbh=rast(args[6])
# the top height raster
top_height=rast(args[7])
# the number of trees raster
NT=rast(args[8])
# the basal area raster
BA=rast(args[9])
# note that NT and BA rasters are required to approximate the spacing between trees
# the ForestGALES soil group classification raster
soil_group=rast(args[10])
# the ForestGALES rooting depth classification raster
rooting=rast(args[11])
# the distance to the forest edge and size of the upwind gap for the four cardinal directions - those are already at the tile level since the calculations to obtain these rasters were done at the same tile level
dist_edge_E=rast(args[12])
gap_size_E=rast(args[13])
dist_edge_W=rast(args[14])
gap_size_W=rast(args[15])
dist_edge_S=rast(args[16])
gap_size_S=rast(args[17])
dist_edge_N=rast(args[18])
gap_size_N=rast(args[19])
# the seasonal scenario considered
season=args[20]
# the directory to save the files
save_dir=args[21]
# a text prefix to save the files
save_prefix=args[22]
# in case of the winter scenario, we add the snow raster
if(length(args)== 23){
snow=rast(args[23])
}

# 1. select ForestGales parameters ####
# Load the ForestGales and the species parameters
fgr_constants <- read.table("fgr_constants_1.txt", header = T, sep = "\t")
sp_params <- read.table("species_parameters_1.txt", sep = "\t", header = T) %>%
column_to_rownames(var="rownames") %>%
filter(reference == "Norge") # we filter for the Norwegian parameters
# if the seasonal scenario is winter, we filter the species parameters to the leafless for birch ("ll") 
if(grepl("winter", season)){
sp_params <- sp_params %>% filter((leafstatus == "ll" | is.na(leafstatus))) 
} else {
sp_params <- sp_params %>% filter((leafstatus == "lo" | is.na(leafstatus))) 
}  
# we set up the factor coding for Breakage and Overturning
damage.factors <- tibble(mode_of_damage = c("Breakage", "Overturning")) %>%
mutate(damage.factor = as.numeric(as.factor(mode_of_damage)))

## We extract the tile and fylke (region) of interest
fylketile <- subset(sr16_tiles, subset = sr16_tiles$SR16_RUTEN == tile_id)
## Check if the tile contains data - some (at least one) doesnt
if(is.empty(ext(fylketile))){
message("This tile has no SR16 data.")
exit()}




## Crop the available rasters to the tile
sp <- crop(species, ext(fylketile))
if(global(sp, fun="notNA") == 0){
    message("There is no forest data in this tile")
	exit()
}
s_g <- crop(soil_group, ext(sp))
root <- crop(rooting, ext(sp))
if(exists("snow")){
snow.r <- crop(snow, ext(sp))
rm(snow)
}
mht <- crop(mean_ht, ext(fylketile))
mdbh <- crop(mean_dbh, ext(fylketile))
tht <- crop(top_height, ext(fylketile))
nt <- crop(NT, ext(fylketile))
ba <- crop(BA, ext(fylketile))

# Note that some forest pixels are at the border of the raster and may thus not have been assigned a fylke so we use a buffer
fylketile <- buffer(fylketile, 100)
# We rasterize to get fylke of each SR16 pixel
fylke <- rasterize(fylketile, sp, field = "FYLKE")

# We remove the full rasters
rm(sr16_tiles, species, soil_group, rooting, mean_ht, mean_dbh, top_height, NT, BA, fylketile)
# We create a stacked forest raster with all the tile rasters currently available
forest.0 <- c(sp, mht, mdbh, tht, nt, ba, s_g, root, fylke)

# We loop through the four cardinal directions and calculate the CWS
for(direction in c("E", "N", "S", "W")){
  d_e <- get(paste0("dist_edge_", direction))
  g_s <- get(paste0("gap_size_", direction))
  # Note that the distance to edge and gap size tiles may be smaller than the other raster tiles because we created them with some constraints which are applicable to the calculation of the CWS
  # a few lines should be added to check this and crop as required.
  message(paste0("Direction - ", direction))
  # prepare output file name
  output_name <- paste0(save_dir, "/", save_prefix, tile_id, "_", direction)
  
  # add the distance to the edge and gap size to the stacked forest raster
  edges <- crop(c(d_e, g_s), ext(forest.0))
    forest.sel <- c(forest.0, edges)
  rm(d_e, g_s) # we clean up and remove the rasters
  # we name the raster layers properly
  names(forest.sel) <- c("species", "mean_ht", "mean_dbh", "top_ht", "N8", "BA5", "soil_group", "rooting", "fylke", "dist_edge", "gap_size")
  if(exists("snow.r")){
    forest.sel <- c(forest.sel, snow.r)
    names(forest.sel) <- c("species", "mean_ht", "mean_dbh", "top_ht", "N8", "BA5", "soil_group", "rooting", "fylke", "dist_edge", "gap_size", "snow_load")
  }
  # we convert the raster into a dataframe with cell identification numbers and coordinates
  forest.sel <- as.data.frame(forest.sel, xy = T, cells = T)
  setDT(forest.sel)
  setkey(forest.sel, cell, x, y)
  message(" -- Forest data preparation --")
  # we remove the the forest pixels not meeting the contraints - pixels with mean height under 10 m (note the units), if mean DBH or BA are null and if number of trees is under 400
  forest.sel <- forest.sel[mean_ht > 100 & !is.na(mean_dbh) & species != 0 & N8 >= 400,][!(N8 == 0 & BA5 == 0),]
  # This filter could remove all forest pixels so we catch this possibility
  if(nrow(forest.sel)== 0){
    message("No data in this tile fits the requirements")
    exit()}
  rm(forest.sel)
  
  # Some forest pixels may have odd values (-999) in the soil_group and rooting depth, so we use a simple gap filling using neighboring values 
  forest.sel[, `:=` (soil_group = nafill(soil_group, type = "locf"), rooting = nafill(rooting, type = "locf")), by = y]
  # We recode the species factor 
  forest.sel[, `:=` (species = fcase(species == 1, "NS", species == 2, "SP", species == 3, "BI"),
                     stand_id = cell)]
  # We split the dataframe into a dataframe with cell and coordinates information (forest.sel.loc) and a dataframe with the key variables 
  forest.sel.loc <- copy(forest.sel)[, c("cell", "x", "y")]
  forest.sel.f <- copy(forest.sel)[, c("stand_id", "species", "mean_dbh", "mean_ht", "spacing", "fylke", "soil_group", "rooting", "dist_edge", "gap_size")]
  if(exists("snow.r")){
    forest.sel.f <- copy(forest.sel2)[, c("stand_id", "species", "mean_dbh", "mean_ht", "spacing", "fylke", "soil_group", "rooting", "dist_edge", "gap_size", "snow_load")]
  }
  # We run the first RCpp function to prepare the input forest datatable into the full datatable with all required variables for calculating the CWS
  forest.dataprep <- fg_rou_dataprep_Norge_cpp(forest.sel.f, fgr_constants = fgr_constants, species_parameters = species_params,
                                               season = season)
  rm(forest.sel) # we remove the input forest datatable
  gc() # we clear the R session memory
  # Calculate the critical wind speed 
  message(" -- Critical Wind Speed Calculation --")
  # we use the second RCpp function
  forest.CWS <- fg_rou_nibio_cpp(forest.dataprep, fgr_constants = fgr_constants_FNFI, species_parameters = species_params, breakage_basecanopy = "no")$results
  setDT(forest.CWS)
  # Convert the character variables to numbers
  forest.CWS <- merge(forest.CWS, damage.factors, by = "mode_of_damage", all.x = T, all.y = F)
  rm(forest.dataprep)
  gc()
  # > Saving files ####
  message(" -- Saving files --")
  # Merge the CWS with the coordinates of the cell
  forest.CWS.f <- merge(forest.CWS[, stand_id := as.numeric(stand_id)], forest.sel.loc, by.x = "stand_id", by.y = "cell", all = T)
  # Get full coverage for the coordinates -- only perhaps for small dataframes - technically the problem arises only for tiles 558_396 & 582_396 with 12 and 3 pixels respectively
  # also applies for two more tiles where only one pixel was found so need to expand a bit to be able to create a raster - tiles 446_228 & 590_380 & 518_316 & 574_404
  if(length(forest.CWS.f$x)<=15){
    xyfull <- setDT(expand.grid(x = seq(min(forest.sel.loc$x), max(forest.sel.loc$x)+16, 16),
                                y = seq(min(forest.sel.loc$y), max(forest.sel.loc$y)+16, 16)))
    forest.CWS.f <- merge(forest.CWS.f, xyfull, by = c("x", "y"), all = T)
  }
  forest.CWS.f <- terra::rast(forest.CWS.f[, c("x", "y", "u_elev_damage", "damage.factor", "u_elev_b", "u_elev_o")], type = "xyz", crs = "epsg:32632")
  terra::writeRaster(forest.CWS.f, filename = paste0(output_name, ".tif"), filetype = "GTiff", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TFW=YES"))
  rm(forest.CWS, forest.CWS.f)
  gc()
}


#Loop through the four directions
for(gap_rasters in list(gap_rasters_E, gap_rasters_W, gap_rasters_N, gap_rasters_S)){
dist_edge <- rast(gap_rasters[1])
gap_size <- rast(gap_rasters[2])
# get direction
direction <- toupper(substr(names(dist_edge), 1, 1))
message(paste0("Direction - ", direction))
# prepare output file name
output_name <- paste0(save_dir, "/", save_prefix, tile_id, "_", direction)
# crop the rasters to the extent of the fylketile and other SR16 rasters
d_e <- crop(dist_edge, ext(sp))
g_s <- crop(gap_size, ext(sp))

rm(dist_edge, gap_size)

# select the mean height, top height and mean DBH for each forest cell based on the dominant species
forest.sel <- c(sp, mht, mdbh, tht, nt, ba, s_g, root, d_e, g_s, fylke) 
names(forest.sel) <- c("species", "mean_ht", "mean_dbh", "top_ht", "N8", "BA5", "soil_group", "rooting", "dist_edge", "gap_size", "fylke")
if(exists("snow.r")){
forest.sel <- c(forest.sel, snow.r)
names(forest.sel) <- c("species", "mean_ht", "mean_dbh", "top_ht", "N8", "BA5", "soil_group", "rooting", "dist_edge", "gap_size", "fylke", "snow_load")
}
forest.sel <- as.data.frame(forest.sel, xy = T, cells = T)
setDT(forest.sel)
setkey(forest.sel, cell, x, y)
message(" -- Forest data preparation --")
# remove the pixels with no trees - clearcuts and pixels with no mean_dbh
forest.sel <- forest.sel[mean_ht > 100 & !is.na(mean_dbh) & species != 0 & N8 >= 400,][!(N8 == 0 & BA5 == 0),]
if(nrow(forest.sel)== 0){
message("No data in this tile fits the requirements in height and/or number of trees")
exit()}

# replace the odd values in soil_group and rooting (-999 which may come from different overlaps between soil raster and SR16 in which pixels might land on land or not)
forest.sel[, `:=` (soil_group = nafill(soil_group, type = "locf"), rooting = nafill(rooting, type = "locf")), by = y]
forest.sel[, `:=` (species = fcase(species == 1, "NS",species == 2, "SP",species == 3, "BI"),
mean_ht = mean_ht/10, top_ht = top_ht/10, stand_id = cell, 
spacing = round(fifelse(N8 == 0, 1/sqrt((BA5/(pi*(mean_dbh/2)^2))/10000), 1/sqrt(N8/10000)), 1))]
forest.sel.loc <- copy(forest.sel)[, c("cell", "x", "y")]
forest.sel.f <- copy(forest.sel)[, c("stand_id", "fylke", "species", "mean_dbh", "top_ht", "mean_ht", "spacing", "soil_group", "rooting", "dist_edge", "gap_size")]
if(exists("snow.r")){
forest.sel.f <- copy(forest.sel)[, c("stand_id", "fylke", "species", "mean_dbh", "top_ht", "mean_ht", "spacing", "soil_group", "rooting", "dist_edge", "gap_size", "snow_load")]
}
forest.dataprep <- fg_rou_dataprep_Norge_cpp(forest.sel.f, fgr_constants = fgr_constants_nibio, species_parameters = sp_params.sr16, season = season)

# Calculate the critical wind speed 
message(" -- Critical Wind Speed Calculation --")
forest.CWS <- fg_rou_nibio_cpp(forest.dataprep, fgr_constants = fgr_constants_nibio, species_parameters = sp_params.sr16, breakage_basecanopy = "no")$results
setDT(forest.CWS)
# Convert the character variables to numbers
forest.CWS <- merge(forest.CWS, damage.factors, by = "mode_of_damage", all.x = T, all.y = F)
 
# > Saving files ####
message(" -- Saving files --")
# Merge the CWS with the coordinates of the cell
forest.CWS.f <- merge(forest.CWS[, stand_id := as.numeric(stand_id)], forest.sel.loc, by.x = "stand_id", by.y = "cell", all = T)
# Get full coverage for the coordinates -- only perhaps for small dataframes - technically the problem arises only for tiles 558_396 & 582_396 with 12 and 3 pixels respectively
# also applies for two more tiles where only one pixel was found so need to expand a bit to be able to create a raster - tiles 446_228 & 590_380 & 518_316 & 574_404
if(length(forest.CWS.f$x)<=15){
xyfull <- setDT(expand.grid(x = seq(min(forest.sel.loc$x), max(forest.sel.loc$x)+16, 16),
y = seq(min(forest.sel.loc$y), max(forest.sel.loc$y)+16, 16)))
forest.CWS.f <- merge(forest.CWS.f, xyfull, by = c("x", "y"), all = T)
}
forest.CWS.f <- terra::rast(forest.CWS.f[, c("x", "y", "u_elev_damage", "damage.factor", "u_elev_b", "u_elev_o")], type = "xyz", crs = "epsg:25833")
terra::writeRaster(forest.CWS.f, filename = paste0(output_name, ".tif"), filetype = "GTiff", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TFW=YES"))
}
rm(list = ls())
gc()
