#!/bin/Rscript
rm(list=ls())
args <- commandArgs(TRUE)

require(pacman)
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## packages and functions
pacman::p_load(purrr, data.table, tidyverse, terra, fgRcpp)
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Path to FutureForests forskergruppe
#working.dir <- commandArgs()[1]
#setwd(file.path("~", working.dir))
path_0 <- args[1]

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Compile the Cpp code for ForestGales
#sourceCpp(paste0(path_0, "/FG_Rcpp/FGdataprep_Norge_cpp.cpp"))
#sourceCpp(paste0(path_0, "/FG_Rcpp/FGfunctions_Norge_cpp.cpp"))
exit <- function() {invokeRestart("abort")}
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# get variables from bash
tile_id=args[2]
sr16_tiles=vect(args[3])
species=rast(args[4])
mean_ht=rast(args[5])["SRRMHOYDE_1"]
mean_dbh=rast(args[6])["SRRDIAMMIDDEL_GE8_1"]
top_height=rast(args[7])["SRROHOYDE_1"]
NT=rast(args[8])["SRRTREANTALL_GE8_1"]
BA=rast(args[9])["SRRGRFLATE_1"]
soil_group=rast(args[10])
rooting=rast(args[11])
dist_edge_E=args[12]
gap_size_E=args[13]
dist_edge_W=args[14]
gap_size_W=args[15]
dist_edge_S=args[16]
gap_size_S=args[17]
dist_edge_N=args[18]
gap_size_N=args[19]
season=args[20]
save_dir=args[21]
save_prefix=args[22]
if(length(args)== 23){
snow=rast(args[23])
}

# group gap rasters
gap_rasters_E <- c(dist_edge_E, gap_size_E)
gap_rasters_W <- c(dist_edge_W, gap_size_W)
gap_rasters_N <- c(dist_edge_N, gap_size_N)
gap_rasters_S <- c(dist_edge_S, gap_size_S)


# 1. select ForestGales parameters ####
# Load the ForestGales and the species parameters
fgr_constants_nibio <- read.table(paste0(path_0, "/FGparams/fgr_constants_nibio.txt"), header = T, sep = "\t")
nibio_params <- read.table(paste0(path_0, "/FGparams/nibio_species_parameters_180624.txt"), sep = "\t", header = T) %>%
column_to_rownames(var="rownames")
if(grepl("winter", season)){
nibio_params.sr16 <- nibio_params %>% filter(reference == "Norge" & (leafstatus == "ll" | is.na(leafstatus))) 
} else {
nibio_params.sr16 <- nibio_params %>% filter(reference == "Norge" & (leafstatus == "lo" | is.na(leafstatus))) 
}  

damage.factors <- tibble(mode_of_damage = c("Breakage", "Overturning")) %>%
mutate(damage.factor = as.numeric(as.factor(mode_of_damage)))

## extract the tile and fylke of interest
fylketile <- subset(sr16_tiles, subset = sr16_tiles$SR16_RUTEN == tile_id)
## Check if the tile contains data - some (at least one) doesnt
if(is.empty(ext(fylketile))){
message("This tile has no SR16 data.")
exit()}


## Crop the rasters to the tile
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


# get buffer around fylketile to include all SR16 pixels (some may be at the very border and thus not get assigned a fylke otherwise)
fylketile <- buffer(fylketile, 100)
# rasterize to get fylke of each SR16 pixel
fylke <- rasterize(fylketile, sp, field = "FYLKE")

rm(sr16_tiles, species, soil_group, rooting, mean_ht, mean_dbh, top_height, NT, BA, fylketile)

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
forest.dataprep <- fg_rou_dataprep_Norge_cpp(forest.sel.f, fgr_constants = fgr_constants_nibio, species_parameters = nibio_params.sr16, season = season)

# Calculate the critical wind speed 
message(" -- Critical Wind Speed Calculation --")
forest.CWS <- fg_rou_nibio_cpp(forest.dataprep, fgr_constants = fgr_constants_nibio, species_parameters = nibio_params.sr16, breakage_basecanopy = "no")$results
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
