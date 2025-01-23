#!/bin/bash
#{ time ./proj_CWS2_sr16tiles_test.sh NO_Management_r1i1p1-remo_rcp26 summer; } > log/20241212_testCWSproj.log 2>&1
# location stuff
# We specify the scenario of interest as an argument to the command
scenario=$1
season=$2
base_dir="/home/nibio/aaa/SR16/CriticalWindSpeed/ImprintProjections"
scfolder_dir="${base_dir}/${scenario}"
mkdir -p "${scfolder_dir}"
tmp_dir="${scfolder_dir}/tmp"
mkdir -p "${tmp_dir}"
# Routine CWS directories
cws_logsdir="${scfolder_dir}/routine/cws"
mkdir -p "${cws_logsdir}"
# Routine Edge directories
edges_logsdir="${scfolder_dir}/routine/edges"
mkdir -p "${edges_logsdir}"


databank_dir="/run/user/1003/gvfs/smb-share:server=int.nibio.no,share=databank"
scenario_dir="${databank_dir}/Prosjekter/51160_IMPRINT/Forest_projections/DataR/RESULTS/Clim_Sens_Projections/${scenario}"

bash_processes="20"

# We first need to convert each of the files into rasters for processing the forest edges using the existing script
# 1. Data preparation ####
# We write a small R script which reads in the rds files, calculates the mean DBH
# and writes rasters for each of the tiles. Then we create a vrt and transform it into single rasters.
#preprocessing_script_for_proj="${base_dir}/scripts/init_process_proj.R"
# Get the list of projection files
#shopt -s nullglob dotglob
# Dont use with the Globiom - there are extra files that we don't need
#projections_id_array=("${scenario_dir}"/*)
#preprocessing_log_file="${log_Rpreproc_dir}/preproccessing_${scenario}.log"

# R scripts
# for getting the species raster
coord_proc_Rscript="${base_dir}/scripts/init_coord.R"
# for formatting the required variables as rasters
#treevar_script_for_proj="${base_dir}/scripts/var_processing.R"


# A. CREATE THE SPECIES RASTER & TILE SHAPEFILE FOR ALL THE PROJECTIONS TIMESTEPS ####
#proj_0=$(basename ${projections_id_array[0]}) # wont work for Globiom files
#proj_0="SR16Pred_000.rds"
# 1. Species raster ####
#Rsp_log_file="${scfolder_dir}/routine/Rsp_log_file.log"
# if it exists remove
#if [ -z "${Rsp_log_file}" ]; then
#	rm "${Rsp_log_file}"
#fi
# create folder
sp_folder_tiles="${scfolder_dir}/species/tile_tifs"
mkdir -p "${sp_folder_tiles}"
# make the tile tifs folder
#Rscript ${coord_proc_Rscript} ${scenario_dir} ${proj_0} ${scfolder_dir} > ${Rsp_log_file} 2>&1
# Create species vrt
#gdalbuildvrt -overwrite "${scfolder_dir}/species/species.vrt" "${scfolder_dir}/species/tile_tifs/species_"*".tif"

# We build the full raster
#gdal_translate -b 1 -co "COMPRESS=DEFLATE" -co "TILED=YES" -co "BIGTIFF=YES" "${scfolder_dir}/species/species.vrt" "${scfolder_dir}/species/species.tif"


# 2. Tile information ####
# The tile information is in "/home/nibio/aaa/SR16/AuxData/Soil_Fylke/NorwayFylkeSR16tiles_UTM32.shp" which was
# created using : S:/Forvaltning/Matching2/NORGE/rutenett/sr16_kartruter80hel.shp and fylke information
# see T:/Aktive/DSU/3737/53280_Future_Forests/Akt2_Mapping_Disturbances/T2_1_Wind/SR16tiles_fylke.R
# We can extract the tiles names from this shapefile
# Note that we need to use two different versions for SR16tiles: the original version "S:/Forvaltning/Matching2/NORGE/rutenett/sr16_kartruter80hel.shp"
# gives us the accurate tile names and extent (most important). The second version with the fylke is for use in the Rscript for calculating CWS

SR16tiles_o="/home/nibio/aaa/SR16/AuxData/Soil_Fylke/sr16_kartruter80hel.shp"
SR16tiles_f="/home/nibio/aaa/SR16/AuxData/Soil_Fylke/NorwayFylkeSR16tiles_UTM32.shp"

# We get the tile array from there
# first, we obtain the attribute list (from fid) in a geojson file from the Tiles shapefile
ogr2ogr -f "GeoJSON" ${scfolder_dir}/routine/sr16tiles.geojson ${SR16tiles_o} -dialect sqlite -sql "SELECT DISTINCT SR16_RUTEN AS fid FROM sr16_kartruter80hel" 
# second, we read the geojson and use jq to create the tile array
tile_id_array=($(jq -r '.features[].properties.fid' ${scfolder_dir}/routine/sr16tiles.geojson  | tr -d '[],"'))
# Check if the tile_id_array is not empty - first if the variable is empty
if [ -z "${tile_id_array// }" ]; then
  echo "ERROR: no tiles found."
  exit 1
fi
# second if the length of the variable is less than one
number_of_tiles=${#tile_id_array[@]}
if (( ${number_of_tiles} < 1 )); then
  echo "ERROR: no tiles found."
  echo "Note: This warning should never occur."
  exit 1
fi
# we remove the temporary geojson file
if [ -f "${scfolder_dir}/routine/sr16tiles.geojson" ]; then
  rm "${scfolder_dir}/routine/sr16tiles.geojson"
fi

# B. PREPARATION OF THE FOREST VARIABLES ####
# 1. Required variables & scripts ####
# Edges
edgescript_for_tile="${base_dir}/scripts/edges/find_gaps_and_distances_for_tile.sh"
find_gaps_and_distances_command="${base_dir}/scripts/edges/find_gaps_and_distances"
# we define a forest gap if there is no forest pixel or the forest pixel is < gap_height = 10m
gap_height="10"
min_treantall_forest="400"
# beyond a certain distance, the distance to the edge/to the forest gap becomes irrelevant in ForestGales (when > than 9xtree_height) so we set up the maximum distance to be 512 m (32 pixels)
max_distance="512"
# same for the gap size
border_length="512"
# prefix for saving the tile rasters
save_prefix="tile"
# EPSG of the rasters
utm_epsg="32632"

# Edges scripts
script_for_tile="${base_dir}/scripts/edges/find_gaps_and_distances_for_tile.sh"
find_gaps_and_distances_command="${base_dir}/scripts/edges/find_gaps_and_distances"
#projections_id_array=("SR16Pred_000.rds" "SR16Pred_001.rds" "SR16Pred_002.rds" "SR16Pred_003.rds" "SR16Pred_004.rds" "SR16Pred_005.rds" "SR16Pred_006.rds" "SR16Pred_007.rds" "SR16Pred_008.rds" "SR16Pred_009.rds" "SR16Pred_010.rds" "SR16Pred_011.rds" "SR16Pred_012.rds" "SR16Pred_013.rds" "SR16Pred_014.rds" "SR16Pred_015.rds")
projections_id_array=("SR16Pred_006.rds" "SR16Pred_007.rds" "SR16Pred_008.rds" "SR16Pred_009.rds" "SR16Pred_010.rds" "SR16Pred_011.rds" "SR16Pred_012.rds" "SR16Pred_013.rds" "SR16Pred_014.rds" "SR16Pred_015.rds")
#projections_id_array=("SR16Pred_005.rds")
# 2. Loop for all the projection timesteps ####
# Get the mean height, mean DBH and spacing variables for all the projection timestamps
for proj_file in "${projections_id_array[@]}"; do
	echo ${proj_file}
	# 1. create the required variable rasters from the RDS files
	proj_file=$(basename "${proj_file}")
	proj_dir=${proj_file::-4} #get the name of the projection timestamp without the file extension
	echo "${proj_dir}"
	# Tree variables processing
	# create folder to save the data
	save_dir="${scfolder_dir}/${proj_dir}"
	mkdir -p "${save_dir}"
	# create the folders for saving the variables
	#height_dir="${save_dir}/height/tile_tifs"
	#mkdir -p "${height_dir}"
	#dbh_dir="${save_dir}/dbh/tile_tifs"
	#mkdir -p "${dbh_dir}"
	#spacing_dir="${save_dir}/spacing/tile_tifs"
	#mkdir -p "${spacing_dir}"
	#treantall_dir="${save_dir}/treantall/tile_tifs"
	#mkdir -p "${treantall_dir}"
	#Rvars_log_file="${scfolder_dir}/routine/Rvars_log_file.log"
	# if it exists, remove
	#if [ -z "${Rvars_log_file}" ]; then
	#	rm "${Rvars_log_file}"
	#fi
	#Rscript ${treevar_script_for_proj} ${scenario_dir} ${proj_file} ${scfolder_dir} ${save_dir} > ${Rvars_log_file} 2>&1
	# Build the vrts
	#gdalbuildvrt -overwrite "${save_dir}/height/meanHt.vrt" "${save_dir}/height/tile_tifs/meanHt_"*".tif"
	#gdalbuildvrt -overwrite "${save_dir}/dbh/meanDBH.vrt" "${save_dir}/dbh/tile_tifs/meanDBH_"*".tif"
	#gdalbuildvrt -overwrite "${save_dir}/spacing/spacing.vrt" "${save_dir}/spacing/tile_tifs/spacing_"*".tif"
	#gdalbuildvrt -overwrite "${save_dir}/treantall/treantall.vrt" "${save_dir}/treantall/tile_tifs/treantall_"*".tif"

	# We need to create actual tifs from the vrt and tiles for the edges calculation - so for height and treantall at least since we need
	# to expand a bit outside of the tiles when finding the forest edges.
	#gdal_translate -co "COMPRESS=DEFLATE" -co "TILED=YES" -co "BIGTIFF=YES" "${save_dir}/treantall/treantall.vrt" "${save_dir}/treantall.tif"
	#gdal_translate -co "COMPRESS=DEFLATE" -co "TILED=YES" -co "BIGTIFF=YES" "${save_dir}/height/meanHt.vrt" "${save_dir}/meanHt.tif"
	#gdal_translate -co "COMPRESS=DEFLATE" -co "TILED=YES" -co "BIGTIFF=YES" "${save_dir}/dbh/meanDBH.vrt" "${save_dir}/meanDBH.tif"
	#gdal_translate -co "COMPRESS=DEFLATE" -co "TILED=YES" -co "BIGTIFF=YES" "${save_dir}/spacing/spacing.vrt" "${save_dir}/spacing.tif"
	# Make the overviews
	#gdaladdo -ro "${save_dir}/meanHt.tif" --config BIGTIFF YES --config COMPRESS_OVERVIEW LZW 4 16 64 256
	#gdaladdo -ro "${save_dir}/meanDBH.tif" --config BIGTIFF YES --config COMPRESS_OVERVIEW LZW 4 16 64 256
        #gdaladdo -ro "${save_dir}/treantall.tif" --config BIGTIFF YES --config COMPRESS_OVERVIEW LZW 4 16 64 256
	#gdaladdo -ro "${save_dir}/spacing.tif" --config BIGTIFF YES --config COMPRESS_OVERVIEW LZW 4 16 64 256
	# Remove the individual tiles and vrt
	#rm -rf "${save_dir}/height/" "${save_dir}/dbh/" "${save_dir}/spacing/" "${save_dir}/treantall/"

	# 2. Calculating the edges
	# Log directories
	#edges_gnu_dir="${edges_logsdir}/${proj_dir}/gnu_parallel"
	#mkdir -p "${edges_gnu_dir}"
	#edges_log_tile_dir="${edges_logsdir}/${proj_dir}/log/log_tile"
	#mkdir -p "${edges_log_tile_dir}"
	#edges_tmp_dir="${edges_logsdir}/${proj_dir}/tmp"
	#mkdir -p "${edges_tmp_dir}"
	# We make a directory for saving the tile rasters
	edge_dir="${save_dir}/edges/tile_tifs"
	mkdir -p "${edge_dir}"
	# Required tifs
	ht_tif="${save_dir}/meanHt.tif"
	treantall_vrt="${save_dir}/treantall.tif"
	# Create the file which records the command for each of the tiles
	edges_gnu_parallel_file="${edges_gnu_dir}/find_gaps_and_distances_gnu_parallel_file.txt"
	edges_gnu_parallel_log_file="${edges_logsdir}/${proj_dir}/log/find_gaps_and_distances_gnu_parallel.log"
	# If the file already exists, we remove it
	#if [ -f "${edges_gnu_parallel_file}" ]; then
	#	rm "${edges_gnu_parallel_file}"
	#fi
	# We run a loop on all the tiles, and write out the bash command to find the distance to the edge and the gap size for each of the tiles.
	#for tile_id in ${tile_id_array[@]}; do
	#	tile_single_log_file="${edges_log_tile_dir}/tile_id_${tile_id}.log"
	#	echo "${script_for_tile} -tile_id ${tile_id} -tmp_dir ${edges_tmp_dir} -sr16proj_tile ${SR16tiles_o} -height_vrt ${ht_tif} -treantall_vrt ${treantall_vrt} -find_gaps_and_distances_command ${find_gaps_and_distances_command} -gap_height ${gap_height} -min_treantall_forest ${min_treantall_forest} -max_distance ${max_distance} -border_length ${border_length} -save_prefix ${save_prefix} -save_dir ${edge_dir} -utm_epsg ${utm_epsg} > ${tile_single_log_file} 2>&1" >> "${edges_gnu_parallel_file}"
	#done
	# We run the commands which are stored in a text file. Each text line is read as a command, and we do so in parallel using the 20 CPUs
	#time parallel --verbose --progress -j $bash_processes :::: ${edges_gnu_parallel_file} > ${edges_gnu_parallel_log_file} 2>&1
	# We build a VRT with all the created tile raster files for the distance to the forest edge and the size of the forest gap in the four cardinal directions using GDAL.
	#for direction in east west south north; do
	#	gdalbuildvrt -overwrite "${save_dir}/edges/${direction}_gap_distance.vrt" "${edge_dir}/tile_"*"${direction}_gap_distance.tif"
	#	gdalbuildvrt -overwrite "${save_dir}/edges/${direction}_gap_size.vrt" "${edge_dir}/tile_"*"${direction}_gap_size.tif"
	#done
	# 3. Calculate the CWS ####
	# Log directories
	cwslog_dir="${cws_logsdir}/${proj_dir}/${season}"
	mkdir -p "${cwslog_dir}"
	cws_gnu_parallel_dir="${cwslog_dir}/gnu_parallel"
	mkdir -p "${cws_gnu_parallel_dir}"
	cws_log_dir="${cwslog_dir}/log"
	mkdir -p "${cws_log_dir}"
	cws_log_tile_dir="${cws_log_dir}/log_tile"
	mkdir -p "${cws_log_tile_dir}"
	cws_tmp_dir="${cwslog_dir}/tmp"
	mkdir -p "${cws_tmp_dir}"
	# We make a directory for saving the tile rasters
	cws_dir="${save_dir}/CWS/${season}/tile_tifs"
	mkdir -p "${cws_dir}"
	# Required tifs - for the edges we have them already as tiles - use SR16 tiles
	species_tif="${scfolder_dir}/species/species.tif"
	dbh_tif="${save_dir}/meanDBH.tif"
	# ht_tif="${save_dir}/meanHt.tif" already defined from the edges
	spacing_tif="${save_dir}/spacing.tif"
	sr16proj_soilgroup="/home/nibio/aaa/SR16/AuxData/Soil_Fylke/FGsoilclass_SR16_UTM32N.tif"
	sr16proj_rooting="/home/nibio/aaa/SR16/AuxData/Soil_Fylke/FGsoildepth_SR16_UTM32N.tif"
	utm_epsg="32632"
	# if the season is winter with snow, then we add a variable and use a slightly different bash script
	if [[ "$season" == "wintersnow" ]]; then
	snow="/home/nibio/aaa/SR16/AuxData/AuxData/Snow/Stot_avgwintermax_20132023_UTM32N_SR16.tif"
		cws_script_for_tile="${base_dir}/scripts/CWS/find_CWS_SR16proj_for_tile_snow.sh"
	else
		cws_script_for_tile="${base_dir}/scripts/CWS/find_CWS_SR16proj_for_tile.sh"
	fi
	find_CWS_command="${base_dir}/scripts/CWS/FGrou_CWS_SR16proj.R"
	cws_save_prefix="CWS_tile"
	# Create the file which records the command for each of the tiles
	cws_gnu_parallel_file="${cws_gnu_parallel_dir}/cws_gnu_parallel_file.txt"
	cws_gnu_parallel_log_file="${cws_log_dir}/cws_gnu_parallel.log"
	# If the file already exists, we remove it
	if [ -f "${cws_gnu_parallel_file}" ]; then
		rm "${cws_gnu_parallel_file}"
	fi
	# We run a loop on all the tiles, and write out the bash command to find the distance to the edge and the gap size for each of the tiles.
	# remove existing files
	if [ -f "${cws_gnu_parallel_file}" ]; then
		rm "${cws_gnu_parallel_file}"
	fi
	for tile_id in ${tile_id_array[@]}; do
		tile_single_log_file="${cws_log_tile_dir}/tile_id_${tile_id}.log"
		if [[ "$season" == "wintersnow" ]]; then
			echo "${cws_script_for_tile} -utm_epsg ${utm_epsg} -tile_id ${tile_id} -dir_0 ${base_dir} -sr16proj_sp ${species_tif} -sr16proj_mht ${ht_tif} -sr16proj_mdbh ${dbh_tif} -sr16proj_spacing ${spacing_tif} -sr16proj_edgedir ${edge_dir} -sr16proj_soilgroup ${sr16proj_soilgroup} -sr16proj_rooting ${sr16proj_rooting} -sr16proj_tile_f ${SR16tiles_f} -sr16proj_tile_o ${SR16tiles_o} -find_CWS_command ${find_CWS_command} -save_prefix ${cws_save_prefix} -save_dir ${cws_dir} -tmp_dir ${cws_tmp_dir} -season ${season} -snow ${snow} > ${tile_single_log_file} 2>&1" >> "${cws_gnu_parallel_file}"
		else
			echo "${cws_script_for_tile} -utm_epsg ${utm_epsg} -tile_id ${tile_id} -dir_0 ${base_dir} -sr16proj_sp ${species_tif} -sr16proj_mht ${ht_tif} -sr16proj_mdbh ${dbh_tif} -sr16proj_spacing ${spacing_tif} -sr16proj_edgedir ${edge_dir} -sr16proj_soilgroup ${sr16proj_soilgroup} -sr16proj_rooting ${sr16proj_rooting} -sr16proj_tile_f ${SR16tiles_f} -sr16proj_tile_o ${SR16tiles_o} -find_CWS_command ${find_CWS_command} -save_prefix ${cws_save_prefix} -save_dir ${cws_dir} -tmp_dir ${cws_tmp_dir} -season ${season} > ${tile_single_log_file} 2>&1" >> "${cws_gnu_parallel_file}"
		fi
	done
	time parallel --verbose --progress -j $bash_processes :::: ${cws_gnu_parallel_file} > ${cws_gnu_parallel_log_file} 2>&1
	# We build a VRT with all the created tile raster files for the distance to the forest edge and the size of the forest gap in the four cardinal directions using GDAL.
	for direction in E W S N; do
		gdalbuildvrt -overwrite  "${save_dir}/CWS/${season}/CWS_${direction}.vrt" "${cws_dir}/${cws_save_prefix}"*"_${direction}.tif"
	done

	# Run the Grass algorithms
	grass /home/nibio/aaa/grass_mappe/GRASS_ENV/PERMANENT --exec sh "${base_dir}/scripts/postprocessing_CWS.sh" "${scenario}" "${season}" "${proj_dir}"

	# Remove the folder containing the tile_tifs - the largest heaviest folder and we can't have that much data stored on the machine or on T:/S:
	#rm -rf "${save_dir}/CWS/${season}/tile_tifs/"
        #rm "${save_dir}/CWS/${season}/"*.vrt

	# Copy files to the prosjekt
	#cp -r -u "${save_dir}/" "/run/user/1003/gvfs/smb-share:server=int.nibio.no,share=databank/Prosjekter/53280_Future_Forests/Projections/${scenario}/"
	#echo "All done !"

done

exit 0

