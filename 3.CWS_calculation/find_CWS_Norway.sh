#!/bin/bash

# we set the variable season from the first argument to the bash command
season=$1

# Set variables for the different locations we are using
# First the base directory where the output files and the scripts are located
base_dir="[your_base_dir]"
# Then we set up the different directories where we will store information while the script is running
# # 1. the GNU parallel directory
gnu_parallel_dir="${base_dir}/gnu_parallel"
mkdir -p "${gnu_parallel_dir}"
# 2. the log directory
log_dir="${base_dir}/log"
mkdir -p "${log_dir}"
# 3. the log directory for each of the forest tiles
log_tile_dir="${base_dir}/log/log_tile"
mkdir -p "${log_tile_dir}"
# 4. a temporary directory where we will store temporary files
tmp_dir="${base_dir}/tmp"
mkdir -p "${tmp_dir}"
# 5. The directory where the files will be stored.
save_dir="${base_dir}/save_dir"
mkdir -p "${save_dir}"
save_tile_dir="${save_dir}/tile_tifs"
mkdir -p "${save_tile_dir}"

# Next are the locations for the input files - fill in with your own.
input_dir="[inputdata_location]"
edges_dir="${base_dir}/edges/Data"

# We indicate the path for the all the required input raster files
# dominant tree species
raster_sp="${input_dir}/species.tif"
# mean height
raster_mht="${input_dir}/mean_height.tif"
# mean DBH
raster_mdbh="${input_dir}/mean_dbh.tif"
# top height
raster_tht="${input_dir}/top_height.tif"
# number of trees
raster_nt="${input_dir}/nbtrees.tif"
# basal area
raster_ba="${input_dir}/basal_area.tif"
# The distance to the forest edge and size of the upwind gap in the four cardinal directions - see scripts from folder NorForestGALES/forest-edge-gaps
# Since those rasters were created with the same tile pattern, we can use those tiles directly, therefore we only write the name of the tile rasters for each of the four cardinal directions
E_gd="east_gap_distance.tif"
E_gs="east_gap_size.tif"
W_gd="west_gap_distance.tif"
W_gs="west_gap_size.tif"
N_gd="north_gap_distance.tif"
N_gs="north_gap_size.tif"
S_gd="south_gap_distance.tif"
S_gs="south_gap_size.tif"
# the ForestGALES soil group and rooting depth classes
raster_soilgroup="${input_dir}/FGsoilclass.tif"
raster_rooting="${input_dir}/FGsoildepth.tif"

# Tiles vector
raster_tiles="${base_dir}/Rastertiles.shp"
# name of the layer of the shapefile we are using
layer="country_tiles"

# We locate the bash script used at the tile level to calculate the CWS and the R script which is called within
script_for_tile="${base_dir}/find_CWS_for_tile.sh"
find_CWS_command="${base_dir}/FGrou_CWS_Norway.R"

# We set the prefix for saving the raster files
save_prefix="CWS_tile"

# We need to extract the identifier of the tiles to be able to write up a command for each of the tiles.
# To do so:
# first, we obtain the attribute list (from fid) in a geojson file from the raster_tiles shapefile
ogr2ogr -f "GeoJSON" ${tmp_dir}/file.geojson ${raster_tiles} -dialect sqlite -sql "SELECT DISTINCT CAST(fid AS INT) AS fid FROM "${layer}"" 
# second, we read the geojson and use jq to create the tile array
tile_id_array=($(jq -r '.features[].properties.fid' ${tmp_dir}/file.geojson  | tr -d '[],"'))

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
if [ -f "${tmp_dir}/file.geojson" ]; then
	rm "${tmp_dir}/file.geojson"
fi

# We set up the file which will hold the full bash command for each of the tiles
gnu_parallel_file="${gnu_parallel_dir}/FGrou_gnu_parallel_file.txt"
# if the file already exists we remove it
if [ -f "${gnu_parallel_file}" ]; then
	rm "${gnu_parallel_file}"
fi

# We run a loop on all the tiles, and write out the bash command to calculate the critical wind speed (CWS) for each of the tiles
# Note that the rasters for the distance to the forest edge and size of the upwind gap are already tiled, therefore we use the tile raster as the input
for tile_id in ${tile_id_array[@]}; do
	tile_enhet_log_file="${log_tile_dir}/tile_id_${tile_id}.log"
	echo "${script_for_tile} -tile_id ${tile_id} -dir_0 ${base_dir} -raster_sp ${raster_sp} -raster_mht ${raster_mht} -raster_mdbh ${raster_mdbh} -raster_tht ${raster_tht} -raster_nt ${raster_nt} -raster_ba ${raster_ba} -raster_distE ${edges_dir}/tile_${tile_id}_${E_gd} -raster_gapE ${edges_dir}/tile_${tile_id}_${E_gs} -raster_distW ${edges_dir}/tile_${tile_id}_${W_gd} -raster_gapW ${edges_dir}/tile_${tile_id}_${W_gs} -raster_distN ${edges_dir}/tile_${tile_id}_${N_gd} -raster_gapN ${edges_dir}/tile_${tile_id}_${N_gs} -raster_distS ${edges_dir}/tile_${tile_id}_${S_gd} -raster_gapS ${edges_dir}/tile_${tile_id}_${S_gs} -raster_soilgroup ${raster_soilgroup} -raster_rooting ${raster_rooting} -raster_tiles ${raster_tiles} -find_CWS_command ${find_CWS_command} -save_prefix ${save_prefix} -save_dir ${save_tile_dir} -season ${season} > ${tile_enhet_log_file} 2>&1" >> "${gnu_parallel_file}"
done

# We set up a log file for following the parallel processes
gnu_parallel_log_file="${log_dir}/find_CWS_gnu_parallel.log"
# We run the commands which are stored in a text file. Each text line is read as a command, and we do so in parallel using the 20 CPUs
time parallel --verbose --progress -j ${bash_processes} :::: ${gnu_parallel_file} > ${gnu_parallel_log_file} 2>&1

# We build a VRT with all the created tile raster files of Critical Wind Speed in the four cardinal directions using GDAL.
for direction in E W S N; do
	gdalbuildvrt -overwrite \
		"${save_dir}/CWS_${direction}.vrt" \
		"${save_tile_dir}/${save_prefix}"*"_${direction}.tif"
done

exit 0

