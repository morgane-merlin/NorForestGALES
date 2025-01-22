#!/bin/bash

# Set date and time
date_str=$(date -u +%y%m%d)
date_seconds=$(date +%s)

date_id="YYYYMMDD"

# Set variables for the different locations we are using
# First the base directory where the output files and the scripts are located
base_dir="[your_dir]"
# Then we set up the different directories where we will store information while the script is running
# 1. the GNU parallel directory
gnu_parallel_dir="${base_dir}/scripts/edges/gnu_parallel"
mkdir -p "${gnu_parallel_dir}"
# 2. the log directory
log_dir="${base_dir}/scripts/edges/log"
mkdir -p "${log_dir}"
# 3. the log directory for each of the forest tiles
log_tile_dir="${base_dir}/scripts/edges/log/log_tile"
mkdir -p "${log_tile_dir}"
# 4. a temporary directory where we will store temporary files
tmp_dir="${base_dir}/scripts/edges/tmp"
mkdir -p "${tmp_dir}"

# The scripts run and save the output in the save directory
save_dir="${base_dir}/edges/save_dir"
mkdir -p "${save_dir}"
# We make a directory for saving the tile rasters
save_tile_dir="${save_dir}/tile_tifs"
mkdir -p "${save_tile_dir}"

# Then we set up the location of the input files - i.e. the forest raster files
raster_dir="[your_raster_location]"
# We need the tree height raster and number of trees rasters
treeht_vrt="${raster_dir}/tree_height.tif"
treenb_vrt="${raster_dir}/nbtrees.tif"
# With nationwide rasters, it is easier to work at the tile level - where the raster is split into tiles.
# We locate the tiles shapefile
raster_tiles="${base_dir}/tiles.shp"
# name of the layer of the shapefile we are using
layer="[tile_layer]"
# Finally we locate the bash script which will be used on each of the tiles to find the distance to the edge and forest gap size
script_for_tile="${base_dir}/find_gaps_and_distances_for_tile.sh"
# The following comman dis created after running through the steps in the forest-edge-gaps/README.md
find_gaps_and_distances_command="${base_dir}/find_gaps_and_distances"

# We set up the parameters
min_nbtree="400"
# we define a forest gap if there is no forest pixel or the forest pixel is < gap_height = 10m or 100 dm
gap_height="100"
# beyond a certain distance, the distance to the edge/to the forest gap becomes irrelevant in ForestGales (when > than 9xtree_height) so we set up the maximum distance to be 512 m (32 16m pixels)
max_distance="512"
# same for the gap size
border_length="512"
# prefix for saving the tile rasters
save_prefix="tile"
# EPSG of the rasters
utm_epsg="32633"

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
gnu_parallel_file="${gnu_parallel_dir}/find_gaps_and_distances_gnu_parallel_file.txt"
# If the file already exists, we remove it
if [ -f "${gnu_parallel_file}" ]; then
  rm "${gnu_parallel_file}"
fi
# We run a loop on all the tiles, and write out the bash command to find the distance to the edge and the gap size for each of the tiles.
for tile_id in ${tile_id_array[@]}; do

  tile_single_log_file="${log_tile_dir}/tile_id_${tile_id}.log"

  echo "${script_for_tile} -tile_id ${tile_id} -tmp_dir ${tmp_dir} -raster_tiles ${raster_tiles} -layer ${layer} -height_vrt ${treeht_vrt} -nbtrees_vrt ${treenb_vrt} -find_gaps_and_distances_command ${find_gaps_and_distances_command} -gap_height ${gap_height} -min_nbtree ${min_nbtree} -max_distance ${max_distance} -border_length ${border_length} -save_prefix ${save_prefix} -save_dir ${save_tile_dir} -utm_epsg ${utm_epsg} > ${tile_single_log_file} 2>&1" >> "${gnu_parallel_file}"

done
