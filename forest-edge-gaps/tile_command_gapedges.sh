#!/bin/bash
# load the geoconda module to access GDAL tools
module load geoconda

# Set date and time
date_str=$(date -u +%y%m%d)
date_seconds=$(date +%s)

date_id="20241111"

# Set variables for the different locations we are using
# First the base directory where the output files and the scripts are located
base_dir="/projappl/project_2011966/WP1/T1-3/Forest_WindSnow"
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

# The scripts run and save the output in the scratch directory - when the outputs are verified, then they can be moved to the projappl and allas directories
# The directory for the edges data is under MOME/AuxData (auxillary data)
scratch_dir="/scratch/project_2011966/MOME/AuxData"
save_dir="${scratch_dir}/edges/save_dir"
mkdir -p "${save_dir}"
# We make a directory for saving the tile rasters
save_tile_dir="${save_dir}/tile_tifs"
mkdir -p "${save_tile_dir}"

# Then we set up the location of the input files - i.e. the MS NFI files
msnfi_dir="/appl/data/geo/luke/vmi/2021"
# We only need the tree height raster
treeht_vrt="${msnfi_dir}/keskipituus_vmi1x_1721.tif"
# there is no raster with tree number information so we discard this constraint (used in Norway) for now
#treenb_vrt
# We locate the tiles information - see the description document for the steps for making this shapefile
msnfi_tiles="${base_dir}/AuxData/tiles/Finland_tiles_100LR_MSNFI.shp"
# name of the layer of the shapefile we are using
layer="Finland_tiles_100LR_MSNFI"

# Finally we locate the bash script which will be used on each of the tiles to find the distance to the edge and forest gap size
script_for_tile="${base_dir}/Scripts/edges/find_gaps_and_distances_for_tile.sh"
find_gaps_and_distances_command="${base_dir}/Scripts/edges/find_gaps_and_distances"

# We set up the parameters
#min_treantall_forest="400"
# NOTE: keskipituus.tif (and gap height) is in dm
# we define a forest gap if there is no forest pixel or the forest pixel is < gap_height = 10m
gap_height="100"
# beyond a certain distance, the distance to the edge/to the forest gap becomes irrelevant in ForestGales (when > than 9xtree_height) so we set up the maximum distance to be 512 m (32 pixels)
max_distance="512"
# same for the gap size
border_length="512"
# prefix for saving the tile rasters
save_prefix="tile"
# EPSG of the rasters
utm_epsg="3067"

# We need to extract the identifier of the tiles to be able to write up a command for each of the tiles.
# To do so:
# first, we obtain the attribute list (from fid) in a geojson file from the msnfi_tiles shapefile
ogr2ogr -f "GeoJSON" ${tmp_dir}/file.geojson ${msnfi_tiles} -dialect sqlite -sql "SELECT DISTINCT CAST(fid AS INT) AS fid FROM "${layer}"" 
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

  echo "${script_for_tile} -tile_id ${tile_id} -tmp_dir ${tmp_dir} -msnfi_tiles ${msnfi_tiles} -layer ${layer} -height_vrt ${treeht_vrt} -find_gaps_and_distances_command ${find_gaps_and_distances_command} -gap_height ${gap_height} -max_distance ${max_distance} -border_length ${border_length} -save_prefix ${save_prefix} -save_dir ${save_tile_dir} -utm_epsg ${utm_epsg} > ${tile_single_log_file} 2>&1" >> "${gnu_parallel_file}"

done
