#!/bin/bash

echo --$(date +%FT%T%z)--
# We set up the date and time
date_str="$(date -u +%Y%m%d)"
date_seconds=$(date +%s)
date_readable=$(date -u +%Y%m%d_%H%M%S%Z)


# We set up the variables which are required
tile_id=""
tmp_dir=""
raster_tiles=""
layer=""
height_vrt=""
nbtrees_vrt=""
# script
find_gaps_and_distances_command=""
# parameters
gap_height=""
min_nbtree=""
max_distance=""
border_length=""
save_prefix=""
save_dir=""
utm_epsg=""

# optional
fixed_height_nodata_value="-1"

# Each of those variables was included as an argument to this bash command, so we get the arguments, and we assign them to the variable
# Note that the names need to match, between the names given to the arguments in the tile_command_gapedges.sh and the names of the variables described just above
args=("$@")

i=0
# bash while loop
while [ $i -lt $# ]; do
  if [ "${args[i]}" = "-tile_id" ]; then
    let i=i+1
    tile_id=${args[i]}

  elif [ "${args[i]}" = "-tmp_dir" ]; then
    let i=i+1
    tmp_dir=${args[i]}

  elif [ "${args[i]}" = "-raster_tiles" ]; then
    let i=i+1
    raster_tiles=${args[i]}

  elif [ "${args[i]}" = "-layer" ]; then
    let i=i+1
    layer=${args[i]}

  elif [ "${args[i]}" = "-height_vrt" ]; then
    let i=i+1
    height_vrt=${args[i]}

  elif [ "${args[i]}" = "-nbtrees_vrt" ]; then
    let i=i+1
    nbtrees_vrt=${args[i]}

  elif [ "${args[i]}" = "-find_gaps_and_distances_command" ]; then
    let i=i+1
    find_gaps_and_distances_command=${args[i]}

  elif [ "${args[i]}" = "-gap_height" ]; then
    let i=i+1
    gap_height=${args[i]}

  elif [ "${args[i]}" = "-min_nbtree" ]; then
    let i=i+1
    min_nbtree=${args[i]}

  elif [ "${args[i]}" = "-max_distance" ]; then
    let i=i+1
    max_distance=${args[i]}

  elif [ "${args[i]}" = "-border_length" ]; then
    let i=i+1
    border_length=${args[i]}

  elif [ "${args[i]}" = "-save_prefix" ]; then
    let i=i+1
    save_prefix=${args[i]}

  elif [ "${args[i]}" = "-save_dir" ]; then
    let i=i+1
    save_dir=${args[i]}

  elif [ "${args[i]}" = "-utm_epsg" ]; then
    let i=i+1
    utm_epsg=${args[i]}

  else
    echo unknown input: ${args[i]}
    exit 1
  fi
  let i=i+1
done

# We check if all variables are present, and if not we print out the list of variables to find which one was missing
if [ -z "$tile_id" ] \
       || [ -z "$tmp_dir" ] \
       || [ -z "${raster_tiles}" ] \
       || [ -z "${layer}" ] \
       || [ -z "${height_vrt}" ] \
       || [ -z "${nbtrees_vrt}" ] \
       || [ -z "$find_gaps_and_distances_command" ] \
       || [ -z "${gap_height}" ] \
       || [ -z "${min_nbtree}" ] \
       || [ -z "${max_distance}" ] \
       || [ -z "${border_length}" ] \
       || [ -z "${save_prefix}" ] \
       || [ -z "${save_dir}" ] \
       || [ -z "${utm_epsg}" ]; then

  echo "ERROR: find_gaps_and_distances_for_tile.sh failed, argument(s) missing"
  echo "tile_id: ${tile_id}"
  echo "tmp_dir: ${tmp_dir}"
  echo "raster_tiles: ${raster_tiles}"
  echo "layer: ${layer}"
  echo "height_vrt: ${height_vrt}"
  echo "nbtrees_vrt: ${nbtrees_vrt}"
  echo "find_gaps_and_distances_command: ${find_gaps_and_distances_command}"
  echo "gap_height: ${gap_height}"
  echo "max_distance: ${max_distance}"
  echo "min_nbtree: ${min_nbtree}"
  echo "border_length: ${border_length}"
  echo "save_prefix: ${save_prefix}"
  echo "save_dir: ${save_dir}"
  echo "utm_epsg: ${utm_epsg}"

  echo ""

  echo "Usage: find_gaps_and_distances_for_tile.sh"
  echo "                   parameters"

  exit 1
fi;

# We get the parameters for the gdal_translate command ready.
# First, we need to get the extent of each tile to be able to crop the full rasters to the tile of interest.
# We extract the minimum bounding box for the supplied geometry of the selected tile (MINX, MINY), (MINX, MAXY), (MAXX, MAXY), (MAXX, MINY), (MINX, MINY) and save it in a geojson file
ogr2ogr -t_srs "EPSG:3067" -s_srs "EPSG:3067" -f GeoJSON -dialect sqlite -sql "select ST_AsText(st_envelope(st_buffer(geometry, 512))) as geometry from "${layer}" where fid="${tile_id}"" ${tmp_dir}/bbox_${tile_id}.json ${msnfi_tiles} 
# Then we extract the coordinates and save them in a string
coord0=($(jq -r '.features[].geometry.coordinates' ${tmp_dir}/bbox_${tile_id}.json | tr -d '[],"'))
# We only need distinct values - the initial order is xmin ymin xmax ymin xmin ymax xmax ymax xmin ymin - final order for gdal_translate needs to be xmin ymax xmax ymin ${ulx_uly_lrx_lry} as a text string
coord1=(${coord0[0]} ${coord0[5]} ${coord0[2]} ${coord0[1]})
ulx_uly_lrx_lry="${coord1[@]}"

# we remove unnecessary files
if [ -f "${tmp_dir}/bbox_${tile_id}.json" ]; then
  rm "${tmp_dir}/bbox_${tile_id}.json"
fi

# We subset the full tree height raster to the tile of interest
echo "creating meanht..."
height_tif="${tmp_dir}/height_r_id_${tile_id}.tif"
# if it already exists for the tile, we remove it
if [ -f "${height_tif}" ]; then
  rm "${height_tif}"
fi

echo gdal_translate -projwin ${ulx_uly_lrx_lry} -b 1 -a_srs "EPSG:${utm_epsg}" -of "GTiff" -co "COMPRESS=LZW" -co "TILED=YES" -co "BIGTIFF=YES" "${height_vrt}" "${height_tif}"
# We crop the input tree height raster to the tile of interest, using its extent that we extracted previously, on the first raster band, using gdal_translate
gdal_translate -projwin ${ulx_uly_lrx_lry} \
               -b 1 \
               -a_srs "EPSG:${utm_epsg}" \
               -of "GTiff" \
               -co "COMPRESS=LZW" \
               -co "TILED=YES" \
               -co "BIGTIFF=YES" \
               "${height_vrt}" "${height_tif}"
               
nbtrees_tif="${tmp_dir}/nbtrees_r_id_${rute_id}.tif"

if [ -f "${nbtrees_tif}" ]; then
  rm "${nbtrees_tif}"
fi

gdal_translate -projwin ${ulx_uly_lrx_lry} \
               -b 1 \
               -a_srs "EPSG:${utm_epsg}" \
               -co "COMPRESS=DEFLATE" \
               -co "TILED=YES" \
               -co "BIGTIFF=YES" \
               "${nbtrees_vrt}" "${nbtrees_tif}"

# This part adds a constraint on the number of trees to decide whether a pixel is considered as a forest or non-forest = gap
fixed_height_tif="${tmp_dir}/fixed_height_r_id_${tile_id}.tif"
if [ -f "${fixed_height_tif}" ]; then
  rm "${fixed_height_tif}"
fi
gdal_calc.py --type='Float32' --quiet \
             -A "${height_tif}" \
             -B "${nbtrees_tif}" \
             --outfile="${fixed_height_tif}" \
             --co "COMPRESS=DEFLATE" \
             --co "TILED=YES" \
             --co "BIGTIFF=YES" \
             --NoDataValue="${fixed_height_nodata_value}" \
             --calc "(A * (B >= ${min_nbtree}))"
if [ -f "${height_tif}" ]; then
  rm "${height_tif}"
fi

tile_save_prefix="${save_prefix}_${tile_id}"

# We run the find_gaps_and_distances algorithm using the tree height tile raster and the other parameters supplied in the input
"${find_gaps_and_distances_command}" \
                       -height_tif "${fixed_height_tif}" \
                       -gap_height "${gap_height}" \
                       -max_distance "${max_distance}" \
                       -border_length "${border_length}" \
                       -save_prefix "${tile_save_prefix}" \
                       -save_dir "${save_dir}"

if [ -f "${fixed_height_tif}" ]; then
  rm "${fixed_height_tif}"
fi

# Finally, we add overviews to the created rasters of distance to the edge (gap_distance) and size of the gap (gap_size) for each tile
for direction in east west south north; do

  gdaladdo -ro "${save_dir}/${tile_save_prefix}_${direction}_gap_distance.tif" \
           --config BIGTIFF YES --config COMPRESS_OVERVIEW LZW \
           4 16 64 256 &

  gdaladdo -ro "${save_dir}/${tile_save_prefix}_${direction}_gap_size.tif" \
           --config BIGTIFF YES --config COMPRESS_OVERVIEW LZW \
           4 16 64 256 &

done

exit 0
