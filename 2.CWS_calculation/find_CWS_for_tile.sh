#!/bin/bash
echo --$(date +%FT%T%z)--
# We set up the date and time
date_str="$(date -u +%Y%m%d)"
date_seconds=$(date +%s)
date_readable=$(date -u +%Y%m%d_%H%M%S%Z)

# We set up the variables which are required
tile_id=""
dir_0=""
raster_sp=""
raster_mht=""
raster_mdbh=""
raster_tht=""
raster_nt=""
raster_ba=""
raster_distE=""
raster_gapE=""
raster_distW=""
raster_gapW=""
raster_distN=""
raster_gapN=""
raster_distS=""
raster_gapS=""
raster_soilgroup=""
raster_rooting=""
raster_tiles=""

# script
find_CWS_command=""

# parameters
save_prefix=""
save_dir=""
season=""

# Each of those variables was included as an argument to this bash command, so we get the arguments, and we assign them to the variable
# Note that the names need to match, between the names given to the arguments in the tile_command_gapedges.sh and the names of the variables described just above
args=("$@")

i=0
# bash while loop
while [ $i -lt $# ]; do
	if [ "${args[i]}" = "-tile_id" ]; then
		let i=i+1
	tile_id=${args[i]}
	elif [ "${args[i]}" = "-dir_0" ]; then
		let i=i+1
	dir_0=${args[i]}
	elif [ "${args[i]}" = "-raster_sp" ]; then
		let i=i+1
	raster_sp=${args[i]}
	elif [ "${args[i]}" = "-raster_mht" ]; then
		let i=i+1
	raster_mht=${args[i]}
	elif [ "${args[i]}" = "-raster_mdbh" ]; then
		let i=i+1
	raster_mdbh=${args[i]}
	elif [ "${args[i]}" = "-raster_tht" ]; then
		let i=i+1
	raster_tht=${args[i]}
	elif [ "${args[i]}" = "-raster_nt" ]; then
		let i=i+1
	raster_nt=${args[i]}
	elif [ "${args[i]}" = "-raster_ba" ]; then
		let i=i+1
	raster_ba=${args[i]}
	elif [ "${args[i]}" = "-raster_distE" ]; then
		let i=i+1
	raster_distE=${args[i]}
	elif [ "${args[i]}" = "-raster_gapE" ]; then
		let i=i+1
	raster_gapE=${args[i]}
	elif [ "${args[i]}" = "-raster_distW" ]; then
		let i=i+1
	raster_distW=${args[i]}
	elif [ "${args[i]}" = "-raster_gapW" ]; then
		let i=i+1
	raster_gapW=${args[i]}
	elif [ "${args[i]}" = "-raster_distN" ]; then
		let i=i+1
	raster_distN=${args[i]}
	elif [ "${args[i]}" = "-raster_gapN" ]; then
		let i=i+1
	raster_gapN=${args[i]}
	elif [ "${args[i]}" = "-raster_distS" ]; then
		let i=i+1
	raster_distS=${args[i]}
	elif [ "${args[i]}" = "-raster_gapS" ]; then
		let i=i+1
	raster_gapS=${args[i]}
	elif [ "${args[i]}" = "-raster_soilgroup" ]; then
		let i=i+1
	raster_soilgroup=${args[i]}
	elif [ "${args[i]}" = "-raster_rooting" ]; then
		let i=i+1
	raster_rooting=${args[i]}
	elif [ "${args[i]}" = "-raster_tiles" ]; then
		let i=i+1
	raster_tiles=${args[i]}
	elif [ "${args[i]}" = "-find_CWS_command" ]; then
		let i=i+1
	find_CWS_command=${args[i]}
	elif [ "${args[i]}" = "-save_prefix" ]; then
		let i=i+1
	save_prefix=${args[i]}
	elif [ "${args[i]}" = "-save_dir" ]; then
		let i=i+1
	save_dir=${args[i]}
	elif [ "${args[i]}" = "-season" ]; then
		let i=i+1
	season=${args[i]}
	else
		echo unknown input: ${args[i]}
		exit 1
	fi
	let i=i+1
done
# We check if all variables are present, and if not we print out the list of variables to find which one was missing
if [ -z "$tile_id" ] \
	|| [ -z "$dir_0" ] \
	|| [ -z "${raster_sp}" ] \
	|| [ -z "${raster_mht}" ] \
	|| [ -z "${raster_mdbh}" ] \
	|| [ -z "${raster_tht}" ] \
	|| [ -z "${raster_nt}" ] \
	|| [ -z "${raster_ba}" ] \
	|| [ -z "${raster_distE}" ] \
	|| [ -z "${raster_gapE}" ] \
	|| [ -z "${raster_distW}" ] \
	|| [ -z "${raster_gapW}" ] \
	|| [ -z "${raster_distN}" ] \
	|| [ -z "${raster_gapN}" ] \
	|| [ -z "${raster_distS}" ] \
	|| [ -z "${raster_gapS}" ] \
	|| [ -z "${raster_soilgroup}" ] \
	|| [ -z "${raster_rooting}" ] \
	|| [ -z "${raster_tiles}" ] \
	|| [ -z "$find_CWS_command" ] \
	|| [ -z "${save_prefix}" ] \
	|| [ -z "${save_dir}" ] \
	|| [ -z "${season}" ] ; then
echo "ERROR: find_CWS_for_tile.sh failed, argument(s) missing"
echo "tile_id: ${tile_id}"
echo "dir_0: ${dir_0}"
echo "raster_sp: ${raster_sp}"
echo "raster_mht: ${raster_mht}"
echo "raster_mdbh: ${raster_mdbh}"
echo "raster_tht: ${raster_tht}"
echo "raster_nt: ${raster_nt}"
echo "raster_ba: ${raster_ba}"
echo "raster_distE: ${raster_distE}"
echo "raster_gapE: ${raster_gapE}"
echo "raster_distW: ${raster_distW}"
echo "raster_gapW: ${raster_gapW}"
echo "raster_distN: ${raster_distN}"
echo "raster_gapN: ${raster_gapN}"
echo "raster_distS: ${raster_distS}"
echo "raster_gapS: ${raster_gapS}"
echo "raster_soilgroup: ${raster_soilgroup}"
echo "raster_rooting: ${raster_rooting}"
echo "raster_tiles: ${raster_tiles}"
echo "find_CWS_command: ${find_CWS_command}"
echo "save_prefix: ${save_prefix}"
echo "save_dir: ${save_dir}"
echo "season: ${season}"

echo ""

echo "Usage: find_CWS_for_tile.sh"
echo "                   parameters"

exit 1
fi;

# We then run the Rscript to calculate the critical wind speed. The R script processes the data in 2 steps.
# First, the data is combined and goes through the ForestGales specific data preparation step to format the variables as ForestGALES expects it, and calculate any missing variable.
# Second, the core ForestGALES function calculate the Critical Wind Speed.
Rscript ${find_CWS_command} ${dir_0} ${tile_id} ${raster_tiles} ${raster_sp} ${raster_mht} ${raster_mdbh} ${raster_tht} ${raster_nt} ${raster_ba} ${raster_soilgroup} ${raster_rooting} ${raster_distE} ${raster_gapE} ${raster_distW} ${raster_gapW} ${raster_distN} ${raster_gapN} ${raster_distS} ${raster_gapS} ${season} ${save_dir} ${save_prefix}

# We then compress the created tif file with gdal
for direction in E W S N; do
	if [ -f "${save_dir}/${save_prefix}_${direction}.tif" ]; then
		gdaladdo -ro "${save_dir}/${save_prefix}_${direction}.tif" \
		--config BIGTIFF YES --config COMPRESS_OVERVIEW LZW \
		4 16 64 256
	fi
done

exit 0
