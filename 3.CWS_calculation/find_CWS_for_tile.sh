#!/bin/bash
echo --$(date +%FT%T%z)--

date_str="$(date -u +%Y%m%d)"
date_seconds=$(date +%s)
date_readable=$(date -u +%Y%m%d_%H%M%S%Z)

tile_id=""
dir_0=""
sr16_sp=""
sr16_mht=""
sr16_mdbh=""
sr16_tht=""
sr16_nt=""
sr16_ba=""
sr16_distE=""
sr16_gapE=""
sr16_distW=""
sr16_gapW=""
sr16_distN=""
sr16_gapN=""
sr16_distS=""
sr16_gapS=""
sr16_soilgroup=""
sr16_rooting=""
sr16_tiles=""

# script
find_CWS_command=""

# parameters
save_prefix=""
save_dir=""
season=""

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
	elif [ "${args[i]}" = "-sr16_sp" ]; then
		let i=i+1
	sr16_sp=${args[i]}
	elif [ "${args[i]}" = "-sr16_mht" ]; then
		let i=i+1
	sr16_mht=${args[i]}
	elif [ "${args[i]}" = "-sr16_mdbh" ]; then
		let i=i+1
	sr16_mdbh=${args[i]}
	elif [ "${args[i]}" = "-sr16_tht" ]; then
		let i=i+1
	sr16_tht=${args[i]}
	elif [ "${args[i]}" = "-sr16_nt" ]; then
		let i=i+1
	sr16_nt=${args[i]}
	elif [ "${args[i]}" = "-sr16_ba" ]; then
		let i=i+1
	sr16_ba=${args[i]}
	elif [ "${args[i]}" = "-sr16_distE" ]; then
		let i=i+1
	sr16_distE=${args[i]}
	elif [ "${args[i]}" = "-sr16_gapE" ]; then
		let i=i+1
	sr16_gapE=${args[i]}
	elif [ "${args[i]}" = "-sr16_distW" ]; then
		let i=i+1
	sr16_distW=${args[i]}
	elif [ "${args[i]}" = "-sr16_gapW" ]; then
		let i=i+1
	sr16_gapW=${args[i]}
	elif [ "${args[i]}" = "-sr16_distN" ]; then
		let i=i+1
	sr16_distN=${args[i]}
	elif [ "${args[i]}" = "-sr16_gapN" ]; then
		let i=i+1
	sr16_gapN=${args[i]}
	elif [ "${args[i]}" = "-sr16_distS" ]; then
		let i=i+1
	sr16_distS=${args[i]}
	elif [ "${args[i]}" = "-sr16_gapS" ]; then
		let i=i+1
	sr16_gapS=${args[i]}
	elif [ "${args[i]}" = "-sr16_soilgroup" ]; then
		let i=i+1
	sr16_soilgroup=${args[i]}
	elif [ "${args[i]}" = "-sr16_rooting" ]; then
		let i=i+1
	sr16_rooting=${args[i]}
	elif [ "${args[i]}" = "-sr16_tiles" ]; then
		let i=i+1
	sr16_tiles=${args[i]}
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
if [ -z "$tile_id" ] \
	|| [ -z "$dir_0" ] \
	|| [ -z "${sr16_sp}" ] \
	|| [ -z "${sr16_mht}" ] \
	|| [ -z "${sr16_mdbh}" ] \
	|| [ -z "${sr16_tht}" ] \
	|| [ -z "${sr16_nt}" ] \
	|| [ -z "${sr16_ba}" ] \
	|| [ -z "${sr16_distE}" ] \
	|| [ -z "${sr16_gapE}" ] \
	|| [ -z "${sr16_distW}" ] \
	|| [ -z "${sr16_gapW}" ] \
	|| [ -z "${sr16_distN}" ] \
	|| [ -z "${sr16_gapN}" ] \
	|| [ -z "${sr16_distS}" ] \
	|| [ -z "${sr16_gapS}" ] \
	|| [ -z "${sr16_soilgroup}" ] \
	|| [ -z "${sr16_rooting}" ] \
	|| [ -z "${sr16_tiles}" ] \
	|| [ -z "$find_CWS_command" ] \
	|| [ -z "${save_prefix}" ] \
	|| [ -z "${save_dir}" ] \
	|| [ -z "${season}" ] ; then
echo "ERROR: find_CWS_for_tile.sh failed, argument(s) missing"
echo "tile_id: ${tile_id}"
echo "dir_0: ${dir_0}"
echo "sr16_sp: ${sr16_sp}"
echo "sr16_mht: ${sr16_mht}"
echo "sr16_mdbh: ${sr16_mdbh}"
echo "sr16_tht: ${sr16_tht}"
echo "sr16_nt: ${sr16_nt}"
echo "sr16_ba: ${sr16_ba}"
echo "sr16_distE: ${sr16_distE}"
echo "sr16_gapE: ${sr16_gapE}"
echo "sr16_distW: ${sr16_distW}"
echo "sr16_gapW: ${sr16_gapW}"
echo "sr16_distN: ${sr16_distN}"
echo "sr16_gapN: ${sr16_gapN}"
echo "sr16_distS: ${sr16_distS}"
echo "sr16_gapS: ${sr16_gapS}"
echo "sr16_soilgroup: ${sr16_soilgroup}"
echo "sr16_rooting: ${sr16_rooting}"
echo "sr16_tiles: ${sr16_tiles}"
echo "find_CWS_command: ${find_CWS_command}"
echo "save_prefix: ${save_prefix}"
echo "save_dir: ${save_dir}"
echo "season: ${season}"

echo ""

echo "Usage: find_CWS_for_tile.sh"
echo "                   parameters"

exit 1
fi;

# Go through the ForestGales specific data preparation step and calculate Critical Wind Speed with a single Rscript
# East
Rscript ${find_CWS_command} ${dir_0} ${tile_id} ${sr16_tiles} ${sr16_sp} ${sr16_mht} ${sr16_mdbh} ${sr16_tht} ${sr16_nt} ${sr16_ba} ${sr16_soilgroup} ${sr16_rooting} ${sr16_distE} ${sr16_gapE} ${sr16_distW} ${sr16_gapW} ${sr16_distN} ${sr16_gapN} ${sr16_distS} ${sr16_gapS} ${season} ${save_dir} ${save_prefix}

# compress the created tif file with gdal
for direction in E W S N; do
	if [ -f "${save_dir}/${save_prefix}_${direction}.tif" ]; then
		gdaladdo -ro "${save_dir}/${save_prefix}_${direction}.tif" \
		--config BIGTIFF YES --config COMPRESS_OVERVIEW LZW \
		4 16 64 256
	fi
done

exit 0
