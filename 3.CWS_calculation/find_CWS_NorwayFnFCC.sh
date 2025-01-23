#!/bin/bash
# { time ./FGrou_SR16_14dec23.sh ; } > log/20231214_FGrou_SR16_14122023.log 2&1
# location stuff
base_dir="/home/nibio/aaa/SR16"
gnu_parallel_dir="${base_dir}/CriticalWindSpeed/Norway_FnFCC/gnu_parallel"
mkdir -p "${gnu_parallel_dir}"
log_dir="${base_dir}/CriticalWindSpeed/Norway_FnFCC/log"
mkdir -p "${log_dir}"
log_tile_dir="${base_dir}/CriticalWindSpeed/Norway_FnFCC/log/log_tile"
mkdir -p "${log_tile_dir}"

save_dir="${base_dir}/CriticalWindSpeed/Norway_FnFCC/save_dir"
mkdir -p "${save_dir}"
save_tile_dir="${save_dir}/tile_tifs"
mkdir -p "${save_tile_dir}"

databank_dir="/run/user/1003/gvfs/smb-share:server=int.nibio.no,share=databank"
sr16_dir="${databank_dir}/Forvaltning/SR16_sluttprodukt/SR16_arsversjon_2023/SR16/UTM_33"
edges_dir="${base_dir}/edges/Data"

# raster files
sr16_sp="${sr16_dir}/SRRTRESLAG.tif"
sr16_mht="${sr16_dir}/SRRMHOYDE.tif"
sr16_mdbh="${sr16_dir}/SRRDIAMMIDDEL_GE8.tif"
sr16_tht="${sr16_dir}/SRROHOYDE.tif"
sr16_nt="${sr16_dir}/SRRTREANTALL_GE8.tif"
sr16_ba="${sr16_dir}/SRRGRFLATE.tif"
sr16_distE="${base_dir}/edges/Data/east_gap_distance.vrt"
sr16_gapE="${base_dir}/edges/Data/east_gap_size.vrt"
sr16_distW="${base_dir}/edges/Data/west_gap_distance.vrt"
sr16_gapW="${base_dir}/edges/Data/west_gap_size.vrt"
sr16_distN="${base_dir}/edges/Data/north_gap_distance.vrt"
sr16_gapN="${base_dir}/edges/Data/north_gap_size.vrt"
sr16_distS="${base_dir}/edges/Data/south_gap_distance.vrt"
sr16_gapS="${base_dir}/edges/Data/south_gap_size.vrt"
sr16_soilgroup="${base_dir}/Soil_Fylke/FGsoilclass_SR16_UTM33N.tif"
sr16_rooting="${base_dir}/Soil_Fylke/FGsoildepth_SR16_UTM33N.tif"

# Tiles vector
sr16_tiles="${base_dir}/Soil_Fylke/NorwayFylkeSR16tiles_UTM33.shp"

#script
script_for_tile="${base_dir}/CriticalWindSpeed/scripts/find_CWS_for_tile.sh"
find_CWS_command="${base_dir}/CriticalWindSpeed/scripts/FGrou_CWS_SR16.R"

#parameters
save_prefix="CWS_tile"
season=$1

# db stuff
db_host="vroom5.int.nibio.no"
dbname="sl"
db_user="sl_mome"
db_pwd="sl_emom"
db_rutenett_schema="org_rutenett"
db_sr16_rutenett_table="sr16_kartruter80_33"

bash_processes="20"

# Find number of tiles
export PGPASSWORD="${db_pwd}"
tile_id_array=($(psql -X -A -t -h "${db_host}" -d "${dbname}" \
		-U "${db_user}" -1 -c "
		SELECT DISTINCT sr16_rutename
			FROM ${db_rutenett_schema}.${db_sr16_rutenett_table};"))

if [ -z "${tile_id_array// }" ]; then
	echo "ERROR: no ruter found."
	exit 1
fi

number_of_tiles=${#tile_id_array[@]}

if (( ${number_of_tiles} < 1 )); then
	echo "ERROR: no tiles found."
	echo "Note: This warning should never occur."
exit 1
fi

gnu_parallel_file="${gnu_parallel_dir}/FGrou_SR16_gnu_parallel_file.txt"

if [ -f "${gnu_parallel_file}" ]; then
	rm "${gnu_parallel_file}"
fi

for tile_id in ${tile_id_array[@]}; do
	tile_enhet_log_file="${log_tile_dir}/tile_id_${tile_id}.log"
	echo "${script_for_tile} -tile_id ${tile_id} -dir_0 ${base_dir} -sr16_sp ${sr16_sp} -sr16_mht ${sr16_mht} -sr16_mdbh ${sr16_mdbh} -sr16_tht ${sr16_tht} -sr16_nt ${sr16_nt} -sr16_ba ${sr16_ba} -sr16_distE ${sr16_distE} -sr16_gapE ${sr16_gapE} -sr16_distW ${sr16_distW} -sr16_gapW ${sr16_gapW} -sr16_distN ${sr16_distN} -sr16_gapN ${sr16_gapN} -sr16_distS ${sr16_distS} -sr16_gapS ${sr16_gapS} -sr16_soilgroup ${sr16_soilgroup} -sr16_rooting ${sr16_rooting} -sr16_tiles ${sr16_tiles} -find_CWS_command ${find_CWS_command} -save_prefix ${save_prefix} -save_dir ${save_tile_dir} -season ${season} > ${tile_enhet_log_file} 2>&1" >> "${gnu_parallel_file}"
done

gnu_parallel_log_file="${log_dir}/find_CWS_gnu_parallel.log"

time parallel --verbose --progress -j ${bash_processes} :::: ${gnu_parallel_file} > ${gnu_parallel_log_file} 2>&1

for direction in E W S N; do
	gdalbuildvrt -overwrite \
		"${save_dir}/CWS_${direction}.vrt" \
		"${save_tile_dir}/${save_prefix}"*"_${direction}.tif"
done

exit 0

