
set -e


pbi_file=/home/rejugi/projets/hydro_mgb/dev/congo/congo_mgb/Input/CHUVAbin.pbi
static_inputs_dir_old_format=/home/rejugi/projets/hydro_mgb/dev/congo/congo_mgb/Input
output_folder=data_congo

mkdir -p ${output_folder}

#convert static inputs
./convert_mgb_static_inputs.py \
    --compression_level=4 --compression_shuffle=1 \
    --congo_region_mode --verbose=1 \
    --climato="default,/home/rejugi/projets/hydro_mgb/dev/congo/congo_mgb/Input/medias2018.cli" \
    --climato="2011,/home/rejugi/projets/hydro_mgb/dev/congo/congo_mgb/Input/medias.cli" \
    --climato="2012,/home/rejugi/projets/hydro_mgb/dev/congo/congo_mgb/Input/medias.cli" \
    --climato="2013,/home/rejugi/projets/hydro_mgb/dev/congo/congo_mgb/Input/medias2013.cli" \
    --climato="2014,/home/rejugi/projets/hydro_mgb/dev/congo/congo_mgb/Input/medias2014.cli" \
    --climato="2015,/home/rejugi/projets/hydro_mgb/dev/congo/congo_mgb/Input/medias2015.cli" \
    --climato="2016,/home/rejugi/projets/hydro_mgb/dev/congo/congo_mgb/Input/medias2016.cli" \
    --climato="2017,/home/rejugi/projets/hydro_mgb/dev/congo/congo_mgb/Input/medias2017.cli" \
    ${static_inputs_dir_old_format} \
    ${output_folder}/static_data.nc

#~ #convert pbi file to netCDF
#~ ./convert_mgb_dynamic_inputs.py \
    #~ -input=${pbi_file} \
    #~ -output=${output_folder}/forcing_data.nc \
    #~ -start_date='2011-09-30T00:00:00' \
    #~ -nc=9220 \
    #~ -dt=1.
#~ #convert netcdf file to forcing database
#~ ./rain_netcdf_to_forcing_database.py ${output_folder}/forcing_data.nc ${output_folder}/forcing_onmesh_db
