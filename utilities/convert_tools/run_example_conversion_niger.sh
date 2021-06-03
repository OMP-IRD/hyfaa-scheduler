
set -e

#change git path to your own !!!
git_path=/home/rejugi/projets/hydro_mgb/dev/cnes-dev_sol_mgb


pbi_file=${git_path}/data/test_configs/data_niger/old/chuvabin_complete.pbi
static_inputs_dir_old_format=${git_path}/data/test_configs/data_niger/old/mgb_iph_static_inputs
output_folder=data

mkdir -p ${output_folder}

#convert static inputs
./convert_mgb_static_inputs.py \
    --compression_level=4 --compression_shuffle=1 \
    --verbose=1 \
    --delta="${static_inputs_dir_old_format}/flag_inland_delta.txt" \
    --delta="${static_inputs_dir_old_format}/flag_northernDelta.txt,500" \
    --delta="${static_inputs_dir_old_format}/flag_southernDelta.txt,30" \
    --specific_outlet="diaka_distribuary,230,600" \
    ${static_inputs_dir_old_format} \
    ${output_folder}/static_data.nc



#~ #convert pbi file to netCDF
#~ ./convert_mgb_dynamic_inputs.py \
    #~ -input=${pbi_file} \
    #~ -output=${output_folder}/forcing_data.nc \
    #~ -start_date='2011-01-01T00:00:00' \
    #~ -nc=11595 \
    #~ -dt=1. \
    #~ --nt=2557
#~ #convert netcdf file to forcing database
#~ ./rain_netcdf_to_forcing_database.py ${output_folder}/forcing_data.nc ${output_folder}/forcing_onmesh_db
