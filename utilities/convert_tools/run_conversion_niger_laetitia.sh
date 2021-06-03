
set -e

#change git path to your own !!!
#~ static_inputs_dir_old_format=/home/rejugi/projets/hydro_mgb/dev/laetitia_reprod/Donnees_Niger_Laetitia
static_inputs_dir_old_format=/media/rejugi/hdd_data_linux/hydro_mgb/work/runs/Niger_laetitia/MGB-Niger_GSMAP-NRT-gauges/Input

output_folder=data

mkdir -p ${output_folder}


#convert static inputs
./convert_mgb_static_inputs.py \
    --compression_level=4 --compression_shuffle=1 \
    --verbose=1 \
    --delta="${static_inputs_dir_old_format}/flag_inland_delta.txt" \
    --delta="${static_inputs_dir_old_format}/flag_northernDelta.txt" \
    --delta="${static_inputs_dir_old_format}/flag_southernDelta.txt" \
    --specific_outlet="diaka_distribuary,230,600" \
    ${static_inputs_dir_old_format} \
    ${output_folder}/static_data_corrected.nc


