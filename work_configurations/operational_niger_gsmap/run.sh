set -e 

#hyfaa_workdir, hyfaa_temp_dir and hyfaa_config_folder variables must be set not only for this script but also because they are used in configuration files (in hyfaa_config_folder)

export hyfaa_workdir=$(pwd)
mkdir -p ${hyfaa_workdir}/databases
if [ -z "$hyfaa_temp_dir" ]
then
    export hyfaa_temp_dir=${hyfaa_workdir}/temp
else
    if [ ! -z "$TMPDIR" ]
    then
        export hyfaa_temp_dir=${TMPDIR}/hyfaa
    fi
fi

hyfaa_preprocessing_forcing.py --input_yaml_file ${hyfaa_workdir}/config/input.yaml
hyfaa_preprocessing_assimilation.py --input_yaml_file ${hyfaa_workdir}/config/input.yaml
hyfaa_processing.py --input_yaml_file ${hyfaa_workdir}/config/input.yaml
hyfaa_postprocessing.py --input_yaml_file ${hyfaa_workdir}/config/input.yaml



