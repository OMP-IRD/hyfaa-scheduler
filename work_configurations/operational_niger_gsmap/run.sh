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

###########
#get forcing and assimilated data : common step for assimilated, ensemblist and mgbstandard solutions
hyfaa_preprocessing_forcing.py --input_yaml_file ${hyfaa_workdir}/config/input_assimilated_solution.yaml
hyfaa_preprocessing_assimilation.py --input_yaml_file ${hyfaa_workdir}/config/input_assimilated_solution.yaml


###########
#process mgbstandard solution
hyfaa_processing.py --input_yaml_file ${hyfaa_workdir}/config/input_mgbstandard_solution.yaml
hyfaa_postprocessing.py --input_yaml_file ${hyfaa_workdir}/config/input_mgbstandard_solution.yaml

#process ensemblist solution
hyfaa_processing.py --input_yaml_file ${hyfaa_workdir}/config/input_ensemblist_solution.yaml
hyfaa_postprocessing.py --input_yaml_file ${hyfaa_workdir}/config/input_ensemblist_solution.yaml

#process assimilated solution
hyfaa_processing.py --input_yaml_file ${hyfaa_workdir}/config/input_assimilated_solution.yaml
hyfaa_postprocessing.py --input_yaml_file ${hyfaa_workdir}/config/input_assimilated_solution.yaml


###########
#create fake previsions using previous years
fake_previsions_using_previous_years_forcing.py --input_yaml_file ${hyfaa_workdir}/config/input_mgbstandard_solution.yaml --ndays 10 --nyearsmax 10
fake_previsions_using_previous_years_forcing.py --input_yaml_file ${hyfaa_workdir}/config/input_ensemblist_solution.yaml --ndays 10 --nyearsmax 10
fake_previsions_using_previous_years_forcing.py --input_yaml_file ${hyfaa_workdir}/config/input_assimilated_solution.yaml --ndays 10 --nyearsmax 10
