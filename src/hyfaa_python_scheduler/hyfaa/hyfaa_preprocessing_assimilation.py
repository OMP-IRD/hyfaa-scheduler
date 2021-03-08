#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################################
#  This code has been developped by Magellium SAS
#
#  Licensing:
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>
#######################################################################


from hyfaa.common.common_functions import *
from multiprocessing import cpu_count
import glob

from hyfaa.common.test.common_test_functions import *
from hyfaa.check_input_parameters.check_input_parameters import *
from hyfaa.common.yaml.yaml_parser import load_yaml
from hyfaa.common.parallel_processing.easy_parallel import simple_parallel_run, Job
from hyfaa.database.assimilation.assimilation_db import Assimilation_Database


from hyfaa.database.assimilation.hysope_download_module import get_hysope_assimilation_data as get_assimilation_data


    

def hyfaa_preprocessing_assimilation(yaml_file_or_dict, verbose=None):
    """main scheduler processing function
    
    :param yaml_file_or_dict: python dict or yaml file
    """
    
    #check if input is a python dict or a yaml file
    valid_obj = False
    if isinstance(yaml_file_or_dict, dict):
        valid_obj = True
        dico = yaml_file_or_dict
    else:
        if hasattr(yaml_file_or_dict, 'replace'): #pythonic way of checking if this is a python3 string, python2 string or python2 unicode
            if os.path.exists(yaml_file_or_dict):
                valid_obj = True
                dico = load_yaml(yaml_file_or_dict, env_vars=True)
    if not valid_obj:
        raise Exception('input must be a python dict or a path to an existing yaml file')
        
    #check parameters
    dico = check_parameters(dico)
    exec_time = dico['exec_time']
    if verbose is None:
        if 'verbose' in dico:
            verbose = dico['verbose']
    if verbose is None:
        verbose = 0
    
    #if no assimilation sources, exit
    if len(dico['assimilation_sources']) == 0:
        return
        
    if verbose >= 1:
        print('Assimilation pre-prococessing launched...')
    
    if verbose > 1:
        print('  Assimilation pre-proc: Parameters and mesh coordinates loaded...')
    os.makedirs(dico['assimilation_database_directory'], exist_ok=True)
    
    min_time = dico['init_conditions']['time_start']
    max_time = dico['exec_time']
    
    if verbose > 1:
        print('  Assimilation pre-proc: Opening assimilation database...')


    #update forcing mesh DB
    with Assimilation_Database(dico['assimilation_database_directory'], mode='w', verbose=0) as db:

        for source_file in dico['assimilation_sources']:
            
            if verbose > 1:
                print('  Assimilation pre-proc: loading DB data')
            
            db_data = db.get_values_between_dates(date_start=min_time, date_end=max_time)
            db_id_tuples = set([(db_data['mesh_id'][ii], db_data['sv_name'][ii], db.str2date(db_data['date_data'][ii]), round(db_data['value'][ii],3), round(db_data['uncertainty'][ii],3)) \
                for ii in range(len(db_data['mesh_id']))])
            db_simpleid_tuples = set([(el[0], el[1], el[2]) for el in db_id_tuples])
            del db_data
            
            if verbose > 1:
                print('  Assimilation pre-proc: loading source data')
            
            if isinstance(source_file, dict):
                dico_source = source_file
            else:
                dico_source = load_yaml(source_file, env_vars=False)
            
            if verbose > 1:
                print('  Assimilation pre-proc: executing get_assimilation_data function')
            
            info_data = get_assimilation_data(dico_source, min_date=min_time, max_date=max_time, verbose=verbose)
            
            if verbose > 1:
                print('  Assimilation pre-proc: adding new data to DB')
            
            info_add = []
            info_remove = []
            for elem in info_data:
                try:
                    id_tuple = (elem['mesh_id'], elem['sv_name'], elem['date_data'], round(elem['value'],3), round(elem['uncertainty'],3))
                except:
                    print(elem)
                    raise
                simpleid_tuple = (id_tuple[0], id_tuple[1], id_tuple[2])
                if simpleid_tuple in db_simpleid_tuples:
                    if id_tuple in db_id_tuples:
                        #already present
                        continue
                    else:
                        #already present but data changed => replacing
                        info_remove.append(simpleid_tuple)
                info_add.append(elem)
            db.remove(info_remove)
            db.add(info_add)

            if verbose >= 1:
                print('Updated database with data requested from conf file %s: %d points added, %d removed'%(source_file, len(info_add), len(info_remove)))
    
    

    
    
def test_main_program():
    
    test_use_web = False
    if 'test_use_web' in os.environ:
        if os.environ['test_use_web'].lower() == 'true':
            test_use_web = internet_available()
    
    if test_use_web:
        print_utest_message('Test full retrieval process for assimilation data: success if runs without errors')
        
        main_test_dir = tempfile.mkdtemp(dir=os.path.abspath('.'), prefix='temppreprocassim_')
        try:
        
            dico_input = get_scheduler_basic_main_parameters_dict()

            dico_input['init_conditions']['time_start'] = datetime(2000,1,1)
            dico_input['exec_time'] = datetime(2000,1,2)
            dico_input['hydrological_states_database_directory'] = '%s/hydrological_states_database_directory'%main_test_dir
            dico_input['forcing_hours'] = [0, 12]
            dico_input['forcing_grid_geo_selection'] = {'lonmin': -12., 'lonmax': 16., 'latmin': 4., 'latmax': 25.}
            dico_input['forcing_grid_database_directory'] = '%s/forcing_grid_database_directory'%main_test_dir
            dico_input['forcing_onmesh_database_directory'] = '%s/forcing_onmesh_database_directory'%main_test_dir
            dico_input['assimilation_database_directory'] = '%s/assimilation_database_directory'%main_test_dir
            dico_input['assimilation_sources'] = [{'function': 'get_hysope_data', \
                'function_info': {'L_baikal':{'mesh_id': 10}, 'R_amz_nap_env_0422_01': {'mesh_id': 12}}}]
            dico_input['temporary_files_directory'] = '%s/temporary_files_directory'%main_test_dir
            dico_input['output_directory']['directory'] = '%s/output_directory'%main_test_dir
            dico_input['assim_params_file'] = 'None'

            hyfaa_preprocessing_assimilation(dico_input, verbose=3)
            success_message_utest()
        except:
            fail_message_utest()
        finally:
            shutil.rmtree(main_test_dir)

    
    
if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="HYFAA pre-processing assimilation chain", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--test", action='store_true', help="launch unit tests")
    parser.add_argument("--input_yaml_file", type=str, help="path to input yaml file") 
    parser.add_argument("--verbose", type=int, help="verbose level overload, default to the one contained in the yaml file")
    args = parser.parse_args()

    if args.test:
        test_main_program()
    else:
        assert args.input_yaml_file is not None
        hyfaa_preprocessing_assimilation(args.input_yaml_file, verbose=args.verbose)


